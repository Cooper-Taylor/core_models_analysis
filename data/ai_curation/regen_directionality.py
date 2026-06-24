"""Regenerate AICurationCacheReactionDirectionality.json with Claude Opus 4.8 (via Argo).

Uses the REAL AICurationUtils.analyze_reaction_directionality method (faithful
prompt / JSON parsing / 'reversed' flip / cache schema). Reactions are rebuilt as
cobra objects from the local ModelSEED database. The reaction set is taken from the
keys of the existing cache file, so we regenerate exactly the same reactions.

Parallelism: the method reloads/saves the cache each call, which isn't safe for
threads, so we point _load_cached_curation at one shared dict and no-op
_save_cached_curation; the main thread checkpoints that dict to disk. The method
also skips ids already present -> the run is resumable.

Usage:
  python regen_directionality.py --limit 6  --workers 3   --out data_out_dirtest
  python regen_directionality.py --limit 0  --workers 10  --out regen        # full
"""

import argparse
import glob
import json
import os
import sys
import threading
import time
import types
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

import httpx

sys.path.insert(0, "/scratch/ctaylor/KBUtils_Local/src")
from kbutillib import AICurationUtils, MSBiochemUtils
from kbutillib.ai_curation_utils import scrub_sensitive_terms
import cobra

WD = "/scratch/ctaylor/ai_curation_argo"
SRC_CACHE = os.path.join(WD, "data", "AICurationCacheReactionDirectionality.json")
MSDB = "/scratch/ctaylor/ModelSEEDDatabase"


def make_driver(**kwargs):
    """Plain AICurationUtils + the two pure MSBiochemUtils methods bound on.

    MSBiochemUtils.__init__ requires modelseedpy (to load the DB), which isn't
    installed; but reaction_to_string / _parse_id are pure (regex + cobra only),
    so we bind them onto an AICurationUtils instance and skip that __init__.
    """
    driver = AICurationUtils(**kwargs)
    driver.reaction_to_string = types.MethodType(MSBiochemUtils.reaction_to_string, driver)
    driver._parse_id = types.MethodType(MSBiochemUtils._parse_id, driver)
    return driver


def load_msdb_index():
    idx = {}
    for f in sorted(glob.glob(os.path.join(MSDB, "Biochemistry", "reaction_*.json"))):
        for r in json.load(open(f)):
            idx[r["id"]] = r
    return idx


def build_reaction(rec):
    """Build a reversible cobra Reaction from an MSDB reaction record."""
    rxn = cobra.Reaction(rec["id"])
    rxn.name = rec.get("name") or rec["id"]
    rxn.bounds = (-1000, 1000)  # reversible -> build_reaction_string uses '<=>'
    mets = {}
    for s in rec["stoichiometry"]:
        mid = f'{s["compound"]}_{s.get("compartment", 0)}'
        if mid not in mets:
            m = cobra.Metabolite(mid, name=s.get("name") or s["compound"],
                                 compartment=str(s.get("compartment", 0)))
            mets[mid] = (m, 0.0)
        m, c = mets[mid]
        mets[mid] = (m, c + float(s["coefficient"]))
    rxn.add_metabolites({m: c for (m, c) in mets.values() if c != 0})
    return rxn


def extract_json_str(txt: str) -> str:
    """Return the first balanced JSON object as a string.

    Handles ```json ... ``` fences AND a JSON object followed by trailing prose
    (the 'Extra data' case), which plain json.loads rejects.
    """
    t = txt.strip()
    if t.startswith("```"):
        nl = t.find("\n")
        if nl != -1:
            t = t[nl + 1:]
        if t.rstrip().endswith("```"):
            t = t.rstrip()[:-3]
        t = t.strip()
    start = t.find("{")
    if start == -1:
        return t
    depth = 0
    in_str = False
    esc = False
    for i in range(start, len(t)):
        c = t[i]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
        else:
            if c == '"':
                in_str = True
            elif c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    return t[start:i + 1]
    return t[start:]


# --------------------------------------------------------------------------
# Fill-gaps: retry only the missing reactions, looping on timeouts.
# Uses the native Anthropic messages endpoint so refusals are detectable
# (stop_reason == "refusal") and cleanly separated from timeouts.
# --------------------------------------------------------------------------
ARGO_MESSAGES_URL = "https://apps.inside.anl.gov/argoapi/v1/messages"

# Verbatim from AICurationUtils.analyze_reaction_directionality so fill-gap
# results use the same prompting as the main run.
DIR_SYSTEM = """
        You are an expert in biochemistry and molecular biology.
        You will receive a biochemical reaction and must evaluate it for stoichiometric
        correctness and biological directionality.

        Respond strictly in valid JSON with **no text outside the JSON**.
        All keys and string values must use double quotes.
        Use only plain ASCII characters.
        """
DIR_SHARED_PROMPT = """Analyze the following reaction for stoichiometric correctness and
        directionality in vivo.

        Return a JSON object in this exact format:

        {
        "errors": ["error 1", "error 2"],
        "directionality": "forward|reverse|reversible|uncertain",
        "other_comments": "Brief general comments about the reaction so I know you understood the input.",
        "confidence": "high|medium|low|none"
        }

        Reaction:
        """


def _fill_attempt(cli, user, model, rid, msdb, driver):
    """Attempt one missing reaction via the native Anthropic messages endpoint.

    Returns ``(rid, klass, result)``. ``klass`` is one of
    ``ok`` | ``refusal`` | ``nonmsdb`` | ``other`` | ``timeout``.
    Only ``timeout`` (network error / 5xx / malformed body) is retryable;
    the rest are terminal.
    """
    rec = msdb.get(rid)
    if rec is None:
        return rid, "nonmsdb", None
    try:
        rxnstring = scrub_sensitive_terms(
            driver.reaction_to_string(build_reaction(rec))["rxnstring"])
    except Exception:
        return rid, "other", None
    body = {"model": model, "max_tokens": 2048, "system": DIR_SYSTEM,
            "messages": [{"role": "user", "content": DIR_SHARED_PROMPT + rxnstring}]}
    hdr = {"x-api-key": user, "anthropic-version": "2023-06-01",
           "content-type": "application/json"}
    try:
        r = cli.post(ARGO_MESSAGES_URL, json=body, headers=hdr)
    except httpx.HTTPError:
        return rid, "timeout", None        # network/timeout -> retry
    if r.status_code >= 500:
        return rid, "timeout", None        # gateway hiccup -> retry
    if r.status_code != 200:
        return rid, "other", None
    try:
        j = r.json()
    except Exception:
        return rid, "timeout", None        # partial/malformed body -> retry
    stop = j.get("stop_reason")
    txt = "".join(b.get("text", "") for b in j.get("content", []) if b.get("type") == "text")
    if stop == "refusal" or not txt.strip():
        return rid, "refusal", None        # terminal
    try:
        obj = json.loads(extract_json_str(txt))
    except Exception:
        return rid, "other", None          # terminal
    if not isinstance(obj, dict) or "directionality" not in obj:
        return rid, "other", None
    return rid, "ok", obj


def fill_gaps(driver, msdb, missing, shared_cache, checkpoint, workers=8,
              max_rounds=20, patience=3, refusal_retries=0):
    """Continually re-run missing reactions until each succeeds or hits a
    terminal error.

    Timeouts (network / 5xx / malformed body) are retried across rounds.
    Refusals are retried up to ``refusal_retries`` extra times (the scrub makes
    refusals nondeterministic, so a borderline reaction may classify on a later
    attempt); after ``refusal_retries + 1`` refusals it is terminal. non-MSDB and
    other parse errors are terminal immediately.

    Stops when nothing is pending, after ``max_rounds``, or after ``patience``
    consecutive rounds in which *every* attempt timed out (Argo down) -- so it
    never spins forever.
    """
    user, model = driver.user, driver.model
    pending = set(missing)
    refusals = Counter()
    terminal = {"refusal": [], "nonmsdb": [], "other": []}
    cli = httpx.Client(timeout=120.0)
    no_progress = 0
    rnd = 0
    try:
        while pending and rnd < max_rounds:
            rnd += 1
            cur = sorted(pending)
            next_pending = set()
            n_ok = n_ref_retry = n_ref_term = n_timeout = n_other = 0
            with ThreadPoolExecutor(max_workers=workers) as ex:
                futs = {ex.submit(_fill_attempt, cli, user, model, rid, msdb, driver): rid
                        for rid in cur}
                for fut in as_completed(futs):
                    rid, klass, result = fut.result()
                    if klass == "ok":
                        shared_cache[rid] = result
                        n_ok += 1
                    elif klass == "timeout":
                        next_pending.add(rid)
                        n_timeout += 1
                    elif klass == "refusal":
                        refusals[rid] += 1
                        if refusals[rid] > refusal_retries:   # 1 + refusal_retries attempts
                            terminal["refusal"].append(rid)
                            n_ref_term += 1
                        else:
                            next_pending.add(rid)
                            n_ref_retry += 1
                    else:  # nonmsdb / other
                        terminal[klass].append(rid)
                        n_other += 1
            checkpoint()
            print(f"[fill] round {rnd}: ok {n_ok} | refusal(retry {n_ref_retry}/give-up {n_ref_term}) "
                  f"| timeout {n_timeout} | other-terminal {n_other} | pending {len(next_pending)} "
                  f"| total {len(shared_cache)}", flush=True)
            # "Progress" = anything other than pure timeout re-queue. Only an
            # all-timeout round counts against patience (detects Argo outage).
            progress = (n_ok + n_ref_retry + n_ref_term + n_other) > 0
            no_progress = 0 if progress else (no_progress + 1)
            pending = next_pending
            if no_progress >= patience:
                print(f"[fill] no progress (all timeouts) for {patience} rounds; stopping "
                      f"with {len(pending)} pending", flush=True)
                break
    finally:
        cli.close()
    checkpoint()
    print(f"[fill] DONE after {rnd} rounds | total now {len(shared_cache)}")
    print(f"[fill] still-pending/timing-out ({len(pending)}): {sorted(pending)}")
    for k in ("refusal", "nonmsdb", "other"):
        print(f"[fill] terminal {k} ({len(terminal[k])}): {sorted(terminal[k])}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=6, help="0 = all reactions")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--out", default="data_out_dirtest", help="output dir under WD")
    ap.add_argument("--checkpoint-every", type=int, default=50)
    ap.add_argument("--source", choices=["cache", "all"], default="cache",
                    help="'cache' = reaction ids from the existing cache file; "
                         "'all' = every buildable reaction in the ModelSEED DB")
    ap.add_argument("--fill-gaps", action="store_true",
                    help="retry ONLY the missing reactions via the native Anthropic "
                         "messages endpoint, looping on timeouts until each succeeds or "
                         "hits a terminal error (refusal / non-MSDB / other)")
    ap.add_argument("--max-rounds", type=int, default=20,
                    help="fill-gaps: max retry rounds for timed-out reactions")
    ap.add_argument("--refusal-retries", type=int, default=0,
                    help="fill-gaps: extra attempts to give a REFUSED reaction before "
                         "treating it as terminal (scrubbed refusals are nondeterministic). "
                         "0 = refusal is terminal on first occurrence (default).")
    args = ap.parse_args()

    out_dir = os.path.join(WD, args.out)
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "AICurationCacheReactionDirectionality.json")

    print("[regen] loading ModelSEED reaction index ...")
    msdb = load_msdb_index()
    print(f"[regen] MSDB reactions indexed: {len(msdb)}")

    if args.source == "all":
        # Every buildable MSDB reaction (has stoichiometry); includes obsolete.
        ids = sorted(rid for rid, rec in msdb.items() if rec.get("stoichiometry"))
    else:
        ids = list(json.load(open(SRC_CACHE)).keys())
    if args.limit:
        ids = ids[: args.limit]
    print(f"[regen] target reactions: {len(ids)} | source: {args.source} "
          f"| workers: {args.workers} | out: {out_file}")

    driver = make_driver(
        name="DirectionalityRegen", backend="argo", proxy_port=None,
        config_file=False, token_file=None, kbase_token_file=None,
    )
    driver.data_directory = out_dir

    # Resume: seed shared cache from any existing output.
    shared_cache = {}
    if os.path.exists(out_file):
        shared_cache = json.load(open(out_file))
        print(f"[regen] resuming; {len(shared_cache)} already done")

    lock = threading.Lock()

    def checkpoint():
        tmp = out_file + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(shared_cache, fh, indent=4, skipkeys=True)
        os.replace(tmp, out_file)

    # ---- fill-gaps mode: retry ONLY missing reactions, looping on timeouts ----
    if args.fill_gaps:
        missing = [r for r in ids if r not in shared_cache]
        print(f"[fill] {len(missing)} missing reactions to attempt "
              f"(max_rounds={args.max_rounds}, workers={args.workers})", flush=True)
        fill_gaps(driver, msdb, missing, shared_cache, checkpoint,
                  workers=args.workers, max_rounds=args.max_rounds,
                  refusal_retries=args.refusal_retries)
        return

    # ---- normal generate mode -------------------------------------------------
    # Make the real method parallel-safe: each call works on its own throwaway
    # dict so worker threads never mutate shared_cache; the method returns the
    # result, which the main thread stores into shared_cache under a lock.
    driver._load_cached_curation = lambda name: {}
    driver._save_cached_curation = lambda name, cache: None
    # Robust JSON: strip markdown fences before the method parses.
    _orig_chat = driver.chat
    driver.chat = lambda prompt, system="": extract_json_str(_orig_chat(prompt, system=system))

    done = [len(shared_cache)]
    failures = []

    def work(rid):
        rec = msdb.get(rid)
        if rec is None:
            return rid, "missing_in_msdb", None
        try:
            rxn = build_reaction(rec)
            result = driver.analyze_reaction_directionality(rxn)
            if result is None:
                return rid, "skipped_prefix", None
            return rid, "ok", result
        except Exception as e:
            return rid, f"ERROR:{type(e).__name__}:{str(e)[:100]}", None

    t0 = time.time()
    todo = [r for r in ids if r not in shared_cache]
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(work, r): r for r in todo}
        for fut in as_completed(futs):
            rid, status, result = fut.result()
            with lock:
                if status == "ok":
                    shared_cache[rid] = result
                    done[0] += 1
                else:
                    failures.append((rid, status))
                    print(f"  [{rid}] {status}", flush=True)
                n = done[0]
                if n and n % args.checkpoint_every == 0:
                    checkpoint()
                    el = time.time() - t0
                    print(f"  ... {n}/{len(ids)} done ({el:.0f}s elapsed)", flush=True)

    checkpoint()
    el = time.time() - t0
    print(f"[regen] FINISHED: {len(shared_cache)} entries in {out_file}")
    print(f"[regen] elapsed: {el:.0f}s | failures: {len(failures)}")
    for rid, st in failures[:20]:
        print("   fail:", rid, st)


if __name__ == "__main__":
    main()
