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
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, "/scratch/ctaylor/KBUtils_Local/src")
from kbutillib import AICurationUtils, MSBiochemUtils
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=6, help="0 = all reactions")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--out", default="data_out_dirtest", help="output dir under WD")
    ap.add_argument("--checkpoint-every", type=int, default=50)
    ap.add_argument("--source", choices=["cache", "all"], default="cache",
                    help="'cache' = reaction ids from the existing cache file; "
                         "'all' = every buildable reaction in the ModelSEED DB")
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

    # Make the real method parallel-safe: each call works on its own throwaway
    # dict so worker threads never mutate shared_cache; the method returns the
    # result, which the main thread stores into shared_cache under a lock.
    driver._load_cached_curation = lambda name: {}
    driver._save_cached_curation = lambda name, cache: None
    # Robust JSON: strip markdown fences before the method parses.
    _orig_chat = driver.chat
    driver.chat = lambda prompt, system="": extract_json_str(_orig_chat(prompt, system=system))

    lock = threading.Lock()
    done = [len(shared_cache)]
    failures = []

    def checkpoint():
        tmp = out_file + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(shared_cache, fh, indent=4, skipkeys=True)
        os.replace(tmp, out_file)

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
