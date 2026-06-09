#!/usr/bin/env python3
"""Tiny http.server-based backend for the heuristics-explorer website — auto-detects static vs live FBA mode.

By default this runs in **static mode**: it serves the static UI and the
JSON snapshots shipped under ``site/data/``. No external dependencies —
stdlib only, no pip install required. Use this for "give someone a fresh
clone and let them browse the variants" workflows.

Pass ``--live`` to enable **live FBA mode**, which exposes
``/api/panel_fba`` and ``/api/reaction_impact`` and runs cobra FBA on
the 100-model panel. Live mode requires the full upstream pipeline:

  * ``cobra`` importable (``pip install cobra``)
  * ``ModelSEEDDatabase`` cloned as a sibling of ``core_models_analysis/``
  * ``notebooks/.kbcache/`` populated (built by notebook 06)
  * ``results/selected_ids.txt``
  * ``site/data/baseline.json`` and ``site/data/panel_rxnsets.json``
    (built by ``scripts/build_site_data.py``)

Pass ``--static`` to force static mode even when those prerequisites
are present (useful for tests / CI).

Routes:
  /             -> static/index.html
  /static/*     -> static files (HTML/JS/CSS)
  /data/*       -> JSON snapshots written by build_site_data.py
  /api/health   -> {"ok": true, "static_mode": bool, "n_variants": int|null}
  /api/rxn/{id} -> per-reaction lookup (panel-first, fallback to compact non-panel index)
  /api/reaction_impact (POST)  [live mode only — 503 in static]
       body: {"rxn_id": "rxnXXXXX",
              "modes": ["off","forward","reverse"],
              "variant": "baseline",
              "models": [optional subset of panel ids]}
       returns: {"baseline": {model_id: {grows, growth_flux}},
                 "by_mode":  {mode: {model_id: {grows, growth_flux, delta_flux}}}}
  /api/panel_fba (POST)  [live mode only — 503 in static]
       body: {"variant": "baseline|<tag>",
              "overrides": {"rxnXXXXX": "off|forward|reverse"} | {},
              "models": [optional subset]}
       returns: [{model_id, grows, growth_flux, status}, ...]

Stdlib only -- no Flask / FastAPI dependency.  multiprocessing-based
panel FBA via growth_heuristics.run_panel (live mode only).
"""

from __future__ import annotations

import argparse
import http.server
import importlib.util
import json
import os
import socket
import socketserver
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse, unquote

SITE_ROOT = Path(__file__).resolve().parent
ANALYSIS_ROOT = SITE_ROOT.parent
SCRIPTS = ANALYSIS_ROOT / "scripts"
DATA_ROOT = SITE_ROOT / "data"
STATIC_ROOT = SITE_ROOT / "static"


# ---------------------------------------------------------------------------
# Mode + lazy FBA runtime
# ---------------------------------------------------------------------------
# STATIC_MODE is set in main() before serve. Never None at request-handling
# time. Default-to-True keeps any out-of-band import (tests, REPL) safe.
STATIC_MODE: Optional[bool] = True
_FBA_READY = False
_FBA_INIT_LOCK = threading.Lock()
gh = None  # filled by _init_fba_runtime in live mode
_GH_BOUNDS_ORIG = None


def _bounds_for_rev_ext(rev, default_bound=1000.0):
    if rev == "off":
        return 0.0, 0.0
    return _GH_BOUNDS_ORIG(rev, default_bound)


_bounds_for_rev_ext._serve_patched = True


def _detect_static_mode() -> "tuple[bool, str]":
    """Return (is_static, reason). True if any FBA prerequisite is missing.

    Predicate order matches _fba_state() reads, so auto-detect can never
    pick live mode and then crash on first FBA call.
    """
    try:
        have_cobra = importlib.util.find_spec("cobra") is not None
    except (ImportError, ValueError):
        have_cobra = False
    if not have_cobra:
        return True, "cobra not importable"
    msdb_libs = ANALYSIS_ROOT.parent / "ModelSEEDDatabase" / "Libs" / "Python"
    if not msdb_libs.is_dir():
        return True, f"MSDB Libs/Python dir missing at {msdb_libs}"
    for path in (
        ANALYSIS_ROOT / "results" / "selected_ids.txt",
        ANALYSIS_ROOT / "notebooks" / ".kbcache",
        DATA_ROOT / "baseline.json",
        DATA_ROOT / "panel_rxnsets.json",
    ):
        if not path.exists():
            return True, f"missing {path}"
    return False, "all FBA prerequisites present"


def _init_fba_runtime():
    """Import growth_heuristics + install the 'off' monkey-patch. Lazy + idempotent."""
    global _FBA_READY, gh, _GH_BOUNDS_ORIG
    if _FBA_READY:
        return
    with _FBA_INIT_LOCK:
        if _FBA_READY:
            return
        if STATIC_MODE:
            raise RuntimeError(
                "FBA disabled in static mode — restart with "
                "`python3 site/serve.py --live` after installing cobra "
                "+ ModelSEEDDatabase"
            )
        sys.path.insert(0, str(SCRIPTS))
        sys.path.insert(0, str(ANALYSIS_ROOT.parent / "ModelSEEDDatabase" / "Libs" / "Python"))
        import growth_heuristics as _gh_mod  # noqa: E402
        gh = _gh_mod
        # Idempotency on re-entry from a separate process / re-import path.
        if getattr(gh._bounds_for_rev, "_serve_patched", False):
            _FBA_READY = True
            return
        _GH_BOUNDS_ORIG = gh._bounds_for_rev
        gh._bounds_for_rev = _bounds_for_rev_ext
        _FBA_READY = True
        print("[fba] runtime initialized", flush=True)


# ---------------------------------------------------------------------------
# Lazy-loaded shared state.
#
# Two namespaces on the shared _STATE dict (one lock):
#   STATIC keys: manifest, reactions_panel, reactions_other
#   FBA keys   : baseline_map, variant_maps, panel_ids, panel_rxnsets
# ---------------------------------------------------------------------------
_STATE_LOCK = threading.Lock()
_STATE: dict = {}


def _static_state() -> dict:
    """Load static-mode data (manifest + reactions_panel). Safe in both modes."""
    if "manifest" in _STATE and "reactions_panel" in _STATE:
        return _STATE
    with _STATE_LOCK:
        if "manifest" not in _STATE:
            _STATE["manifest"] = json.loads((DATA_ROOT / "manifest.json").read_text())
        if "reactions_panel" not in _STATE:
            t = time.time()
            _STATE["reactions_panel"] = json.loads(
                (DATA_ROOT / "reactions_panel.json").read_text())
            print(f"[static] loaded reactions_panel "
                  f"({len(_STATE['reactions_panel'])} entries) "
                  f"in {time.time()-t:.1f}s", flush=True)
    return _STATE


def _fba_state() -> dict:
    """Load FBA-mode data. Raises RuntimeError if STATIC_MODE."""
    _init_fba_runtime()
    if all(k in _STATE for k in ("baseline_map", "variant_maps", "panel_ids", "panel_rxnsets")):
        return _STATE
    with _STATE_LOCK:
        if "baseline_map" not in _STATE:
            print("[fba] loading baseline + per-variant maps...", flush=True)
            t = time.time()
            baseline = json.loads((DATA_ROOT / "baseline.json").read_text())
            _STATE["baseline_map"] = baseline["map"]
            variant_maps = {"baseline": dict(baseline["map"])}
            for vfile in (DATA_ROOT / "variants").glob("*.json"):
                tag = vfile.stem
                payload = json.loads(vfile.read_text())
                vmap = dict(baseline["map"])
                for d in payload["diffs"]:
                    vmap[d["rxn"]] = d["new"]
                variant_maps[tag] = vmap
            _STATE["variant_maps"] = variant_maps
            _STATE["panel_ids"] = (ANALYSIS_ROOT / "results" /
                                   "selected_ids.txt").read_text().split()
            _STATE["panel_rxnsets"] = json.loads(
                (DATA_ROOT / "panel_rxnsets.json").read_text())
            print(f"[fba] loaded {len(variant_maps)} variants, "
                  f"{len(_STATE['panel_ids'])} panel models "
                  f"in {time.time()-t:.1f}s", flush=True)
    return _STATE


def get_reactions_other() -> dict:
    """Loaded on first request only -- 4MB lazy load. Static-safe."""
    if "reactions_other" not in _STATE:
        with _STATE_LOCK:
            if "reactions_other" not in _STATE:
                t = time.time()
                _STATE["reactions_other"] = json.loads(
                    (DATA_ROOT / "reactions_other.json").read_text())
                print(f"[static] loaded reactions_other ({len(_STATE['reactions_other'])} "
                      f"entries) in {time.time()-t:.1f}s", flush=True)
    return _STATE["reactions_other"]


# ---------------------------------------------------------------------------
# Override modes -- map a 1-letter "mode" to the bounds the model gets.
# ---------------------------------------------------------------------------
MODE_TO_REV = {
    # Use the cascade's rev for this rxn.  Special: handled by leaving the map
    # entry as-is (no override).
    "as_is":   None,
    # Force forward / reverse / reversible / off.
    "forward": ">",
    "reverse": "<",
    "free":    "=",
    "off":     "off",
}


def _apply_overrides(rev_map: dict, overrides: dict) -> dict:
    """Return a copy of ``rev_map`` with per-rxn overrides applied.

    Bounds for ``"off"`` get coded as the new sentinel ``"X"``; we handle that
    in the gh.run_panel call by post-processing the cobra model below.
    """
    out = dict(rev_map)
    for rxn, mode in overrides.items():
        if mode == "as_is":
            continue
        rev = MODE_TO_REV.get(mode)
        if rev is None:
            raise ValueError(f"unknown mode for {rxn}: {mode}")
        out[rxn] = rev
    return out


def run_panel_fba(variant: str, overrides: dict,
                  model_ids: Optional[list] = None,
                  n_workers: int = 4) -> list:
    """Drive gh.run_panel with the requested variant + per-rxn overrides."""
    _init_fba_runtime()
    state = _fba_state()
    base_map = state["variant_maps"].get(variant)
    if base_map is None:
        raise KeyError(f"unknown variant: {variant}")
    eff_map = _apply_overrides(base_map, overrides or {})
    mids = model_ids or state["panel_ids"]
    return gh.run_panel(mids, reversibility_map=eff_map,
                        baseline_map=None, n_workers=n_workers)


# ---------------------------------------------------------------------------
# HTTP handler
# ---------------------------------------------------------------------------
class Handler(http.server.SimpleHTTPRequestHandler):
    # Quieter logs.
    def log_message(self, fmt, *args):
        sys.stderr.write("%s - - %s\n" % (self.address_string(), fmt % args))

    # --- helpers -----------------------------------------------------------
    def _send_json(self, payload, code=200):
        body = json.dumps(payload, default=str).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, path: Path, content_type: Optional[str] = None):
        if not path.exists() or not path.is_file():
            self._send_json({"error": "not found", "path": str(path)}, code=404)
            return
        body = path.read_bytes()
        self.send_response(200)
        if content_type:
            self.send_header("Content-Type", content_type)
        elif path.suffix == ".html":
            self.send_header("Content-Type", "text/html; charset=utf-8")
        elif path.suffix == ".js":
            self.send_header("Content-Type", "application/javascript; charset=utf-8")
        elif path.suffix == ".css":
            self.send_header("Content-Type", "text/css; charset=utf-8")
        elif path.suffix == ".json":
            self.send_header("Content-Type", "application/json")
        elif path.suffix == ".svg":
            self.send_header("Content-Type", "image/svg+xml")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_json_body(self) -> dict:
        n = int(self.headers.get("Content-Length", "0"))
        if not n:
            return {}
        raw = self.rfile.read(n)
        try:
            return json.loads(raw.decode())
        except Exception as exc:
            raise ValueError(f"invalid JSON body: {exc}")

    # --- GET ---------------------------------------------------------------
    def do_GET(self):
        url = urlparse(self.path)
        path = unquote(url.path)
        try:
            if path == "/" or path == "/index.html":
                self._send_file(STATIC_ROOT / "index.html")
                return
            if path.startswith("/static/"):
                self._send_file(STATIC_ROOT / path[len("/static/"):])
                return
            if path.startswith("/data/"):
                self._send_file(DATA_ROOT / path[len("/data/"):])
                return
            if path == "/api/health":
                n_var = None
                try:
                    m = json.loads((DATA_ROOT / "manifest.json").read_text())
                    n_var = len(m.get("variants", []))
                except Exception:
                    pass
                self._send_json({
                    "ok": True,
                    "static_mode": bool(STATIC_MODE),
                    "n_variants": n_var,
                })
                return
            if path.startswith("/api/rxn/"):
                rxn_id = path[len("/api/rxn/"):]
                state = _static_state()
                r = state["reactions_panel"].get(rxn_id)
                if r is None:
                    r = get_reactions_other().get(rxn_id)
                if r is None:
                    self._send_json({"error": "not found", "rxn_id": rxn_id}, code=404)
                else:
                    self._send_json(r)
                return
            self._send_json({"error": "not found", "path": path}, code=404)
        except Exception as exc:
            traceback.print_exc()
            self._send_json({"error": str(exc), "type": type(exc).__name__}, code=500)

    # --- POST --------------------------------------------------------------
    def do_POST(self):
        url = urlparse(self.path)
        path = unquote(url.path)
        try:
            if STATIC_MODE and path in ("/api/panel_fba", "/api/reaction_impact"):
                self._send_json({
                    "error": "FBA disabled in static mode — restart with "
                             "`python3 site/serve.py --live` after installing "
                             "cobra + ModelSEEDDatabase",
                    "static_mode": True,
                }, code=503)
                return
            body = self._read_json_body()
            if path == "/api/panel_fba":
                variant = body.get("variant", "baseline")
                overrides = body.get("overrides") or {}
                models = body.get("models") or None
                workers = int(body.get("n_workers", 4))
                t = time.time()
                results = run_panel_fba(variant, overrides, models, workers)
                self._send_json({
                    "variant": variant,
                    "n_overrides": len(overrides),
                    "n_models": len(results),
                    "elapsed_s": round(time.time() - t, 2),
                    "results": results,
                })
                return
            if path == "/api/reaction_impact":
                rxn_id = body["rxn_id"]
                modes = body.get("modes") or ["off", "forward", "reverse"]
                variant = body.get("variant", "baseline")
                models = body.get("models") or None
                state = _fba_state()
                if models is None:
                    # Only models that actually contain rxn_id are interesting.
                    panel_rxnsets = state["panel_rxnsets"]
                    interesting = [mid for mid in state["panel_ids"]
                                   if rxn_id in panel_rxnsets.get(mid, [])]
                    models = interesting or state["panel_ids"]
                t = time.time()
                # Baseline (no override) under the chosen variant.
                base = run_panel_fba(variant, {}, models, n_workers=4)
                base_by_id = {r["model_id"]: r for r in base}
                out = {
                    "rxn_id": rxn_id,
                    "variant": variant,
                    "n_models": len(models),
                    "baseline": {r["model_id"]: {
                        "grows": r["grows"], "growth_flux": r["growth_flux"],
                        "status": r["status"],
                    } for r in base},
                    "by_mode": {},
                }
                for mode in modes:
                    res = run_panel_fba(variant, {rxn_id: mode}, models, n_workers=4)
                    by_id = {}
                    for r in res:
                        b = base_by_id.get(r["model_id"], {})
                        by_id[r["model_id"]] = {
                            "grows": r["grows"], "growth_flux": r["growth_flux"],
                            "status": r["status"],
                            "delta_flux": r["growth_flux"] - float(b.get("growth_flux", 0.0)),
                            "grew_before": bool(b.get("grows", False)),
                        }
                    out["by_mode"][mode] = by_id
                out["elapsed_s"] = round(time.time() - t, 2)
                self._send_json(out)
                return
            self._send_json({"error": "not found", "path": path}, code=404)
        except Exception as exc:
            traceback.print_exc()
            self._send_json({"error": str(exc), "type": type(exc).__name__}, code=500)


# Threaded so the FBA endpoint doesn't block static requests.
class ThreadedHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True
    allow_reuse_address = True


def _bind(host: str, port_pref: int, explicit: bool):
    """Bind ThreadedHTTPServer. If port_pref is default, try +0..+10. If
    explicit, fail fast. Returns (srv, bound_port). Raises OSError on
    exhaustion (caller decides the exit message)."""
    tries = [port_pref] if explicit else list(range(port_pref, port_pref + 11))
    last_exc = None
    for p in tries:
        try:
            srv = ThreadedHTTPServer((host, p), Handler)
            return srv, p
        except OSError as e:
            last_exc = e
            continue
    raise last_exc


# ---------------------------------------------------------------------------
def main(argv: Optional[list] = None) -> None:
    if sys.version_info < (3, 9):
        sys.exit("error: site/serve.py requires Python 3.9+ "
                 "(got %d.%d)" % sys.version_info[:2])

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765,
                    help="default 8765; if taken (and not explicitly set), tries 8766..8775")
    ap.add_argument("--preload", action="store_true",
                    help="load all state up-front instead of on first request")
    mode_group = ap.add_mutually_exclusive_group()
    mode_group.add_argument("--static", action="store_true",
                            help="force static mode (no FBA endpoints)")
    mode_group.add_argument("--live", action="store_true",
                            help="force live FBA mode (requires cobra + MSDB + data pipeline)")
    args = ap.parse_args(argv)

    # Mode resolution: default is static even if all prereqs present.
    global STATIC_MODE
    is_static_auto, reason = _detect_static_mode()
    if args.static:
        STATIC_MODE = True
        mode_reason = "forced via --static"
    elif args.live:
        if is_static_auto:
            sys.exit("error: --live requires all FBA prerequisites — "
                     "auto-detect says: " + reason)
        STATIC_MODE = False
        mode_reason = "forced via --live"
    else:
        STATIC_MODE = True
        mode_reason = ("auto: " + reason) if is_static_auto else \
            "default (prereqs present, but use --live to opt into FBA)"
    mode_label = "static (FBA disabled)" if STATIC_MODE else "live (FBA enabled)"
    print(f"[mode] {mode_label}: {mode_reason}", flush=True)

    if args.preload:
        _static_state()
        get_reactions_other()
        if not STATIC_MODE:
            _fba_state()

    explicit_port = (args.port != 8765)
    try:
        srv, bound_port = _bind(args.host, args.port, explicit=explicit_port)
    except OSError as e:
        sys.exit(f"error: could not bind {args.host}:{args.port}"
                 f"{'' if explicit_port else ' (and 10 fallbacks)'}: {e}")

    actual_host = args.host if args.host != "0.0.0.0" else socket.gethostname()
    print(f"\n  Heuristics Explorer running at: http://{actual_host}:{bound_port}/")
    print(f"  Mode        : {mode_label}")
    print(f"  Repo root   : {ANALYSIS_ROOT}")
    print(f"  Static root : {STATIC_ROOT}")
    print(f"  Data root   : {DATA_ROOT}")
    if STATIC_MODE:
        print(f"  Note        : FBA endpoints return 503; pass --live to enable")
    else:
        print(f"  Note        : first /api/panel_fba may take 5-30s")
    print(f"  Ctrl-C to stop.\n", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        srv.shutdown()


if __name__ == "__main__":
    main()
