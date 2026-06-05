#!/usr/bin/env python3
"""Tiny http.server-based backend for the heuristics-explorer website.

Serves:
  /             -> static/index.html
  /static/*     -> static files (HTML/JS/CSS)
  /data/*       -> JSON snapshots written by build_site_data.py
  /api/health   -> {"ok": true}
  /api/rxn/{id} -> per-reaction lookup (panel-first, fallback to compact non-panel index)
  /api/reaction_impact (POST)
       body: {"rxn_id": "rxnXXXXX",
              "modes": ["off","forward","reverse"],
              "variant": "baseline",
              "models": [optional subset of panel ids]}
       returns: {"baseline": {model_id: {grows, growth_flux}},
                 "by_mode":  {mode: {model_id: {grows, growth_flux, delta_flux}}}}

  /api/panel_fba (POST)
       body: {"variant": "baseline|<tag>",
              "overrides": {"rxnXXXXX": "off|forward|reverse"} | {},
              "models": [optional subset]}
       returns: [{model_id, grows, growth_flux, status}, ...]

Stdlib only -- no Flask / FastAPI dependency.  multiprocessing-based
panel FBA via growth_heuristics.run_panel.
"""

from __future__ import annotations

import argparse
import http.server
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

sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(ANALYSIS_ROOT.parent / "ModelSEEDDatabase" / "Libs" / "Python"))

import growth_heuristics as gh  # noqa: E402


# ---------------------------------------------------------------------------
# Lazy-loaded shared state
# ---------------------------------------------------------------------------
_STATE_LOCK = threading.Lock()
_STATE: dict = {}


def _state() -> dict:
    if not _STATE:
        with _STATE_LOCK:
            if not _STATE:
                print("[state] loading site data + per-variant maps...", flush=True)
                t = time.time()
                baseline = json.loads((DATA_ROOT / "baseline.json").read_text())
                _STATE["baseline_map"] = baseline["map"]
                variant_maps = {"baseline": dict(baseline["map"])}
                for vfile in (DATA_ROOT / "variants").glob("*.json"):
                    tag = vfile.stem
                    payload = json.loads(vfile.read_text())
                    # Build variant_map from baseline + per-rxn diffs
                    vmap = dict(baseline["map"])
                    for d in payload["diffs"]:
                        vmap[d["rxn"]] = d["new"]
                    variant_maps[tag] = vmap
                _STATE["variant_maps"] = variant_maps
                _STATE["panel_ids"] = (ANALYSIS_ROOT / "results" /
                                       "selected_ids.txt").read_text().split()
                _STATE["reactions_panel"] = json.loads(
                    (DATA_ROOT / "reactions_panel.json").read_text())
                _STATE["manifest"] = json.loads(
                    (DATA_ROOT / "manifest.json").read_text())
                print(f"[state] loaded {len(variant_maps)} variants, "
                      f"{len(_STATE['panel_ids'])} panel models, "
                      f"{len(_STATE['reactions_panel'])} panel reactions "
                      f"in {time.time()-t:.1f}s", flush=True)
    return _STATE


def get_reactions_other() -> dict:
    """Loaded on first request only -- 4MB lazy load."""
    if "reactions_other" not in _STATE:
        with _STATE_LOCK:
            if "reactions_other" not in _STATE:
                t = time.time()
                _STATE["reactions_other"] = json.loads(
                    (DATA_ROOT / "reactions_other.json").read_text())
                print(f"[state] loaded reactions_other ({len(_STATE['reactions_other'])} "
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


# Custom bounds helper: extend gh._bounds_for_rev to support "off".
# Keep a reference to the original to avoid self-recursion after monkey-patch.
_GH_BOUNDS_ORIG = gh._bounds_for_rev


def _bounds_for_rev_ext(rev: str, default_bound: float = 1000.0):
    if rev == "off":
        return 0.0, 0.0
    return _GH_BOUNDS_ORIG(rev, default_bound)


# Monkey-patch growth_heuristics so the rev "off" works inside worker procs.
gh._bounds_for_rev = _bounds_for_rev_ext


def run_panel_fba(variant: str, overrides: dict,
                  model_ids: Optional[list] = None,
                  n_workers: int = 4) -> list:
    """Drive gh.run_panel with the requested variant + per-rxn overrides."""
    state = _state()
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
                self._send_json({"ok": True, "n_variants": len(_state()["variant_maps"])})
                return
            if path.startswith("/api/rxn/"):
                rxn_id = path[len("/api/rxn/"):]
                state = _state()
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
                state = _state()
                if models is None:
                    # Only models that actually contain rxn_id are interesting.
                    panel_rxnsets = json.loads(
                        (DATA_ROOT / "panel_rxnsets.json").read_text())
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


# ---------------------------------------------------------------------------
def main(argv: Optional[list] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--preload", action="store_true",
                    help="load all state up-front instead of on first request")
    args = ap.parse_args(argv)

    if args.preload:
        _state()
        get_reactions_other()

    srv = ThreadedHTTPServer((args.host, args.port), Handler)
    actual_host = args.host if args.host != "0.0.0.0" else socket.gethostname()
    print(f"\n  Heuristics Explorer running at: http://{actual_host}:{args.port}/")
    print(f"  Static root : {STATIC_ROOT}")
    print(f"  Data root   : {DATA_ROOT}")
    print(f"  Ctrl-C to stop.\n", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        srv.shutdown()


if __name__ == "__main__":
    main()
