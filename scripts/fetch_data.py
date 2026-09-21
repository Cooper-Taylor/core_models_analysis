#!/usr/bin/env python3
"""Fetch the external data this repo needs but deliberately does not vendor.

    python3 scripts/fetch_data.py --list          # what is needed, what you have
    python3 scripts/fetch_data.py msdb            # the ModelSEEDDatabase clone
    python3 scripts/fetch_data.py gs-models       # ModelSEED v2 genome-scale models
    python3 scripts/fetch_data.py --all --into ../

Three external inputs, in rough order of how much you need them:

``msdb``       ModelSEEDDatabase, a git clone. Required by almost everything:
               the biochemistry, the media, the per-source thermodynamics and
               the release tags the direction sets are built from.
``gs-models``  The 5,420 ModelSEED v2 genome-scale models, from two public
               KBase workspaces. No auth token needed. About 2 GB converted.
``core-models`` The 5,683 core model JSONs. Not redistributable from here; the
               ModelSEED v2 supplement publishes them as `core_models.tar.gz`.

Everything lands next to the repo by default, which is one of the locations
`cma.paths` searches, so no configuration is needed afterwards.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from cma import paths  # noqa: E402

MSDB_URL = "https://github.com/ModelSEED/ModelSEEDDatabase"
CORE_MODELS_URL = "https://bioseed.mcs.anl.gov/~fliu/modelseed2/core_models.tar.gz"

ITEMS = {
    "msdb": ("the ModelSEEDDatabase clone", "msdb"),
    "gs-models": ("the ModelSEED v2 genome-scale models", "gs_models"),
    "core-models": ("the 5,683 core model JSONs", "core_models"),
}


def have(key: str) -> Path | None:
    from cma.paths import _resolve  # noqa: PLC2701 -- same package

    try:
        return _resolve(key, required=False)
    except Exception:  # noqa: BLE001
        return None


def cmd_list() -> int:
    print("external data\n")
    rc = 0
    for name, (what, key) in ITEMS.items():
        p = have(key)
        print(f"  {name:<13} {'ok     ' + str(p) if p else 'MISSING -- ' + what}")
        if not p:
            rc = 1
    print("\nFetch with:  python3 scripts/fetch_data.py --all")
    return rc


def run(cmd: list[str], **kw) -> int:
    print("  $ " + " ".join(cmd), flush=True)
    return subprocess.call(cmd, **kw)


def fetch_msdb(into: Path, ref: str) -> int:
    dest = into / "ModelSEEDDatabase"
    if dest.exists():
        print(f"  {dest} already exists; fetching tags instead of cloning")
        return run(["git", "-C", str(dest), "fetch", "--tags", "--all"])
    into.mkdir(parents=True, exist_ok=True)
    rc = run(["git", "clone", MSDB_URL, str(dest)])
    if rc == 0 and ref:
        # a release tag is what the direction sets are built from
        run(["git", "-C", str(dest), "fetch", "--tags"])
    return rc


def fetch_gs_models(into: Path, jobs: int) -> int:
    """Delegate to the downloader that ships with the model set, if present."""
    existing = have("gs_models")
    root = existing or (into / "modelseed2_gs_models")
    dl = root / "scripts" / "download_ms2_models.py"
    if not dl.exists():
        print(f"  the downloader is not at {dl}.")
        print("  It ships alongside the model set. If you do not have it, the models")
        print("  are in the public KBase workspaces 155807 (glucose minimal) and")
        print("  155808 (auxotrophy) of https://narrative.kbase.us/#/orgs/ms2-manuscript")
        print("  and need no auth token.")
        return 1
    rc = run([paths.python_exe(), str(dl), "--workers", str(jobs)])
    if rc == 0:
        rc = run([paths.python_exe(), str(root / "scripts" / "convert_to_cobra.py"),
                  "--jobs", str(jobs)])
    return rc


def fetch_core_models(into: Path) -> int:
    dest = into / "core_models_kegg2"
    if dest.exists():
        print(f"  {dest} already exists")
        return 0
    dest.mkdir(parents=True, exist_ok=True)
    tar = dest / "core_models.tar.gz"
    print(f"  downloading {CORE_MODELS_URL}")
    rc = run(["curl", "-fL", "--retry", "3", "-o", str(tar), CORE_MODELS_URL])
    if rc != 0:
        print("  download failed. The tarball is Supplementary Data S2 of the")
        print("  ModelSEED v2 paper: https://bioseed.mcs.anl.gov/~fliu/modelseed2/")
        return rc
    return run(["tar", "xzf", str(tar), "-C", str(dest), "--strip-components=1"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("items", nargs="*", choices=[*ITEMS, []], default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--into", type=Path, default=None,
                    help="where to put things (default: beside the repo)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--ref", default="v2.0.1", help="MSDB release tag to ensure is fetched")
    args = ap.parse_args()

    if args.list or (not args.items and not args.all):
        return cmd_list()

    into = (args.into or paths.analysis_dir().parent).resolve()
    wanted = list(ITEMS) if args.all else args.items
    print(f"fetching into {into}\n")
    rc = 0
    for name in wanted:
        print(f"--- {name} ---")
        if name == "msdb":
            rc |= fetch_msdb(into, args.ref)
        elif name == "gs-models":
            rc |= fetch_gs_models(into, args.jobs)
        elif name == "core-models":
            rc |= fetch_core_models(into)
    print()
    return cmd_list() if rc == 0 else rc


if __name__ == "__main__":
    sys.exit(main())
