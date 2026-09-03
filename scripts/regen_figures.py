#!/usr/bin/env python3
"""Regenerate figure sets from scripts/figures.tsv, in parallel.

The registry exists because "which script draws this figure?" was costing more
time than the edit itself. Every figure set has one row: output directory, the
script, the exact argv, and tags.

    python3 scripts/regen_figures.py --list
    python3 scripts/regen_figures.py --all
    python3 scripts/regen_figures.py --tag transition     # the coloured scatters
    python3 scripts/regen_figures.py kegg_fix_impact thermo_source_dg_scatter
    python3 scripts/regen_figures.py --all --exclude-tag slow

Each job runs as its own process; logs go to reports/thermoComparison/figures/
_regen_logs/<name>.log and only the tail is echoed on failure. Exit status is
non-zero if any job failed.
"""
from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "scripts" / "figures.tsv"
LOG_DIR = ROOT / "reports" / "thermoComparison" / "figures" / "_regen_logs"
PYTHON = os.environ.get(
    "CMA_PYTHON",
    "/mnt/homes/ctaylor/conda/miniforge3/envs/core_models_analysis/bin/python")


def load_registry() -> list[dict]:
    rows = []
    with open(REGISTRY) as fh:
        lines = [ln for ln in fh if not ln.startswith(">")]
    for row in csv.DictReader(lines, delimiter="\t"):
        if not row.get("name"):
            continue
        row["tags"] = [t.strip() for t in (row.get("tags") or "").split(",") if t.strip()]
        rows.append(row)
    return rows


def run_one(job: dict, dry: bool) -> tuple[str, int, float, Path]:
    cmd = [PYTHON, str(ROOT / "scripts" / job["script"])]
    cmd += shlex.split(job.get("args") or "")
    log = LOG_DIR / f"{job['name']}.log"
    if dry:
        print(f"  would run: {shlex.join(cmd)}")
        return job["name"], 0, 0.0, log
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write(shlex.join(cmd) + "\n\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    return job["name"], rc, time.time() - t0, log


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="*", help="figure set names (default: none)")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--tag", action="append", default=[],
                    help="include every set carrying this tag (repeatable)")
    ap.add_argument("--exclude-tag", action="append", default=[])
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-j", "--jobs", type=int, default=4)
    args = ap.parse_args()

    reg = load_registry()
    by_name = {r["name"]: r for r in reg}

    if args.list:
        w = max(len(r["name"]) for r in reg)
        print(f"{'name'.ljust(w)}  {'script':46s}  tags")
        for r in reg:
            print(f"{r['name'].ljust(w)}  {r['script']:46s}  {','.join(r['tags'])}")
        return 0

    selected: list[dict] = []
    if args.all:
        selected = list(reg)
    else:
        for n in args.names:
            if n not in by_name:
                print(f"unknown figure set: {n}  (try --list)", file=sys.stderr)
                return 2
            selected.append(by_name[n])
        for t in args.tag:
            selected += [r for r in reg if t in r["tags"]]
    for t in args.exclude_tag:
        selected = [r for r in selected if t not in r["tags"]]
    # de-duplicate, preserve order
    seen, jobs = set(), []
    for r in selected:
        if r["name"] not in seen:
            seen.add(r["name"])
            jobs.append(r)
    if not jobs:
        ap.print_help()
        return 2

    print(f"regenerating {len(jobs)} figure set(s) with {args.jobs} worker(s):")
    for j in jobs:
        print(f"  {j['name']}")
    print()

    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for name, rc, dt, log in ex.map(lambda j: run_one(j, args.dry_run), jobs):
            if args.dry_run:
                continue
            status = "ok  " if rc == 0 else "FAIL"
            print(f"  [{status}] {name:38s} {dt:6.1f}s   {log.relative_to(ROOT)}")
            if rc != 0:
                failures.append((name, log))
    for name, log in failures:
        print(f"\n----- tail of {name} -----")
        print("\n".join(log.read_text().splitlines()[-15:]))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
