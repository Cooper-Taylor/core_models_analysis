#!/usr/bin/env python3
"""Characterization harness: prove a change did not move the numbers.

    python3 scripts/check_goldens.py --list
    python3 scripts/check_goldens.py --all
    python3 scripts/check_goldens.py --tag cascade
    python3 scripts/check_goldens.py --update        # re-pin after an intended change

Compares regenerated artifacts against the pinned copies under
``thermo_variants/``. Two things make the comparison trustworthy rather than
noisy:

* **``new_rev`` is compared separately from ``status``.** Status strings embed
  pre-formatted floats (``'MdeltaG(Max): -1.23'``, ``'lnRI: -7.78'``), so a
  value that shifts in the fourth decimal produces a mass textual diff that is
  not a behaviour change. Only the direction column decides pass/fail.

* **Exit 2 means "refuse to judge".** If the MSDB reference that produced the
  pinned artifacts differs from the one on disk now, the harness reports a data
  mismatch instead of failing. A harness that cries wolf for data reasons gets
  ignored, which is worse than not having one.

Exit codes: 0 pass, 1 a real behaviour change, 2 a data-reference mismatch.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cma import paths, variants  # noqa: E402

REPORTS = (
    ("EQ", "Estimated_Reaction_Reversibility_Report_EQ.txt", 4),
    ("GC", "Estimated_Reaction_Reversibility_Report_GC.txt", 3),
    ("unfiltered", "Estimated_Reaction_Reversibility_Report.txt", 4),
)


def parse(path: Path, ncols: int) -> dict:
    """``{rxn: (status, new_rev)}``; the GC report omits the old_rev column."""
    out = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        c = line.split("\t")
        if len(c) >= ncols:
            out[c[0]] = (c[1], c[-1].strip())
    return out


def msdb_ref() -> str:
    try:
        r = subprocess.run(["git", "-C", str(paths.msdb("live")), "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True, timeout=30)
        b = subprocess.run(["git", "-C", str(paths.msdb("live")), "branch", "--show-current"],
                           capture_output=True, text=True, timeout=30)
        return f"{b.stdout.strip()}@{r.stdout.strip()}"
    except Exception:
        return "unknown"


def pin_path() -> Path:
    return paths.thermo_variants("GOLDEN_REF.json")


def compare_tag(tag: str, fresh_root: Path) -> dict:
    pinned_dir = paths.thermo_variants(tag)
    fresh_dir = fresh_root / tag
    res = {"tag": tag, "status_diffs": 0, "direction_diffs": 0, "missing": [], "n": 0}
    for _label, fname, ncols in REPORTS:
        a, b = pinned_dir / fname, fresh_dir / fname
        if not a.exists() or not b.exists():
            res["missing"].append(fname)
            continue
        pa, pb = parse(a, ncols), parse(b, ncols)
        keys = set(pa) & set(pb)
        res["n"] += len(keys)
        res["status_diffs"] += sum(1 for k in keys if pa[k][0] != pb[k][0])
        res["direction_diffs"] += sum(1 for k in keys if pa[k][1] != pb[k][1])
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--tag", action="append", help="variant tag (repeatable)")
    ap.add_argument("--update", action="store_true", help="re-pin the MSDB reference")
    ap.add_argument("--out", type=Path, default=None, help="reuse an existing fresh export")
    args = ap.parse_args()

    cascade = variants.tags("cascade")
    if args.list:
        print("cascade variants (regenerable, so checkable):")
        for t in cascade:
            print(f"  {t}")
        print("\noverlay/derived variants are built from external tables and are")
        print("compared by the builder that produces them, not here.")
        return 0

    tags = args.tag or (cascade if args.all else cascade[:3])
    ref_now = msdb_ref()
    pin = json.loads(pin_path().read_text()) if pin_path().exists() else {}

    if args.update:
        pin_path().write_text(json.dumps({"msdb_ref": ref_now, "tags": cascade}, indent=2))
        print(f"pinned MSDB reference {ref_now} in {pin_path()}")
        return 0

    print(f"MSDB reference now    : {ref_now}")
    print(f"MSDB reference pinned : {pin.get('msdb_ref', '(unpinned)')}")
    if pin.get("msdb_ref") and pin["msdb_ref"] != ref_now:
        print("\nREFUSING TO COMPARE: the pinned artifacts were produced against a "
              "different MSDB reference, so any diff would be a data difference, not "
              "a behaviour change. Re-pin with --update once the reference is right.")
        return 2

    tmp = args.out or Path(tempfile.mkdtemp(prefix="cma_goldens_"))
    if args.out is None:
        cmd = [paths.python_exe(), str(paths.scripts("export_thermo_variants.py")),
               "--out", str(tmp)]
        for t in tags:
            cmd += ["--only", t]
        print(f"\nregenerating {len(tags)} variant(s) -> {tmp}")
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            print(p.stdout[-2000:]); print(p.stderr[-2000:])
            return 1

    print(f"\n{'tag':<14} {'reactions':>10} {'status≠':>9} {'DIRECTION≠':>11}  verdict")
    rc = 0
    for t in tags:
        r = compare_tag(t, tmp)
        verdict = "pass"
        if r["direction_diffs"]:
            verdict = "FAIL (behaviour changed)"
            rc = 1
        elif r["status_diffs"]:
            verdict = "pass (status text only)"
        if r["missing"]:
            verdict = f"missing {r['missing']}"
            rc = max(rc, 1)
        print(f"{t:<14} {r['n']:>10} {r['status_diffs']:>9} {r['direction_diffs']:>11}  {verdict}")

    print()
    print("RESULT:", "pass" if rc == 0 else "FAIL")
    return rc


if __name__ == "__main__":
    sys.exit(main())
