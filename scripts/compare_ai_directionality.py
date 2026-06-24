#!/usr/bin/env python3
"""Compare AI reaction-directionality predictions: new (Claude Opus 4.8) vs old.

Both inputs are AICurationCacheReactionDirectionality.json files keyed by base
reaction id, value ``{"directionality": "forward|reverse|reversible|uncertain", ...}``.

Reports how many shared reactions changed direction, the old->new transition
matrix, and the distribution shift; writes a figure and a changed-reactions CSV.

Usage (defaults compare curated_8848 new vs its .original.json old):
  python scripts/compare_ai_directionality.py
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CUR = ROOT / "data" / "ai_curation" / "curated_8848"
ORDER = ["forward", "reverse", "reversible", "uncertain"]


def load_dir(path: Path) -> dict:
    data = json.loads(Path(path).read_text())
    return {k: (v.get("directionality") if isinstance(v, dict) else None)
            for k, v in data.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--new", default=str(CUR / "AICurationCacheReactionDirectionality.json"),
                    help="new predictions (Claude Opus 4.8)")
    ap.add_argument("--old", default=str(CUR / "AICurationCacheReactionDirectionality.original.json"),
                    help="old/original predictions")
    ap.add_argument("--fig", default=str(ROOT / "reports" / "figures" /
                                         "ai_directionality_opus48_vs_old.png"))
    ap.add_argument("--csv", default=str(ROOT / "data" / "ai_curation" /
                                         "ai_directionality_opus48_vs_old_changes.csv"))
    args = ap.parse_args(argv)

    new = load_dir(Path(args.new))
    old = load_dir(Path(args.old))
    shared = sorted(set(new) & set(old))
    n = len(shared)

    changed = [(r, old[r], new[r]) for r in shared if old[r] != new[r]]
    n_changed = len(changed)
    n_same = n - n_changed

    # transition matrix old(rows) -> new(cols)
    idx = {d: i for i, d in enumerate(ORDER)}
    M = np.zeros((len(ORDER), len(ORDER)), dtype=int)
    for r in shared:
        o, nw = old[r], new[r]
        if o in idx and nw in idx:
            M[idx[o], idx[nw]] += 1
    old_counts = Counter(old[r] for r in shared)
    new_counts = Counter(new[r] for r in shared)

    # ---- text summary -----------------------------------------------------
    print(f"shared reactions compared : {n}")
    print(f"  unchanged direction     : {n_same} ({100*n_same/n:.1f}%)")
    print(f"  CHANGED direction       : {n_changed} ({100*n_changed/n:.1f}%)")
    print(f"  (old-only, not in new)  : {len(set(old) - set(new))}")
    print("\ndistribution (old -> new):")
    for d in ORDER:
        print(f"  {d:11s}: {old_counts.get(d,0):5d} -> {new_counts.get(d,0):5d} "
              f"({new_counts.get(d,0)-old_counts.get(d,0):+d})")
    print("\ntop transitions (old -> new):")
    trans = [((ORDER[i], ORDER[j]), int(M[i, j]))
             for i in range(len(ORDER)) for j in range(len(ORDER)) if i != j and M[i, j]]
    for (o, nw), c in sorted(trans, key=lambda x: -x[1])[:8]:
        print(f"  {o:11s} -> {nw:11s}: {c}")

    # ---- changed-reactions CSV -------------------------------------------
    Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rxn_id", "old_directionality", "new_directionality"])
        for r, o, nw in sorted(changed):
            w.writerow([r, o, nw])
    print(f"\nwrote {args.csv} ({n_changed} rows)")

    # ---- figure -----------------------------------------------------------
    fig, (axM, axB) = plt.subplots(1, 2, figsize=(13, 5.4))

    # Panel A: transition heatmap (off-diagonal = changes)
    disp = M.astype(float).copy()
    for i in range(len(ORDER)):
        disp[i, i] = np.nan  # grey out the unchanged diagonal so changes pop
    im = axM.imshow(disp, cmap="Reds")
    axM.set_xticks(range(len(ORDER))); axM.set_xticklabels(ORDER, rotation=30, ha="right")
    axM.set_yticks(range(len(ORDER))); axM.set_yticklabels(ORDER)
    axM.set_xlabel("new (Claude Opus 4.8)"); axM.set_ylabel("old (original)")
    axM.set_title(f"Direction transitions  (off-diagonal = changed)\n"
                  f"{n_changed:,} of {n:,} changed ({100*n_changed/n:.1f}%)")
    for i in range(len(ORDER)):
        for j in range(len(ORDER)):
            if M[i, j] == 0:
                continue
            on_diag = (i == j)
            axM.text(j, i, f"{M[i,j]:,}", ha="center", va="center",
                     color=("0.5" if on_diag else
                            ("white" if disp[i, j] > np.nanmax(disp) * 0.6 else "black")),
                     fontsize=9, fontweight=("normal" if on_diag else "bold"))
    fig.colorbar(im, ax=axM, fraction=0.046, pad=0.04, label="reactions changed")

    # Panel B: distribution old vs new
    x = np.arange(len(ORDER)); w = 0.38
    axB.bar(x - w/2, [old_counts.get(d, 0) for d in ORDER], w, label="old", color="#9aa7b1")
    axB.bar(x + w/2, [new_counts.get(d, 0) for d in ORDER], w, label="new (Opus 4.8)", color="#c0392b")
    axB.set_xticks(x); axB.set_xticklabels(ORDER, rotation=30, ha="right")
    axB.set_ylabel("reactions"); axB.set_title("Directionality distribution")
    for i, d in enumerate(ORDER):
        axB.text(i - w/2, old_counts.get(d, 0), str(old_counts.get(d, 0)),
                 ha="center", va="bottom", fontsize=8, color="#5b6770")
        axB.text(i + w/2, new_counts.get(d, 0), str(new_counts.get(d, 0)),
                 ha="center", va="bottom", fontsize=8, color="#c0392b")
    axB.legend()

    fig.suptitle("AI reaction-directionality: curated Claude Opus 4.8 vs old data "
                 f"({n:,} shared reactions)", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    Path(args.fig).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.fig, dpi=150)
    print(f"wrote {args.fig}")


if __name__ == "__main__":
    main()
