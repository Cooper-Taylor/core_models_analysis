#!/usr/bin/env bash
# Rebuild the pre-2022-05-20 ModelSEEDDatabase layout that ModelSEEDpy's
# FullThermoPkg still expects.
#
# MSDB commit 8bf2281b (2022-05-20, "Massive re-orientation towards breaking
# down biochemistry into chunks") deleted Biochemistry/compounds.tsv and
# reactions.tsv in favour of 50 shards, and dropped the `direction` column.
# modelseedpy's from_local() dispatches on whether compounds.tsv exists: with
# it you get ModelSEEDBiochem (which has get_seed_compound and .deltag, the
# API FullThermoPkg was written against); without it you get ModelSEEDDatabase
# and the package raises AttributeError.
#
# This reassembles that layout from the CURRENT sharded data, so the stock
# package runs against today's numbers. ModelSEEDDatabase itself is only read.
#
# See reports/THERMODYNAMIC_FBA.md.
set -euo pipefail

SRC="${MSDB_ROOT:-/scratch/ctaylor/ModelSEEDDatabase}"
DST="${1:-${TFBA_LEGACY_MSDB:-$(git rev-parse --show-toplevel)/results/tfba/legacy_msdb}}"

[ -d "$SRC/Biochemistry" ] || { echo "no Biochemistry/ under $SRC" >&2; exit 1; }

rm -rf "$DST"
mkdir -p "$DST/Biochemistry"
ln -s "$SRC/Biochemistry/Aliases"    "$DST/Biochemistry/Aliases"
ln -s "$SRC/Biochemistry/Structures" "$DST/Biochemistry/Structures"

for kind in compound reaction; do
    out="$DST/Biochemistry/${kind}s.tsv"
    first=1
    for f in $(ls "$SRC/Biochemistry/${kind}_"[0-9][0-9].tsv | sort); do
        if [ $first -eq 1 ]; then cat "$f" > "$out"; first=0
        else tail -n +2 "$f" >> "$out"; fi
    done
    echo "$(basename "$out"): $(( $(wc -l < "$out") - 1 )) rows"
done

# The legacy loader reads a `direction` column (old col 10) that the sharded
# schema no longer carries; re-derive it from `reversibility`.
python3 - "$DST/Biochemistry/reactions.tsv" <<'PY'
import csv, sys
src = sys.argv[1]
rows = list(csv.DictReader(open(src), delimiter="\t"))
cols = list(rows[0].keys())
if "direction" not in cols:
    cols.insert(cols.index("reversibility") + 1, "direction")
    with open(src, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            r["direction"] = r.get("reversibility") or "="
            w.writerow(r)
    print(f"added `direction` column ({len(rows)} rows)")
PY

echo "legacy layout ready at $DST"
