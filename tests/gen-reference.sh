#!/bin/bash
# Build tests/data/reference-polynomials.txt by picking real rows from
# AllKnownAdvanpix. This avoids transcription errors.

set -e
SRC=/home/advanpix/Development/PSMM/AllKnownAdvanpix
OUT=/home/advanpix/Development/PSMM/tests/data/reference-polynomials.txt

header='# Curated reference set for PSMM unit/regression tests.
#
# Extracted verbatim from AllKnownAdvanpix (Pavel Holoborodko, Advanpix).
# Format matches AllKnownAdvanpix:
#
#   N M NNZ H L K U Q R c_0 c_1 ... c_{N/2}
#
# Selection:
#   - lowest-M polynomial at each even degree N in {10, 12, ..., 30}
#     (baseline correctness)
#   - several near-threshold polynomials (M in [1.29, 1.30]) at mixed
#     degrees (critical for A1 regression: none may be falsely rejected
#     when threshold T = 1.3)
#
# Regeneration: tests/gen-reference.sh (if ever needed).
'
echo "$header" > "$OUT"

# Baseline: lowest-M per degree for small N.
echo "# --- Lowest-M baseline per degree ---" >> "$OUT"
for N in 10 12 14 16 18 20 22 24 26 28 30; do
    # First row after "# degree = N," header.
    awk -v deg="$N" '
        $0 ~ "^# degree = "deg"," { in_block = 1; next }
        /^$/ { in_block = 0 }
        in_block && $1 ~ /^[0-9]+$/ {
            print $0
            in_block = 0
        }
    ' "$SRC" >> "$OUT"
done

echo "" >> "$OUT"
echo "# --- Near-threshold (M in [1.29, 1.30]); A1 regression set ---" >> "$OUT"
# Pick at most one per degree so we spread across N.
awk '
    $1 ~ /^[0-9]+$/ && $2+0 >= 1.29 && $2+0 <= 1.30 {
        if (!(($1) in seen)) {
            print $0
            seen[$1] = 1
        }
    }
' "$SRC" | head -10 >> "$OUT"
