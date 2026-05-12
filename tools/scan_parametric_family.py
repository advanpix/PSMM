#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
scan_parametric_family.py — sweep the parametric family

  P_m(x) = (x+1)(x^m + 1) - x^((m-1)/2) * Phi_3(x),   m odd

for m in a configurable range. For each m we report:
  - degree (m+1)
  - irreducibility (over Z)
  - factor structure (cyclotomic decomposition + smallest non-cyclotomic
    irreducible factor F_m)
  - Mahler measure M(P_m) and M(F_m)
  - whether F_m matches any polynomial already in AllKnownAdvanpix
    (by 13-digit Mahler-measure key + degree)

Output: CSV at the path given by --output (default: doc/parametric-family.csv).
Run in background; the script writes incrementally so progress is visible.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"


def load_db_keys(path: Path) -> dict:
    """Map (degree, Mahler-prefix-13-digits) -> polynomial line."""
    keys = {}
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 2:
                continue
            try:
                deg = int(toks[0])
            except ValueError:
                continue
            M_prefix = toks[1][:15]  # 13 fractional digits + "1."
            keys[(deg, M_prefix)] = line
    return keys


def gp_script_for_m(m: int) -> str:
    """Build a gp script that prints CSV rows for a single m."""
    return f"""
default(realprecision, 30);
m = {m};
mid = (m - 1)/2;
P = (x+1)*(x^m + 1) - x^mid * (x^2 + x + 1);
N = poldegree(P);
\\\\ Total Mahler measure of P_m
rts = polroots(P);
Mp = 1.0;
for (i = 1, #rts, if (abs(rts[i]) > 1, Mp *= abs(rts[i])));
\\\\ Factor over Z
fac = factor(P);
nf = #fac~;
\\\\ Find smallest non-cyclotomic factor (M > 1.001 to skip cyclo).
best_deg = 0;
best_M   = 100.0;
for (i = 1, nf, F = fac[i,1]; rs = polroots(F); MF = 1.0; \
    for (j = 1, #rs, if (abs(rs[j]) > 1, MF *= abs(rs[j]))); \
    if (MF > 1.001 && MF < best_M, best_M = MF; best_deg = poldegree(F)));
\\\\ Output: CSV row
print(m, ",", N, ",", polisirreducible(P), ",", nf, ",", Mp, ",", best_deg, ",", best_M);
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m-min", type=int, default=5)
    ap.add_argument("--m-max", type=int, default=1001)
    ap.add_argument("--output", default=str(REPO / "doc" / "parametric-family.csv"))
    args = ap.parse_args()

    # Load DB once for cross-checks.
    if DB_PATH.exists():
        print(f"loading {DB_PATH} ...", file=sys.stderr)
        db_keys = load_db_keys(DB_PATH)
        print(f"  {len(db_keys)} entries indexed by (degree, M-prefix)", file=sys.stderr)
    else:
        db_keys = {}

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", buffering=1) as out:
        writer = csv.writer(out)
        writer.writerow([
            "m", "deg_Pm", "Pm_irreducible", "n_factors",
            "M_Pm", "deg_smallest_factor", "M_smallest_factor",
            "factor_in_database",
        ])
        for m in range(args.m_min, args.m_max + 1, 2):  # odd m only
            if m % 2 == 0:
                continue
            try:
                proc = subprocess.run(
                    ["gp", "-q", "--default", "parisize=1000000000"],
                    input=gp_script_for_m(m),
                    capture_output=True, text=True,
                    timeout=600,
                )
            except subprocess.TimeoutExpired:
                print(f"m={m}: timeout", file=sys.stderr)
                writer.writerow([m, "", "", "", "", "", "", "timeout"])
                continue
            if proc.returncode != 0:
                print(f"m={m}: gp failure: {proc.stderr[:200]}", file=sys.stderr)
                continue
            # The gp output is a single CSV line "m, N, irr, nf, Mp, deg_F, M_F"
            lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
            if not lines:
                continue
            row = lines[-1].split(",")
            row = [c.strip() for c in row]
            try:
                deg_F = int(row[5])
                M_F   = row[6]  # high-precision string
                key = (deg_F, M_F[:15])
                in_db = key in db_keys
            except (ValueError, IndexError):
                in_db = False
            writer.writerow(row + [str(in_db)])
            print(f"m={m:5d} deg={row[1]:>4} nf={row[3]:>2} "
                  f"M(P_m)={row[4][:8]} "
                  f"smallest: deg={row[5]:>4} M={row[6][:8]} in_db={in_db}",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
