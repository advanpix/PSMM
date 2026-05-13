#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
sweep_ad.py — empirical (a, d, m, sign) sweep for the cyclotomic-perturbation
family with VARYING BACKGROUND:

    P_{a,d,m,s}(x) = Phi_a(x) (x^m + 1) + s * x^k * Phi_d(x)
    k = (phi(a) + m - phi(d)) / 2          (requires phi(a) + m + phi(d) even)

For each (a, d, sign, m) with phi(d) >= phi(a) (so k >= (m-1)/2 sits in the
polynomial interior, and the bivariate lift is a polynomial not a Laurent):
  - construct P, factor over Z
  - report Mahler measure of smallest non-cyclotomic factor (the "Salem")

Writes incremental CSV to doc/ad_sweep.csv.

Goal: find the global minimum across all (a, d, sign, m) in our range.
The current best is (a=2, d=3, m=21) with M=1.18837 (the second-smallest
known Salem polynomial).
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent


CYCLOTOMICS = {
    2:  ("(x+1)", 1),
    3:  ("(x^2+x+1)", 2),
    4:  ("(x^2+1)", 2),
    5:  ("(x^4+x^3+x^2+x+1)", 4),
    6:  ("(x^2-x+1)", 2),
    7:  ("(x^6+x^5+x^4+x^3+x^2+x+1)", 6),
    8:  ("(x^4+1)", 4),
    9:  ("(x^6+x^3+1)", 6),
    10: ("(x^4-x^3+x^2-x+1)", 4),
    12: ("(x^4-x^2+1)", 4),
    15: ("(x^8-x^7+x^5-x^4+x^3-x+1)", 8),
}


def gp_script(a: int, d: int, m: int, sign: int) -> str:
    phi_a, pa = CYCLOTOMICS[a]
    phi_d, pd = CYCLOTOMICS[d]
    k = (pa + m - pd) // 2
    s_str = "+" if sign > 0 else "-"
    return f"""
default(realprecision, 25);
P = {phi_a} * (x^{m} + 1) {s_str} x^{k} * {phi_d};
N = poldegree(P);
rts = polroots(P);
Mp = 1.0;
for (i = 1, #rts, if (abs(rts[i]) > 1, Mp *= abs(rts[i])));
fac = factor(P);
nf = #fac~;
best_deg = 0;
best_M   = 100.0;
for (i = 1, nf, F = fac[i,1]; rs = polroots(F); MF = 1.0; \
    for (j = 1, #rs, if (abs(rs[j]) > 1, MF *= abs(rs[j]))); \
    if (MF > 1.001 && MF < best_M, best_M = MF; best_deg = poldegree(F)));
print({a}, ",", {d}, ",", {sign}, ",", {m}, ",", N, ",", nf, ",", Mp, ",", best_deg, ",", best_M);
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m-min", type=int, default=5)
    ap.add_argument("--m-max", type=int, default=201)
    ap.add_argument("--ads",
                    default="2,3;2,5;2,7;2,8;2,9;2,10;2,12;3,5;3,7;3,8;3,9;3,10;3,12;4,5;4,7;4,8;4,9;4,10;4,12;6,5;6,7;6,8;6,9;6,10;6,12",
                    help="semicolon-separated 'a,d' pairs")
    ap.add_argument("--signs", default="-1,+1")
    ap.add_argument("--output", default=str(REPO / "doc" / "ad_sweep.csv"))
    args = ap.parse_args()

    ads = [tuple(int(z) for z in pair.split(","))
           for pair in args.ads.split(";")]
    signs = [int(s) for s in args.signs.split(",")]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", buffering=1) as out:
        w = csv.writer(out)
        w.writerow(["a", "d", "sign", "m", "deg_P", "n_factors",
                    "M_P", "deg_smallest_factor", "M_smallest_factor"])
        best_overall = (None, 100.0)
        for (a, d) in ads:
            pa = CYCLOTOMICS[a][1]
            pd = CYCLOTOMICS[d][1]
            if pd < pa:
                continue
            # parity: pa + m + pd must be even, so m parity = (pa + pd) mod 2
            m_parity = (pa + pd) % 2
            for s in signs:
                # Smallest-M-in-family tracker
                best_in_family = (None, 100.0)
                for m in range(args.m_min, args.m_max + 1):
                    if m % 2 != m_parity:
                        continue
                    if m < 1:
                        continue
                    try:
                        out_str = subprocess.run(
                            ["gp", "-q", "--default", "parisize=1000000000"],
                            input=gp_script(a, d, m, s),
                            capture_output=True, text=True, timeout=120,
                        ).stdout
                    except subprocess.TimeoutExpired:
                        continue
                    lines = [ln.strip() for ln in out_str.splitlines()
                             if ln.strip() and not ln.startswith("***")]
                    if not lines:
                        continue
                    row = [c.strip() for c in lines[-1].split(",")]
                    w.writerow(row)
                    try:
                        m_sf = float(row[8])
                    except (ValueError, IndexError):
                        continue
                    if m_sf < best_in_family[1]:
                        best_in_family = ((a, d, s, m, row), m_sf)
                    if m_sf < best_overall[1]:
                        best_overall = ((a, d, s, m), m_sf)
                key, val = best_in_family
                if key is not None:
                    a_, d_, s_, m_, _ = key
                    print(f"family (a={a_}, d={d_}, sign={s_}): "
                          f"best m={m_}, M={val:.6f}",
                          file=sys.stderr)
        print(f"\nGLOBAL BEST: {best_overall}", file=sys.stderr)


if __name__ == "__main__":
    main()
