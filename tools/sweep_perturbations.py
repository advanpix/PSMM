#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
sweep_perturbations.py — generalised sweep over the cyclotomic-perturbation
family

    P_{m,d,s}(x) = (x+1)(x^m + 1) + s * x^((m-1)/2) * Phi_d(x),
                                                      s in {-1, +1}, m odd.

For each (d, s) and each odd m in [m_min, m_max], the script:
  1. Computes the polynomial in PARI.
  2. Factors it over Z.
  3. Records the Mahler measure of the smallest non-cyclotomic factor
     (the "best" Salem-like irreducible factor).
  4. Optionally computes the Boyd-Lawton analytical limit
     L_{d,s} = exp(m(F_{d,s})) where F is the bivariate lift
        F_{d,s}(x, u) = x(x+1)*u^2 + s * Phi_d(x)*u + (x+1).

Writes incremental CSV at the path given by --output.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from math import gcd
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent


# Cyclotomic polynomials Phi_d(x) for small d, written as PARI expressions.
PHI = {
    3:  "(x^2 + x + 1)",
    4:  "(x^2 + 1)",
    5:  "(x^4 + x^3 + x^2 + x + 1)",
    6:  "(x^2 - x + 1)",
    7:  "(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)",
    8:  "(x^4 + 1)",
    9:  "(x^6 + x^3 + 1)",
    10: "(x^4 - x^3 + x^2 - x + 1)",
    12: "(x^4 - x^2 + 1)",
    14: "(x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)",
    15: "(x^8 - x^7 + x^5 - x^4 + x^3 - x + 1)",
    20: "(x^8 - x^6 + x^4 - x^2 + 1)",
}


def phi_degree(d: int) -> int:
    """Degree of Phi_d(x), the d-th cyclotomic polynomial (Euler totient)."""
    return sum(1 for k in range(1, d + 1) if gcd(k, d) == 1)


def gp_script(m: int, d: int, sign: int) -> str:
    s_str = "+" if sign > 0 else "-"
    return f"""
default(realprecision, 30);
P = (x + 1) * (x^{m} + 1) {s_str} x^{(m-1)//2} * {PHI[d]};
N = poldegree(P);
\\\\ Total Mahler measure of P
rts = polroots(P);
Mp = 1.0;
for (i = 1, #rts, if (abs(rts[i]) > 1, Mp *= abs(rts[i])));
\\\\ Factor over Z, find smallest non-cyclotomic factor.
fac = factor(P);
nf = #fac~;
best_deg = 0;
best_M   = 100.0;
for (i = 1, nf, F = fac[i,1]; rs = polroots(F); MF = 1.0; \
    for (j = 1, #rs, if (abs(rs[j]) > 1, MF *= abs(rs[j]))); \
    if (MF > 1.001 && MF < best_M, best_M = MF; best_deg = poldegree(F)));
print({m}, ",", {d}, ",", {sign}, ",", N, ",", nf, ",", Mp, ",", best_deg, ",", best_M);
"""


def gp_boyd_lawton(d: int, sign: int) -> str:
    """gp script: compute Boyd-Lawton limit for the family (d, sign).

    F(x, u) = x(x+1)*u^2 + s*Phi_d(x)*u + (x+1)
    For each x on the unit circle this is quadratic in u with roots u_{1,2}.
    By Jensen, m(F(x,.)) = log|a(x)| + max(log|u_1|,0) + max(log|u_2|,0).
    Then m(F) = (1/2pi) integral over [0, 2pi] of m(F(x,.)) dtheta.
    """
    s_str = "+" if sign > 0 else "-"
    return f"""
default(realprecision, 50);
mF_integrand(theta) = {{
    my(x, a, b, c, disc, u1, u2, lam);
    x = exp(I * theta);
    a = x * (x + 1);
    b = {s_str} {PHI[d]};
    c = x + 1;
    disc = b^2 - 4 * a * c;
    u1 = (-b + sqrt(disc)) / (2 * a);
    u2 = (-b - sqrt(disc)) / (2 * a);
    lam = log(abs(a)) + max(log(abs(u1)), 0.0) + max(log(abs(u2)), 0.0);
    return(lam);
}}
mF = intnum(theta = 0, 2*Pi, mF_integrand(theta)) / (2 * Pi);
print({d}, ",", {sign}, ",", mF, ",", exp(mF));
"""


def run_pari(script: str, parisize: str = "1000000000",
             timeout: int = 600) -> str:
    proc = subprocess.run(
        ["gp", "-q", "--default", f"parisize={parisize}"],
        input=script, capture_output=True, text=True, timeout=timeout,
    )
    return proc.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m-min", type=int, default=5)
    ap.add_argument("--m-max", type=int, default=201)
    ap.add_argument("--ds", type=str, default="3,4,5,6,7,8,9,10,12,14,15,20",
                    help="comma-separated list of d (cyclotomic indices)")
    ap.add_argument("--signs", type=str, default="-1,+1")
    ap.add_argument("--output", default=str(REPO / "doc" / "perturbation-sweep.csv"))
    ap.add_argument("--limits-output", default=str(REPO / "doc" / "boyd-lawton-limits.csv"))
    args = ap.parse_args()

    ds    = [int(x) for x in args.ds.split(",")]
    signs = [int(x) for x in args.signs.split(",")]

    # 1) Compute Boyd-Lawton limits per (d, sign) — small, fast, once.
    limit_path = Path(args.limits_output)
    limit_path.parent.mkdir(parents=True, exist_ok=True)
    print("Computing Boyd-Lawton limits...", file=sys.stderr)
    with limit_path.open("w", buffering=1) as lf:
        lw = csv.writer(lf)
        lw.writerow(["d", "sign", "log_mF", "limit_M"])
        for d in ds:
            for s in signs:
                try:
                    out = run_pari(gp_boyd_lawton(d, s), timeout=300)
                    line = [ln.strip() for ln in out.splitlines()
                            if ln.strip() and not ln.startswith("***")][-1]
                    lw.writerow([c.strip() for c in line.split(",")])
                    print(f"  d={d:3d} s={s:+d}: {line}", file=sys.stderr)
                except Exception as e:
                    print(f"  d={d} s={s}: FAILED ({e})", file=sys.stderr)

    # 2) Sweep (m, d, sign).
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSweeping m={args.m_min}..{args.m_max}, ds={ds}, signs={signs}",
          file=sys.stderr)
    with out_path.open("w", buffering=1) as out:
        w = csv.writer(out)
        w.writerow(["m", "d", "sign", "deg_P", "n_factors",
                    "M_P", "deg_smallest_factor", "M_smallest_factor"])
        for d in ds:
            for s in signs:
                pd = phi_degree(d)
                for m in range(args.m_min, args.m_max + 1, 2):
                    if m % 2 == 0:
                        continue
                    # Require deg(x^((m-1)/2) * Phi_d) <= m+1
                    if (m - 1) // 2 + pd > m + 1:
                        continue
                    try:
                        out_str = run_pari(gp_script(m, d, s), timeout=120)
                    except subprocess.TimeoutExpired:
                        print(f"m={m} d={d} s={s}: timeout", file=sys.stderr)
                        continue
                    lines = [ln.strip() for ln in out_str.splitlines()
                             if ln.strip() and not ln.startswith("***")]
                    if not lines:
                        continue
                    row = [c.strip() for c in lines[-1].split(",")]
                    w.writerow(row)
                    if m % 20 == 1 or m == args.m_max:
                        try:
                            print(f"d={d:3d} s={s:+d} m={m:4d} deg={row[3]:>4} "
                                  f"nf={row[4]:>2} M(P)={row[5][:8]} "
                                  f"smallest: deg={row[6]:>3} M={row[7][:9]}",
                                  file=sys.stderr)
                        except IndexError:
                            pass


if __name__ == "__main__":
    main()
