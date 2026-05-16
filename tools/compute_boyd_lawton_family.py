#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
compute_boyd_lawton_family.py — compute the Boyd-Lawton limit
exp(m(F_{a,d})) for any (a, d) cyclotomic-perturbation family

    P_{a,d,m,s}(x) = Phi_a(x)(x^m + 1) + s * x^{(phi(a)+m-phi(d))/2} * Phi_d(x)

with phi(d) >= phi(a). The bivariate companion is

    F_{a,d}(x, u) = Phi_a(x) * (x^{phi(d)-phi(a)} * u^2 + 1)
                  + s * Phi_d(x) * u
                  = Phi_a(x) * x^{phi(d)-phi(a)} * u^2 + s * Phi_d(x) * u + Phi_a(x).

The Mahler measure is invariant under s -> -s (u -> -u substitution), so
we omit s and compute m(F_{a,d}) once per (a, d) pair.

Reduction (Jensen's lemma): at fixed x = e^{i t}, F is a quadratic in u
with leading A(x) = Phi_a(x) * x^{phi(d)-phi(a)}, middle Phi_d(x), constant
Phi_a(x). The product of roots is |C/A| = 1 on |x| = 1, so one root sits
outside the unit circle (or both on it), and

    m(F_{a,d}) = (1/(2 pi)) integral_0^{2 pi} log+ max(|u_+|, |u_-|) dt
                + (1/(2 pi)) integral log|A(x)| dt
              = (1/(2 pi)) integral log+ max dt + (phi(d)-phi(a)) * m(Phi_a)
              = (1/(2 pi)) integral log+ max dt + 0     (cyclotomic Phi_a has m = 0)

so the constant + cyclotomic part contributes zero, just as for (2, 3).

Usage:
    python3 tools/compute_boyd_lawton_family.py             # survey all 25 pairs
    python3 tools/compute_boyd_lawton_family.py 2 3         # just (a=2, d=3)
    python3 tools/compute_boyd_lawton_family.py --dps 80    # higher precision
"""

from __future__ import annotations

import argparse
import sys

import mpmath as mp


def cyclotomic(n: int, x):
    """Evaluate Phi_n(x) using sympy-style approach but inline for small n.
       For our purposes n in {2, 3, 4, 5, 6, 7, 8, 9, 10, 12}."""
    table = {
        2:  lambda x: x + 1,
        3:  lambda x: x*x + x + 1,
        4:  lambda x: x*x + 1,
        5:  lambda x: x**4 + x**3 + x**2 + x + 1,
        6:  lambda x: x*x - x + 1,
        7:  lambda x: x**6 + x**5 + x**4 + x**3 + x**2 + x + 1,
        8:  lambda x: x**4 + 1,
        9:  lambda x: x**6 + x**3 + 1,
        10: lambda x: x**4 - x**3 + x**2 - x + 1,
        12: lambda x: x**4 - x**2 + 1,
    }
    return table[n](x)


# Euler totient — small n only
PHI = {2:1, 3:2, 4:2, 5:4, 6:2, 7:6, 8:4, 9:6, 10:4, 12:4}


def integrand(a: int, d: int):
    """Return f(t) = log+ max(|u_+|, |u_-|) for the (a, d) bivariate companion."""
    pa = PHI[a]
    pd = PHI[d]
    diff = pd - pa
    def f(t):
        x = mp.expj(t)
        Phi_a = cyclotomic(a, x)
        Phi_d = cyclotomic(d, x)
        if abs(Phi_a) < mp.mpf("1e-40"):
            return mp.mpf(0)   # degenerate at zero of Phi_a
        A = Phi_a * x**diff
        B = Phi_d
        C = Phi_a
        disc = B*B - 4 * A * C
        sq = mp.sqrt(disc)
        u_plus  = (-B + sq) / (2 * A)
        u_minus = (-B - sq) / (2 * A)
        M = max(abs(u_plus), abs(u_minus))
        return max(mp.mpf(0), mp.log(M))
    return f


def cyclotomic_angles(a: int):
    """Angles t in [0, 2*pi) where Phi_a(e^{i t}) = 0 (primitive a-th roots
    of unity). These are log singularities of the integrand: include as
    explicit split points so mpmath.quad handles them correctly."""
    import math
    return [2 * math.pi * k / a for k in range(1, a) if math.gcd(k, a) == 1]


def find_split_points(f, a: int, samples: int = 400):
    """Build integration split points: known cyclotomic singularities of
    Phi_a, plus heuristically-detected |u_+|=1 transitions where the
    integrand crosses 0."""
    import math
    pts = {0.0, 2 * math.pi}
    for ang in cyclotomic_angles(a):
        pts.add(ang)
    # Heuristic sampling for transitions of (|u_+| > 1)
    prev = None
    prev_t = None
    for i in range(samples + 1):
        t = i * 2 * math.pi / samples
        try:
            v = f(t)
        except Exception:
            v = None
        if v is None:
            continue
        if prev is not None:
            crossed = (v > 1e-30) != (prev > 1e-30)
            if crossed:
                pts.add(prev_t)
                pts.add(t)
        prev = v
        prev_t = t
    return sorted(pts)


def compute_for_pair(a: int, d: int, dps: int = 60):
    mp.mp.dps = dps
    f = integrand(a, d)
    splits = find_split_points(f, a)
    splits_mp = [mp.mpf(p) for p in splits]
    val = mp.quad(f, splits_mp) / (2 * mp.pi)
    return val, mp.exp(val)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a", nargs="?", type=int, help="(omit to survey all 25 pairs)")
    ap.add_argument("d", nargs="?", type=int)
    ap.add_argument("--dps", type=int, default=60, help="precision (digits)")
    args = ap.parse_args()

    if args.a is not None and args.d is not None:
        m_F, M_F = compute_for_pair(args.a, args.d, args.dps)
        print(f"(a, d) = ({args.a}, {args.d}):")
        print(f"  m(F) = {mp.nstr(m_F, 20)}")
        print(f"  M(F) = exp(m(F)) = {mp.nstr(M_F, 20)}")
        if float(M_F) < 1.3:
            print("  ** L < 1.3 -> family CAN produce sub-1.3 polynomials **")
        else:
            print("  L > 1.3 -> family CANNOT produce sub-1.3 polynomials")
        return

    # Survey all 25 (a, d) pairs from sweep_ad.py's default
    pairs = [
        (2,3), (2,5), (2,7), (2,8), (2,9), (2,10), (2,12),
        (3,5), (3,7), (3,8), (3,9), (3,10), (3,12),
        (4,5), (4,7), (4,8), (4,9), (4,10), (4,12),
        (6,5), (6,7), (6,8), (6,9), (6,10), (6,12),
    ]
    print(f"Surveying {len(pairs)} (a, d) pairs at dps={args.dps}:")
    print(f"{'a':>3} {'d':>3} {'M(F) = exp(m(F))':>22}  flag")
    sub13 = []
    for (a, d) in pairs:
        try:
            _, M_F = compute_for_pair(a, d, args.dps)
            M_f = float(M_F)
            flag = "*** L<1.3 ***" if M_f < 1.3 else ""
            print(f"{a:>3} {d:>3}  {mp.nstr(M_F, 16):>22}  {flag}",
                  flush=True)
            if M_f < 1.3:
                sub13.append((a, d, M_F))
        except Exception as e:
            print(f"{a:>3} {d:>3}  ERROR: {str(e)[:60]}", flush=True)

    print()
    print(f"L < 1.3 families (can produce sub-1.3 entries at large m):")
    for (a, d, M_F) in sub13:
        print(f"  (a={a}, d={d}):  L = {mp.nstr(M_F, 18)}")
    if not sub13:
        print("  (none)")


if __name__ == "__main__":
    main()
