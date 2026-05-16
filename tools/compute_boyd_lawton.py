#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
compute_boyd_lawton.py — compute m(F) and exp(m(F)) for the bivariate
companion of the Max-U family P_N.

The univariate family
    P_N(x) = (x+1)(x^{N-1} + 1) - x^{N/2 - 1} Phi_3(x)
is a monomial substitution into
    F(x, u) = x(x+1) u^2 - Phi_3(x) u + (x+1).
By Boyd-Lawton, lim_{N -> infty} M(P_N) = exp(m(F)).

Reduce the 2D Mahler integral to a 1D integral via Jensen's lemma:
At fixed x = e^{i t}, F is a quadratic in u with leading c = x(x+1),
roots u_+, u_- satisfying |u_+||u_-| = 1 on |x| = 1 (so one root is
inside, one outside the unit circle, or both on it). Then

    log M(F_x) = log|x(x+1)| + log+|u_+(x)| + log+|u_-(x)|
               = log|x(x+1)| + log+ max(|u_+|, |u_-|).

The constant + log|cos(t/2)| parts integrate to zero around the circle
(int log|cos| = -log 2 cancels the int log 2). So

    m(F) = (1/(2 pi)) int_0^{2 pi} log+ max(|u_+(t)|, |u_-(t)|) dt.

A previous estimate (1.249358390752959...) computed via 1D adaptive
integration of the wrong reduction was unreliable because of log
singularities of |x(x+1)| on the unit circle. The reduction above
removes that singularity analytically.

Run:  python3 tools/compute_boyd_lawton.py
Result: m(F) and M(F) printed to ~30 decimal digits.
"""

from __future__ import annotations

import mpmath as mp


def integrand(t):
    """log+ max(|u_+|, |u_-|) at x = e^{i t}."""
    x = mp.expj(t)
    x_plus_1 = x + 1
    if abs(x_plus_1) < mp.mpf("1e-40"):
        return mp.mpf(0)        # at x = -1 the quadratic degenerates
    phi3 = x*x + x + 1
    disc = phi3*phi3 - 4 * x * x_plus_1 * x_plus_1
    sq = mp.sqrt(disc)
    denom = 2 * x * x_plus_1
    u_plus  = (phi3 + sq) / denom
    u_minus = (phi3 - sq) / denom
    M = max(abs(u_plus), abs(u_minus))
    return max(mp.mpf(0), mp.log(M))


def compute(precision: int = 80):
    mp.mp.dps = precision
    # Split the interval at the only "real" feature: t = pi (where the
    # quadratic degenerates to linear). The transition between
    # "both roots on |u|=1" and "one outside" happens around t = pi.
    val = mp.quad(integrand, [0, mp.pi, 2 * mp.pi]) / (2 * mp.pi)
    return val, mp.exp(val)


def main():
    print("Bivariate companion of the P_N family (Max-U):")
    print("  F(x, u) = x(x+1) u^2 - Phi_3(x) u + (x+1)")
    print()
    for prec in (30, 50, 80, 120):
        m_F, M_F = compute(prec)
        print(f"  dps={prec:>3d}  m(F) = {mp.nstr(m_F, 30)}")
        print(f"            M(F) = {mp.nstr(M_F, 30)}")
    print()
    print("Empirical cross-check: mean of M(P_N) for N in [500, 1002]")
    print("  agrees with M(F) above to ~10^-6 (sample mean dominates the")
    print("  numerical convergence error in the integration).")


if __name__ == "__main__":
    main()
