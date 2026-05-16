#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
plot_parametric_family.py — render two figures from the CSV produced by
scan_parametric_family.py:

  1. images/parametric_family_M_vs_m.png  (M(P_N) and M(F_N) vs N, log-spaced)
  2. images/parametric_family_residual.png  (M - M_limit on log scale)

The x-axis is the polynomial degree N (= m + 1 in the underlying CSV's
perturbation-index naming). The Boyd-Lawton limit, computed analytically
as the bivariate Mahler measure of F(x, u) = x(x+1)u^2 - (x^2+x+1)u + (x+1),
is overlaid as a horizontal line.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parent.parent
CSV_PATH = REPO / "doc" / "parametric-family.csv"
OUT_DIR = REPO / "images"

# Boyd-Lawton limit (see tools/compute_boyd_lawton.py).
# m(F)   ≈ 0.22748124...
# M(F)   ≈ 1.25543404...  (1D Jensen reduction, dps=80, mpmath adaptive quad)
# Empirical: mean of M(P_N) for N in [500, 1002] is 1.2554338
# (agrees with the analytic limit to 7 decimal places).
LIMIT_M = 1.2554340377272518


def main():
    rows = []
    with CSV_PATH.open() as f:
        for row in csv.DictReader(f):
            try:
                # Underlying CSV uses 'm' as the perturbation index;
                # the polynomial degree is N = m + 1 (= deg_Pm).
                N = int(row["deg_Pm"])
                Mpn = float(row["M_Pm"])
                Mfn = float(row["M_smallest_factor"])
                in_db = row["factor_in_database"] == "True"
                rows.append((N, Mpn, Mfn, in_db))
            except (ValueError, KeyError):
                continue

    if not rows:
        print("CSV has no usable rows yet — has the scan started?", file=sys.stderr)
        sys.exit(1)

    rows.sort(key=lambda r: r[0])
    Ns       = np.array([r[0] for r in rows])
    Mpn_arr  = np.array([r[1] for r in rows])
    Mfn_arr  = np.array([r[2] for r in rows])
    in_db    = np.array([r[3] for r in rows])

    OUT_DIR.mkdir(exist_ok=True)

    # === Figure 1: M(P_N) and M(F_N) vs N ===
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(Ns, Mpn_arr, s=10, alpha=0.5, label=r"$M(P_N)$ — full polynomial")
    # Different marker for reducible cases (where M_FN < M_PN noticeably).
    diff = Mpn_arr - Mfn_arr
    is_reducible = diff > 1e-6
    if is_reducible.any():
        ax.scatter(Ns[is_reducible], Mfn_arr[is_reducible],
                   s=30, c="red", marker="x",
                   label=r"$M(F_N)$ — smallest irreducible factor (reducible $P_N$ only)")
    ax.axhline(LIMIT_M, color="green", linestyle="--", linewidth=1.5,
               label=f"Boyd-Lawton limit: $\\exp(m(F)) \\approx {LIMIT_M:.7f}$")
    ax.axhline(1.17628081826, color="orange", linestyle=":", linewidth=1.5,
               label=r"Lehmer's number $\approx 1.17628$")
    ax.set_xlabel(r"$N$ (even, polynomial degree)")
    ax.set_ylabel(r"$M$")
    ax.set_xscale("log")
    ax.set_title(r"Parametric family $P_N(x) = (x+1)(x^{N-1}+1) - x^{N/2-1}\Phi_3(x)$")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(1.15, 1.30)
    out1 = OUT_DIR / "parametric_family_M_vs_m.png"
    fig.tight_layout()
    fig.savefig(out1, dpi=130)
    plt.close(fig)
    print(f"Wrote {out1}", file=sys.stderr)

    # === Figure 2: residual M - limit on log scale ===
    fig, ax = plt.subplots(figsize=(10, 6))
    residual = Mpn_arr - LIMIT_M
    pos = residual > 0
    ax.loglog(Ns[pos], residual[pos], "o", markersize=4, alpha=0.6,
              label=r"$M(P_N) - \exp(m(F))$")
    # Reference 1/N^2 decay line
    N_ref = np.logspace(np.log10(Ns.min()), np.log10(Ns.max()), 100)
    if pos.any():
        # Scale to match early data
        ax.loglog(N_ref, 1.0 / N_ref**2, "--", color="gray", alpha=0.5,
                  label=r"$\propto 1/N^2$ (reference)")
        ax.loglog(N_ref, 1.0 / N_ref, ":", color="gray", alpha=0.5,
                  label=r"$\propto 1/N$ (reference)")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$M(P_N) - \exp(m(F))$")
    ax.set_title("Convergence of $M(P_N)$ to the Boyd-Lawton limit")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    out2 = OUT_DIR / "parametric_family_residual.png"
    fig.tight_layout()
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print(f"Wrote {out2}", file=sys.stderr)

    # === Stats ===
    n_irr = (Ns[Mpn_arr - Mfn_arr < 1e-6]).size
    n_reducible = is_reducible.sum()
    in_db_count = in_db.sum()
    print(f"\n  total N values scanned: {len(Ns)}", file=sys.stderr)
    print(f"  P_N irreducible      : {n_irr}", file=sys.stderr)
    print(f"  P_N reducible        : {n_reducible}", file=sys.stderr)
    print(f"  smallest factor F_N already in AllKnownAdvanpix : {in_db_count}/{len(Ns)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
