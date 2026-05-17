#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
plot_two_family_convergence.py — overlay M(P_N) for two parametric
families against N, with each family's Boyd-Lawton limit drawn as a
dashed horizontal line. Shows that distinct (a, d) families converge
to distinct L.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent


def read_23_minus(path: Path):
    """Parse doc/parametric-family.csv (scan_parametric_family.py output).
    Columns: m, deg_Pm, Pm_irreducible, n_factors, M_Pm, ..."""
    xs, ys = [], []
    with path.open() as f:
        for row in csv.DictReader(f):
            try:
                N = int(row["deg_Pm"])
                M = float(row["M_Pm"])
            except (ValueError, KeyError):
                continue
            if 1.001 < M < 1.5:  # drop M=1 cyclotomic-only and M>>1 outliers
                xs.append(N)
                ys.append(M)
    return xs, ys


def read_convergence_csv(path: Path):
    """Parse scan_pn_convergence.py --output (CSV with N, M_PN, ...)."""
    xs, ys = [], []
    with path.open() as f:
        for row in csv.DictReader(f):
            try:
                N = int(row["N"])
                M = float(row["M_PN"])
            except (ValueError, KeyError):
                continue
            if 1.001 < M < 1.5:
                xs.append(N)
                ys.append(M)
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv-a", type=Path,
                    default=REPO / "doc" / "parametric-family.csv",
                    help="family A CSV (default: scan_parametric_family format "
                         "with M_Pm column)")
    ap.add_argument("--csv-a-format",
                    choices=["scan_parametric", "scan_pn"],
                    default="scan_parametric")
    ap.add_argument("--label-a", default="(2,3,-1) Max-U family")
    ap.add_argument("--color-a", default="#cc1010")
    ap.add_argument("--limit-a", type=float, default=1.2554340377272518)
    ap.add_argument("--csv-b", type=Path,
                    default=REPO / "doc" / "parametric-family-convergence-2_12.csv",
                    help="family B CSV (default: scan_pn_convergence format "
                         "with N,M_PN columns)")
    ap.add_argument("--csv-b-format",
                    choices=["scan_parametric", "scan_pn"],
                    default="scan_pn")
    ap.add_argument("--label-b", default="(2,12,-1) family")
    ap.add_argument("--color-b", default="#1060d8")
    ap.add_argument("--limit-b", type=float, default=1.3091)
    ap.add_argument("--output", type=Path,
                    default=REPO / "images" / "two-family-convergence.png")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    reader = {"scan_parametric": read_23_minus,
              "scan_pn":         read_convergence_csv}
    xs_a, ys_a = reader[args.csv_a_format](args.csv_a)
    xs_b, ys_b = reader[args.csv_b_format](args.csv_b)
    print(f"family A: {len(xs_a):,} points")
    print(f"family B: {len(xs_b):,} points")

    fig, ax = plt.subplots(figsize=(13, 7.5), dpi=args.dpi)
    ax.scatter(xs_a, ys_a, s=10, c=args.color_a, alpha=0.6,
               linewidths=0,
               label=f"{args.label_a}  ({len(xs_a):,} N values)",
               zorder=3)
    ax.scatter(xs_b, ys_b, s=10, c=args.color_b, alpha=0.6,
               linewidths=0,
               label=f"{args.label_b}  ({len(xs_b):,} N values)",
               zorder=3)
    ax.axhline(args.limit_a, color=args.color_a, linestyle="--",
               linewidth=1.3, alpha=0.7,
               label=f"$L_{{(2,3)}}$ = {args.limit_a:.4f}",
               zorder=2)
    ax.axhline(args.limit_b, color=args.color_b, linestyle="--",
               linewidth=1.3, alpha=0.7,
               label=f"$L_{{(2,12)}}$ = {args.limit_b:.4f}",
               zorder=2)
    ax.set_xlabel(r"$N$  (polynomial degree)", fontsize=14)
    ax.set_ylabel(r"$M(P_N)$", fontsize=14)
    ax.set_title("Mahler-measure convergence for two parametric families",
                 fontsize=16)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.25, linewidth=0.4, zorder=0)
    ax.legend(loc="upper right", fontsize="medium", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(args.output)
    plt.close()
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
