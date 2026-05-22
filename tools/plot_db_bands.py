#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
plot_db_bands.py — M(P) vs N scatter, configurable for README inclusion.

Default: NNZ-stratified scatter of AllKnownAdvanpix, truncated to
N <= 456 (the brute-force boundary; bands beyond that become
needle-sparse), M <= 1.30 (the two-loci story window). L values
annotated on the right side of each band.

Common variants for tonight's iteration:

  --variant two-bands
      Default. N in [6, 456], M in [1.17, 1.30], L lines at 1.2554
      and 1.2857.

  --variant three-bands
      Same but extended to M in [1.17, 1.32] to include the (2,12,-)
      family band at L = 1.30911 (useful for the Future Work section).

  --variant full
      No truncation: N up to DB max, M up to 1.5 or whatever's in DB.
      Cluttered but complete.

  --variant log-x   /  --variant linear-x
      Swap the x-axis scale.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parent.parent

# Boyd-Lawton limits we know about, with their family labels.
KNOWN_L = [
    (1.2554340377, r"$L_{(2,3)} = 1.2554$"),
    (1.2857144758, r"$L_{S5[\Phi_3]} = 1.2857$"),
    (1.3091112,    r"$L_{(2,12)} = 1.3091$"),
]


def load_db():
    Ns, Ms, NNZs = [], [], []
    with (REPO / "AllKnownAdvanpix").open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            t = s.split()
            try:
                N = int(t[0]); M = float(t[1]); NNZ = int(t[2])
            except (ValueError, IndexError):
                continue
            Ns.append(N); Ms.append(M); NNZs.append(NNZ)
    return np.array(Ns), np.array(Ms), np.array(NNZs)


def render(variant: str, n_max: int, m_min: float, m_max: float,
           xscale: str, out_path: Path,
           dpi: int = 150, figsize=(12, 7), show_lehmer: bool = True):
    Ns, Ms, NNZs = load_db()
    mask = (Ns <= n_max) & (Ms >= m_min) & (Ms <= m_max)
    Ns_p = Ns[mask]; Ms_p = Ms[mask]; NNZs_p = NNZs[mask]

    lehmer_mask = (Ns_p == 10) & (Ms_p < 1.18)
    sparse_mask = (NNZs_p <= 3) & (~lehmer_mask)
    dense_mask  = (NNZs_p > 3) & (~lehmer_mask)
    n_sparse = sparse_mask.sum()
    n_dense  = dense_mask.sum()

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.grid(True, alpha=0.20, linewidth=0.4, zorder=0)

    # L lines + right-side labels showing the limiting M value each
    # band converges to. Position inside the plot's right edge, above
    # the band line; numeric value only (family naming deferred to
    # prose since multiple polynomial families share each L).
    label_x = n_max * 0.96
    for L, _ in KNOWN_L:
        if m_min <= L <= m_max:
            ax.axhline(L, color="black", linestyle="--",
                       linewidth=0.7, alpha=0.30, zorder=1)
            ax.text(label_x, L + 0.0015, f"{L:.4f}",
                    fontsize=11, va="bottom", ha="left", color="#333")

    # Sparse layer (blue), dense layer (red on top).
    ax.scatter(Ns_p[sparse_mask], Ms_p[sparse_mask],
               s=3.5, c="#1060d8", alpha=0.50, linewidths=0, zorder=2,
               label=f"sparse  (NNZ ≤ 3,  {n_sparse:,})")
    ax.scatter(Ns_p[dense_mask], Ms_p[dense_mask],
               s=3.5, c="#cc1010", alpha=0.55, linewidths=0, zorder=3,
               label=f"dense   (NNZ > 3,  {n_dense:,})")

    # Lehmer marker (only when his M is in range). Label sits to the
    # right of the star, with a small absolute N offset that's
    # readable on both linear and log axes.
    if show_lehmer and lehmer_mask.any():
        ax.scatter(Ns_p[lehmer_mask], Ms_p[lehmer_mask], s=120,
                   c="#000", marker="*", linewidths=0, zorder=5)
        x_off = (7 if xscale == "linear"
                 else Ns_p[lehmer_mask][0] * 0.35)
        ax.text(Ns_p[lehmer_mask][0] + x_off if xscale == "linear"
                else Ns_p[lehmer_mask][0] * 1.35,
                Ms_p[lehmer_mask][0],
                "Lehmer", fontsize=11, va="center", ha="left",
                color="#000")

    ax.set_xlabel(rf"$N$ (polynomial degree, ${xscale}$ scale)"
                  if xscale == "log" else r"$N$ (polynomial degree)",
                  fontsize=13)
    ax.set_ylabel(r"$M(P)$", fontsize=13)
    if xscale == "log":
        ax.set_xscale("log")

    # Small right extension so the L-value labels (positioned near the
    # right edge, just above each band) fit entirely inside the plot.
    ax.set_xlim(max(6, Ns_p.min() - 1), n_max + 15)
    ax.set_ylim(m_min - 0.002, m_max + 0.002)
    ax.tick_params(labelsize=11)

    title_main = f"{n_sparse + n_dense + lehmer_mask.sum():,} entries"
    title_sub  = (f"$N \\leq {n_max}$,  $M \\in [{m_min:.2f}, {m_max:.2f}]$")
    ax.set_title(f"AllKnownAdvanpix: $M(P)$ vs $N$ — {title_main}\n{title_sub}",
                 fontsize=14, pad=10)

    leg = ax.legend(loc="lower right", fontsize=11, framealpha=0.92,
                    markerscale=3)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close()
    print(f"wrote {out_path} (sparse={n_sparse:,}, dense={n_dense:,})",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variant", choices=["two-bands", "three-bands",
                                          "full"],
                    default="two-bands")
    ap.add_argument("--xscale", choices=["log", "linear"], default="linear")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--n-max", type=int, default=None,
                    help="override the variant's N truncation")
    ap.add_argument("--m-min", type=float, default=None)
    ap.add_argument("--m-max", type=float, default=None)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    presets = {
        "two-bands":   dict(n_max=456,  m_min=1.17, m_max=1.30),
        "three-bands": dict(n_max=456,  m_min=1.17, m_max=1.32),
        "full":        dict(n_max=1500, m_min=1.17, m_max=1.50),
    }
    p = presets[args.variant]
    n_max = args.n_max if args.n_max is not None else p["n_max"]
    m_min = args.m_min if args.m_min is not None else p["m_min"]
    m_max = args.m_max if args.m_max is not None else p["m_max"]

    out = args.output
    if out is None:
        suffix = f"-{args.xscale}" if args.xscale != "linear" else ""
        out = REPO / "work" / f"db-bands-{args.variant}{suffix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    render(args.variant, n_max, m_min, m_max, args.xscale, out,
           dpi=args.dpi)


if __name__ == "__main__":
    main()
