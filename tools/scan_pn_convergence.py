#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
scan_pn_convergence.py — convergence study for the Max-U family

    P_N(x) = (x+1)(x^{N-1} + 1) - x^{N/2 - 1} Phi_3(x),  N even, N >= 6.

For each even N in a configurable range, computes M(P_N) DIRECTLY via
PARI's polroots (no factoring). The point of the family is that
M(P_N) -> L = 1.249358390752959362866...  (Boyd-Lawton limit), so this
script is meant to be pushed to large N to characterise the convergence
rate empirically. Skipping factoring is what makes large N tractable.

Output: CSV at the path given by --output (default
doc/parametric-family-convergence.csv) with columns

    N, M_PN, abs_diff_from_limit, gap_from_previous

where gap_from_previous is |M(P_N) - M(P_{N-2})| (blank for the first row).

The script writes incrementally so progress is visible, and prints a
status line per row to stderr.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent

BOYD_LAWTON_LIMIT = "1.249358390752959362866"   # high-precision string


def gp_script(N: int, precision: int) -> str:
    """gp script that prints just M(P_N) — no factoring."""
    half = N // 2 - 1
    return (
        f"default(realprecision, {precision});\n"
        f"P = (x+1)*(x^{N-1} + 1) - x^{half} * (x^2 + x + 1);\n"
        "rts = polroots(P);\n"
        "M = 1.0;\n"
        f"for (i = 1, #rts, if (abs(rts[i]) > 1, M *= abs(rts[i])));\n"
        "print(M);\n"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-min", type=int, default=6,
                    help="minimum even N (default 6)")
    ap.add_argument("--n-max", type=int, default=2002,
                    help="maximum N (default 2002, i.e. m up to 2001)")
    ap.add_argument("--n-step", type=int, default=2,
                    help="step in N (default 2 — every even value)")
    ap.add_argument("--precision", type=int, default=40,
                    help="PARI realprecision in decimal digits (default 40)")
    ap.add_argument("--timeout", type=int, default=1800,
                    help="per-call gp timeout in seconds (default 1800)")
    ap.add_argument("--output",
                    default=str(REPO / "doc" / "parametric-family-convergence.csv"))
    args = ap.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Ensure even N values (the family is only defined for even N)
    n_min = args.n_min if args.n_min % 2 == 0 else args.n_min + 1
    if n_min < 6:
        n_min = 6
    n_max = args.n_max
    step = args.n_step
    if step % 2 != 0:
        step += 1  # round up to even step so we stay on even N

    started = time.time()
    n_done = 0
    last_M = None  # for gap_from_previous

    with out_path.open("w", buffering=1) as out:
        w = csv.writer(out)
        w.writerow(["N", "M_PN", "abs_diff_from_limit", "gap_from_previous"])

        # Reference limit as Python Decimal-like string for printing only;
        # we compare via PARI for the abs_diff, but for plotting we just
        # store the absolute decimal-string difference using string-arith.
        # Simpler: compute the difference inline in PARI.
        for N in range(n_min, n_max + 1, step):
            try:
                proc = subprocess.run(
                    ["gp", "-q", "--default", "parisize=4000000000"],
                    input=gp_script(N, args.precision),
                    capture_output=True, text=True,
                    timeout=args.timeout,
                )
            except subprocess.TimeoutExpired:
                print(f"N={N}: timeout", file=sys.stderr, flush=True)
                w.writerow([N, "timeout", "", ""])
                last_M = None
                continue
            if proc.returncode != 0:
                print(f"N={N}: gp rc={proc.returncode}: "
                      f"{proc.stderr[:200]}",
                      file=sys.stderr, flush=True)
                last_M = None
                continue
            lines = [ln.strip() for ln in proc.stdout.splitlines()
                     if ln.strip() and not ln.startswith("***")]
            if not lines:
                print(f"N={N}: no output", file=sys.stderr, flush=True)
                last_M = None
                continue
            M_str = lines[-1]

            # Numerical comparisons via Python float (good to ~15 digits;
            # for high-precision arithmetic, do it in the analysis stage).
            try:
                M_f = float(M_str)
                limit_f = float(BOYD_LAWTON_LIMIT)
                abs_diff = abs(M_f - limit_f)
                if last_M is not None:
                    gap_f = abs(M_f - last_M)
                    gap_str = f"{gap_f:.6e}"
                else:
                    gap_str = ""
                last_M = M_f
                w.writerow([N, M_str, f"{abs_diff:.6e}", gap_str])
            except ValueError:
                w.writerow([N, M_str, "", ""])
                last_M = None

            n_done += 1
            elapsed = time.time() - started
            rate = n_done / elapsed if elapsed > 0 else 0
            print(f"N={N:6d} M={M_str[:14]} "
                  f"|M-L|={abs_diff:.3e}  "
                  f"({n_done} done, {rate:.2f}/s, elapsed {elapsed/60:.1f} min)",
                  file=sys.stderr, flush=True)

    elapsed = time.time() - started
    print(f"\ndone: {n_done} values, elapsed {elapsed/60:.1f} min",
          file=sys.stderr)


if __name__ == "__main__":
    main()
