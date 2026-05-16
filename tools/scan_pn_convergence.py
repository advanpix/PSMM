#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
scan_pn_convergence.py — convergence study and DB extraction for the Max-U family

    P_N(x) = (x+1)(x^{N-1} + 1) - x^{N/2 - 1} Phi_3(x),  N even, N >= 6.

For each even N in a configurable range:
  - Compute M(P_N) via PARI polroots (no factoring; factoring is the
    bottleneck past N ~ 1000).
  - Classify roots into K (outside |z|=1), U (on the circle, within eps),
    Q (off-circle complex), R (off-circle real).
  - Sanity-check 2K + U = N and Q + R = 2K.

Two outputs are written incrementally:

  --output (default doc/parametric-family-convergence.csv):
      N, M_PN, abs_diff_from_limit, gap_from_previous, K, U, Q, R
  (convergence-rate study)

  --db-output (default doc/new_finds_d_only_pn_extended.txt):
      AllKnownAdvanpix-format entries (one per line):
      N M_PN NNZ H L K U Q R c_0 c_1 ... c_{N/2}
  (ready for `psmm -merge` into AllKnownAdvanpix)

The half-coefficient vector for P_N is sparse with a known shape:
  c_0 = c_N = 1, c_1 = c_{N-1} = 1, c_{N/2 +/- 1} = -1, c_{N/2} = -1, rest 0.
So NNZ=3, H=1, L=7 for every member. No analytic uncertainty there; the
PARI step is only needed for K, U, Q, R (and the high-precision M).

For all N >= ~100 in the family, P_N is empirically irreducible. The
script does not run polisirreducible() (too slow at large N); a separate
verification pass via tools/bulk_verify.py is the safety net.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent

# Corrected Boyd-Lawton limit for the Max-U family P_N
# (see tools/compute_boyd_lawton.py for the 1D Jensen-reduction proof).
BOYD_LAWTON_LIMIT = "1.2554340377272518"


def gp_script(N: int, precision: int) -> str:
    """gp script: compute M(P_N), K, U, Q, R via polroots — no factoring."""
    half = N // 2 - 1
    classify = (
        "r = abs(rts[i]); "
        "if (abs(r - 1) < eps, U += 1, "
        "    if (r > 1, M *= r; K += 1); "
        "    if (abs(imag(rts[i])) < eps, R += 1, Q += 1) )"
    )
    return (
        f"default(realprecision, {precision});\n"
        f"P = (x+1)*(x^{N-1} + 1) - x^{half} * (x^2 + x + 1);\n"
        "rts = polroots(P);\n"
        "M = 1.0; K = 0; U = 0; Q = 0; R = 0;\n"
        "eps = 1e-30;\n"
        f"for (i = 1, #rts, {classify});\n"
        "print(M, \" \", K, \" \", U, \" \", Q, \" \", R);\n"
    )


def db_line_for(N: int, M_str: str, K: int, U: int, Q: int, R: int) -> str:
    """Build an AllKnownAdvanpix-format line for P_N (sparse 5-term)."""
    half_len = N // 2 + 1
    # Initialise to zero; then set the known non-zero half-coefficients.
    half = [0] * half_len
    half[0] = 1                     # c_0 = leading = 1
    half[1] = 1                     # c_1 = 1
    half[N // 2 - 1] = -1           # c_{N/2 - 1} = -1 (off-middle)
    half[N // 2] = -1               # c_{N/2}     = -1 (middle)
    # NNZ, H, L from the known sparse structure
    NNZ = 3
    H   = 1
    L   = 7
    coeffs = " ".join(str(c) for c in half)
    return f"{N} {M_str} {NNZ} {H} {L} {K} {U} {Q} {R} {coeffs}"


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
                    default=str(REPO / "doc" / "parametric-family-convergence.csv"),
                    help="convergence-study CSV path")
    ap.add_argument("--db-output",
                    default=str(REPO / "doc" / "new_finds_d_only_pn_extended.txt"),
                    help="AllKnownAdvanpix-format entries path "
                         "(set to empty string to disable DB output)")
    args = ap.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    db_path = Path(args.db_output) if args.db_output else None
    if db_path:
        db_path.parent.mkdir(parents=True, exist_ok=True)

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
    n_db_emitted = 0
    n_sanity_failed = 0
    last_M = None  # for gap_from_previous

    db_f = db_path.open("w", buffering=1) if db_path else None
    try:
        with out_path.open("w", buffering=1) as out:
            w = csv.writer(out)
            w.writerow(["N", "M_PN", "abs_diff_from_limit",
                        "gap_from_previous", "K", "U", "Q", "R"])

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
                    w.writerow([N, "timeout", "", "", "", "", "", ""])
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
                # gp output is "M K U Q R" on a single line
                parts = lines[-1].split()
                if len(parts) < 5:
                    print(f"N={N}: malformed gp output: {lines[-1][:80]}",
                          file=sys.stderr, flush=True)
                    last_M = None
                    continue
                M_str = parts[0]
                try:
                    K = int(parts[1])
                    U = int(parts[2])
                    Q = int(parts[3])
                    R = int(parts[4])
                except ValueError:
                    print(f"N={N}: KUQR parse error", file=sys.stderr, flush=True)
                    last_M = None
                    continue

                # Sanity: 2K + U = N and Q + R = 2K
                ok = (2*K + U == N) and (Q + R == 2*K)
                if not ok:
                    n_sanity_failed += 1
                    print(f"N={N}: SANITY FAIL  2K+U={2*K+U} (want {N}), "
                          f"Q+R={Q+R} (want {2*K}) "
                          f"— likely reducible; skipping DB emit",
                          file=sys.stderr, flush=True)

                # Convergence CSV row
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
                    w.writerow([N, M_str, f"{abs_diff:.6e}", gap_str,
                                K, U, Q, R])
                except ValueError:
                    w.writerow([N, M_str, "", "", K, U, Q, R])
                    last_M = None

                # DB entry (only if sanity passed AND within DB scope:
                #   1.001 < M < 1.3, i.e. has a non-cyclotomic factor and
                #   sits inside AllKnownAdvanpix's M < 1.3 coverage.)
                if db_f and ok and 1.001 < M_f < 1.3:
                    db_f.write(db_line_for(N, M_str, K, U, Q, R) + "\n")
                    n_db_emitted += 1

                n_done += 1
                elapsed = time.time() - started
                rate = n_done / elapsed if elapsed > 0 else 0
                print(f"N={N:6d} M={M_str[:14]} K={K:>3d} U={U:>4d} "
                      f"Q={Q:>3d} R={R:>2d}  |M-L|={abs_diff:.3e}  "
                      f"({n_done} done, {rate:.2f}/s, "
                      f"elapsed {elapsed/60:.1f} min)",
                      file=sys.stderr, flush=True)
    finally:
        if db_f:
            db_f.close()

    elapsed = time.time() - started
    print(f"\ndone: {n_done} N values processed in {elapsed/60:.1f} min",
          file=sys.stderr)
    print(f"  DB entries emitted: {n_db_emitted}", file=sys.stderr)
    print(f"  sanity check failures: {n_sanity_failed}", file=sys.stderr)


if __name__ == "__main__":
    main()
