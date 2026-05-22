#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
factor_2_12_minus_reducibles.py — recover (2, 12, -1) family entries
in the M >= 1.3 region missed by the original scan.

The (2, 12, -1) family polynomial is

    P_N(x) = (x + 1) * (x^(N-1) + 1) - x^((N-4)/2) * Phi_12(x),
             Phi_12(x) = x^4 - x^2 + 1

with even N >= 6. The Boyd-Lawton limit is L_{(2,12)} = 1.30911..., so
the bulk of the family has M > 1.3 and was filtered out by the
original scan_pn_convergence.py (--m-max-db was effectively hardcoded
at 1.3 at the time). Their roots ARE in work/roots-2-12-minus.txt
(roots were emitted regardless of M), so we have everything we need
without redoing polroots(P).

For each N in the configurable M range:
  - factor P_N over Z
  - for each non-cyclotomic palindromic factor F:
      compute polroots(F) and the classification, emit DB-format line
      + canonical roots-block.

Mirrors tools/factor_5term_reducibles.py in structure; family-specific
piece is only the polynomial construction.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import subprocess
import sys
import time
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_MAHLER_DIGITS = 72


def round_to_db_format(M_str: str, digits: int = DB_MAHLER_DIGITS) -> str:
    if "." not in M_str:
        return M_str
    int_part, frac = M_str.split(".", 1)
    if len(frac) <= digits:
        return int_part + "." + frac.ljust(digits, "0")
    getcontext().prec = max(len(int_part) + len(frac) + 4, digits + 8)
    d = Decimal(M_str)
    quant = Decimal("1." + "0" * digits)
    rounded = d.quantize(quant, rounding=ROUND_HALF_EVEN)
    s = str(rounded)
    if "." in s:
        ip, fp = s.split(".", 1)
        return ip + "." + fp.ljust(digits, "0")
    return s + "." + "0" * digits


def gp_factor_script(N: int, m_max_log: str, precision: int) -> str:
    """gp script: build P_N for (a=2, d=12, sign=-1), factor, emit each
    non-cyclotomic factor's degree / M / K / U / Q / R / coefficients
    / roots.

    P_N(x) = (x+1)*(x^(N-1)+1) - x^((N-4)/2) * (x^4 - x^2 + 1).
    """
    shift = (N - 4) // 2
    poly_expr = (
        f"(x + 1) * (x^{N - 1} + 1) "
        f"- x^{shift} * (x^4 - x^2 + 1)"
    )
    return (
        f"default(realprecision, {precision});\n"
        f"P = {poly_expr};\n"
        "fac = factor(P); nf = #fac~;\n"
        "for (i = 1, nf, "
        "  F = fac[i, 1]; deg_F = poldegree(F); "
        "  if (deg_F < 2, next); "
        "  rts = polroots(F); "
        "  M = 1.0; K = 0; U = 0; Q = 0; R = 0; eps = 1e-30; "
        "  for (j = 1, #rts, "
        "    r = abs(rts[j]); "
        "    if (abs(r - 1) < eps, U += 1, "
        "      if (r > 1, M *= r; K += 1); "
        "      if (abs(imag(rts[j])) < eps, R += 1, Q += 1) ) ); "
        f"  if (M <= 1.001 || M >= {m_max_log}, next); "
        "  print(\"FACTOR_HEADER \", deg_F, \" \", M, \" \", K, "
        "        \" \", U, \" \", Q, \" \", R); "
        "  for (k = 0, deg_F, print(\"COEF \", polcoef(F, k))); "
        "  for (j = 1, #rts, "
        "    print(\"ROOT \", real(rts[j]), \" \", imag(rts[j]))) "
        ");\n"
        "print(\"END\");\n"
    )


def parse_gp_output(stdout: str):
    factors = []
    current = None
    for raw in stdout.splitlines():
        line = raw.strip()
        if not line or line.startswith("***"):
            continue
        if line.startswith("FACTOR_HEADER"):
            if current is not None:
                factors.append(current)
            toks = line.split()
            try:
                current = {
                    "deg": int(toks[1]),
                    "M_str": toks[2],
                    "K": int(toks[3]),
                    "U": int(toks[4]),
                    "Q": int(toks[5]),
                    "R": int(toks[6]),
                    "coefs": [],
                    "roots": [],
                }
            except (ValueError, IndexError):
                current = None
        elif line.startswith("COEF") and current is not None:
            try:
                current["coefs"].append(int(line.split(maxsplit=1)[1]))
            except (ValueError, IndexError):
                pass
        elif line.startswith("ROOT") and current is not None:
            rest = line.split(maxsplit=1)[1].strip()
            parts = rest.split()
            if len(parts) == 2:
                current["roots"].append((parts[0], parts[1]))
        elif line == "END":
            break
    if current is not None:
        factors.append(current)
    return factors


def process_N(args_tuple):
    N, m_max_log, precision, timeout = args_tuple
    script = gp_factor_script(N, m_max_log, precision)
    try:
        proc = subprocess.run(
            ["gp", "-q", "--default", "parisize=4000000000"],
            input=script, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return (N, "timeout", "gp timed out")
    if proc.returncode != 0:
        return (N, "error",
                f"gp rc={proc.returncode}: {proc.stderr[:200]}")
    return (N, "ok", parse_gp_output(proc.stdout))


def is_palindromic(coefs: list[int]) -> bool:
    n = len(coefs) - 1
    for i in range(n // 2 + 1):
        if coefs[i] != coefs[n - i]:
            return False
    return True


def db_line_for_factor(factor: dict) -> str | None:
    coefs = factor["coefs"]
    deg = factor["deg"]
    if len(coefs) != deg + 1 or deg % 2 != 0:
        return None
    if not is_palindromic(coefs):
        return None
    K = factor["K"]; U = factor["U"]; Q = factor["Q"]; R = factor["R"]
    if 2 * K + U != deg or Q + R != 2 * K:
        return None
    half = coefs[:deg // 2 + 1]
    NNZ = sum(1 for c in half[1:] if c != 0)
    H = max(abs(c) for c in coefs) if coefs else 0
    L = sum(abs(c) for c in coefs)
    M_db = round_to_db_format(factor["M_str"])
    coef_str = " ".join(str(c) for c in half)
    return f"{deg} {M_db} {NNZ} {H} {L} {K} {U} {Q} {R} {coef_str}"


def collect_N_targets(roots_file: Path, m_min: float, m_max_log: float):
    """Walk the existing (2,12,-1) roots file headers and return the
    list of N values whose P_N has M in [m_min, m_max_log). These are
    the candidates that the original scan computed but failed to
    emit to the DB."""
    targets = []
    with roots_file.open() as f:
        for line in f:
            if not line.startswith("# "):
                continue
            toks = line[2:].strip().split()
            if len(toks) < 10:
                continue
            try:
                N = int(toks[0])
                M = float(toks[1])
            except (ValueError, IndexError):
                continue
            if m_min <= M < m_max_log:
                targets.append((N, toks[1]))
    return targets


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--roots-input", type=Path,
                    default=REPO / "work" / "roots-2-12-minus.txt",
                    help="existing (2,12,-1) scan roots file used to "
                         "identify which N values fall in the recovery "
                         "M range (default work/roots-2-12-minus.txt)")
    ap.add_argument("--m-min-recover", type=float, default=1.3,
                    help="recover entries with M >= this. Default 1.3 "
                         "(the historical DB threshold the original "
                         "scan filtered against).")
    ap.add_argument("--m-max-db", type=float, default=1.3,
                    help="DB-format emission cutoff (default 1.3 matches "
                         "the historical DB threshold; raise to widen "
                         "DB merge eligibility)")
    ap.add_argument("--m-max-log", type=str, default="1.5",
                    help="logging cutoff (default 1.5; covers the "
                         "(2,12,-1) family which converges to 1.30911)")
    ap.add_argument("--m-min", type=float, default=1.001,
                    help="lower M cutoff (default 1.001)")
    ap.add_argument("--precision", type=int, default=120)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--factor-log", type=Path,
                    default=REPO / "work" / "factor_2_12_minus" / "factor-log.csv")
    ap.add_argument("--db-output", type=Path,
                    default=REPO / "work" / "factor_2_12_minus" / "new_finds_2_12_minus_factors.txt")
    ap.add_argument("--roots-output", type=Path,
                    default=REPO / "work" / "factor_2_12_minus" / "roots-2-12-minus-factors.txt")
    ap.add_argument("--limit", type=int, default=0,
                    help="for smoke-testing: process at most this many N values")
    args = ap.parse_args()

    m_max_log_str = args.m_max_log.strip() or "100.0"
    m_max_log_num = 100.0 if m_max_log_str in ("0", "0.0", "none", "inf") else float(m_max_log_str)

    args.factor_log.parent.mkdir(parents=True, exist_ok=True)
    args.db_output.parent.mkdir(parents=True, exist_ok=True)
    args.roots_output.parent.mkdir(parents=True, exist_ok=True)

    targets = collect_N_targets(args.roots_input, args.m_min_recover,
                                m_max_log_num)
    if args.limit > 0:
        targets = targets[:args.limit]
    print(f"recovery targets (N values with {args.m_min_recover} <= M < "
          f"{m_max_log_num}): {len(targets)}", file=sys.stderr)

    if not targets:
        print("nothing to recover.", file=sys.stderr)
        return

    factor_log_f = args.factor_log.open("w", buffering=1, newline="")
    factor_log_w = csv.writer(factor_log_f)
    factor_log_w.writerow([
        "N_parent", "factor_deg", "factor_M", "factor_K", "factor_U",
        "factor_Q", "factor_R", "palindromic", "in_db", "in_log", "coefs",
    ])
    db_f = args.db_output.open("w", buffering=1)
    roots_f = args.roots_output.open("w", buffering=1)

    started = time.time()
    n_done = 0
    n_db_emitted = 0
    n_log_emitted = 0
    n_errors = 0
    n_skip_nonpali = 0
    n_total = len(targets)
    tasks = [(N, m_max_log_str, args.precision, args.timeout) for N, _ in targets]

    def handle(result):
        nonlocal n_done, n_db_emitted, n_log_emitted, n_errors, n_skip_nonpali
        N, status, payload = result
        n_done += 1
        if status != "ok":
            n_errors += 1
            factor_log_w.writerow([N, "", "", "", "", "", "", "", "", "",
                                   f"<{status}: {payload}>"])
            return
        for fac in payload:
            try:
                M_f = float(fac["M_str"])
            except ValueError:
                M_f = float("nan")
            pal = is_palindromic(fac["coefs"])
            reciprocity_ok = (2 * fac["K"] + fac["U"] == fac["deg"]
                              and fac["Q"] + fac["R"] == 2 * fac["K"])
            in_db = (pal and reciprocity_ok
                     and args.m_min < M_f < args.m_max_db
                     and fac["deg"] % 2 == 0)
            in_log = pal and reciprocity_ok and args.m_min < M_f < m_max_log_num

            factor_log_w.writerow([
                N, fac["deg"], fac["M_str"], fac["K"], fac["U"],
                fac["Q"], fac["R"], int(pal), int(in_db), int(in_log),
                " ".join(str(c) for c in fac["coefs"]),
            ])

            if not pal:
                n_skip_nonpali += 1
                continue
            if not reciprocity_ok:
                continue

            db_line = db_line_for_factor(fac)
            if db_line is None:
                continue
            if in_db:
                db_f.write(db_line + "\n")
                n_db_emitted += 1
            if in_log:
                roots_f.write("# " + db_line + "\n")
                for re_s, im_s in fac["roots"]:
                    roots_f.write(f"{re_s} {im_s}\n")
                roots_f.write("\n")
                n_log_emitted += 1

        if n_done % 25 == 0 or n_done == n_total:
            elapsed = time.time() - started
            rate = n_done / elapsed if elapsed > 0 else 0
            eta = (n_total - n_done) / rate if rate > 0 else 0
            print(f"  [{n_done}/{n_total} = {100*n_done/n_total:.1f}%, "
                  f"{elapsed/60:.1f} min, ETA {eta/60:.1f} min] "
                  f"db_emit={n_db_emitted} log_emit={n_log_emitted} "
                  f"nonpali={n_skip_nonpali} err={n_errors}",
                  file=sys.stderr, flush=True)

    if args.workers > 1:
        with mp.Pool(args.workers) as pool:
            for result in pool.imap_unordered(process_N, tasks):
                handle(result)
    else:
        for task in tasks:
            handle(process_N(task))

    elapsed = time.time() - started
    print(f"\ndone in {elapsed/60:.1f} min: "
          f"{n_done} targets processed, "
          f"{n_db_emitted} DB-format emissions, "
          f"{n_log_emitted} roots-block emissions, "
          f"{n_skip_nonpali} non-palindromic skipped, "
          f"{n_errors} errors", file=sys.stderr)
    factor_log_f.close()
    db_f.close()
    roots_f.close()


if __name__ == "__main__":
    main()
