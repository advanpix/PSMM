#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
factor_5term_reducibles.py — backfill the dense factor-extract band
missed by tools/scan_5term_salem.py.

Reads work/scan-5term-all-combos.csv, identifies reducible candidates
(irr=0) with M in a configurable range, factors each
  P_{N,a,s1,s2}(x) = 1 + s1*x^a + s2*x^{N/2} + s1*x^{N-a} + x^N
over Z, and for every non-cyclotomic reciprocal factor F emits:

  --factor-log   one CSV row per (parent, factor) pair, ALWAYS (no
                 threshold). Full audit trail.
  --db-output    DB-format coefficient line for the factor, gated by
                 1.001 < M(F) < --m-max-db (default 1.3).
  --roots-output canonical roots-block file (header = the DB-format
                 line, body = one 'real imag' per root), gated by
                 1.001 < M(F) < --m-max-log (default 1.5). Looser
                 than the DB gate so future relaxation of the cutoff
                 doesn't require rerunning this scan.

Parallel via multiprocessing.Pool (one candidate per task). Tasks are
deterministic; safe to resume by leaving the partial output files in
place and re-running (the tool re-reads them on startup and skips
candidates already processed).

Per-candidate compute is fast (factor + polroots of each factor) since
the factors are smaller than the parent P, so the slow polroots(P) of
the original scan is NOT redone. Expected runtime on 6,838 candidates
with 16 workers: tens of minutes to a few hours, family-dependent.
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
    """Round M to `digits` fractional digits, ROUND_HALF_EVEN. Mirrors
    scan_5term_salem.round_to_db_format / scan_pn_convergence's logic."""
    if "." not in M_str:
        return M_str
    int_part, frac = M_str.split(".", 1)
    if len(frac) <= digits:
        if len(frac) < digits:
            frac = frac + "0" * (digits - len(frac))
        return int_part + "." + frac
    getcontext().prec = max(len(int_part) + len(frac) + 4, digits + 8)
    d = Decimal(M_str)
    quant = Decimal("1." + "0" * digits)
    rounded = d.quantize(quant, rounding=ROUND_HALF_EVEN)
    s = str(rounded)
    if "." in s:
        ip, fp = s.split(".", 1)
        if len(fp) < digits:
            fp = fp + "0" * (digits - len(fp))
        return ip + "." + fp
    return s + "." + "0" * digits


def gp_factor_script(N: int, a: int, s1: int, s2: int,
                     m_max_log: str, precision: int) -> str:
    """gp script: build P_{N,a,s1,s2}, factor over Z, emit each
    non-cyclotomic factor's degree / M / K / U / Q / R / coefficients
    / roots. Threshold to drop trivially-cyclotomic factors is
    1.001 < M < m_max_log; the DB-cutoff filter is applied later in
    Python so we don't have to rerun if the cutoff changes.

    Output is a sequence of blocks separated by FACTOR_HEADER lines:

      FACTOR_HEADER <deg> <M_str> <K> <U> <Q> <R>
      COEF c_0
      COEF c_1
      ...
      COEF c_deg
      ROOT re im
      ROOT re im
      ...
      ROOT re im
      [next FACTOR_HEADER ...]
      END
    """
    s1c = "+" if s1 > 0 else "-"
    s2c = "+" if s2 > 0 else "-"
    half = N // 2
    co = N - a
    poly = f"1 {s1c} x^{a} {s2c} x^{half} {s1c} x^{co} + x^{N}"
    return (
        f"default(realprecision, {precision});\n"
        f"P = {poly};\n"
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
    """Parse the gp output into a list of factor records. Each record:
       dict(deg, M_str, K, U, Q, R, coefs: list[int], roots: list[(re,im)])"""
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
                current["coefs"].append(int(line.split(maxsplit=1)[1].strip()))
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


def process_candidate(args_tuple):
    """Worker: factor one candidate P. Returns
       (N, a, s1, s2, status, factors_list_or_errmsg)."""
    N, a, s1, s2, m_max_log, precision, timeout = args_tuple
    script = gp_factor_script(N, a, s1, s2, m_max_log, precision)
    try:
        proc = subprocess.run(
            ["gp", "-q", "--default", "parisize=4000000000"],
            input=script, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return (N, a, s1, s2, "timeout", "gp timed out")
    if proc.returncode != 0:
        return (N, a, s1, s2, "error",
                f"gp rc={proc.returncode}: {proc.stderr[:200]}")
    factors = parse_gp_output(proc.stdout)
    return (N, a, s1, s2, "ok", factors)


def is_palindromic(coefs: list[int]) -> bool:
    n = len(coefs) - 1
    for i in range(n // 2 + 1):
        if coefs[i] != coefs[n - i]:
            return False
    return True


def db_line_for_factor(factor: dict) -> str | None:
    """Build an AllKnownAdvanpix-format line for a factor. Returns None
    if the factor isn't palindromic, has wrong degree parity, or fails
    reciprocity. Caller is expected to have already gated on M."""
    coefs = factor["coefs"]
    deg = factor["deg"]
    if len(coefs) != deg + 1:
        return None
    if deg % 2 != 0:
        return None
    if not is_palindromic(coefs):
        return None
    K = factor["K"]; U = factor["U"]; Q = factor["Q"]; R = factor["R"]
    if 2 * K + U != deg or Q + R != 2 * K:
        return None  # reciprocity sanity failed
    half = coefs[:deg // 2 + 1]
    NNZ = sum(1 for c in half[1:] if c != 0)
    H = max(abs(c) for c in coefs) if coefs else 0
    L = sum(abs(c) for c in coefs)
    M_db = round_to_db_format(factor["M_str"])
    coef_str = " ".join(str(c) for c in half)
    return f"{deg} {M_db} {NNZ} {H} {L} {K} {U} {Q} {R} {coef_str}"


def load_already_processed(factor_log_path: Path):
    """Return set of (N, a, s1, s2) tuples already in the factor log."""
    seen = set()
    if not factor_log_path.exists():
        return seen
    with factor_log_path.open() as f:
        for row in csv.reader(f):
            if len(row) < 4:
                continue
            try:
                seen.add((int(row[0]), int(row[1]),
                          int(row[2]), int(row[3])))
            except ValueError:
                continue
    return seen


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path,
                    default=REPO / "work" / "scan-5term-all-combos.csv",
                    help="scan CSV with the irr column (default: "
                         "work/scan-5term-all-combos.csv)")
    ap.add_argument("--m-max-db", type=float, default=1.3,
                    help="DB-format emission cutoff (default 1.3, "
                         "matches the project DB threshold)")
    ap.add_argument("--m-max-log", type=str, default="1.5",
                    help="logging cutoff: factors with M < this go to "
                         "the factor-log and roots files even if they "
                         "exceed --m-max-db (default 1.5; '0' or '' "
                         "for no upper cap)")
    ap.add_argument("--m-min", type=float, default=1.001,
                    help="lower cutoff (default 1.001 to skip "
                         "cyclotomic-only factors)")
    ap.add_argument("--precision", type=int, default=120)
    ap.add_argument("--timeout", type=int, default=900,
                    help="per-candidate gp timeout (default 900s)")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--factor-log", type=Path,
                    default=REPO / "work" / "factor_5term" / "factor-log.csv",
                    help="per-factor CSV audit log (no threshold)")
    ap.add_argument("--db-output", type=Path,
                    default=REPO / "work" / "factor_5term" / "new_finds_5term_factors.txt",
                    help="DB-format file for merge into AllKnownAdvanpix "
                         "(gated by --m-max-db)")
    ap.add_argument("--roots-output", type=Path,
                    default=REPO / "work" / "factor_5term" / "roots-5term-factors.txt",
                    help="canonical roots-block file (gated by "
                         "--m-max-log)")
    ap.add_argument("--limit", type=int, default=0,
                    help="for smoke-testing: process at most this many "
                         "candidates (0 = all)")
    args = ap.parse_args()

    # Normalise the m_max_log string for gp: 'auto'/'0'/'' means no cap.
    m_max_log_val = args.m_max_log.strip() if args.m_max_log else ""
    if m_max_log_val in ("", "0", "0.0", "none", "inf"):
        m_max_log_str = "100.0"
        m_max_log_num = 100.0
    else:
        m_max_log_str = m_max_log_val
        m_max_log_num = float(m_max_log_val)

    # Ensure output dirs exist.
    args.factor_log.parent.mkdir(parents=True, exist_ok=True)
    args.db_output.parent.mkdir(parents=True, exist_ok=True)
    args.roots_output.parent.mkdir(parents=True, exist_ok=True)

    # Resume: load already-processed candidates.
    seen = load_already_processed(args.factor_log)
    print(f"already processed candidates: {len(seen)}", file=sys.stderr)

    # Build the candidate list from the CSV.
    candidates = []
    with args.csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                N = int(row["N"])
                a = int(row["a"])
                s1 = int(row["s1"])
                s2 = int(row["s2"])
                M_str = row.get("M", "")
                irr = row.get("irr", "").strip()
            except (KeyError, ValueError):
                continue
            if irr != "0":
                continue  # only reducible candidates
            try:
                M_f = float(M_str)
            except ValueError:
                continue
            if not (args.m_min < M_f < m_max_log_num):
                continue
            if (N, a, s1, s2) in seen:
                continue
            candidates.append((N, a, s1, s2))
    if args.limit > 0:
        candidates = candidates[:args.limit]
    print(f"candidates to process: {len(candidates)}", file=sys.stderr)

    if not candidates:
        print("nothing to do.", file=sys.stderr)
        return

    # Open output files in append mode (so resumes don't clobber).
    factor_log_f = args.factor_log.open("a", buffering=1, newline="")
    factor_log_w = csv.writer(factor_log_f)
    if args.factor_log.stat().st_size == 0:
        factor_log_w.writerow([
            "N", "a", "s1", "s2",
            "factor_deg", "factor_M", "factor_K", "factor_U",
            "factor_Q", "factor_R", "palindromic",
            "in_db", "in_log",
            "coefs",
        ])
    db_f = args.db_output.open("a", buffering=1)
    roots_f = args.roots_output.open("a", buffering=1)

    started = time.time()
    n_done = 0
    n_total = len(candidates)
    n_db_emitted = 0
    n_log_emitted = 0
    n_skip_nonpali = 0
    n_skip_reciprocity = 0
    n_errors = 0

    tasks = [
        (N, a, s1, s2, m_max_log_str, args.precision, args.timeout)
        for (N, a, s1, s2) in candidates
    ]

    def handle_result(result):
        nonlocal n_done, n_db_emitted, n_log_emitted, n_skip_nonpali
        nonlocal n_skip_reciprocity, n_errors
        N, a, s1, s2, status, payload = result
        n_done += 1
        if status != "ok":
            n_errors += 1
            factor_log_w.writerow([N, a, s1, s2, "", "", "", "", "", "",
                                   "", "", "", f"<{status}: {payload}>"])
            return
        factors = payload
        for fac in factors:
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
                N, a, s1, s2,
                fac["deg"], fac["M_str"], fac["K"], fac["U"],
                fac["Q"], fac["R"], int(pal),
                int(in_db), int(in_log),
                " ".join(str(c) for c in fac["coefs"]),
            ])

            if not pal:
                n_skip_nonpali += 1
                continue
            if not reciprocity_ok:
                n_skip_reciprocity += 1
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
        if n_done % 50 == 0 or n_done == n_total:
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
            for result in pool.imap_unordered(process_candidate, tasks):
                handle_result(result)
    else:
        for task in tasks:
            handle_result(process_candidate(task))

    elapsed = time.time() - started
    print(f"\ndone in {elapsed/60:.1f} min: "
          f"{n_done} candidates processed, "
          f"{n_db_emitted} DB-format emissions, "
          f"{n_log_emitted} roots-block emissions, "
          f"{n_skip_nonpali} non-palindromic skipped, "
          f"{n_skip_reciprocity} reciprocity-fail skipped, "
          f"{n_errors} errors",
          file=sys.stderr)
    factor_log_f.close()
    db_f.close()
    roots_f.close()


if __name__ == "__main__":
    main()
