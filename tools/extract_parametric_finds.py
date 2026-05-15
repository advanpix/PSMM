#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
extract_parametric_finds.py — convert parametric-sweep CSV rows into
AllKnownAdvanpix-format polynomial entries, ready for `psmm -merge`.

For each row of the input CSV (--src) where:
  - the smallest non-cyclotomic factor exists,
  - its Mahler measure is < 1.3,
  - it is NOT already in AllKnownAdvanpix (in_db / factor_in_database == False),
re-construct the polynomial via PARI at high precision, factor over Z,
isolate the smallest non-cyclotomic factor F, and emit one DB-format line:

    N M NNZ H L K U Q R c_0 c_1 ... c_{N/2}

where c_0 is the leading coefficient.

Two input CSV layouts are supported (auto-detected from the header):
  - d-only sweep (`scan_parametric_family.py`): m, deg_Pm, ..., factor_in_database
  - (a,d) sweep (`sweep_ad.py`): a, d, sign, m, ..., in_db
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"


CYCLOTOMICS = {
    2:  ("(x+1)", 1),
    3:  ("(x^2+x+1)", 2),
    4:  ("(x^2+1)", 2),
    5:  ("(x^4+x^3+x^2+x+1)", 4),
    6:  ("(x^2-x+1)", 2),
    7:  ("(x^6+x^5+x^4+x^3+x^2+x+1)", 6),
    8:  ("(x^4+1)", 4),
    9:  ("(x^6+x^3+1)", 6),
    10: ("(x^4-x^3+x^2-x+1)", 4),
    12: ("(x^4-x^2+1)", 4),
    15: ("(x^8-x^7+x^5-x^4+x^3-x+1)", 8),
}


def construct_d_only(m: int) -> str:
    """P_m(x) = (x+1)(x^m + 1) - x^((m-1)/2) (x^2 + x + 1)"""
    return f"(x+1)*(x^{m} + 1) - x^{(m - 1) // 2} * (x^2 + x + 1)"


def construct_ad(a: int, d: int, m: int, sign: int) -> str:
    """P_{a,d,m,s}(x) = Phi_a(x)(x^m + 1) + s x^k Phi_d(x)
       k = (phi(a) + m - phi(d)) / 2"""
    phi_a, pa = CYCLOTOMICS[a]
    phi_d, pd = CYCLOTOMICS[d]
    k = (pa + m - pd) // 2
    s_str = "+" if sign > 0 else "-"
    return f"{phi_a} * (x^{m} + 1) {s_str} x^{k} * {phi_d}"


_CLASSIFY_BODY = (
    "r = abs(rts[i]); "
    "if (abs(r - 1) < eps, U += 1, "
    "    d_off = abs(r - 1); if (d_off > maxoff, maxoff = d_off); "
    "    if (r > 1, M *= r; K += 1); "
    "    if (abs(imag(rts[i])) < eps, R += 1, Q += 1) )"
)

_FACTOR_BODY = (
    "F = fac[i,1]; rs = polroots(F); MF = 1.0; "
    "for (j = 1, #rs, if (abs(rs[j]) > 1, MF *= abs(rs[j]))); "
    "if (MF > 1.001 && MF < best_M, best_M = MF; best_idx = i)"
)

GP_EXTRACT = (
    "default(realprecision, %(prec)d);\n"
    "P = %(poly)s;\n"
    "fac = factor(P); nf = #fac~;\n"
    "best_idx = 0; best_M = 1e10;\n"
    f"for (i = 1, nf, {_FACTOR_BODY});\n"
    "if (best_idx == 0, print(\"SKIP no_factor\"); quit);\n"
    "if (best_M >= 1.3, print(\"SKIP M>=1.3\"); quit);\n"
    "F = fac[best_idx, 1]; N = poldegree(F);\n"
    "rts = polroots(F);\n"
    "M = 1.0; K = 0; U = 0; Q = 0; R = 0; eps = 1e-30; maxoff = 0.0;\n"
    f"for (i = 1, #rts, {_CLASSIFY_BODY});\n"
    "if (2*K + U != N, print(\"SKIP KUQR_failed_2K+U=\", 2*K+U, \"_N=\", N); quit);\n"
    "if (Q + R != 2*K, print(\"SKIP KUQR_failed_Q+R=\", Q+R, \"_2K=\", 2*K); quit);\n"
    "if (polcoef(F, N) != polcoef(F, 0), print(\"SKIP non_palindromic\"); quit);\n"
    "half = vector(N\\2 + 1, i, polcoef(F, N - (i-1)));\n"
    "NNZ = 0; for (i = 2, #half, if (half[i] != 0, NNZ += 1));\n"
    "H = 0;   for (i = 1, #half, if (abs(half[i]) > H, H = abs(half[i])));\n"
    "L = 0;   for (i = 0, N, L += abs(polcoef(F, i)));\n"
    "out = Str(N, \" \", M, \" \", NNZ, \" \", H, \" \", L, \" \", K, \" \", U, \" \", Q, \" \", R);\n"
    "for (i = 1, #half, out = Str(out, \" \", half[i]));\n"
    "print(\"OK \", out);\n"
)


def run_gp(poly: str, precision: int = 120, timeout: int = 600):
    proc = subprocess.run(
        ["gp", "-q", "--default", "parisize=2000000000"],
        input=GP_EXTRACT % {"prec": precision, "poly": poly},
        capture_output=True, text=True,
        timeout=timeout,
    )
    out = [ln.strip() for ln in proc.stdout.splitlines()
           if ln.strip() and not ln.startswith("***")]
    err = proc.stderr.strip()
    return out, err, proc.returncode


def detect_layout(header):
    if "factor_in_database" in header:
        return "d_only"
    if "in_db" in header:
        return "ad"
    raise SystemExit(f"unknown CSV layout (header: {header})")


def candidate_filter(row, layout):
    """Return (ok, poly_str, label) or (False, None, reason)."""
    if layout == "d_only":
        m = int(row["m"])
        in_db = row["factor_in_database"].strip()
        M_str = row["M_smallest_factor"]
        deg_F = row["deg_smallest_factor"]
        label = f"d-only m={m}"
        poly = construct_d_only(m)
    else:
        a = int(row["a"]); d = int(row["d"])
        sign = int(row["sign"]); m = int(row["m"])
        in_db = row["in_db"].strip()
        M_str = row["M_smallest_factor"]
        deg_F = row["deg_smallest_factor"]
        label = f"(a={a}, d={d}, sign={sign:+d}, m={m})"
        poly = construct_ad(a, d, m, sign)
    if in_db == "True":
        return False, None, None, "already in DB"
    try:
        M = float(M_str)
    except ValueError:
        return False, None, None, "M parse error"
    if M >= 1.3:
        return False, None, None, "M >= 1.3"
    if M <= 1.001:
        return False, None, None, "no non-cyclo factor"
    if not deg_F or int(deg_F) == 0:
        return False, None, None, "deg = 0"
    return True, poly, label, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True,
                    help="parametric sweep CSV (annotated with in_db)")
    ap.add_argument("--output", required=True,
                    help="output: AllKnownAdvanpix-format polynomial entries")
    ap.add_argument("--log", default=None,
                    help="optional per-row log file")
    ap.add_argument("--precision", type=int, default=120,
                    help="PARI realprecision in decimal digits (default 120)")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after this many candidates (0 = no limit)")
    args = ap.parse_args()

    src = Path(args.src)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    log_f = open(args.log, "w", buffering=1) if args.log else None

    def logmsg(msg):
        print(msg, file=sys.stderr, flush=True)
        if log_f:
            log_f.write(msg + "\n")

    with src.open(newline="") as f:
        reader = csv.DictReader(f)
        layout = detect_layout(reader.fieldnames)
        logmsg(f"layout: {layout}")
        rows = list(reader)
    logmsg(f"{len(rows)} rows in {src}")

    started = time.time()
    n_processed = 0
    n_emitted = 0
    n_skipped = 0
    n_errors = 0
    seen_keys = set()  # within-file dedup by (N, M-prefix-15)

    with out.open("w", buffering=1) as out_f:
        for i, row in enumerate(rows):
            ok, poly, label, reason = candidate_filter(row, layout)
            if not ok:
                continue
            n_processed += 1
            if args.limit and n_processed > args.limit:
                logmsg(f"hit --limit {args.limit}, stopping")
                break
            try:
                out_lines, err, rc = run_gp(poly, args.precision, args.timeout)
            except subprocess.TimeoutExpired:
                n_errors += 1
                logmsg(f"[{label}] TIMEOUT")
                continue
            if rc != 0:
                n_errors += 1
                logmsg(f"[{label}] gp rc={rc}: {err[:100]}")
                continue
            ok_line = next((ln for ln in out_lines if ln.startswith("OK ")),
                           None)
            if not ok_line:
                skip_line = next((ln for ln in out_lines
                                  if ln.startswith("SKIP")), "(no output)")
                n_skipped += 1
                logmsg(f"[{label}] {skip_line}")
                continue
            db_line = ok_line[3:]  # strip "OK "
            toks = db_line.split()
            try:
                N = int(toks[0])
                key = (N, toks[1][:15])
            except (ValueError, IndexError):
                n_errors += 1
                logmsg(f"[{label}] malformed output: {db_line[:80]}")
                continue
            if key in seen_keys:
                logmsg(f"[{label}] dedup (intra-file)")
                continue
            seen_keys.add(key)
            out_f.write(db_line + "\n")
            n_emitted += 1
            if n_emitted % 10 == 0:
                elapsed = time.time() - started
                rate = n_emitted / elapsed if elapsed > 0 else 0
                logmsg(f"  emitted {n_emitted} ({rate:.2f}/s)")

    elapsed = time.time() - started
    logmsg("")
    logmsg(f"done in {elapsed:.1f}s")
    logmsg(f"  candidates processed: {n_processed}")
    logmsg(f"  emitted entries:      {n_emitted}")
    logmsg(f"  skipped (gp SKIP):    {n_skipped}")
    logmsg(f"  errors:               {n_errors}")
    logmsg(f"  output:               {out}")
    if log_f:
        log_f.close()


if __name__ == "__main__":
    main()
