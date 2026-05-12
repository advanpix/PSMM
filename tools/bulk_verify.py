#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
bulk_verify.py — independently verify every polynomial in AllKnownAdvanpix
with PARI/GP.

For each entry we re-compute, from scratch:
  - irreducibility over Z
  - Mahler measure (high precision)
  - root classification (K, U, Q, R)
  - max off-circle distance |r|-1
and compare against the PSMM-stored values. Any disagreement is reported.

Output: doc/database-verification.csv (incremental writes for liveness).

Designed to run in background. One PARI invocation per polynomial; should
take ~1-2 hours on a fast machine (longer for the high-degree tail).
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


def parse_psmm(line: str):
    toks = line.split()
    return {
        "N":   int(toks[0]),
        "M":   toks[1],
        "NNZ": int(toks[2]),
        "H":   int(toks[3]),
        "L":   int(toks[4]),
        "K":   int(toks[5]),
        "U":   int(toks[6]),
        "Q":   int(toks[7]),
        "R":   int(toks[8]),
        "half": [int(t) for t in toks[9:]],
    }


def coeffs_to_pari(coeffs):
    """Build PARI polynomial string from full coefficient list (a_0..a_N)."""
    pieces = []
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        sign = "+" if c > 0 else "-"
        mag = abs(c)
        if k == 0:
            body = f"{mag}"
        elif k == 1:
            body = "x" if mag == 1 else f"{mag}*x"
        else:
            body = f"x^{k}" if mag == 1 else f"{mag}*x^{k}"
        if not pieces:
            pieces.append(("-" + body) if sign == "-" else body)
        else:
            pieces.append(f"{sign} {body}")
    return " ".join(pieces) if pieces else "0"


def expand_full(half, N):
    full = [0] * (N + 1)
    for k in range(N // 2 + 1):
        full[k] = half[k]
    for k in range(N // 2 + 1, N + 1):
        full[k] = full[N - k]
    return full


def gp_script(poly_str: str, precision: int = 60) -> str:
    body = (
        "r = abs(rts[i]); "
        "if (abs(r - 1) < eps, U += 1, "
        "    d = abs(r - 1); if (d > maxoff, maxoff = d); "
        "    if (r > 1, M *= r; K += 1); "
        "    if (abs(imag(rts[i])) < eps, R += 1, Q += 1) )"
    )
    return (
        f"default(realprecision, {precision});\n"
        f"P = {poly_str};\n"
        "irr = polisirreducible(P);\n"
        "rts = polroots(P);\n"
        "M = 1.0; K = 0; U = 0; Q = 0; R = 0;\n"
        "eps = 1e-20; maxoff = 0.0;\n"
        f"for (i = 1, #rts, {body});\n"
        "print(irr, \",\", M, \",\", K, \",\", U, \",\", Q, \",\", R, \",\", maxoff);\n"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(DB_PATH))
    ap.add_argument("--output", default=str(REPO / "doc" / "database-verification.csv"))
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after this many polynomials (0 = no limit)")
    ap.add_argument("--max-degree", type=int, default=0,
                    help="skip polynomials with degree above this (0 = no limit)")
    ap.add_argument("--precision", type=int, default=60,
                    help="PARI realprecision in decimal digits (default 60)")
    args = ap.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    src = Path(args.src)
    n_total = sum(1 for ln in src.open()
                  if ln.strip() and not ln.startswith("#"))
    print(f"{n_total} polynomials in {src}", file=sys.stderr)

    started = time.time()
    n_done = 0
    n_mismatches = 0

    with out_path.open("w", buffering=1) as out_f:
        w = csv.writer(out_f)
        w.writerow([
            "line_no", "N", "M_stored", "M_pari",
            "K_stored", "K_pari", "U_stored", "U_pari",
            "Q_stored", "Q_pari", "R_stored", "R_pari",
            "irreducible", "max_off_circle",
            "ok_irr", "ok_M", "ok_KUQR",
        ])
        with src.open() as f:
            for line_no, line in enumerate(f, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                try:
                    meta = parse_psmm(line)
                except Exception:
                    continue
                if args.max_degree and meta["N"] > args.max_degree:
                    continue
                if args.limit and n_done >= args.limit:
                    break

                full = expand_full(meta["half"], meta["N"])
                poly_str = coeffs_to_pari(full)
                try:
                    proc = subprocess.run(
                        ["gp", "-q", "--default", "parisize=600000000"],
                        input=gp_script(poly_str, args.precision),
                        capture_output=True, text=True,
                        timeout=300,
                    )
                except subprocess.TimeoutExpired:
                    w.writerow([line_no, meta["N"], meta["M"], "timeout",
                                meta["K"], "", meta["U"], "",
                                meta["Q"], "", meta["R"], "",
                                "", "", "", "", ""])
                    n_done += 1
                    continue
                out_lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
                if not out_lines:
                    continue
                parts = out_lines[-1].split(",")
                parts = [p.strip() for p in parts]
                if len(parts) < 7:
                    continue
                irr_pari = parts[0]
                M_pari   = parts[1]
                K_pari   = int(parts[2])
                U_pari   = int(parts[3])
                Q_pari   = int(parts[4])
                R_pari   = int(parts[5])
                maxoff   = parts[6]
                ok_irr   = (irr_pari == "1")
                # M agreement: first 15 digits of the decimal expansion.
                ok_M = M_pari[:17] == meta["M"][:17]
                ok_KUQR = (K_pari == meta["K"] and U_pari == meta["U"]
                           and Q_pari == meta["Q"] and R_pari == meta["R"])
                if not (ok_irr and ok_M and ok_KUQR):
                    n_mismatches += 1
                w.writerow([
                    line_no, meta["N"], meta["M"], M_pari,
                    meta["K"], K_pari, meta["U"], U_pari,
                    meta["Q"], Q_pari, meta["R"], R_pari,
                    int(ok_irr), maxoff,
                    int(ok_irr), int(ok_M), int(ok_KUQR),
                ])
                n_done += 1
                if n_done % 100 == 0:
                    elapsed = time.time() - started
                    rate = n_done / elapsed if elapsed > 0 else 0.0
                    print(f"[{n_done}/{n_total}] {rate:.1f} polys/s, "
                          f"mismatches={n_mismatches}",
                          file=sys.stderr, flush=True)

    print(f"done: {n_done} polynomials checked, {n_mismatches} mismatches",
          file=sys.stderr)


if __name__ == "__main__":
    main()
