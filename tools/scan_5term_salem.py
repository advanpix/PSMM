#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
scan_5term_salem.py — sweep the 5-term reciprocal Salem family

    P_{N, a, s1, s2}(x) = 1 + s1 x^a + s2 x^(N/2) + s1 x^(N-a) + x^N

over (N, a) and compute M(P), classification (K, U, Q, R), and
irreducibility per polynomial. Output CSV mirrors the format of
tools/scan_pn_convergence.py.

Defaults: sweep s1 = s2 = -1 (the densest cluster of M ~ 1.286 in
AllKnownAdvanpix) over N in [10, 1000] in steps of 2, and over
a in [1, N/2 - 1].  Sub-select with --a-list to restrict to a few
slices for the Boyd-Lawton convergence picture.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

DB_MAHLER_DIGITS = 72


def round_to_db_format(M_str: str, digits: int = DB_MAHLER_DIGITS) -> str:
    """Round to exactly `digits` fractional digits, ROUND_HALF_EVEN.
    Mirrors tools/scan_pn_convergence.py."""
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


def gp_script(N: int, a: int, s1: int, s2: int, precision: int,
              emit_roots: bool = False) -> str:
    """gp script: build P, run polroots, classify, optionally print roots."""
    s1c = "+" if s1 > 0 else "-"
    s2c = "+" if s2 > 0 else "-"
    half = N // 2
    co = N - a
    poly = f"1 {s1c} x^{a} {s2c} x^{half} {s1c} x^{co} + x^{N}"
    classify = (
        "r = abs(rts[i]); "
        "if (abs(r - 1) < eps, U += 1, "
        "    if (r > 1, M *= r; K += 1); "
        "    if (abs(imag(rts[i])) < eps, R += 1, Q += 1) )"
    )
    root_emit = (
        "print(\"ROOT \", real(rts[i]), \" \", imag(rts[i])); "
        if emit_roots else ""
    )
    return (
        f"default(realprecision, {precision});\n"
        f"P = {poly};\n"
        "rts = polroots(P);\n"
        "M = 1.0; K = 0; U = 0; Q = 0; R = 0;\n"
        "eps = 1e-30;\n"
        f"for (i = 1, #rts, {root_emit}{classify});\n"
        "irr = polisirreducible(P);\n"
        "print(M, \" \", K, \" \", U, \" \", Q, \" \", R, \" \", irr);\n"
    )


def parse_combos(spec: str) -> list[tuple[int, int]]:
    """Parse 'MM,PP,MP,PM' into [(-1,-1), (+1,+1), (-1,+1), (+1,-1)].
    Each combo is a 2-letter code: M = -1, P = +1, in order (s1, s2).
    Combos are comma-separated.  We use letters so the spec is free of
    shell metacharacters (', ;, |, etc.)."""
    out = []
    code_map = {"M": -1, "P": 1}
    for token in spec.split(","):
        token = token.strip().upper()
        if not token:
            continue
        if len(token) != 2 or token[0] not in code_map or token[1] not in code_map:
            raise SystemExit(f"bad combo code {token!r}; use MM/PP/MP/PM")
        out.append((code_map[token[0]], code_map[token[1]]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1", type=int, choices=[-1, 1], default=-1,
                    help="legacy single-combo flag (ignored when --combos set)")
    ap.add_argument("--s2", type=int, choices=[-1, 1], default=-1,
                    help="legacy single-combo flag (ignored when --combos set)")
    ap.add_argument("--combos", default="",
                    help="comma-separated list of sign-combo codes; each "
                         "code is two letters from {M, P} mapping to "
                         "(s1, s2). Examples: 'MM' = (-1, -1); "
                         "'MM,PP,MP,PM' = all four. Codes used so the "
                         "spec is free of shell metacharacters. When "
                         "empty, falls back to --s1/--s2.")
    ap.add_argument("--n-min", type=int, default=10)
    ap.add_argument("--n-max", type=int, default=1000)
    ap.add_argument("--n-step", type=int, default=2)
    ap.add_argument("--a-list", default="",
                    help="comma-separated list of a values to sweep "
                         "(default: all a in [1, N/2-1])")
    ap.add_argument("--a-max", type=int, default=0,
                    help="when --a-list is empty, cap a at this value "
                         "(default 0 = all)")
    ap.add_argument("--precision", type=int, default=120)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--output", type=Path,
                    default=REPO / "work" / "scan-5term-salem.csv")
    ap.add_argument("--db-output", type=Path, default=None,
                    help="optional AllKnownAdvanpix-format dump of "
                         "irreducible entries with 1.001 < M < 1.3")
    ap.add_argument("--roots-output", type=Path, default=None)
    args = ap.parse_args()

    if args.combos:
        combos = parse_combos(args.combos)
    else:
        combos = [(args.s1, args.s2)]
    print(f"sweeping {len(combos)} sign combo(s): "
          + ", ".join(f"({s1:+d},{s2:+d})" for s1, s2 in combos),
          file=sys.stderr)

    a_fixed = [int(x) for x in args.a_list.split(",") if x.strip()]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    db_f = args.db_output.open("w", buffering=1) if args.db_output else None
    roots_f = (args.roots_output.open("w", buffering=1)
               if args.roots_output else None)

    out = args.output.open("w", buffering=1)
    w = csv.writer(out)
    w.writerow(["N", "a", "s1", "s2", "M", "K", "U", "Q", "R", "irr"])

    started = time.time()
    n_done = 0
    n_db_emitted = 0

    n_step = args.n_step if args.n_step % 2 == 0 else args.n_step + 1
    try:
        for N in range(args.n_min if args.n_min % 2 == 0 else args.n_min + 1,
                       args.n_max + 1, n_step):
            if a_fixed:
                a_list = [a for a in a_fixed if 1 <= a < N // 2]
            else:
                a_max = args.a_max if args.a_max > 0 else N // 2 - 1
                a_list = list(range(1, min(a_max + 1, N // 2)))
            for a in a_list:
                for (s1, s2) in combos:
                    emit_roots = roots_f is not None
                    try:
                        proc = subprocess.run(
                            ["gp", "-q", "--default",
                             "parisize=4000000000"],
                            input=gp_script(N, a, s1, s2,
                                            args.precision,
                                            emit_roots=emit_roots),
                            capture_output=True, text=True,
                            timeout=args.timeout,
                        )
                    except subprocess.TimeoutExpired:
                        print(f"N={N} a={a} ({s1:+d},{s2:+d}): timeout",
                              file=sys.stderr)
                        w.writerow([N, a, s1, s2,
                                    "timeout", "", "", "", "", ""])
                        continue
                    if proc.returncode != 0:
                        print(f"N={N} a={a} ({s1:+d},{s2:+d}): "
                              f"gp rc={proc.returncode}",
                              file=sys.stderr)
                        continue
                    lines = [ln.strip()
                             for ln in proc.stdout.splitlines()
                             if ln.strip() and not ln.startswith("***")]
                    root_lines = [ln for ln in lines
                                  if ln.startswith("ROOT ")]
                    summary_lines = [ln for ln in lines
                                     if not ln.startswith("ROOT ")]
                    if not summary_lines:
                        continue
                    parts = summary_lines[-1].split()
                    if len(parts) < 6:
                        continue
                    M_str = parts[0]
                    try:
                        K = int(parts[1]); U = int(parts[2])
                        Q = int(parts[3]); R = int(parts[4])
                        irr = int(parts[5])
                    except ValueError:
                        continue
                    w.writerow([N, a, s1, s2, M_str,
                                K, U, Q, R, irr])

                    try:
                        M_f = float(M_str)
                    except ValueError:
                        M_f = float("nan")
                    if (db_f and irr and 1.001 < M_f < 1.3
                            and 2 * K + U == N and Q + R == 2 * K):
                        half = [0] * (N // 2 + 1)
                        half[0] = 1
                        half[a] = s1
                        half[N // 2] = s2
                        H = 1
                        L_full = 5
                        NNZ = 2
                        M_db = round_to_db_format(M_str)
                        coefs = " ".join(str(c) for c in half)
                        db_f.write(f"{N} {M_db} {NNZ} {H} {L_full} "
                                   f"{K} {U} {Q} {R} {coefs}\n")
                        n_db_emitted += 1
                    if roots_f is not None:
                        roots_f.write(
                            f"# N={N} a={a} s1={s1} s2={s2} "
                            f"M={M_str} K={K} U={U} Q={Q} R={R}\n"
                        )
                        for rl in root_lines:
                            toks = rl.split(maxsplit=2)
                            if len(toks) == 3:
                                roots_f.write(f"{toks[1]} {toks[2]}\n")
                        roots_f.write("\n")

                    n_done += 1
                    if n_done % 100 == 0:
                        elapsed = time.time() - started
                        rate = n_done / elapsed if elapsed else 0
                        print(f"N={N} a={a} ({s1:+d},{s2:+d}) "
                              f"M={M_str[:14]} K={K} U={U}  "
                              f"({n_done} done, {rate:.1f}/s, "
                              f"{elapsed/60:.1f} min)",
                              file=sys.stderr, flush=True)
    finally:
        out.close()
        if db_f:
            db_f.close()
        if roots_f:
            roots_f.close()

    elapsed = time.time() - started
    print(f"\ndone: {n_done} polynomials in {elapsed/60:.1f} min",
          file=sys.stderr)
    print(f"  DB-mergeable entries: {n_db_emitted}", file=sys.stderr)


if __name__ == "__main__":
    main()
