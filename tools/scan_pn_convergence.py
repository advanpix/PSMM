#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
scan_pn_convergence.py — convergence study and DB extraction for the
cyclotomic-perturbation families

    P_{a,d,sign,N}(x) = Phi_a(x) * (x^m + 1)
                      + sign * x^k * Phi_d(x),
    where N = phi(a) + m  (the polynomial degree, even),
          k = (phi(a) + m - phi(d)) / 2.

By default (a=2, d=3, sign=-1) this is the Max-U family P_N studied first.

For each even N in a configurable range:
  - Compute M(P_{a,d,sign,N}) via PARI polroots (no factoring; factoring
    is the bottleneck past N ~ 1000).
  - Classify roots into K (outside |z|=1), U (on the circle within eps),
    Q (off-circle complex), R (off-circle real).
  - Sanity-check 2K + U = N and Q + R = 2K.

Two outputs are written incrementally:

  --output  (convergence CSV):
      N, M_PN, abs_diff_from_limit, gap_from_previous, K, U, Q, R
  (for the convergence-rate study)

  --db-output  (AllKnownAdvanpix-format entries, one line per N):
      N M_PN NNZ H L K U Q R c_0 c_1 ... c_{N/2}
  (ready for `psmm -merge`)  Emitted only if 1.001 < M < 1.3 AND the
  2K+U=N / Q+R=2K sanity passes.

Half-coefficients are constructed analytically from the (a, d, sign)
polynomial form (no PARI for the coefficients).

For all sufficiently large N in a given (a, d, sign) family, the
polynomial P_{a,d,sign,N} is empirically irreducible (a separate
verification pass via tools/bulk_verify.py is the safety net).
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

# AllKnownAdvanpix stores M as 1. + exactly 72 fractional digits, rounded
# ROUND_HALF_EVEN. See doc/DB-FORMAT.md. Every M written into a DB-format
# line by this tool must pass through round_to_db_format() first; we round
# at the source so the merge pipeline is a single step with no follow-up
# precision-patch.
DB_MAHLER_DIGITS = 72


def round_to_db_format(M_str: str, digits: int = DB_MAHLER_DIGITS) -> str:
    """Round a high-precision M string to `digits` fractional digits with
    ROUND_HALF_EVEN, matching AllKnownAdvanpix format. Pads with trailing
    zeros if the input is shorter than `digits` (the caller is responsible
    for using a precision high enough that this padding doesn't happen)."""
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

# Corrected Boyd-Lawton limit for the Max-U family P_N (a=2,d=3,sign=-1)
# (see tools/compute_boyd_lawton.py for the 1D Jensen-reduction proof).
# For other (a, d) pairs, use compute_boyd_lawton_family.py to get the
# right limit; the value below is just the default for the convergence-CSV
# diff column. Surveys: only (a=2, d=3) has L < 1.3.
BOYD_LAWTON_LIMIT_DEFAULT = "1.2554340377272518"


# Cyclotomic Phi_n: coefficient vector low-degree-first
# (Phi_2 = x+1 -> [1, 1] meaning 1 + 1*x).
CYCLOTOMICS = {
    2:  [1, 1],
    3:  [1, 1, 1],
    4:  [1, 0, 1],
    5:  [1, 1, 1, 1, 1],
    6:  [1, -1, 1],
    7:  [1, 1, 1, 1, 1, 1, 1],
    8:  [1, 0, 0, 0, 1],
    9:  [1, 0, 0, 1, 0, 0, 1],
    10: [1, -1, 1, -1, 1],
    12: [1, 0, -1, 0, 1],
}


def phi_n(n: int) -> int:
    return len(CYCLOTOMICS[n]) - 1


def pari_poly_expr(coeffs: list[int], var: str = "x") -> str:
    """Build a PARI polynomial expression from a coefficient vector
       coeffs[i] = coefficient of x^i. Output like '1 + x + x^2'."""
    terms = []
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        sign = "+" if c > 0 else "-"
        mag = abs(c)
        if i == 0:
            body = str(mag)
        elif i == 1:
            body = var if mag == 1 else f"{mag}*{var}"
        else:
            body = f"{var}^{i}" if mag == 1 else f"{mag}*{var}^{i}"
        if not terms:
            terms.append(("-" + body) if sign == "-" else body)
        else:
            terms.append(f"{sign} {body}")
    return " ".join(terms) if terms else "0"


def construct_full_coeffs(a: int, d: int, sign: int, N: int) -> list[int]:
    """Return the full coefficient list (length N+1) of P_{a,d,sign,N}.
       P = Phi_a(x) * (x^m + 1) + sign * x^k * Phi_d(x),
       with m = N - phi(a) and k = (N - phi(d)) / 2.
       Caller is responsible for verifying that k is a non-negative integer
       and that the polynomial is palindromic."""
    pa = phi_n(a)
    pd = phi_n(d)
    m = N - pa
    k_num = N - pd
    assert k_num >= 0 and k_num % 2 == 0, \
        f"k = (N-phi(d))/2 not a non-negative integer: N={N}, d={d}"
    k = k_num // 2
    coeffs = [0] * (N + 1)
    phi_a = CYCLOTOMICS[a]
    phi_d = CYCLOTOMICS[d]
    # Phi_a(x) * (x^m + 1) -> add phi_a at positions 0..pa AND at positions m..m+pa
    for i, c in enumerate(phi_a):
        coeffs[i] += c
        coeffs[i + m] += c
    # sign * x^k * Phi_d(x)
    for i, c in enumerate(phi_d):
        coeffs[i + k] += sign * c
    return coeffs


_CLASSIFY_BODY = (
    "r = abs(rts[j]); "
    "if (abs(r - 1) < eps, U += 1, "
    "    if (r > 1, M *= r; K += 1); "
    "    if (abs(imag(rts[j])) < eps, R += 1, Q += 1) )"
)


def _build_p_expr(a: int, d: int, sign: int, N: int) -> str:
    """Return the PARI expression for P_{a,d,sign,N}."""
    pa = phi_n(a)
    pd = phi_n(d)
    m = N - pa
    k = (N - pd) // 2
    s_str = "+" if sign > 0 else "-"
    phi_a_expr = pari_poly_expr(CYCLOTOMICS[a])
    phi_d_expr = pari_poly_expr(CYCLOTOMICS[d])
    return (
        f"({phi_a_expr}) * (x^{m} + 1) {s_str} x^{k} * ({phi_d_expr})"
    )


def gp_script(a: int, d: int, sign: int, N: int, precision: int,
              emit_roots: bool = False) -> str:
    """Cheap pass: compute M(P), K, U, Q, R, irreducibility (no factoring).
    With emit_roots=True, also prints one "ROOT re im" line per root of
    P inside the same polroots call, so a downstream --roots-output path
    can collect them without a second polroots invocation."""
    p_expr = _build_p_expr(a, d, sign, N)
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
        f"P = {p_expr};\n"
        "rts = polroots(P);\n"
        "M = 1.0; K = 0; U = 0; Q = 0; R = 0;\n"
        "eps = 1e-30;\n"
        f"for (i = 1, #rts, {root_emit}{classify});\n"
        "irr = polisirreducible(P);\n"
        "print(M, \" \", K, \" \", U, \" \", Q, \" \", R, \" \", irr);\n"
    )


def gp_script_factor(a: int, d: int, sign: int, N: int, precision: int) -> str:
    """Reducible-case pass: factor P over Z, emit each non-cyclotomic factor
       F with 1.001 < M(F) < 1.3 as an AllKnownAdvanpix-format line.
       Mirrors PSMM's brute-force behaviour of factoring candidates and
       collecting all sub-threshold irreducible Salem-type factors."""
    p_expr = _build_p_expr(a, d, sign, N)
    return (
        f"default(realprecision, {precision});\n"
        f"P = {p_expr};\n"
        "fac = factor(P); nf = #fac~;\n"
        "for (i = 1, nf, "
        "  F = fac[i, 1]; deg_F = poldegree(F); "
        "  if (deg_F < 2, next); "
        "  rts = polroots(F); "
        "  M = 1.0; K = 0; U = 0; Q = 0; R = 0; eps = 1e-30; "
        f"  for (j = 1, #rts, {_CLASSIFY_BODY}); "
        "  if (M <= 1.001 || M >= 1.3, next); "
        "  if (2*K + U != deg_F || Q + R != 2*K, "
        "      print(\"SKIP_KUQR \", deg_F, \" \", K, \" \", U, \" \", Q, \" \", R); next); "
        "  if (polcoef(F, deg_F) != polcoef(F, 0), "
        "      print(\"SKIP_NONPALI \", deg_F); next); "
        "  half = vector(deg_F\\2 + 1, k, polcoef(F, deg_F - (k-1))); "
        "  NNZ = 0; for (k = 2, #half, if (half[k] != 0, NNZ += 1)); "
        "  Hh = 0; for (k = 1, #half, if (abs(half[k]) > Hh, Hh = abs(half[k]))); "
        "  L = 0; for (k = 0, deg_F, L += abs(polcoef(F, k))); "
        "  out = Str(deg_F, \" \", M, \" \", NNZ, \" \", Hh, \" \", L, "
        "            \" \", K, \" \", U, \" \", Q, \" \", R); "
        "  for (k = 1, #half, out = Str(out, \" \", half[k])); "
        "  print(\"FACTOR \", out));\n"
        "print(\"END\");\n"
    )


def db_line_for(a: int, d: int, sign: int, N: int, M_str: str,
                K: int, U: int, Q: int, R: int):
    """Build an AllKnownAdvanpix-format line. Returns the line string or
       None if the polynomial isn't palindromic (so doesn't belong in
       AllKnownAdvanpix as a single half-coefficient block)."""
    full = construct_full_coeffs(a, d, sign, N)
    # Palindromic check: c_i == c_{N-i}
    for i in range(N // 2 + 1):
        if full[i] != full[N - i]:
            return None
    half = full[N // 2:][::-1]   # c_{N/2}, c_{N/2+1}, ..., c_N  reversed
    # PSMM convention: half[0] = leading = c_N, half[N/2] = middle = c_{N/2}
    half = [full[N - i] for i in range(N // 2 + 1)]
    NNZ = sum(1 for c in half[1:] if c != 0)
    H = max(abs(c) for c in half)
    L = sum(abs(c) for c in full)
    coeffs = " ".join(str(c) for c in half)
    return f"{N} {M_str} {NNZ} {H} {L} {K} {U} {Q} {R} {coeffs}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=int, default=2,
                    help="cyclotomic background index Phi_a (default 2)")
    ap.add_argument("--d", type=int, default=3,
                    help="cyclotomic perturbation index Phi_d (default 3)")
    ap.add_argument("--sign", type=int, choices=[-1, 1], default=-1,
                    help="sign of the perturbation (default -1)")
    ap.add_argument("--n-min", type=int, default=6,
                    help="minimum even N (default 6)")
    ap.add_argument("--n-max", type=int, default=2002,
                    help="maximum N (default 2002)")
    ap.add_argument("--n-step", type=int, default=2,
                    help="step in N (default 2 — every even value)")
    ap.add_argument("--precision", type=int, default=120,
                    help="PARI realprecision in decimal digits (default 120 "
                         "to fully populate the DB's 72-digit Mahler measure "
                         "with safety margin)")
    ap.add_argument("--timeout", type=int, default=1800,
                    help="per-call gp timeout in seconds (default 1800)")
    ap.add_argument("--limit-m", type=str, default=BOYD_LAWTON_LIMIT_DEFAULT,
                    help="Boyd-Lawton limit for the convergence CSV "
                         "abs_diff column (default matches a=2,d=3)")
    ap.add_argument("--output",
                    default=str(REPO / "doc" / "parametric-family-convergence.csv"),
                    help="convergence-study CSV path")
    ap.add_argument("--db-output",
                    default=str(REPO / "doc" / "new_finds_d_only_pn_extended.txt"),
                    help="AllKnownAdvanpix-format entries path "
                         "(set to empty string to disable DB output)")
    ap.add_argument("--roots-output",
                    default="",
                    help="path to save root coordinates per N in "
                         "compute_roots block format (header is the "
                         "DB-style line for P_N, body is one 'real imag' "
                         "per root). NOT filtered by M<1.3 -- all N values "
                         "are emitted. Default: empty = no roots output.")
    args = ap.parse_args()
    a, d, sign = args.a, args.d, args.sign
    if phi_n(d) < phi_n(a):
        raise SystemExit(f"need phi(d) >= phi(a); got phi({d})={phi_n(d)} "
                         f"< phi({a})={phi_n(a)}")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    db_path = Path(args.db_output) if args.db_output else None
    if db_path:
        db_path.parent.mkdir(parents=True, exist_ok=True)
    roots_path = Path(args.roots_output) if args.roots_output else None
    if roots_path:
        roots_path.parent.mkdir(parents=True, exist_ok=True)
    emit_roots = roots_path is not None

    # N parity is determined by phi(a) + phi(d): m = N - phi(a) must give
    # integer k = (phi(a) + m - phi(d))/2 = (N - phi(d))/2. So (N - phi(d))
    # must be even, i.e. N must have the same parity as phi(d).
    n_parity = phi_n(d) % 2
    n_min_floor = max(args.n_min, phi_n(a) + 1)  # need m >= 1
    n_min = n_min_floor
    while n_min % 2 != n_parity:
        n_min += 1
    n_max = args.n_max
    step = args.n_step
    if step % 2 != 0:
        step += 1  # keep parity

    started = time.time()
    n_done = 0
    n_db_emitted = 0
    n_sanity_failed = 0
    last_M = None  # for gap_from_previous

    db_f = db_path.open("w", buffering=1) if db_path else None
    roots_f = roots_path.open("w", buffering=1) if roots_path else None
    try:
        with out_path.open("w", buffering=1) as out:
            w = csv.writer(out)
            w.writerow(["N", "M_PN", "abs_diff_from_limit",
                        "gap_from_previous", "K", "U", "Q", "R"])

            for N in range(n_min, n_max + 1, step):
                try:
                    proc = subprocess.run(
                        ["gp", "-q", "--default", "parisize=4000000000"],
                        input=gp_script(a, d, sign, N, args.precision,
                                        emit_roots=emit_roots),
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
                # Separate ROOT lines (when emit_roots is on) from the
                # single summary line. Summary is "M K U Q R irr".
                root_lines = [ln for ln in lines if ln.startswith("ROOT ")]
                summary_lines = [ln for ln in lines
                                 if not ln.startswith("ROOT ")]
                if not summary_lines:
                    print(f"N={N}: no summary line", file=sys.stderr,
                          flush=True)
                    last_M = None
                    continue
                parts = summary_lines[-1].split()
                if len(parts) < 6:
                    print(f"N={N}: malformed gp output: "
                          f"{summary_lines[-1][:80]}",
                          file=sys.stderr, flush=True)
                    last_M = None
                    continue
                M_str = parts[0]
                try:
                    K = int(parts[1])
                    U = int(parts[2])
                    Q = int(parts[3])
                    R = int(parts[4])
                    irr = int(parts[5])
                except ValueError:
                    print(f"N={N}: KUQR/irr parse error", file=sys.stderr, flush=True)
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
                    limit_f = float(args.limit_m)
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

                # Roots emit (unconditional w.r.t. M cutoff). Block format
                # matches roots/deg-NNNN.txt so plot_root_heatmap can read
                # it directly. Header is the DB-style line for P_N (with
                # M rounded to 72-digit DB format for consistency) when
                # palindromic, else a simpler "# N=... M=..." string.
                if roots_f is not None:
                    M_db = round_to_db_format(M_str)
                    header = db_line_for(a, d, sign, N, M_db, K, U, Q, R)
                    if header is None:
                        header = (f"N={N} M={M_db} K={K} U={U} Q={Q} R={R} "
                                  "(non-palindromic)")
                    roots_f.write(f"# {header}\n")
                    for rl in root_lines:
                        toks = rl.split(maxsplit=2)
                        if len(toks) == 3:
                            roots_f.write(f"{toks[1]} {toks[2]}\n")
                    roots_f.write("\n")

                # DB-emit logic — mirrors PSMM brute-force search:
                #   - if irreducible AND 1.001 < M < 1.3: emit P directly
                #   - if reducible AND M < 1.3: factor P, extract each
                #     non-cyclotomic factor F with 1.001 < M(F) < 1.3
                if db_f and ok and 1.001 < M_f < 1.3:
                    # Round M to the DB's 72-digit format AT THE EMIT STEP,
                    # so the merge pipeline doesn't need a precision-patch
                    # follow-up.
                    M_db = round_to_db_format(M_str)
                    if irr:
                        db_line = db_line_for(a, d, sign, N, M_db,
                                              K, U, Q, R)
                        if db_line is None:
                            print(f"N={N}: non-palindromic, skipping DB emit",
                                  file=sys.stderr, flush=True)
                        else:
                            db_f.write(db_line + "\n")
                            n_db_emitted += 1
                    else:
                        # Reducible: factor and emit each non-cyclo factor.
                        try:
                            proc2 = subprocess.run(
                                ["gp", "-q", "--default",
                                 "parisize=4000000000"],
                                input=gp_script_factor(a, d, sign, N,
                                                       args.precision),
                                capture_output=True, text=True,
                                timeout=args.timeout,
                            )
                        except subprocess.TimeoutExpired:
                            print(f"N={N}: factor timeout",
                                  file=sys.stderr, flush=True)
                            proc2 = None
                        if proc2 and proc2.returncode == 0:
                            n_factors_emitted = 0
                            for ln in proc2.stdout.splitlines():
                                ln = ln.strip()
                                if ln.startswith("FACTOR "):
                                    # gp built a DB-format line whose M
                                    # field is at PARI realprecision. Round
                                    # it to 72 digits before emitting.
                                    body = ln[len("FACTOR "):]
                                    toks = body.split(maxsplit=2)
                                    if len(toks) >= 3:
                                        toks[1] = round_to_db_format(toks[1])
                                        body = " ".join(toks)
                                    db_f.write(body + "\n")
                                    n_factors_emitted += 1
                                    n_db_emitted += 1
                                elif ln.startswith("SKIP_"):
                                    print(f"N={N}: factor skipped: {ln}",
                                          file=sys.stderr, flush=True)
                            print(f"N={N}: reducible -> "
                                  f"emitted {n_factors_emitted} factor(s)",
                                  file=sys.stderr, flush=True)

                n_done += 1
                elapsed = time.time() - started
                rate = n_done / elapsed if elapsed > 0 else 0
                print(f"(a={a},d={d},s={sign:+d}) N={N:6d} "
                      f"M={M_str[:14]} K={K:>3d} U={U:>4d} "
                      f"Q={Q:>3d} R={R:>2d}  |M-L|={abs_diff:.3e}  "
                      f"({n_done} done, {rate:.2f}/s, "
                      f"elapsed {elapsed/60:.1f} min)",
                      file=sys.stderr, flush=True)
    finally:
        if db_f:
            db_f.close()
        if roots_f:
            roots_f.close()

    elapsed = time.time() - started
    print(f"\ndone: {n_done} N values processed in {elapsed/60:.1f} min",
          file=sys.stderr)
    print(f"  DB entries emitted: {n_db_emitted}", file=sys.stderr)
    print(f"  sanity check failures: {n_sanity_failed}", file=sys.stderr)


if __name__ == "__main__":
    main()
