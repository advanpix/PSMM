#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
recompute_pari_mahler.py — recompute the Mahler measure of selected
AllKnownAdvanpix entries with PARI at high precision.

Targets entries whose stored Mahler measure string ends in a long run
of zeros: that signals the value was originally computed at low PARI
precision and then padded out to 72 digits. The actual digits beyond
the low-precision cut are unknown.

For each targeted line:
  - reconstruct the polynomial from its half-coefficients
  - compute M = product of |root| for |root|>1 via polroots at
    --precision digits (default 120)
  - splice the new high-precision M back into the line (degree, NNZ,
    H, L, K, U, Q, R, and coefficients are unchanged)

Rewrites the database file in place; writes a timestamped backup first.

Usage:
    python3 tools/recompute_pari_mahler.py
    python3 tools/recompute_pari_mahler.py --src AllKnownAdvanpix
    python3 tools/recompute_pari_mahler.py --min-trailing-zeros 25
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"


def gp_script(coeffs_full: list[int], precision: int) -> str:
    """gp script: compute M(P) where P has the given coefficient vector
       (low-to-high)."""
    poly_terms = []
    for k, c in enumerate(coeffs_full):
        if c == 0:
            continue
        sign = "+" if c > 0 else "-"
        mag = abs(c)
        if k == 0:
            body = str(mag)
        elif k == 1:
            body = "x" if mag == 1 else f"{mag}*x"
        else:
            body = f"x^{k}" if mag == 1 else f"{mag}*x^{k}"
        if not poly_terms:
            poly_terms.append(("-" + body) if sign == "-" else body)
        else:
            poly_terms.append(f"{sign} {body}")
    poly_expr = " ".join(poly_terms) if poly_terms else "0"
    return (
        f"default(realprecision, {precision});\n"
        f"P = {poly_expr};\n"
        "rts = polroots(P);\n"
        "M = 1.0;\n"
        "for (i = 1, #rts, if (abs(rts[i]) > 1, M *= abs(rts[i])));\n"
        "print(M);\n"
    )


def expand_full(half: list[int], N: int) -> list[int]:
    full = [0] * (N + 1)
    for k in range(N // 2 + 1):
        full[k] = half[k]
    for k in range(N // 2 + 1, N + 1):
        full[k] = full[N - k]
    return full


def needs_recompute(M_str: str, min_trailing_zeros: int) -> bool:
    # Strip "1." prefix if present, then count trailing zeros.
    if not M_str:
        return False
    m = re.match(r"^\d+\.(\d+)$", M_str)
    if not m:
        return False
    fractional = m.group(1)
    n_trailing = len(fractional) - len(fractional.rstrip("0"))
    return n_trailing >= min_trailing_zeros


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(DB_PATH),
                    help="path to AllKnownAdvanpix-format file")
    ap.add_argument("--precision", type=int, default=120,
                    help="PARI realprecision in decimal digits (default 120)")
    ap.add_argument("--timeout", type=int, default=1800,
                    help="per-call gp timeout in seconds (default 1800)")
    ap.add_argument("--min-trailing-zeros", type=int, default=25,
                    help="recompute entries whose Mahler measure has at "
                         "least this many trailing zeros (default 25)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change without writing")
    args = ap.parse_args()

    src = Path(args.src)
    if not src.exists():
        raise SystemExit(f"not found: {src}")

    with src.open() as f:
        raw_lines = f.readlines()

    # Identify candidate lines
    candidates = []   # list of (line_index, N, M_str, full_coeffs)
    for i, raw in enumerate(raw_lines):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        toks = line.split()
        if len(toks) < 10:
            continue
        try:
            N = int(toks[0])
            M_str = toks[1]
            half = [int(t) for t in toks[9:]]
        except ValueError:
            continue
        if len(half) != N // 2 + 1:
            continue
        if not needs_recompute(M_str, args.min_trailing_zeros):
            continue
        full = expand_full(half, N)
        candidates.append((i, N, M_str, full))

    print(f"{src}: {len(raw_lines)} lines total, "
          f"{len(candidates)} candidates for recompute "
          f"(M has >= {args.min_trailing_zeros} trailing zeros).",
          file=sys.stderr)

    if not candidates:
        print("nothing to do.", file=sys.stderr)
        return

    if args.dry_run:
        for i, N, M_str, _ in candidates[:5]:
            print(f"  line {i+1}: deg {N}, M={M_str[:40]}...", file=sys.stderr)
        if len(candidates) > 5:
            print(f"  ... ({len(candidates)} total)", file=sys.stderr)
        return

    # Backup
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    backup = src.with_name(f"{src.name}.bak.{ts}")
    backup.write_bytes(src.read_bytes())
    print(f"backup: {backup}", file=sys.stderr)

    # Recompute and patch
    n_done = 0
    n_errors = 0
    started = time.time()
    for (line_idx, N, old_M, full_coeffs) in candidates:
        try:
            proc = subprocess.run(
                ["gp", "-q", "--default", "parisize=4000000000"],
                input=gp_script(full_coeffs, args.precision),
                capture_output=True, text=True,
                timeout=args.timeout,
            )
        except subprocess.TimeoutExpired:
            print(f"  deg {N} line {line_idx+1}: TIMEOUT", file=sys.stderr)
            n_errors += 1
            continue
        if proc.returncode != 0:
            print(f"  deg {N} line {line_idx+1}: gp rc={proc.returncode}: "
                  f"{proc.stderr[:200]}", file=sys.stderr)
            n_errors += 1
            continue
        out_lines = [ln.strip() for ln in proc.stdout.splitlines()
                     if ln.strip() and not ln.startswith("***")]
        if not out_lines:
            print(f"  deg {N} line {line_idx+1}: no output", file=sys.stderr)
            n_errors += 1
            continue
        new_M = out_lines[-1]

        # Sanity: first 30 digits should match the old value
        common = min(len(old_M), len(new_M), 30)
        if old_M[:common] != new_M[:common]:
            print(f"  deg {N} line {line_idx+1}: PREFIX MISMATCH",
                  file=sys.stderr)
            print(f"    old: {old_M[:60]}", file=sys.stderr)
            print(f"    new: {new_M[:60]}", file=sys.stderr)
            n_errors += 1
            continue

        # Replace M in the raw line (preserve trailing newline, separators).
        # IMPORTANT: re.match's default `$` doesn't capture the trailing
        # '\n', so we strip it, do the replacement, then put it back.
        # Forgetting this merges patched lines with the next line.
        raw = raw_lines[line_idx]
        trailing = ""
        if raw.endswith("\r\n"):
            body, trailing = raw[:-2], "\r\n"
        elif raw.endswith("\n"):
            body, trailing = raw[:-1], "\n"
        else:
            body = raw
        m = re.match(r"^(\s*\d+\s+)(\S+)(\s.*)?$", body)
        if m is None:
            print(f"  deg {N} line {line_idx+1}: parse error",
                  file=sys.stderr)
            n_errors += 1
            continue
        prefix = m.group(1)
        suffix = m.group(3) or ""
        raw_lines[line_idx] = prefix + new_M + suffix + trailing

        n_done += 1
        if n_done % 10 == 0 or n_done == len(candidates):
            elapsed = time.time() - started
            rate = n_done / elapsed if elapsed else 0
            print(f"  recomputed {n_done}/{len(candidates)} "
                  f"({rate:.2f}/s, {elapsed/60:.1f} min)",
                  file=sys.stderr, flush=True)

    # Sanity: line count of the spliced content must equal the original.
    # An older version of the splicer ate the trailing newline of each
    # patched line and silently merged it with the next line, destroying
    # one polynomial per patch (commit 10c1f2b). Refuse to write if the
    # count doesn't match.
    out_text = "".join(raw_lines)
    out_line_count = out_text.count("\n") + (0 if out_text.endswith("\n") else 1)
    src_text = src.read_text()
    src_line_count = src_text.count("\n") + (0 if src_text.endswith("\n") else 1)
    if out_line_count != src_line_count:
        raise SystemExit(
            f"ABORT: line count would change from {src_line_count} to "
            f"{out_line_count}. Refusing to write a corrupted file. "
            f"Original retained at {backup}."
        )

    # Write back
    src.write_text(out_text)
    elapsed = time.time() - started
    print(f"\ndone: {n_done} entries patched, {n_errors} errors, "
          f"{elapsed/60:.1f} min wall.", file=sys.stderr)
    print(f"backup retained at {backup}", file=sys.stderr)


if __name__ == "__main__":
    main()
