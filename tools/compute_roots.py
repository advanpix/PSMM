#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
compute_roots.py — precompute all roots of every polynomial in
AllKnownAdvanpix and write them to per-degree text files in roots/.

Output layout: one file per even degree, roots/deg-NNNN.txt. Each
polynomial in the degree-N section becomes a block whose header line
is the FULL corresponding line from AllKnownAdvanpix (degree, M, NNZ,
H, L, K, U, Q, R, and half-coefficients). The block body has one root
per line as two whitespace-separated columns, "re im", at the chosen
precision. Blocks are separated by a single blank line.

Example block:

    # 10 1.176280818259917506544070338474035050693415806564695259830106347029688377 4 1 9 1 8 0 2 1 1 0 -1 -1 -1
    1.176280818259917506544070338474... 0.0
    0.850138998415552677470302437834878237793... 0.0
    ...

Roots are computed via PARI/GP polroots at --precision (default 80)
decimal digits. Use --workers to compute multiple degrees in parallel.

Resumability: a degree whose file already exists is skipped. Delete the
file to force a recomputation.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"
ROOTS_DIR = REPO / "roots"


def parse_db(path: Path):
    """Yield (raw_line, N, half_coeffs) for each data row of AllKnownAdvanpix."""
    with path.open() as f:
        for raw in f:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            toks = stripped.split()
            if len(toks) < 10:
                continue
            try:
                N = int(toks[0])
                half = [int(t) for t in toks[9:]]
            except ValueError:
                continue
            if len(half) != N // 2 + 1:
                continue
            yield stripped, N, half


def expand_full(half: list[int], N: int) -> list[int]:
    full = [0] * (N + 1)
    for k in range(N // 2 + 1):
        full[k] = half[k]
    for k in range(N // 2 + 1, N + 1):
        full[k] = full[N - k]
    return full


def pari_poly_expr(coeffs: list[int]) -> str:
    pieces = []
    for k, c in enumerate(coeffs):
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
        if not pieces:
            pieces.append(("-" + body) if sign == "-" else body)
        else:
            pieces.append(f"{sign} {body}")
    return " ".join(pieces) if pieces else "0"


def gp_script_for_roots(coeffs_full: list[int], precision: int) -> str:
    poly = pari_poly_expr(coeffs_full)
    return (
        f"default(realprecision, {precision});\n"
        f"P = {poly};\n"
        "rts = polroots(P);\n"
        "for (i = 1, #rts, print(real(rts[i]), \" \", imag(rts[i])));\n"
    )


def compute_one_polynomial(args_tuple):
    """Worker function: compute roots for a single polynomial.
       Returns (raw_db_line, roots_text or None on error, error_msg)."""
    raw_line, full_coeffs, precision, timeout = args_tuple
    try:
        proc = subprocess.run(
            ["gp", "-q", "--default", "parisize=4000000000"],
            input=gp_script_for_roots(full_coeffs, precision),
            capture_output=True, text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return raw_line, None, "timeout"
    if proc.returncode != 0:
        return raw_line, None, f"gp rc={proc.returncode}: {proc.stderr[:200]}"
    out_lines = [ln.strip() for ln in proc.stdout.splitlines()
                 if ln.strip() and not ln.startswith("***")]
    return raw_line, "\n".join(out_lines), None


def compute_one_degree(N: int, entries: list, precision: int, timeout: int,
                       out_dir: Path):
    """Compute and write all polynomials at degree N to one deg-file."""
    out_path = out_dir / f"deg-{N:04d}.txt"
    blocks = []
    n_ok = 0
    n_err = 0
    for raw_line, half in entries:
        full = expand_full(half, N)
        _, roots_text, err = compute_one_polynomial(
            (raw_line, full, precision, timeout)
        )
        if roots_text is None:
            blocks.append(f"# {raw_line}\n# ERROR: {err}\n")
            n_err += 1
        else:
            blocks.append(f"# {raw_line}\n{roots_text}\n")
            n_ok += 1
    out_path.write_text("\n".join(blocks))
    return (N, len(entries), n_ok, n_err, out_path)


def _worker(payload):
    N, entries, precision, timeout, out_dir = payload
    return compute_one_degree(N, entries, precision, timeout, out_dir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(DB_PATH))
    ap.add_argument("--out-dir", default=str(ROOTS_DIR))
    ap.add_argument("--precision", type=int, default=80,
                    help="PARI realprecision in decimal digits (default 80)")
    ap.add_argument("--timeout", type=int, default=1800,
                    help="per-polynomial gp timeout in seconds (default 1800)")
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel workers, one degree per task (default 1)")
    ap.add_argument("--min-degree", type=int, default=0,
                    help="skip degrees below this")
    ap.add_argument("--max-degree", type=int, default=10000,
                    help="skip degrees above this")
    args = ap.parse_args()

    src = Path(args.src)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Group entries by degree (preserve order within each degree)
    by_deg = defaultdict(list)
    for raw, N, half in parse_db(src):
        if N < args.min_degree or N > args.max_degree:
            continue
        by_deg[N].append((raw, half))

    # Filter out degrees whose file already exists (resumability)
    tasks = []
    skipped = []
    for N in sorted(by_deg):
        out_path = out_dir / f"deg-{N:04d}.txt"
        if out_path.exists():
            skipped.append((N, len(by_deg[N])))
            continue
        tasks.append((N, by_deg[N], args.precision, args.timeout, out_dir))

    print(f"degrees to process: {len(tasks)}  (skipped {len(skipped)} "
          f"already present)", file=sys.stderr)
    print(f"workers: {args.workers}", file=sys.stderr)
    if not tasks:
        return

    started = time.time()
    n_polys_done = 0
    n_polys_total = sum(len(entries) for _, entries, *_ in tasks)

    if args.workers > 1:
        with mp.Pool(args.workers) as pool:
            for (N, ntotal, n_ok, n_err, path) in pool.imap_unordered(
                _worker, tasks
            ):
                n_polys_done += ntotal
                elapsed = time.time() - started
                pct = 100 * n_polys_done / n_polys_total if n_polys_total else 0
                print(f"  deg {N:>4}: {n_ok} ok, {n_err} err -> "
                      f"{path.name} "
                      f"[{n_polys_done}/{n_polys_total} = {pct:.1f}%, "
                      f"{elapsed/60:.1f} min]",
                      file=sys.stderr, flush=True)
    else:
        for task in tasks:
            (N, ntotal, n_ok, n_err, path) = _worker(task)
            n_polys_done += ntotal
            elapsed = time.time() - started
            pct = 100 * n_polys_done / n_polys_total if n_polys_total else 0
            print(f"  deg {N:>4}: {n_ok} ok, {n_err} err -> {path.name} "
                  f"[{n_polys_done}/{n_polys_total} = {pct:.1f}%, "
                  f"{elapsed/60:.1f} min]",
                  file=sys.stderr, flush=True)

    elapsed = time.time() - started
    print(f"\ndone: {n_polys_done} polynomials in {elapsed/60:.1f} min",
          file=sys.stderr)


if __name__ == "__main__":
    main()
