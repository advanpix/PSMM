#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
scan_ad_below_14.py — sweep all (a, d, s) cyclotomic-perturbation
cells whose Boyd-Lawton limit sits below 1.5, by running
tools/scan_pn_convergence.py once per cell.

Targets the 27 cells with L < 1.5 from the (a, d) L-table that have
not yet been scanned. The (2, 3, -/+) and (2, 12, -) cells are
excluded because their data is already in AllKnownAdvanpix.

(Originally named *below_14* when the scope was L < 1.4; broadened to
L < 1.5 on 2026-05-24 to match the unified DB's M-coverage policy in
doc/SCOPE.md.)

Per the scan-completeness discipline in
~/.claude/.../feedback_scan_completeness.md the per-cell run:
  - factors every reducible P_{a,d,k,s} and emits each non-cyclotomic
    palindromic factor (scan_pn_convergence handles this);
  - logs polynomial coefficients for every emitted entry
    (--db-output);
  - logs roots for every evaluated polynomial regardless of M
    (--roots-output, NOT filtered by threshold);
  - uses configurable thresholds (--m-max-db default 1.3 for the DB
    merge gate; the convergence CSV always logs every N).

Cells are run in parallel via a multiprocessing pool. Each cell uses
one CPU (scan_pn_convergence's gp calls are single-threaded), so the
default --workers 8 keeps half the CPU free for other work. Resumable
because each cell skips itself if its output files exist (caller
controls via --force).

Outputs in work/scan_ad/:
  <a>_<d>_<sign>.csv         convergence trace (always written)
  finds_<a>_<d>_<sign>.txt   DB-format finds, gated by --m-max-db
  roots_<a>_<d>_<sign>.txt   canonical roots blocks (no threshold)
  <a>_<d>_<sign>.log         stdout+stderr from scan_pn_convergence
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
SCAN_PN = REPO / "tools" / "scan_pn_convergence.py"
OUT_DIR = REPO / "work" / "scan_ad"


# (a, d, sign, L) tuples for cells with L < 1.5 in the (a, d) L-table
# from the README's "Generalising via cyclotomic perturbation" section.
# Already-scanned cells excluded: (2, 3, -) and (2, 3, +) at L=1.2554,
# (2, 12, -) at L=1.30911.
CELLS = [
    # (a, d, sign, L_predicted)
    (3, 5, -1,   1.3157),
    (6, 10, -1,  1.3157),
    (2, 5, -1,   1.3321),
    (4, 3, -1,   1.3405),
    (4, 6, -1,   1.3405),
    (8, 12, -1,  1.3405),
    (3, 7, -1,   1.3500),
    (6, 14, -1,  1.3500),
    (5, 7, -1,   1.3602),
    (10, 14, -1, 1.3602),
    (7, 13, -1,  1.3645),
    (5, 11, -1,  1.3647),
    (7, 11, -1,  1.3689),
    (11, 13, -1, 1.3776),
    (3, 13, -1,  1.3794),
    (5, 13, -1,  1.3814),
    (2, 7, -1,   1.3883),
    # added 2026-05-24 to extend scope from L<1.4 to L<1.5
    (2, 8, -1,   1.4098),
    (4, 5, -1,   1.4227),
    (4, 10, -1,  1.4227),
    (8, 9, -1,   1.4299),
    (2, 9, -1,   1.4497),
    (3, 4, -1,   1.4751),
    (6, 4, -1,   1.4751),
    (12, 8, -1,  1.4751),
    (5, 8, -1,   1.4819),
    (10, 8, -1,  1.4819),
]


def cell_id(a: int, d: int, sign: int) -> str:
    s = "p" if sign > 0 else "m"
    return f"{a}_{d}_{s}"


def cell_outputs(a: int, d: int, sign: int, out_dir: Path):
    cid = cell_id(a, d, sign)
    return {
        "csv":   out_dir / f"{cid}.csv",
        "db":    out_dir / f"finds_{cid}.txt",
        "roots": out_dir / f"roots_{cid}.txt",
        "log":   out_dir / f"{cid}.log",
    }


def run_one_cell(cell_info):
    """Worker: invoke scan_pn_convergence.py for a single (a, d, sign).
    Returns (cell_id, status, wall_seconds)."""
    a, d, sign, L_pred, n_min, n_max, m_max_db, precision, timeout, \
        out_dir, force = cell_info
    cid = cell_id(a, d, sign)
    outs = cell_outputs(a, d, sign, out_dir)

    if not force and all(p.exists() for p in (outs["csv"], outs["db"],
                                              outs["roots"])):
        return cid, "skipped (outputs exist)", 0.0

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(SCAN_PN),
        f"--a={a}", f"--d={d}", f"--sign={sign}",
        f"--n-min={n_min}", f"--n-max={n_max}",
        f"--precision={precision}",
        f"--timeout={timeout}",
        f"--output={outs['csv']}",
        f"--db-output={outs['db']}",
        f"--roots-output={outs['roots']}",
        f"--m-min=1.001",
        f"--m-max-db={m_max_db}",
        # Predicted L is informational; pass it through for the CSV
        # convergence trace's abs_diff_from_limit column.
        f"--limit-m={L_pred:.6f}",
    ]
    start = time.time()
    with outs["log"].open("w") as logf:
        logf.write(f"# scan_ad_below_14.py: cell ({a}, {d}, {sign:+d}), L={L_pred}\n")
        logf.write(f"# command: {' '.join(cmd)}\n")
        logf.write(f"# started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        logf.flush()
        # No outer-process timeout: scans run to completion. The
        # per-gp-call timeout inside scan_pn_convergence (passed via
        # --timeout) handles individual hung polroots calls; the
        # outer subprocess is allowed to take as long as it needs.
        proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
        status = f"rc={proc.returncode}"
    wall = time.time() - start
    return cid, status, wall


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR,
                    help="output directory (default work/scan_ad/)")
    ap.add_argument("--n-min", type=int, default=6,
                    help="lower N bound passed to scan_pn_convergence (default 6)")
    ap.add_argument("--n-max", type=int, default=500,
                    help="upper N bound passed to scan_pn_convergence "
                         "(default 500; the (2,12,-) baseline scan used "
                         "n-max=1000, but Phase E only needs enough N to "
                         "establish convergence)")
    ap.add_argument("--m-max-db", type=float, default=1.5,
                    help="DB-format emission cutoff for each cell "
                         "(default 1.5 — matches the unified DB's "
                         "M-coverage policy in doc/SCOPE.md. Lower this "
                         "to restrict to sub-1.3 emission only.)")
    ap.add_argument("--precision", type=int, default=120,
                    help="PARI realprecision in decimal digits (default 120)")
    ap.add_argument("--timeout", type=int, default=1800,
                    help="per-call gp timeout in seconds (default 1800)")
    ap.add_argument("--workers", type=int, default=8,
                    help="parallel cells (default 8; each uses 1 CPU "
                         "for its gp subprocesses)")
    ap.add_argument("--force", action="store_true",
                    help="re-run cells even if their output files exist")
    ap.add_argument("--cells", default="",
                    help="comma-separated list of cell ids (e.g., "
                         "'3_5_m,5_7_m') to restrict the run; empty = all")
    args = ap.parse_args()

    selected = CELLS
    if args.cells.strip():
        wanted = {c.strip() for c in args.cells.split(",") if c.strip()}
        selected = [c for c in CELLS if cell_id(c[0], c[1], c[2]) in wanted]
        if not selected:
            raise SystemExit(f"no cells matched {args.cells!r}; available: "
                             + ", ".join(cell_id(*c[:3]) for c in CELLS))

    print(f"scan_ad_below_14: {len(selected)} cell(s), workers={args.workers}",
          file=sys.stderr)
    print(f"  out_dir:   {args.out_dir}", file=sys.stderr)
    print(f"  N range:   [{args.n_min}, {args.n_max}]", file=sys.stderr)
    print(f"  m-max-db:  {args.m_max_db}", file=sys.stderr)
    print(f"  precision: {args.precision}", file=sys.stderr)
    print(f"  timeout:   {args.timeout}s per N", file=sys.stderr)
    print(f"  force:     {args.force}", file=sys.stderr)
    print("", file=sys.stderr)
    print("Cells to run:", file=sys.stderr)
    for (a, d, sign, L) in selected:
        sl = "+" if sign > 0 else "-"
        print(f"  ({a:>2}, {d:>2}, {sl})   L = {L:.4f}   "
              f"-> {cell_id(a, d, sign)}", file=sys.stderr)
    print("", file=sys.stderr)

    tasks = [(a, d, sign, L, args.n_min, args.n_max, args.m_max_db,
              args.precision, args.timeout, args.out_dir, args.force)
             for (a, d, sign, L) in selected]

    started = time.time()
    results = []
    if args.workers > 1:
        with mp.Pool(args.workers) as pool:
            for cid, status, wall in pool.imap_unordered(run_one_cell, tasks):
                results.append((cid, status, wall))
                elapsed = time.time() - started
                print(f"  [{len(results)}/{len(tasks)}] {cid:>12} "
                      f"{status:<24} wall={wall/60:.1f} min   "
                      f"(elapsed {elapsed/60:.1f} min)",
                      file=sys.stderr, flush=True)
    else:
        for task in tasks:
            r = run_one_cell(task)
            results.append(r)
            elapsed = time.time() - started
            print(f"  [{len(results)}/{len(tasks)}] {r[0]:>12} "
                  f"{r[1]:<24} wall={r[2]/60:.1f} min   "
                  f"(elapsed {elapsed/60:.1f} min)",
                  file=sys.stderr, flush=True)

    total = time.time() - started
    print(f"\ndone in {total/60:.1f} min ({total/3600:.1f} hours).",
          file=sys.stderr)
    print(f"\nPer-cell summary:", file=sys.stderr)
    for cid, status, wall in results:
        print(f"  {cid:>12}  {status:<24}  {wall/60:.1f} min", file=sys.stderr)


if __name__ == "__main__":
    main()
