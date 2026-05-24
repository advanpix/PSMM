#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
calibrate_psmm.py — benchmark PSMM brute-force search for a given
alphabet/NNZ grid to produce a realistic ETA table.

Runs ./build/psmm at fixed degree N for each (alphabet, NNZ) cell,
captures wall time, polynomial count, find count, and computes
polys/sec. Emits a Markdown-formatted summary at the end.

Used to decide Phase D budget before committing to large-N runs.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from math import comb
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
PSMM = REPO / "build" / "psmm"
KNOWN = REPO / "AllKnownAdvanpix"


def Q(N: int, k: int, b: int) -> int:
    """Search-space size: choose k non-zero positions from N/2, then
    assign each a value from the b-element alphabet."""
    return comb(N // 2, k) * (b ** k)


def parse_psmm_output(stdout: str):
    """Extract find count, eval count, and any summary stats from
    psmm stdout. PSMM emits '***'/'+++'/'---' lines for finds, plus
    periodic progress lines. Returns dict with counts."""
    stars = 0
    pluses = 0
    dashes = 0
    last_pps = None  # polynomials/sec from last progress line
    last_done = None
    for line in stdout.splitlines():
        if line.startswith("*** "):
            stars += 1
        elif line.startswith("+++ "):
            pluses += 1
        elif line.startswith("--- "):
            dashes += 1
        # PSMM progress lines look like:
        #   "tested 12345 of 67890 (18.2%) at 1234567 polys/s, ETA 0:00:23"
        m = re.search(r"tested\s+(\d+)\s+of\s+(\d+).*?at\s+(\d[\d,]*)\s+polys", line)
        if m:
            last_done = int(m.group(1))
            last_pps = int(m.group(3).replace(",", ""))
    return {
        "new_finds": stars,
        "session_dups": pluses,
        "known_dups": dashes,
        "progress_done": last_done,
        "progress_pps": last_pps,
    }


def run_cell(N: int, coeffs: str, nnz: int, threads: int, threshold: float,
             period: int, addto: Path | None, log_dir: Path):
    """Run one psmm cell, return (wall_seconds, parsed_stats, log_path).

    No outer-process timeout — scans run to completion. If a run needs
    to be aborted, kill the PID externally. (A previous version of this
    tool imposed a default 1h cap which silently truncated long N=40
    NNZ>=7 runs; that is the wrong policy for an exhaustive search
    tool.)"""
    coeffs_safe = coeffs.replace(",", "_").replace("-", "m")
    log_path = log_dir / f"calib_N{N:03d}_nnz{nnz}_{coeffs_safe}.log"
    cmd = [
        str(PSMM),
        f"-degree={N}",
        f"-coeffs={coeffs}",
        f"-nnz={nnz}",
        f"-threshold={threshold}",
        f"-threads={threads}",
        f"-period={period}",
        f"-known={KNOWN}",
    ]
    if addto is not None:
        cmd.append(f"-addto={addto}")
    start = time.time()
    with log_path.open("w") as logf:
        proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT,
                              text=True)
        rc = proc.returncode
    wall = time.time() - start
    stats = parse_psmm_output(log_path.read_text())
    stats["timed_out"] = False
    stats["returncode"] = rc
    return wall, stats, log_path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--N", type=int, default=20,
                    help="degree to calibrate at (default 20)")
    ap.add_argument("--nnz-list", default="2,3,4,5",
                    help="comma-separated NNZ values (default 2,3,4,5)")
    ap.add_argument("--coeffs", default="-2,-1,1,2",
                    help="comma-separated alphabet for the cell (default "
                         "'-2,-1,1,2' for H<=2)")
    ap.add_argument("--threshold", type=float, default=1.3)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--period", type=int, default=60,
                    help="psmm -period (progress report interval) in "
                         "seconds (default 60)")
    ap.add_argument("--addto", type=Path, default=None,
                    help="if set, pass -addto=FILE to psmm so new finds "
                         "are written to FILE. Default: do not save "
                         "finds (calibration mode).")
    ap.add_argument("--log-dir", type=Path,
                    default=REPO / "work" / "calibration")
    args = ap.parse_args()

    args.log_dir.mkdir(parents=True, exist_ok=True)
    nnz_list = [int(x) for x in args.nnz_list.split(",") if x.strip()]
    b = len([c for c in args.coeffs.split(",") if c.strip()])

    print(f"PSMM calibration", file=sys.stderr)
    print(f"  binary:      {PSMM}", file=sys.stderr)
    print(f"  N:           {args.N}", file=sys.stderr)
    print(f"  coeffs:      {args.coeffs}  (b={b})", file=sys.stderr)
    print(f"  nnz_list:    {nnz_list}", file=sys.stderr)
    print(f"  threshold:   {args.threshold}", file=sys.stderr)
    print(f"  threads:     {args.threads}", file=sys.stderr)
    print(f"  log_dir:     {args.log_dir}", file=sys.stderr)
    print("", file=sys.stderr)

    results = []
    for nnz in nnz_list:
        q = Q(args.N, nnz, b)
        print(f"=== N={args.N}  NNZ={nnz}  Q={q:,} ===", file=sys.stderr)
        wall, stats, log_path = run_cell(
            args.N, args.coeffs, nnz, args.threads, args.threshold,
            args.period, args.addto, args.log_dir,
        )
        rate = q / wall if wall > 0 else 0
        prog_rate = stats.get("progress_pps") or 0
        finds = stats["new_finds"]
        timed_out = stats["timed_out"]
        results.append({
            "N": args.N, "nnz": nnz, "Q": q,
            "wall_sec": wall, "rate_overall": rate,
            "rate_progress": prog_rate,
            "new_finds": finds,
            "known_dups": stats["known_dups"],
            "session_dups": stats["session_dups"],
            "timed_out": timed_out, "rc": stats["returncode"],
            "log": log_path,
        })
        marker = " (TIMEOUT)" if timed_out else ""
        print(f"  wall: {wall:.2f}s   rate: {rate:>11,.0f} polys/s   "
              f"finds: {finds}   known-dups: {stats['known_dups']}{marker}",
              file=sys.stderr)
        print("", file=sys.stderr)

    # ----- Markdown summary -----
    print("\n## Calibration summary\n")
    print(f"Degree: N = {args.N}  •  Alphabet: `{args.coeffs}`  "
          f"(b = {b})  •  Threshold: M < {args.threshold}  •  "
          f"Threads: {args.threads}\n")
    print("| NNZ | $Q(N, k, b)$ | wall time | polys/sec | new finds | "
          "known dups |")
    print("|---:|---:|---:|---:|---:|---:|")
    for r in results:
        wall_str = (f"{r['wall_sec']:.1f}s" if r['wall_sec'] < 60
                    else f"{r['wall_sec']/60:.1f} min")
        rate_str = f"{r['rate_overall']:,.0f}"
        flag = "  ⏱" if r['timed_out'] else ""
        print(f"| {r['nnz']} | {r['Q']:,} | {wall_str}{flag} | "
              f"{rate_str} | {r['new_finds']} | {r['known_dups']} |")
    print()

    # ETA projections for higher N
    print("### Extrapolation to higher N (same NNZ/alphabet, "
          "assuming polys/sec stays constant)\n")
    print("| NNZ | rate (polys/sec) | "
          + " | ".join(f"$N = {n}$" for n in (30, 40, 50)) + " |")
    print("|---:" * (2 + 3) + "|")
    for r in results:
        rate = r['rate_overall']
        if rate <= 0:
            continue
        cells = [r['nnz'], f"{rate:,.0f}"]
        for n in (30, 40, 50):
            qn = Q(n, r['nnz'], b)
            sec = qn / rate
            if sec < 60:
                cells.append(f"{sec:.1f}s")
            elif sec < 3600:
                cells.append(f"{sec/60:.1f} min")
            elif sec < 86400:
                cells.append(f"{sec/3600:.1f} hours")
            else:
                cells.append(f"{sec/86400:.1f} days")
        print("| " + " | ".join(str(c) for c in cells) + " |")


if __name__ == "__main__":
    main()
