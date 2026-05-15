#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
verify_db_uniqueness.py — assert that no two entries in an AllKnownAdvanpix
file share the same coefficient vector.

This catches the failure mode psmm's `-merge` cannot catch on its own:
two distinct polynomials of the same degree whose Mahler measures happen to
agree to within the merge tolerance (~10^-69) get conflated. Coefficient
comparison is exact (integer-valued) so it sidesteps any precision concern.

Also reports:
  - any pair where (degree, half-coeffs) match exactly (true duplicates)
  - any pair where M-prefix-15 collides but coefficients differ (the
    dangerous case the merge would silently merge)
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"


def parse_db(path: Path):
    """Yield (line_no, N, M_str, half_tuple, raw_line) per data row."""
    with path.open() as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 10:
                continue
            try:
                N = int(toks[0])
            except ValueError:
                continue
            M_str = toks[1]
            try:
                half = tuple(int(t) for t in toks[9:])
            except ValueError:
                continue
            yield line_no, N, M_str, half, raw.rstrip("\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(DB_PATH),
                    help="DB-format file to check (default: AllKnownAdvanpix)")
    ap.add_argument("--show-collisions", type=int, default=20,
                    help="show up to N M-prefix collisions in detail")
    args = ap.parse_args()

    src = Path(args.src)
    print(f"checking {src}")
    by_coeffs = defaultdict(list)   # (N, half) -> [line_no...]
    by_m_prefix = defaultdict(list) # (N, M[:15]) -> [(line_no, half)]
    n = 0

    for line_no, N, M_str, half, _raw in parse_db(src):
        n += 1
        by_coeffs[(N, half)].append(line_no)
        by_m_prefix[(N, M_str[:15])].append((line_no, half))

    print(f"  {n} data rows parsed")

    # 1. Exact coefficient duplicates
    coeff_dupes = {k: v for k, v in by_coeffs.items() if len(v) > 1}
    print("")
    print(f"exact coefficient duplicates: {len(coeff_dupes)} group(s)")
    if coeff_dupes:
        for (N, half), lines in list(coeff_dupes.items())[:10]:
            print(f"  deg={N}, lines {lines}, "
                  f"first 10 half-coeffs: {list(half[:10])}{'...' if len(half) > 10 else ''}")

    # 2. M-prefix collisions where coefficients differ
    #    (this is the failure mode of M-tolerance-based dedup)
    dangerous = []
    for (N, prefix), entries in by_m_prefix.items():
        if len(entries) < 2:
            continue
        # group entries by coefficients
        coeff_groups = defaultdict(list)
        for line_no, half in entries:
            coeff_groups[half].append(line_no)
        if len(coeff_groups) > 1:
            dangerous.append((N, prefix, coeff_groups))

    print("")
    print(f"M-prefix-15 collisions with DIFFERENT coefficients: {len(dangerous)} group(s)")
    if dangerous:
        for N, prefix, groups in dangerous[:args.show_collisions]:
            print(f"  deg={N}, M-prefix={prefix}, "
                  f"{len(groups)} distinct coefficient vectors:")
            for half, lines in groups.items():
                print(f"    lines {lines}: half[:10]={list(half[:10])}"
                      f"{'...' if len(half) > 10 else ''}")

    # Exit non-zero if anything looks wrong.
    if coeff_dupes or dangerous:
        print("")
        print("PROBLEMS FOUND -- review above.")
        sys.exit(1)
    print("")
    print("OK -- no exact duplicates, no M-prefix collisions with differing coeffs.")


if __name__ == "__main__":
    main()
