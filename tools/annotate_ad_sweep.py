#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
annotate_ad_sweep.py — add an `in_db` column to an existing ad_sweep.csv
that was produced by an older sweep_ad.py (before the cross-check was added).

Reads the input CSV, looks up each row's (deg_smallest_factor, M_smallest_factor)
in AllKnownAdvanpix using the same 15-char Mahler-prefix key as
scan_parametric_family.py, and rewrites the CSV with `in_db` appended.

Idempotent: if the input already has an `in_db` column, it is refreshed.

Also prints a summary to stderr:
  - total rows
  - rows where smallest non-cyclotomic factor exists (M_smallest_factor < 100)
  - of those: how many are in the DB vs not
  - top non-DB rows by smallest M (candidate new finds)
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"


def load_db_keys(path: Path) -> dict:
    keys = {}
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 2:
                continue
            try:
                deg = int(toks[0])
            except ValueError:
                continue
            keys[(deg, toks[1][:15])] = line
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(REPO / "doc" / "ad_sweep.csv"))
    ap.add_argument("--output", default=None,
                    help="default: overwrite input")
    ap.add_argument("--top", type=int, default=20,
                    help="show top N candidate non-DB rows by M_smallest_factor")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path

    print(f"loading {DB_PATH} ...", file=sys.stderr)
    db_keys = load_db_keys(DB_PATH)
    print(f"  {len(db_keys)} entries indexed", file=sys.stderr)

    with in_path.open() as f:
        reader = csv.reader(f)
        header = next(reader)
        had_in_db = header and header[-1] == "in_db"
        base_header = header[:-1] if had_in_db else header
        rows = []
        for row in reader:
            row = row[:-1] if had_in_db and len(row) == len(header) else row
            rows.append(row)

    n_total = len(rows)
    n_with_factor = 0
    n_in_db = 0
    n_not_in_db = 0
    candidates = []  # (M_float, row + [in_db_str])

    annotated = []
    for row in rows:
        try:
            deg_F = int(row[7])
            M_F = row[8]
        except (ValueError, IndexError):
            annotated.append(row + ["False"])
            continue
        # 100.000... means no non-cyclotomic factor was found
        is_real_factor = not M_F.startswith("100.")
        if is_real_factor:
            n_with_factor += 1
        in_db = (deg_F, M_F[:15]) in db_keys
        if is_real_factor:
            if in_db:
                n_in_db += 1
            else:
                n_not_in_db += 1
                try:
                    candidates.append((float(M_F), row + ["False"]))
                except ValueError:
                    pass
        annotated.append(row + [str(in_db)])

    out_header = base_header + ["in_db"]
    with out_path.open("w") as f:
        w = csv.writer(f)
        w.writerow(out_header)
        w.writerows(annotated)

    print("", file=sys.stderr)
    print(f"wrote {out_path} ({len(annotated)} rows)", file=sys.stderr)
    print("", file=sys.stderr)
    print(f"total rows:                       {n_total}", file=sys.stderr)
    print(f"rows with non-cyclotomic factor:  {n_with_factor}", file=sys.stderr)
    print(f"  of those, already in DB:        {n_in_db}", file=sys.stderr)
    print(f"  of those, NOT in DB (novel):    {n_not_in_db}", file=sys.stderr)

    if candidates:
        candidates.sort(key=lambda t: t[0])
        print("", file=sys.stderr)
        print(f"top {min(args.top, len(candidates))} non-DB rows by M_smallest_factor:",
              file=sys.stderr)
        print(f"{'a':>3} {'d':>3} {'sign':>5} {'m':>5} "
              f"{'deg_F':>6} {'M_smallest_factor':>26}",
              file=sys.stderr)
        for _, row in candidates[:args.top]:
            print(f"{row[0]:>3} {row[1]:>3} {row[2]:>5} {row[3]:>5} "
                  f"{row[7]:>6} {row[8]:>26}", file=sys.stderr)


if __name__ == "__main__":
    main()
