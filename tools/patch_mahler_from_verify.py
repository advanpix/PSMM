#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
patch_mahler_from_verify.py — patch Mahler measures in AllKnownAdvanpix
from a bulk_verify.py output CSV.

When entries were emitted at a PARI precision below the database's 72
stored digits, the stored Mahler measure has trailing zeros past the
emit precision. A subsequent bulk_verify pass at higher precision
recomputes M and writes it to its CSV as M_pari. This script reads
that CSV and splices each row's M_pari back into the matching
AllKnownAdvanpix entry.

Matching rule: same degree N AND the stored M agrees with M_pari to the
first 30 significant digits. This is unique by construction (in any
realistic dataset two distinct polynomials of the same degree can't
share 30 digits of Mahler measure) and tolerates the trailing-zero pad.

Backups the file before editing. Idempotent: a row whose stored M
already equals M_pari is skipped.

Usage:
    python3 tools/patch_mahler_from_verify.py \
        --verify-csv doc/.verify_part1_inflight.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(DB_PATH),
                    help="AllKnownAdvanpix-format file to patch")
    ap.add_argument("--verify-csv", required=True,
                    help="bulk_verify.py output with M_pari column")
    ap.add_argument("--match-digits", type=int, default=30,
                    help="number of leading digits of M used to match "
                         "DB entry to verify row (default 30)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # Build a map (N, M_prefix) -> M_pari from the verify CSV
    target = {}
    with open(args.verify_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                N = int(row["N"])
                M_stored = row["M_stored"].strip()
                M_pari = row["M_pari"].strip()
            except (KeyError, ValueError):
                continue
            if not M_stored or not M_pari or M_pari == "timeout":
                continue
            key = (N, M_stored[:args.match_digits])
            target[key] = M_pari

    print(f"{args.verify_csv}: {len(target)} usable patch rows",
          file=sys.stderr)
    if not target:
        return

    src = Path(args.src)
    with src.open() as f:
        raw_lines = f.readlines()

    patched = 0
    already_ok = 0
    no_match = 0
    for i, raw in enumerate(raw_lines):
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
        M_stored = toks[1]
        # Match by first match-digits digits
        key = (N, M_stored[:args.match_digits])
        if key not in target:
            continue
        M_new = target[key]
        if M_new == M_stored:
            already_ok += 1
            continue
        # Splice M into the raw line, preserving everything else
        m = re.match(r"^(\s*\d+\s+)(\S+)(\s+.*)$", raw)
        if m is None:
            no_match += 1
            continue
        prefix, _old, suffix = m.group(1), m.group(2), m.group(3)
        raw_lines[i] = prefix + M_new + suffix
        patched += 1

    print(f"{src}: {patched} entries patched, "
          f"{already_ok} already had matching M, "
          f"{no_match} parse failures", file=sys.stderr)

    if args.dry_run:
        print("dry-run; not writing", file=sys.stderr)
        return

    if patched == 0:
        print("no changes; not writing", file=sys.stderr)
        return

    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    backup = src.with_name(f"{src.name}.bak.{ts}")
    backup.write_bytes(src.read_bytes())
    print(f"backup: {backup}", file=sys.stderr)
    src.write_text("".join(raw_lines))
    print(f"wrote {src}", file=sys.stderr)


if __name__ == "__main__":
    main()
