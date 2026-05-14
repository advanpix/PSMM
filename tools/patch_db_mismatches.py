#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
patch_db_mismatches.py — apply PARI-verified corrections to AllKnownAdvanpix
for the 4 entries whose K/U/Q root counts violated the reciprocity identity
2K + U = N. The Mahler measure and coefficients are correct as stored;
only K, U, Q need to be rewritten.

Computed corrections (from doc/database-verification.csv, all detected
by tools/bulk_verify.py):

    line 31272, N=360:  K  61->60, U 239->240, Q 121->120
    line 31875, N=364:  K  53->52, U 259->260, Q 105->104
    line 41843, N=420:  K  56->56, U 307->308, Q 113->112
    line 45319, N=438:  K  51->49, U 337->340, Q  99-> 96

Idempotent: running twice does nothing.
"""

from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DB = REPO / "AllKnownAdvanpix"

# line_no (1-based) -> (new_K, new_U, new_Q)
PATCHES = {
    31272: (60, 240, 120),
    31875: (52, 260, 104),
    41843: (56, 308, 112),
    45319: (49, 340,  96),
}


def main():
    lines = DB.read_bytes().splitlines(keepends=True)
    for line_no, (new_K, new_U, new_Q) in PATCHES.items():
        idx = line_no - 1
        raw = lines[idx]
        # Detect line ending
        if raw.endswith(b"\r\n"):
            eol = b"\r\n"
        elif raw.endswith(b"\n"):
            eol = b"\n"
        else:
            eol = b""
        content = raw[:-len(eol)] if eol else raw
        toks = content.decode("ascii").split()
        # Format: N M NNZ H L K U Q R c_0 ... c_{N/2}
        N = int(toks[0])
        old_K, old_U, old_Q = int(toks[5]), int(toks[6]), int(toks[7])
        if (old_K, old_U, old_Q) == (new_K, new_U, new_Q):
            print(f"  line {line_no} already patched (skip)")
            continue
        toks[5] = str(new_K)
        toks[6] = str(new_U)
        toks[7] = str(new_Q)
        # Reproduce the printp() format: "%3zu " for degree, then rest joined by spaces, trailing " "
        body = " ".join(toks[1:])
        new_line = (f"{N:3d} {body} ".encode("ascii")) + eol
        # Sanity check: 2K + U should equal N
        assert 2 * new_K + new_U == N, f"line {line_no}: 2K+U != N after patch"
        lines[idx] = new_line
        print(f"  line {line_no}, N={N}: K {old_K}->{new_K}, "
              f"U {old_U}->{new_U}, Q {old_Q}->{new_Q}")
    DB.write_bytes(b"".join(lines))
    print(f"\nWrote {DB}")


if __name__ == "__main__":
    main()
