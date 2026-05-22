#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
inject_scan_roots.py — merge scan-output roots files into the canonical
roots/ per-degree store without recomputing polroots.

Two input formats supported via --format:

  db-line   Block header is already in canonical AllKnownAdvanpix line
            format (degree, M, NNZ, H, L, K, U, Q, R, half-coefficients).
            Used for work/roots-2-3-plus.txt and any scan-output file
            produced by the fixed scan_pn_convergence.py --roots-output
            emission.

  s5-params Block header is the parametric form used by scan_5term_salem:
              # N=... a=... s1=... s2=... M=... K=... U=... Q=... R=...
            Half-coefficients are reconstructed from (N, a, s1, s2) per
            the S5 family form
              P_{N,a,s1,s2}(x) = 1 + s1*x^a + s2*x^{N/2}
                                    + s1*x^{N-a} + x^N
            i.e. half[0]=1, half[a]=s1, half[N/2]=s2, rest 0.

For each scan block, the (N, half_coeffs) key is looked up in
AllKnownAdvanpix. If present, the canonical DB line is used as the
output block header (the scan-header's full-precision M is replaced
with the rounded DB M). If absent (e.g. a polynomial that did not make
it into the DB because of the M < 1.3 cutoff or because of a duplicate
under another parametrisation), the block is skipped.

Operation is idempotent: blocks whose (N, half_coeffs) already appear
in the target roots/deg-NNNN.txt are skipped. Re-running the tool is
safe and produces no duplicates.

Usage examples:
  python3 tools/inject_scan_roots.py --format=db-line \\
      work/roots-2-3-plus.txt
  python3 tools/inject_scan_roots.py --format=s5-params \\
      work/roots-5term-all-combos.txt
  python3 tools/inject_scan_roots.py --format=db-line --dry-run \\
      work/roots-2-3-plus.txt
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "AllKnownAdvanpix"
ROOTS_DIR = REPO / "roots"


def parse_db_line(line: str):
    """Return (N, half_tuple, canonical_line_stripped) or None."""
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None
    toks = stripped.split()
    if len(toks) < 10:
        return None
    try:
        N = int(toks[0])
        half = tuple(int(t) for t in toks[9:])
    except ValueError:
        return None
    if len(half) != N // 2 + 1:
        return None
    return N, half, stripped


def load_db_index(db_path: Path) -> dict:
    """Return dict (N, half_tuple) -> canonical_db_line."""
    idx = {}
    with db_path.open() as f:
        for line in f:
            parsed = parse_db_line(line)
            if parsed is None:
                continue
            N, half, raw = parsed
            idx[(N, half)] = raw
    return idx


def load_existing_keys(roots_dir: Path) -> dict:
    """Return dict N -> set of half_tuples already in roots/deg-NNNN.txt."""
    seen = defaultdict(set)
    for path in roots_dir.glob("deg-*.txt"):
        with path.open() as f:
            for line in f:
                if not line.startswith("# "):
                    continue
                parsed = parse_db_line(line[2:])
                if parsed is None:
                    continue
                N, half, _ = parsed
                seen[N].add(half)
    return seen


def iter_blocks(path: Path):
    """Yield (header_line, body_lines) for each '# ...' block."""
    header = None
    body: list = []
    with path.open() as f:
        for line in f:
            stripped_nl = line.rstrip("\n")
            if line.startswith("#"):
                if header is not None:
                    yield header, body
                header = stripped_nl
                body = []
            elif line.strip() == "":
                continue
            else:
                if header is not None:
                    body.append(stripped_nl)
        if header is not None:
            yield header, body


def parse_dbline_header(header: str):
    """Header form: '# <DB-line>'. Return (N, half_tuple) or None."""
    # Strip leading '# ' or '#'
    if header.startswith("# "):
        body = header[2:]
    elif header.startswith("#"):
        body = header[1:].lstrip()
    else:
        return None
    parsed = parse_db_line(body)
    if parsed is None:
        return None
    N, half, _ = parsed
    return N, half


def parse_s5_header(header: str):
    """Header form: '# N=.. a=.. s1=.. s2=.. ...'. Return (N, half_tuple)."""
    # Drop the leading '#'
    if not header.startswith("#"):
        return None
    body = header[1:].strip()
    kv = {}
    for tok in body.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            kv[k] = v
    try:
        N = int(kv["N"])
        a = int(kv["a"])
        s1 = int(kv["s1"])
        s2 = int(kv["s2"])
    except (KeyError, ValueError):
        return None
    if N % 2 != 0 or N < 4:
        return None
    if not (1 <= a < N // 2):
        return None
    if s1 not in (-1, 1) or s2 not in (-1, 1):
        return None
    half = [0] * (N // 2 + 1)
    half[0] = 1
    half[a] = s1
    half[N // 2] = s2
    return N, tuple(half)


def parse_header(header: str, fmt: str):
    if fmt == "db-line":
        return parse_dbline_header(header)
    if fmt == "s5-params":
        return parse_s5_header(header)
    raise ValueError(f"unknown format: {fmt}")


def append_blocks_to_file(path: Path, blocks: list, dry_run: bool):
    """Append `blocks` (list of (canonical_header_with_hash, body)) to path.

    Uses a single blank line as block separator. Idempotent dedupe is
    expected to be done by the caller before invocation.
    """
    if dry_run:
        return

    if not path.exists():
        # Sanity-shout: we shouldn't typically be CREATING a new degree
        # file via injection alone. Degrees missing from roots/ usually
        # indicate the original compute_roots.py run skipped them.
        new_chunks = []
        for header, body in blocks:
            new_chunks.append(header)
            new_chunks.extend(body)
        # Use newline join across all lines; trailing newline at EOF.
        content = "\n".join(new_chunks) + "\n"
        # Insert blank lines between blocks (already correct because
        # the join above puts the next "# ..." right after the last
        # body line — but we need a real blank line). Easier: build
        # block-by-block and join with "\n\n" then "\n" at end.
        # Recompute cleanly:
        formatted_blocks = []
        for header, body in blocks:
            block_text = header + "\n" + "\n".join(body)
            formatted_blocks.append(block_text)
        content = "\n\n".join(formatted_blocks) + "\n"
        path.write_text(content)
        return

    existing = path.read_text()
    # Ensure existing ends with exactly two newlines (blank-line
    # separator) before appending the next block, so the new block's
    # first character is the '#' of its header right after a blank line.
    if not existing.endswith("\n"):
        existing += "\n"
    if not existing.endswith("\n\n"):
        existing += "\n"

    formatted_blocks = []
    for header, body in blocks:
        block_text = header + "\n" + "\n".join(body)
        formatted_blocks.append(block_text)
    new_content = "\n\n".join(formatted_blocks) + "\n"

    path.write_text(existing + new_content)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path,
                    help="scan-output roots file to inject")
    ap.add_argument("--format", required=True,
                    choices=["db-line", "s5-params"])
    ap.add_argument("--db", default=str(DB_PATH), type=Path)
    ap.add_argument("--roots-dir", default=str(ROOTS_DIR), type=Path)
    ap.add_argument("--dry-run", action="store_true",
                    help="Compute and report counts without writing")
    args = ap.parse_args()

    print(f"[inject_scan_roots] format={args.format}", file=sys.stderr)
    print(f"[inject_scan_roots] input={args.input}", file=sys.stderr)
    print(f"[inject_scan_roots] db={args.db}", file=sys.stderr)
    print(f"[inject_scan_roots] roots_dir={args.roots_dir}", file=sys.stderr)
    if args.dry_run:
        print("[inject_scan_roots] DRY RUN — no files will be modified",
              file=sys.stderr)

    print("Loading DB index...", file=sys.stderr)
    db_idx = load_db_index(args.db)
    print(f"  {len(db_idx)} DB entries indexed.", file=sys.stderr)

    print("Scanning existing roots/ keys...", file=sys.stderr)
    existing = load_existing_keys(args.roots_dir)
    n_existing = sum(len(s) for s in existing.values())
    print(f"  {n_existing} blocks already in roots/.", file=sys.stderr)

    print("Iterating input blocks...", file=sys.stderr)
    appended_per_degree: dict = defaultdict(list)
    n_appended = 0
    n_already_present = 0
    n_not_in_db = 0
    n_parse_fail = 0
    n_empty_body = 0
    n_total_blocks = 0

    for header, body in iter_blocks(args.input):
        n_total_blocks += 1
        parsed = parse_header(header, args.format)
        if parsed is None:
            n_parse_fail += 1
            continue
        N, half = parsed
        key = (N, half)
        if key not in db_idx:
            n_not_in_db += 1
            continue
        if half in existing[N]:
            n_already_present += 1
            continue
        # Refuse to inject a block with no root rows (e.g. from a failed
        # compute_roots.py invocation that emitted only an "# ERROR: ..."
        # line). Better to leave the entry missing-from-roots/ than to
        # write a zero-root block that hides the gap.
        if not body:
            n_empty_body += 1
            continue
        canonical_header = "# " + db_idx[key]
        appended_per_degree[N].append((canonical_header, body))
        existing[N].add(half)  # avoid double-counting within input
        n_appended += 1

    print("", file=sys.stderr)
    print(f"Input blocks read:     {n_total_blocks}", file=sys.stderr)
    print(f"  appended:            {n_appended}", file=sys.stderr)
    print(f"  already in roots/:   {n_already_present}", file=sys.stderr)
    print(f"  not in DB:           {n_not_in_db}", file=sys.stderr)
    print(f"  header parse fail:   {n_parse_fail}", file=sys.stderr)
    print(f"  empty body (skip):   {n_empty_body}", file=sys.stderr)
    print(f"  degrees touched:     {len(appended_per_degree)}",
          file=sys.stderr)

    if args.dry_run:
        print("(dry-run: no files written)", file=sys.stderr)
        return

    print("Writing per-degree appends...", file=sys.stderr)
    for N in sorted(appended_per_degree):
        path = args.roots_dir / f"deg-{N:04d}.txt"
        blocks = appended_per_degree[N]
        if not path.exists():
            print(f"  WARN: creating new {path.name} (no prior file)",
                  file=sys.stderr)
        append_blocks_to_file(path, blocks, dry_run=False)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
