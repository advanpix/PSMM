#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
Render the dense champion polynomials from AllKnownAdvanpix as PNG images.

For each champion we:
  1. Parse the half-coefficient line from AllKnownAdvanpix.
  2. Expand to the full polynomial via reciprocity a_k = a_{N-k}.
  3. Emit a LaTeX document with the full polynomial typeset in an `align*`
     environment broken every TERMS_PER_LINE non-zero terms.
  4. Compile with pdflatex.
  5. Crop / convert to PNG with pdftoppm.

Output goes to /home/advanpix/Development/PSMM/images/champion_*.png
"""

import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "AllKnownAdvanpix"
OUT_DIR = REPO / "images"
WORK = Path("/tmp/champion_render")

TERMS_PER_LINE = 8  # how many polynomial terms per line of LaTeX

# Each champion: (output stem, line in AllKnownAdvanpix (1-indexed), title)
CHAMPIONS = [
    ("champion_nnz212", 44151, r"NNZ\,=\,212 / L\,=\,4173 / H\,=\,23 champion"),
    ("champion_h29",    29032, r"Height H\,=\,29 champion"),
    ("champion_k76",    48028, r"K\,=\,76 / Q\,=\,152 champion"),
]


def parse_line(line):
    """Parse one AllKnownAdvanpix line into (N, M_str, fields_dict, half_coeffs)."""
    toks = line.split()
    N = int(toks[0])
    M_str = toks[1]
    fields = {
        "NNZ": int(toks[2]),
        "H":   int(toks[3]),
        "L":   int(toks[4]),
        "K":   int(toks[5]),
        "U":   int(toks[6]),
        "Q":   int(toks[7]),
        "R":   int(toks[8]),
    }
    half = [int(t) for t in toks[9:]]
    assert len(half) == N // 2 + 1, f"expected {N//2+1} coeffs, got {len(half)}"
    return N, M_str, fields, half


def expand_full(half, N):
    """Reconstruct the full coefficient vector a_0..a_N via reciprocity."""
    full = [0] * (N + 1)
    for k in range(N // 2 + 1):
        full[k] = half[k]
    for k in range(N // 2 + 1, N + 1):
        full[k] = full[N - k]
    return full


def term_latex(c, k):
    """LaTeX for a single signed term `c * x^k`, with sign separator."""
    sign = "+" if c > 0 else "-"
    mag = abs(c)
    if k == 0:
        body = f"{mag}"
    elif k == 1:
        body = "x" if mag == 1 else f"{mag}\\,x"
    else:
        body = f"x^{{{k}}}" if mag == 1 else f"{mag}\\,x^{{{k}}}"
    return sign, body


def polynomial_to_align(full):
    """Render the full polynomial as an align* body string.

    Each line holds TERMS_PER_LINE non-zero terms, aligned on the sign column.
    """
    # Iterate from highest degree to lowest (descending), include only non-zero.
    nz = [(k, c) for k, c in enumerate(full) if c != 0]
    nz.reverse()  # descending degree
    lines = []
    cur = []
    first = True
    for k, c in nz:
        sign, body = term_latex(c, k)
        if first:
            # Leading term: no leading "+", inline negative sign onto body.
            if sign == "-":
                term_str = f"-\\,{body}"
            else:
                term_str = body
            cur.append((None, term_str))  # use & alignment after first column
            first = False
        else:
            cur.append((sign, body))
        if len(cur) >= TERMS_PER_LINE:
            lines.append(cur)
            cur = []
    if cur:
        lines.append(cur)
    # Build align* body
    parts = []
    for i, line in enumerate(lines):
        if i == 0:
            # First line begins with "P(x) =" then the first term (no separator)
            leading = line[0][1]
            rest = line[1:]
            row = f"P(x) &= {leading}"
            for sign, body in rest:
                row += f" \\, &\\; {sign}\\,{body}"
            # actually keep it simple: don't worry about column alignment beyond LHS
        else:
            row = "&"
            for j, (sign, body) in enumerate(line):
                if j == 0:
                    row += f"\\; {sign}\\,{body}"
                else:
                    row += f" \\,\\; {sign}\\,{body}"
        parts.append(row)
    return " \\\\\n".join(parts)


def polynomial_simple(full):
    """Simpler: produce one long inline expression broken at term boundaries
    via a manual `multline*` with explicit \\ every TERMS_PER_LINE terms.
    """
    nz = [(k, c) for k, c in enumerate(full) if c != 0]
    nz.reverse()
    tokens = []
    for i, (k, c) in enumerate(nz):
        sign, body = term_latex(c, k)
        if i == 0:
            tokens.append(f"-\\,{body}" if sign == "-" else body)
        else:
            tokens.append(f"{sign}\\,{body}")
    # Group into lines
    lines = []
    for i in range(0, len(tokens), TERMS_PER_LINE):
        lines.append(" ".join(tokens[i:i+TERMS_PER_LINE]))
    return " \\\\\n".join(lines)


def build_tex(N, M_str, fields, full, title):
    body = polynomial_simple(full)
    # Truncate M to ~25 digits for header readability
    m_short = M_str[:25] + r"\ldots"
    return rf"""\documentclass[11pt]{{article}}
\usepackage[paperwidth=10in,paperheight=20in,margin=0.6in]{{geometry}}
\usepackage{{amsmath,amssymb}}
\usepackage{{microtype}}
\pagestyle{{empty}}
\setlength{{\parindent}}{{0pt}}
\allowdisplaybreaks
\begin{{document}}
\large
\noindent\textbf{{{title}}}\par
\medskip
\noindent$\deg P = {N}$\quad
$M(P) \approx {m_short}$\quad
$\mathrm{{NNZ}} = {fields['NNZ']}$\quad
$H = {fields['H']}$\quad
$L = {fields['L']}$\par
\noindent$K = {fields['K']}$ roots outside $|z|=1$;\quad
$U = {fields['U']}$ roots on $|z|=1$ (pairs);\quad
$Q = {fields['Q']}$ Salem quadruplets;\quad
$R = {fields['R']}$ real reciprocal pairs.
\bigskip

\small
\begin{{align*}}
{body}
\end{{align*}}

\end{{document}}
"""


def compile_and_convert(stem, tex_src):
    WORK.mkdir(exist_ok=True)
    tex_file = WORK / f"{stem}.tex"
    tex_file.write_text(tex_src)
    # Compile
    r = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_file.name],
        cwd=WORK, capture_output=True, text=True,
    )
    if r.returncode != 0:
        print(f"pdflatex failed for {stem}:")
        print(r.stdout[-3000:])
        sys.exit(1)
    pdf_file = WORK / f"{stem}.pdf"
    if not pdf_file.exists():
        print(f"No PDF produced for {stem}")
        sys.exit(1)
    # Convert PDF -> PNG (single page, ~150 DPI)
    png_prefix = WORK / stem
    subprocess.run(
        ["pdftoppm", "-png", "-r", "150", str(pdf_file), str(png_prefix)],
        check=True,
    )
    # pdftoppm names file as stem-1.png (or just stem.png for single page with -singlefile)
    # use -singlefile to get exactly stem.png
    subprocess.run(
        ["pdftoppm", "-png", "-r", "150", "-singlefile", str(pdf_file), str(png_prefix)],
        check=True,
    )
    out_path = OUT_DIR / f"{stem}.png"
    final = WORK / f"{stem}.png"
    if not final.exists():
        # pdftoppm may have produced -1.png
        alt = WORK / f"{stem}-1.png"
        if alt.exists():
            final = alt
    import shutil
    shutil.move(str(final), str(out_path))
    print(f"Wrote {out_path}  ({out_path.stat().st_size//1024} KB)")


def main():
    OUT_DIR.mkdir(exist_ok=True)
    with open(SRC) as f:
        lines = f.readlines()
    for stem, lineno, title in CHAMPIONS:
        line = lines[lineno - 1].rstrip("\n").rstrip("\r")
        N, M_str, fields, half = parse_line(line)
        full = expand_full(half, N)
        tex = build_tex(N, M_str, fields, full, title)
        compile_and_convert(stem, tex)


if __name__ == "__main__":
    main()
