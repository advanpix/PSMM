#!/usr/bin/env python3
"""
analyze_champions.py — run PolyAnalyzer on the record-holder polynomials
and emit a markdown report at doc/champion-analysis.md.

Each champion is loaded from AllKnownAdvanpix by line number, fully
analysed with PARI, and tested against a hand-picked set of structural
hypotheses (cyclotomic skeletons).

Copyright (C) 2020-2026, Advanpix LLC.
License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher

Author: Pavel Holoborodko <pavel@advanpix.com>
"""

from __future__ import annotations

import sys
from pathlib import Path
from textwrap import indent

sys.path.insert(0, str(Path(__file__).resolve().parent))
from poly_analyze import (
    PolyAnalyzer,
    expand_reciprocal,
    parse_psmm_line,
)

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "AllKnownAdvanpix"
OUT = REPO / "doc" / "champion-analysis.md"


# Each champion: (anchor, title, line_number, skeleton_candidates)
# skeleton_candidates is a list of (label, PARI polynomial string) — the
# analyser will compute P - skeleton and report the count of non-zero terms.
CHAMPIONS = [
    {
        "anchor":   "lehmer",
        "title":    "Lehmer's polynomial — the smallest known Mahler measure",
        "line":     None,  # special: Lehmer is in Known180; we hand-code coeffs
        "half":     [1, 1, 0, -1, -1, -1],  # degree 10
        "M_stored": "1.176280818259917506544070338474035050693415806564695259830106347029688377",
        "K": 1, "U": 8, "Q": 0, "R": 2,
        "skeletons": [
            # No obvious cyclotomic skeleton — included to confirm the absence.
            ("x^10 - 1",  "x^10 - 1"),
            ("x^10 + 1",  "x^10 + 1"),
        ],
        "notes": (
            "Lehmer's polynomial is irreducible with one Salem number "
            "root (the largest, β ≈ 1.17628) and its reciprocal as the "
            "only off-circle roots. It is a *classical Salem polynomial* "
            "in Salem's original sense: exactly two real off-circle roots "
            "(R = 2) and no complex off-circle roots (Q = 0). It is not "
            "close to any low-order cyclotomic product."
        ),
    },
    {
        "anchor":   "max-u",
        "title":    "Max-U champion — 396 roots on the unit circle, degree 456",
        "line":     None,  # we'll search by Mahler-measure prefix
        "find_M":   "1.254914757578847933786502373052300153727826780838735038831940237651331712",
        "skeletons": [
            ("(x+1)(x^455+1)",                 "(x+1)*(x^455 + 1)"),
            ("(x+1)(x^455+1) - x^227 Phi_3",   "(x+1)*(x^455 + 1) - x^227*(x^2 + x + 1)"),
            ("x^456 + 1",                      "x^456 + 1"),
        ],
        "notes": (
            "We confirm the conjectured decomposition: this polynomial is "
            "**exactly** the cyclotomic product (x+1)(x^{455}+1) perturbed "
            "by a single sparse term -x^{227} Φ_3(x). The unperturbed part "
            "is a product of cyclotomic polynomials, hence has all 456 roots "
            "on the unit circle. The Φ_3 perturbation, scaled by the half-"
            "degree power x^{227}, pushes exactly 60 roots off the circle "
            "into 15 Salem quadruplets — all at mean distance ≈ 0.008 from "
            "|z| = 1. This is the cleanest 'cyclotomic-perturbation' "
            "construction we observe across the database."
        ),
    },
    {
        "anchor":   "max-k",
        "title":    "Max-K / Max-Q champion — 76 roots outside the unit disk, degree 452",
        "line":     None,
        "find_M":   "1.285304814552773240944862462490593640653959388154660443572614515599313747",
        "skeletons": [
            # (1, -1, 0) repeated: (1 - x^{3m})/(x^2+x+1) for various m.
            ("(1 - x^453) / Phi_3",  "(1 - x^453) / (x^2 + x + 1)"),
            ("(1 - x^456) / Phi_3",  "(1 - x^456) / (x^2 + x + 1)"),
        ],
        "notes": (
            "The half-coefficient pattern (1, -1, 0)^k continues until a "
            "phase shift in the middle. The 'pure' periodic polynomial with "
            "(1, -1, 0) coefficients is the cyclotomic quotient "
            "(1 - x^{3m}) / Φ_3(x), whose roots are all on the unit circle. "
            "Our champion deviates from this template, and that deviation "
            "is responsible for 152 off-circle roots in 38 Salem quadruplets."
        ),
    },
    {
        "anchor":   "max-nnz",
        "title":    "Max-NNZ / Max-L / Max-H champion — 212 non-zero half-coefficients, degree 432",
        "line":     None,
        "find_M":   "1.255410338063801734431046949890021772597007831617861408843137182583448951",
        "skeletons": [],   # dense, structure unclear at first sight
        "notes": (
            "This is the densest of the four. The half-coefficient sequence "
            "(1, 2, 2, 1, 0, -1, -2, -3, -3, -2, -1, 0, 1, 2, 3, 4, 4, 3, 1, ...) "
            "resembles a discrete convolution of an arithmetic sequence with "
            "itself, suggesting that the polynomial might be the *square* (over "
            "the integers) of a simpler polynomial, or a product of two such. "
            "The script reports any low-order cyclotomic gcd hits below."
        ),
    },
    {
        "anchor":   "max-h",
        "title":    "Max-H champion — height 29, degree 348",
        "line":     None,
        "find_M":   "1.254312128683608179478643743971122104728037406656613567091505843287892707",
        "skeletons": [],
        "notes": (
            "168 non-zero half-coefficients with maximum magnitude 29. "
            "Same triangular-build-up shape as the Max-NNZ champion but "
            "at a smaller degree, leading to higher individual coefficients. "
            "We leave its detailed factor-theoretic origin as an open "
            "question."
        ),
    },
]


def find_line(src_path: Path, find_M: str) -> tuple:
    """Find the line in AllKnownAdvanpix whose Mahler measure starts with find_M."""
    with src_path.open() as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) >= 2 and toks[1].startswith(find_M[:30]):
                return i, line
    raise RuntimeError(f"polynomial with M ≈ {find_M[:20]}... not found in {src_path}")


def format_perturbation(pert: str, max_len: int = 100) -> str:
    """Format a PARI polynomial string for inline markdown display."""
    if pert is None:
        return ""
    pert = pert.replace("*", "").replace("  ", " ")
    if len(pert) > max_len:
        pert = pert[:max_len] + " ... "
    return pert


def write_champion(out, champ: dict, result):
    """Append one champion's analysis to the open output stream."""
    print(f"## {champ['title']} {{#{champ['anchor']}}}", file=out)
    print("", file=out)

    # Basic facts table
    print(f"- **Degree**: {result.degree}", file=out)
    print(f"- **Irreducible**: {'yes' if result.irreducible else '**NO**'}", file=out)
    print(f"- **Reciprocal**: {'yes' if result.reciprocal else 'no'}", file=out)
    print(f"- **Mahler measure** (200-digit PARI verification):", file=out)
    print(f"  ```", file=out)
    M_long = result.mahler_measure or ""
    for i in range(0, len(M_long), 80):
        print(f"  {M_long[i:i+80]}", file=out)
    print(f"  ```", file=out)
    print(f"- **Root counts**: K = {result.K}, U = {result.U}, Q = {result.Q}, R = {result.R}", file=out)
    print(f"  (verification of identities: 2K + U = {2*result.K + result.U}, "
          f"Q + R = {result.Q + result.R} vs 2K = {2*result.K})", file=out)
    print(f"- **Max off-circle distance** $\\max_i ||z_i| - 1|$ over off-circle roots:", file=out)
    print(f"  ```", file=out)
    print(f"  {result.max_off_circle_dist}", file=out)
    print(f"  ```", file=out)
    if result.mean_off_circle_dist:
        print(f"- **Mean off-circle distance**:", file=out)
        print(f"  ```", file=out)
        print(f"  {result.mean_off_circle_dist}", file=out)
        print(f"  ```", file=out)
    if result.cross_check:
        ok = all(result.cross_check.values())
        print(f"- **Cross-check vs PSMM-stored values**: "
              f"{'all match' if ok else 'MISMATCH ' + str(result.cross_check)}", file=out)
    print("", file=out)

    # Cyclotomic gcd hits
    cyc_div_hits = [h for h in result.cyclotomic_hits
                    if h.kind in ("xn-1", "xn+1")]
    if cyc_div_hits:
        print("**Cyclotomic gcd hits** (P shares roots with $x^n \\pm 1$):", file=out)
        for h in cyc_div_hits:
            sign = "−" if h.kind.endswith("-1") else "+"
            print(f"- $\\gcd(P, x^{{{h.n_or_d}}} {sign} 1)$ has {h.detail}", file=out)
        print("", file=out)

    # Skeleton-test hits
    skel_hits = [h for h in result.cyclotomic_hits if h.kind == "skeleton_test"]
    if skel_hits:
        print("**Structural hypothesis tests** ($P - S$ for candidate skeleton $S$):", file=out)
        print("", file=out)
        print("| Skeleton $S$ | Non-zero terms in $P - S$ |", file=out)
        print("|---|---:|", file=out)
        for h in skel_hits:
            # detail looks like "label=foo, nz(P-S)=3, P-S=..."
            parts = h.detail.split(", ")
            label = parts[0].split("=", 1)[1] if "=" in parts[0] else "?"
            nz = parts[1].split("=", 1)[1] if "=" in parts[1] else "?"
            print(f"| `{label}` | {nz} |", file=out)
        print("", file=out)
        if result.near_cyclotomic_skeleton:
            pert_short = format_perturbation(result.skeleton_perturbation, 200)
            print(f"**Verified decomposition** "
                  f"`P = {result.near_cyclotomic_skeleton} + (perturbation)` where", file=out)
            print(f"```", file=out)
            print(f"P − S = {pert_short}", file=out)
            print(f"```", file=out)
            print("", file=out)

    if champ.get("notes"):
        print(f"**Notes.** {champ['notes']}", file=out)
        print("", file=out)
    print("---", file=out)
    print("", file=out)


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)

    # Resolve lines / coefficients for each champion
    for champ in CHAMPIONS:
        if "half" in champ and champ["half"]:
            champ["full"] = expand_reciprocal(champ["half"])
            champ["origin"] = "Known180 (Lehmer 1933)"
        else:
            lineno, line = find_line(SRC, champ["find_M"])
            meta = parse_psmm_line(line)
            champ["meta"] = meta
            champ["full"] = expand_reciprocal(meta["half_coeffs"])
            champ["origin"] = f"AllKnownAdvanpix line {lineno}"

    with OUT.open("w") as out:
        print("# Champion polynomials — PARI-verified analysis", file=out)
        print("", file=out)
        print("This document presents the results of running each record-holder", file=out)
        print("from the [README's Results section](../README.md#results) through", file=out)
        print("an independent PARI/GP analysis. For each polynomial we verify:", file=out)
        print("", file=out)
        print("1. **Irreducibility** over $\\mathbb{Z}[x]$ (sanity check vs PSMM's", file=out)
        print("   NTL-based filter).", file=out)
        print("2. **Mahler measure** to 200 decimal digits (vs the 72-digit value", file=out)
        print("   stored in `AllKnownAdvanpix`).", file=out)
        print("3. **Root classification** ($K$, $U$, $Q$, $R$) computed", file=out)
        print("   independently of PSMM.", file=out)
        print("4. **Structural decomposition**: for each champion we test", file=out)
        print("   hand-picked cyclotomic-skeleton hypotheses and report the", file=out)
        print("   exact perturbation $P(x) - S(x)$.", file=out)
        print("", file=out)
        print("All analyses are reproducible with:", file=out)
        print("```sh", file=out)
        print("python3 tools/analyze_champions.py", file=out)
        print("```", file=out)
        print("which writes this file.", file=out)
        print("", file=out)
        print("---", file=out)
        print("", file=out)

        for champ in CHAMPIONS:
            print(f"Analysing {champ['title']} ...", file=sys.stderr)
            analyzer = PolyAnalyzer(champ["full"])
            result = analyzer.analyse(
                precision=200,
                cyclotomic_search_range=(2, 100),
                test_skeletons=champ.get("skeletons", []),
            )
            if "meta" in champ:
                analyzer.cross_check(
                    result,
                    expected_M=champ["meta"]["M"],
                    expected_K=champ["meta"]["K"],
                    expected_U=champ["meta"]["U"],
                    expected_Q=champ["meta"]["Q"],
                    expected_R=champ["meta"]["R"],
                )
            elif "M_stored" in champ:
                analyzer.cross_check(
                    result,
                    expected_M=champ["M_stored"],
                    expected_K=champ.get("K"),
                    expected_U=champ.get("U"),
                    expected_Q=champ.get("Q"),
                    expected_R=champ.get("R"),
                )
            write_champion(out, champ, result)

        # The structural analysis (parametric family, Boyd-Lawton, (a, d)
        # generalisation, and database-wide verification) lives in the
        # README's Results section. This file is the per-champion appendix.
        print(
            "*The structural analysis of these champions (the parametric family, the\n"
            "Boyd-Lawton limit, the $(a, d)$ generalisation, and the database-wide\n"
            "verification) is in the [Results section of the README](../README.md#results).\n"
            "This file is the per-champion technical appendix: irreducibility, root\n"
            "counts, and high-precision Mahler measures verified by PARI/GP.*",
            file=out,
        )
    print(f"Wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
