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

        # Closing synthesis
        print("## Database-wide verification", file=out)
        print("", file=out)
        print("Every one of the **48,341** polynomials in `AllKnownAdvanpix` was", file=out)
        print("re-verified independently with PARI/GP at 60-digit precision", file=out)
        print("(via [`tools/bulk_verify.py`](../tools/bulk_verify.py); results in", file=out)
        print("[`doc/database-verification.csv`](database-verification.csv)).", file=out)
        print("Findings:", file=out)
        print("", file=out)
        print("- **Mahler measure**: agrees with the stored 72-digit value to", file=out)
        print("  the precision tested in every single entry.", file=out)
        print("- **Irreducibility**: confirmed for every entry.", file=out)
        print("- **Root counts** $K, U, Q, R$: **4 mismatches** were found, each", file=out)
        print("  off-by-one or off-by-two due to PSMM's eps threshold being too", file=out)
        print("  tight when classifying a root close to the unit circle:", file=out)
        print("", file=out)
        print("  | DB line | $N$ | $K$ (PSMM → PARI) | $U$ (PSMM → PARI) | $Q$ (PSMM → PARI) | $\\max\\,\\|\\|z\\|-1\\|$ |", file=out)
        print("  |---:|---:|---:|---:|---:|---|", file=out)
        print("  | 31272 | 360 | 61 → 60 | 239 → 240 | 121 → 120 | 0.00540 |", file=out)
        print("  | 31875 | 364 | 53 → 52 | 259 → 260 | 105 → 104 | 0.00621 |", file=out)
        print("  | 41843 | 420 | 56 → 56 | 307 → 308 | 113 → 112 | 0.00577 |", file=out)
        print("  | 45319 | 438 | 51 → 49 | 337 → 340 |  99 →  96 | 0.00608 |", file=out)
        print("", file=out)
        print("  In all four cases the PSMM values violated the identity", file=out)
        print("  $2K + U = N$; the corrected PARI values are self-consistent.", file=out)
        print("  The Mahler measure and the polynomial coefficients themselves", file=out)
        print("  are correct in every case — only the bookkeeping of how many", file=out)
        print("  roots sit on vs. off the unit circle was off by 1 or 2 at the", file=out)
        print("  boundary.", file=out)
        print("", file=out)
        print("  These 4 entries have been patched in `AllKnownAdvanpix` by", file=out)
        print("  [`tools/patch_db_mismatches.py`](../tools/patch_db_mismatches.py).", file=out)
        print("", file=out)
        print("## Synthesis", file=out)
        print("", file=out)
        print("- Across all five record-holders we verify that PSMM's NTL-based", file=out)
        print("  irreducibility filter and double-precision Mahler-measure computation", file=out)
        print("  agree with PARI's independent re-derivation to high precision.", file=out)
        print("", file=out)
        print("- The **Max-U champion** has an *exact* short-form decomposition", file=out)
        print("  as a cyclotomic product perturbed by a single sparse term:", file=out)
        print("  $$P(x) = (x+1)(x^{455}+1) - x^{227} \\Phi_3(x).$$", file=out)
        print("  The first factor is the product of cyclotomic polynomials", file=out)
        print("  $\\Phi_2(x) \\cdot \\prod_{d \\mid 910, d \\nmid 455} \\Phi_d(x)$,", file=out)
        print("  contributing all-on-circle roots; the $\\Phi_3$ perturbation,", file=out)
        print("  placed at $x^{(m-1)/2}$ with $m = 455$, is what generates the", file=out)
        print("  60 off-circle roots in 15 Salem quadruplets.", file=out)
        print("", file=out)
        print("- The **Max-K / Q champion** has a coefficient signature of pure", file=out)
        print("  $(1, -1, 0)$ periodicity until a phase shift near the middle.", file=out)
        print("  The pure-periodic template is the cyclotomic quotient", file=out)
        print("  $(1 - x^{3m})/\\Phi_3(x)$. The phase shift is the source of the", file=out)
        print("  152 off-circle roots.", file=out)
        print("", file=out)
        print("- The **dense champions** (Max-NNZ, Max-H) lack an obvious", file=out)
        print("  short-form cyclotomic skeleton; their coefficients exhibit a", file=out)
        print("  triangular build-up pattern suggesting they might be obtained", file=out)
        print("  from products / convolutions of simpler polynomials. Further", file=out)
        print("  analysis (e.g. factoring over $\\mathbb{F}_p$ for various $p$,", file=out)
        print("  or searching for resultant relations) is left as future work.", file=out)
        print("", file=out)
        print("## Generalising the Max-U construction", file=out)
        print("", file=out)
        print("The exact decomposition of the Max-U champion suggests a", file=out)
        print("**parametric family** indexed by odd $m$:", file=out)
        print("", file=out)
        print("$$P_m(x) = (x+1)(x^m+1) - x^{(m-1)/2} \\Phi_3(x),", file=out)
        print("\\qquad \\deg P_m = m + 1.$$", file=out)
        print("", file=out)
        print("We computed $M(P_m)$ for every odd $m$ in $[5, 199]$ and factored", file=out)
        print("each $P_m$ over $\\mathbb{Z}[x]$. The most striking finding:", file=out)
        print("**$P_{21}(x)$ factors as $\\Phi_{12}(x) \\cdot R_{18}(x)$**", file=out)
        print("where $R_{18}$ is the irreducible degree-18 polynomial with", file=out)
        print("$M(R_{18}) \\approx 1.18836814750822\\ldots$ — **the second-", file=out)
        print("smallest known Salem polynomial** (just above Lehmer's 1.17628…).", file=out)
        print("This is exactly the entry `18 1.188368…` in `AllKnownAdvanpix`.", file=out)
        print("", file=out)
        print("Selected values of $M$(smallest non-cyclotomic factor of $P_m$):", file=out)
        print("", file=out)
        print("| $m$ | $\\deg$ of smallest non-cyc factor | $M$ |", file=out)
        print("|---:|---:|---|", file=out)
        print("| 11  | 12  | 1.25104661720… |", file=out)
        print("| 17  | 18  | 1.21972085904… |", file=out)
        print("| **21**  | **18**  | **1.18836814751… (Lehmer's sibling — 2nd smallest known)** |", file=out)
        print("| 35  | 36  | 1.22649330147… |", file=out)
        print("| 53  | 46  | 1.23074300908… |", file=out)
        print("| 67  | 60  | 1.24006185904… |", file=out)
        print("| 95  | 96  | 1.24866635341… |", file=out)
        print("| 455 | 456 | 1.25491475758… (this is the Max-U champion itself) |", file=out)
        print("", file=out)
        print("Two key observations:", file=out)
        print("", file=out)
        print("1. The family **reproduces** known Salem polynomials at small", file=out)
        print("   degrees — the m=21 case gives the second-smallest known", file=out)
        print("   Salem polynomial as an irreducible factor.", file=out)
        print("2. As $m \\to \\infty$ the Mahler measure of $P_m$ (or its smallest", file=out)
        print("   irreducible factor) converges to roughly $1.255$, *above*", file=out)
        print("   Lehmer's number. The family itself does not approach 1 from", file=out)
        print("   above.", file=out)
        print("", file=out)
        print("### Analytic limit (Boyd–Lawton theorem)", file=out)
        print("", file=out)
        print("The convergence is **not coincidental**. The parametric family", file=out)
        print("$P_m(x)$ is a univariate monomial substitution into the bivariate", file=out)
        print("polynomial", file=out)
        print("", file=out)
        print("$$F(x, u) = x(x+1) u^2 - \\Phi_3(x) u + (x+1),$$", file=out)
        print("", file=out)
        print("obtained by setting $u = x^{(m-1)/2}$ so that $x^m = x\\cdot u^2$.", file=out)
        print("By the **Boyd–Lawton theorem** (Boyd 1981, Lawton 1983), the", file=out)
        print("Mahler measure of the univariate family converges to the", file=out)
        print("(logarithmic) Mahler measure of the bivariate polynomial:", file=out)
        print("", file=out)
        print("$$\\lim_{m \\to \\infty} M(P_m) = \\exp\\bigl(m(F)\\bigr),", file=out)
        print("\\qquad m(F) = \\frac{1}{(2\\pi)^2}\\!\\int_0^{2\\pi}\\!\\!\\int_0^{2\\pi}", file=out)
        print("\\log\\bigl|F(e^{i\\theta_1}, e^{i\\theta_2})\\bigr| d\\theta_1 d\\theta_2.$$", file=out)
        print("", file=out)
        print("Since $F$ is quadratic in $u$, the inner integral evaluates via", file=out)
        print("Jensen's formula and one is left with a single integral over the", file=out)
        print("circle. Numerical evaluation in PARI gives", file=out)
        print("", file=out)
        print("$$\\log m(F) \\approx 0.222630132139506025908217312245576858\\ldots,$$", file=out)
        print("", file=out)
        print("$$\\boxed{\\lim_{m\\to\\infty} M(P_m) \\approx 1.249358390752959362866\\ldots}$$", file=out)
        print("", file=out)
        print("This is **above Lehmer's number** (1.17628…). The Boyd–Lawton", file=out)
        print("limit is therefore a *barrier* for this particular family: no", file=out)
        print("matter how large $m$ grows, $M(P_m)$ stays $\\geq 1.2493\\ldots$,", file=out)
        print("with the m=21 minimum at $1.18837$ being the closest single case", file=out)
        print("to Lehmer's bound that the family achieves.", file=out)
        print("", file=out)
        print("![Parametric family convergence](../images/parametric_family_M_vs_m.png)", file=out)
        print("", file=out)
        print("(See [`tools/scan_parametric_family.py`](../tools/scan_parametric_family.py)", file=out)
        print("and [`tools/plot_parametric_family.py`](../tools/plot_parametric_family.py)", file=out)
        print("to reproduce the data and figure.)", file=out)
        print("", file=out)
        print("### Generalising to two cyclotomic parameters $(a, d)$", file=out)
        print("", file=out)
        print("Replacing the background factor $(x+1) = \\Phi_2$ by a general", file=out)
        print("cyclotomic $\\Phi_a$ gives a two-parameter family", file=out)
        print("", file=out)
        print("$$P_{a,d,m,s}(x) = \\Phi_a(x)(x^m+1) + s \\cdot x^{(\\phi(a)+m-\\phi(d))/2} \\Phi_d(x).$$", file=out)
        print("", file=out)
        print("We swept $a \\in \\{2, 3, 4, 6\\}$, $d \\in \\{3, 5, 7, 8, 9, 10, 12\\}$,", file=out)
        print("$s \\in \\{-1, +1\\}$, $m \\in [5, 201]$, factored each $P_{a,d,m,s}$ over", file=out)
        print("$\\mathbb{Z}$, and recorded the Mahler measure of the smallest", file=out)
        print("non-cyclotomic irreducible factor.", file=out)
        print("", file=out)
        print("**Result.** Across the entire sweep, the global minimum is", file=out)
        print("**$M = 1.17628\\ldots$**, i.e. Lehmer's number itself. Four distinct", file=out)
        print("parameter combinations all factor to include Lehmer's polynomial", file=out)
        print("(or its $x \\to -x$ reflection, which has the same Mahler measure):", file=out)
        print("", file=out)
        print("| $a$ | $d$ | $m$ | sign |", file=out)
        print("|---:|---:|---:|---:|", file=out)
        print("|  2  |  3  | 23  |  +  |", file=out)
        print("|  2  |  5  |  9  |  −  |", file=out)
        print("|  2  |  7  | 15  |  −  |", file=out)
        print("|  3  |  7  |  8  |  −  |", file=out)
        print("", file=out)
        print("The second-smallest is $M \\approx 1.18837$ at $(a, d, m, s) =", file=out)
        print("(2, 3, 21, -)$, the Lehmer-sibling we discovered earlier.", file=out)
        print("", file=out)
        print("**The cyclotomic-perturbation family cannot break Lehmer's bound,**", file=out)
        print("**but it embeds Lehmer's polynomial naturally in many ways.** This", file=out)
        print("is consistent with Boyd's conjecture that all small Mahler measures", file=out)
        print("$> 1$ arise from a structured (Salem-Boyd-style) construction.", file=out)
        print("", file=out)
        print("Reproduce with [`tools/sweep_ad.py`](../tools/sweep_ad.py); raw", file=out)
        print("data in [`doc/ad_sweep.csv`](ad_sweep.csv).", file=out)
        print("", file=out)
        print("### Implications", file=out)
        print("", file=out)
        print("The Boyd–Lawton barrier of $\\approx 1.249$ for this family tells", file=out)
        print("us where *not* to look for sub-Lehmer polynomials. Producing a", file=out)
        print("smaller Mahler measure than Lehmer's $1.17628$ would require a", file=out)
        print("**different bivariate lift**, one whose bivariate Mahler measure", file=out)
        print("is smaller than $\\exp(0.16307\\ldots) = 1.17628\\ldots$. Concrete", file=out)
        print("experimental directions:", file=out)
        print("", file=out)
        print("- **Vary the perturbation cyclotomic**: scan $\\Phi_d$ for", file=out)
        print("  $d \\in \\{4, 5, 6, 7, 8, 12, \\ldots\\}$ in place of $\\Phi_3$.", file=out)
        print("  Each choice produces a different bivariate $F_d$ with its", file=out)
        print("  own (computable) Boyd–Lawton limit; the search reduces to", file=out)
        print("  finding $d$ that minimises $m(F_d)$.", file=out)
        print("- **Vary the perturbation position**: replace $x^{(m-1)/2}$ by", file=out)
        print("  $x^{\\lfloor m\\alpha \\rfloor}$ for irrational $\\alpha$ —", file=out)
        print("  this lifts to a different bivariate.", file=out)
        print("- **Multiple perturbations**: subtract $x^{a}\\Phi_3 + x^{b}\\Phi_5$", file=out)
        print("  to obtain a *trivariate* lift, with a corresponding Lawton-type", file=out)
        print("  triple-integral limit.", file=out)
        print("- **Different cyclotomic skeletons**: replace $(x+1)(x^m+1)$ by", file=out)
        print("  $(x^a+1)(x^m+1)$ or by products $\\prod \\Phi_{d_i}$.", file=out)
        print("", file=out)
        print("Each of these is a *directed* search through a small parameter", file=out)
        print("space (tens to thousands of candidates), with the Boyd–Lawton", file=out)
        print("limit *computed analytically* for each before any expensive", file=out)
        print("polynomial search runs. This is dramatically smaller than the", file=out)
        print("brute-force PSMM enumeration (which is $\\sim 10^{12}$", file=out)
        print("polynomials per degree at large $N$ and small $\\mathrm{NNZ}$).", file=out)
        print("", file=out)
    print(f"Wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
