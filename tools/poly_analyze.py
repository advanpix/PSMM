#!/usr/bin/env python3
"""
poly_analyze.py — flexible polynomial analysis using PARI/GP.

This module provides a class `PolyAnalyzer` that runs PARI to analyse
any integer polynomial. Designed to be reusable: feed it any
polynomial as a coefficient list (full or PSMM-style half) and get a
structured result.

CLI usage:
    python3 poly_analyze.py --half '1 1 0 -1 -1 -1' --degree 10
    python3 poly_analyze.py --full '1 1 0 -1 -1 -1 -1 -1 0 1 1'
    python3 poly_analyze.py --psmm-line '10 1.1762... 4 1 9 1 8 0 2 1 1 0 -1 -1 -1'

Each PARI invocation costs ~50 ms of startup. For bulk analysis of
many polynomials, batch them with `PolyAnalyzer.analyse_many`.

This file is part of "Polynomials with Small Mahler Measure" (PSMM).
Copyright (C) 2020-2026, Advanpix LLC.
License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher

Author: Pavel Holoborodko <pavel@advanpix.com>
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


# -----------------------------------------------------------------------------
# Polynomial representation helpers
# -----------------------------------------------------------------------------

def expand_reciprocal(half: Sequence[int]) -> List[int]:
    """Expand PSMM half-coefficients (a_0..a_{N/2}) to full (a_0..a_N).

    Uses reciprocity a_k = a_{N-k}.
    """
    N = 2 * (len(half) - 1)
    full = [0] * (N + 1)
    for k in range(N // 2 + 1):
        full[k] = half[k]
    for k in range(N // 2 + 1, N + 1):
        full[k] = full[N - k]
    return full


def coeffs_to_pari(coeffs: Sequence[int]) -> str:
    """Render full coefficient list as a PARI polynomial in `x`.

    e.g. [1, -2, 0, 3] -> "1 - 2*x + 3*x^3".
    """
    pieces: List[str] = []
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        sign = "+" if c > 0 else "-"
        mag = abs(c)
        if k == 0:
            body = f"{mag}"
        elif k == 1:
            body = "x" if mag == 1 else f"{mag}*x"
        else:
            body = f"x^{k}" if mag == 1 else f"{mag}*x^{k}"
        if not pieces:
            pieces.append(("-" + body) if sign == "-" else body)
        else:
            pieces.append(f"{sign} {body}")
    return " ".join(pieces) if pieces else "0"


def is_reciprocal(coeffs: Sequence[int]) -> bool:
    """Test whether a_k == a_{N-k} for all k."""
    n = len(coeffs) - 1
    return all(coeffs[k] == coeffs[n - k] for k in range(n // 2 + 1))


# -----------------------------------------------------------------------------
# PARI driver
# -----------------------------------------------------------------------------

PARI_TIMEOUT_SEC = 600


def run_pari(script: str, parisize: int = 200_000_000) -> str:
    """Execute a gp script and return its stdout.

    Raises RuntimeError on failure.
    """
    proc = subprocess.run(
        ["gp", "-q", "--default", f"parisize={parisize}"],
        input=script,
        capture_output=True,
        text=True,
        timeout=PARI_TIMEOUT_SEC,
    )
    if proc.returncode != 0 or "***" in proc.stderr:
        raise RuntimeError(
            f"PARI failed (rc={proc.returncode}):\n"
            f"  stderr: {proc.stderr[:500]}\n"
            f"  stdout: {proc.stdout[:500]}"
        )
    return proc.stdout


# -----------------------------------------------------------------------------
# Result types
# -----------------------------------------------------------------------------

@dataclass
class CyclotomicHit:
    """A cyclotomic structural relation involving Phi_d(x)."""
    kind: str            # "x^n-1", "x^n+1", "Phi_d", or "perturbation_match"
    n_or_d: int          # the n (for x^n±1) or d (for Phi_d)
    detail: str = ""     # extra info, e.g. "exact divisor", "P - S = ..."


@dataclass
class AnalysisResult:
    # Input description
    degree: int
    coeffs_pari: str
    reciprocal: bool

    # Verification
    irreducible: Optional[bool] = None
    mahler_measure: Optional[str] = None      # high-precision string
    K: Optional[int] = None
    U: Optional[int] = None
    Q: Optional[int] = None
    R: Optional[int] = None

    # Root-magnitude statistics on off-circle roots
    max_off_circle_dist: Optional[str] = None  # max |z|-1 over outside roots
    mean_off_circle_dist: Optional[str] = None

    # Structural findings
    cyclotomic_hits: List[CyclotomicHit] = field(default_factory=list)
    near_cyclotomic_skeleton: Optional[str] = None   # description if found
    skeleton_perturbation: Optional[str] = None      # P - S as PARI string

    # Sanity cross-checks (vs externally supplied values)
    cross_check: dict = field(default_factory=dict)

    # Raw script output for debugging
    raw_output: str = ""


# -----------------------------------------------------------------------------
# Analyzer
# -----------------------------------------------------------------------------

class PolyAnalyzer:
    """Run PARI-backed analyses on a single integer polynomial."""

    def __init__(self, full_coeffs: Sequence[int]):
        self.coeffs: List[int] = list(full_coeffs)
        self.N: int = len(self.coeffs) - 1
        if self.N < 1:
            raise ValueError("polynomial degree must be >= 1")
        self.poly_str: str = coeffs_to_pari(self.coeffs)
        self.reciprocal: bool = is_reciprocal(self.coeffs)

    # ----- High-level entry point -----

    def analyse(self, precision: int = 200,
                cyclotomic_search_range: Tuple[int, int] = (1, 100),
                test_skeletons: Optional[List[Tuple[str, str]]] = None
                ) -> AnalysisResult:
        """Run the full bundle of analyses and return a structured result.

        `test_skeletons` is an optional list of (label, PARI-polynomial-string)
        candidates for the "cyclotomic skeleton". For each we'll compute
        the perturbation P(x) - skeleton(x) and report its non-zero count.
        """
        result = AnalysisResult(
            degree=self.N,
            coeffs_pari=self.poly_str,
            reciprocal=self.reciprocal,
        )

        script = self._build_script(precision, cyclotomic_search_range,
                                    test_skeletons or [])
        output = run_pari(script)
        result.raw_output = output

        # Parse the line-oriented output
        for line in output.splitlines():
            line = line.strip()
            if not line or ":" not in line:
                continue
            key, _, val = line.partition(":")
            key = key.strip()
            val = val.strip()
            if key == "IRR":
                result.irreducible = (val == "1")
            elif key == "M":
                result.mahler_measure = val
            elif key in ("K", "U", "Q", "R"):
                setattr(result, key, int(val))
            elif key == "MAX_OFF":
                result.max_off_circle_dist = val
            elif key == "MEAN_OFF":
                result.mean_off_circle_dist = val
            elif key.startswith("CYC_DIV"):
                # CYC_DIV xn-1 n=12 deg=4
                mtok = re.match(r"(xn[+-]1)\s+n=(\d+)\s+deg=(\d+)", val)
                if mtok:
                    result.cyclotomic_hits.append(CyclotomicHit(
                        kind=mtok.group(1),
                        n_or_d=int(mtok.group(2)),
                        detail=f"gcd has degree {mtok.group(3)}",
                    ))
            elif key.startswith("SKEL"):
                # SKEL label=foo nz=3 perturbation=...
                mtok = re.match(r"label=(\S+)\s+nz=(\d+)\s+pert=(.+)", val)
                if mtok:
                    label, nz, pert = mtok.group(1), mtok.group(2), mtok.group(3)
                    if int(nz) <= 12:  # report only when perturbation is small
                        result.near_cyclotomic_skeleton = label
                        result.skeleton_perturbation = pert
                    result.cyclotomic_hits.append(CyclotomicHit(
                        kind="skeleton_test",
                        n_or_d=0,
                        detail=f"label={label}, nz(P-S)={nz}, P-S={pert}",
                    ))

        return result

    # ----- Script generation -----

    def _build_script(self, precision: int,
                      cyc_range: Tuple[int, int],
                      skeletons: List[Tuple[str, str]]) -> str:
        """Build a gp script. PARI requires each for(...) body on one line.

        We compose multi-statement loop bodies as semicolon-separated expressions.
        """
        n_lo, n_hi = cyc_range
        # Inner body of the root-classification loop, all one PARI expression.
        # Note: in PARI, `if(c, then, else)` returns a value; we use parens to
        # sequence side-effecting assignments inside `then`/`else`.
        root_body = (
            "r = abs(rts[i]); "
            "if (abs(r - 1) < eps_circle, "
            "    U += 1, "
            "    d = abs(r - 1); "
            "    if (d > maxoff, maxoff = d); "
            "    sumoff += d; "
            "    if (r > 1, M *= r; K += 1); "
            "    if (abs(imag(rts[i])) < eps_real, R += 1, Q += 1) "
            ")"
        )
        cyc_body = (
            "g = gcd(P, x^n - 1); "
            'if (poldegree(g) > 0, print("CYC_DIV1: xn-1 n=", n, " deg=", poldegree(g))); '
            "g = gcd(P, x^n + 1); "
            'if (poldegree(g) > 0, print("CYC_DIV2: xn+1 n=", n, " deg=", poldegree(g)))'
        )

        lines = [
            f"default(realprecision, {precision});",
            f"P = {self.poly_str};",
            "N = poldegree(P);",
            'print("IRR: ", polisirreducible(P));',
            "rts = polroots(P);",
            "M = 1.0; K = 0; U = 0; Q = 0; R = 0;",
            "eps_circle = 1e-30; eps_real = 1e-30;",
            "maxoff = 0.0; sumoff = 0.0;",
            f"for (i = 1, #rts, {root_body});",
            'print("M: ", M);',
            'print("K: ", K);',
            'print("U: ", U);',
            'print("Q: ", Q);',
            'print("R: ", R);',
            'print("MAX_OFF: ", maxoff);',
            'if (Q + R > 0, print("MEAN_OFF: ", sumoff / (Q + R)));',
            f"for (n = {n_lo}, {n_hi}, {cyc_body});",
        ]

        # Skeleton tests: for each candidate skeleton S, compute P - S and count
        # non-zero coefficients. Small NZ means S is a good "skeleton" of P.
        for i, (label, sk) in enumerate(skeletons):
            lines.append(f"S{i} = {sk};")
            lines.append(f"D{i} = P - S{i};")
            lines.append(
                f"nz{i} = 0; "
                f"for (k = 0, poldegree(D{i}), "
                f"if (polcoeff(D{i}, k) != 0, nz{i} += 1)); "
                f'print("SKEL: label={label} nz=", nz{i}, " pert=", D{i});'
            )

        return "\n".join(lines) + "\n"

    # ----- Convenience -----

    def cross_check(self, result: AnalysisResult, *,
                    expected_M: Optional[str] = None,
                    expected_K: Optional[int] = None,
                    expected_U: Optional[int] = None,
                    expected_Q: Optional[int] = None,
                    expected_R: Optional[int] = None) -> dict:
        """Compare PARI-computed values against externally supplied expectations."""
        cc: dict = {}
        if expected_M is not None and result.mahler_measure is not None:
            # Compare first 15 digits (decimal portion)
            cc["M_match"] = result.mahler_measure[:17] == expected_M[:17]
        if expected_K is not None:
            cc["K_match"] = (result.K == expected_K)
        if expected_U is not None:
            cc["U_match"] = (result.U == expected_U)
        if expected_Q is not None:
            cc["Q_match"] = (result.Q == expected_Q)
        if expected_R is not None:
            cc["R_match"] = (result.R == expected_R)
        result.cross_check = cc
        return cc


# -----------------------------------------------------------------------------
# PSMM line parsing
# -----------------------------------------------------------------------------

def parse_psmm_line(line: str) -> dict:
    """Parse a single line from AllKnownAdvanpix.

    Format: N M NNZ H L K U Q R c_0 c_1 ... c_{N/2}
    Returns a dict with keys: N, M, NNZ, H, L, K, U, Q, R, half_coeffs.
    """
    toks = line.split()
    N = int(toks[0])
    return {
        "N":   N,
        "M":   toks[1],
        "NNZ": int(toks[2]),
        "H":   int(toks[3]),
        "L":   int(toks[4]),
        "K":   int(toks[5]),
        "U":   int(toks[6]),
        "Q":   int(toks[7]),
        "R":   int(toks[8]),
        "half_coeffs": [int(t) for t in toks[9:]],
    }


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def _parse_intlist(s: str) -> List[int]:
    s = s.replace(",", " ")
    return [int(t) for t in s.split() if t.strip()]


def main() -> int:
    ap = argparse.ArgumentParser(description="PARI-backed polynomial analysis.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--full",  help="full coefficients a_0..a_N (whitespace or comma separated)")
    g.add_argument("--half",  help="half coefficients a_0..a_{N/2} (reciprocal expansion assumed)")
    g.add_argument("--psmm-line", help="a single line in AllKnownAdvanpix format")
    g.add_argument("--from-file", help="path to a file containing AllKnownAdvanpix-format lines (analyse all)")
    ap.add_argument("--precision", type=int, default=200,
                    help="PARI real precision (decimal digits, default 200)")
    ap.add_argument("--cyc-range", default="1,100",
                    help="lo,hi for cyclotomic gcd scan (default 1,100)")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of text")
    args = ap.parse_args()

    cyc_lo, cyc_hi = (int(x) for x in args.cyc_range.split(","))

    def _print(result: AnalysisResult, meta: Optional[dict] = None):
        if args.json:
            d = asdict(result)
            d.pop("raw_output", None)
            if meta:
                d["meta"] = meta
            print(json.dumps(d, indent=2, default=str))
            return
        print(f"--- degree {result.degree}{(' (' + meta['origin'] + ')') if meta and 'origin' in meta else ''} ---")
        print(f"  irreducible : {result.irreducible}")
        print(f"  reciprocal  : {result.reciprocal}")
        print(f"  M(P)        : {result.mahler_measure}")
        print(f"  K U Q R     : {result.K} {result.U} {result.Q} {result.R}")
        if result.max_off_circle_dist is not None:
            print(f"  max |r|-1   : {result.max_off_circle_dist}")
        if result.cross_check:
            print(f"  cross-check : {result.cross_check}")
        for hit in result.cyclotomic_hits:
            print(f"  hit         : {hit.kind} {hit.n_or_d}  ({hit.detail})")

    if args.full:
        coeffs = _parse_intlist(args.full)
        a = PolyAnalyzer(coeffs)
        _print(a.analyse(precision=args.precision, cyclotomic_search_range=(cyc_lo, cyc_hi)))
    elif args.half:
        half = _parse_intlist(args.half)
        coeffs = expand_reciprocal(half)
        a = PolyAnalyzer(coeffs)
        _print(a.analyse(precision=args.precision, cyclotomic_search_range=(cyc_lo, cyc_hi)))
    elif args.psmm_line:
        meta = parse_psmm_line(args.psmm_line)
        coeffs = expand_reciprocal(meta["half_coeffs"])
        a = PolyAnalyzer(coeffs)
        result = a.analyse(precision=args.precision, cyclotomic_search_range=(cyc_lo, cyc_hi))
        a.cross_check(result,
                      expected_M=meta["M"],
                      expected_K=meta["K"],
                      expected_U=meta["U"],
                      expected_Q=meta["Q"],
                      expected_R=meta["R"])
        _print(result, meta)
    elif args.from_file:
        path = Path(args.from_file)
        with path.open() as f:
            for ln in f:
                if not ln.strip() or ln.lstrip().startswith("#"):
                    continue
                try:
                    meta = parse_psmm_line(ln)
                except Exception as e:
                    print(f"skip: {e}", file=sys.stderr)
                    continue
                coeffs = expand_reciprocal(meta["half_coeffs"])
                a = PolyAnalyzer(coeffs)
                try:
                    res = a.analyse(precision=args.precision,
                                    cyclotomic_search_range=(cyc_lo, cyc_hi))
                    a.cross_check(res, expected_M=meta["M"],
                                  expected_K=meta["K"], expected_U=meta["U"],
                                  expected_Q=meta["Q"], expected_R=meta["R"])
                    _print(res, meta)
                except Exception as e:
                    print(f"failed degree {meta['N']}: {e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
