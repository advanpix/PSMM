[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# PSMM — Polynomials with Small Mahler Measure

Exhaustive search for primitive, irreducible, reciprocal integer polynomials
whose [Mahler measure](https://en.wikipedia.org/wiki/Mahler_measure) is below
a given threshold (typically M(p) < 1.3). The smallest known Mahler measure
greater than 1 is **Lehmer's number**, achieved by a degree-10 polynomial
discovered by D. H. Lehmer in 1933:

$$P(x) = x^{10} + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1, \qquad M(P) \approx 1.17628081825991750654\ldots$$

Whether any smaller value $M(P) > 1$ exists is **Lehmer's problem** —
one of the oldest open questions in number theory.

<table>
  <tr>
    <td align="center">
      <a href="images/lehmer.png">
        <img src="images/lehmer_1k.png" width="100%" alt="Newton convergence basins for Lehmer's polynomial"/>
      </a>
      <br/>
      <sub><b>Newton convergence basins</b> for Lehmer's polynomial.<br/>
      The 10 roots are marked with black crosses; each colored region is the<br/>
      basin of attraction of one root under Newton's iteration.</sub>
    </td>
    <td align="center">
      <a href="images/lehmer_0.25.png">
        <img src="images/lehmer_0.25_1k.png" width="100%" alt="Zoom-in showing fractal basin boundaries"/>
      </a>
      <br/>
      <sub><b>Zoom-in</b> on a 0.25&times;0.25 region of the basin boundary,<br/>
      revealing the fractal self-similar structure.<br/>
      Click any image for the full 4096&times;4096 version.</sub>
    </td>
  </tr>
</table>

## Results

The `AllKnownAdvanpix` database contains **48,341 primitive irreducible
reciprocal polynomials** with Mahler measure $M(P) < 1.3$, collected by
exhaustive search up to degree 456. It is a strict superset of
Mossinghoff's historical `Known180` list (~8,500 polynomials up to
degree 180) — **all polynomials of degree > 180 in our database are new
discoveries made with this search program**.

### Definitions

For a reciprocal polynomial
$P(x) = a_N x^N + a_{N-1} x^{N-1} + \cdots + a_1 x + a_0$
with $a_k = a_{N-k}$ (so $a_0 = a_N = 1$):

| Symbol | Definition |
|---|---|
| $\deg P = N$ | Polynomial degree (even). |
| $M(P) = \|a_N\| \cdot \prod_{\|z_i\|>1} \|z_i\|$ | **Mahler measure** — product of $\|a_N\|$ and root magnitudes outside the unit disk. |
| $\mathrm{NNZ}$ | Number of non-zero coefficients among $a_1, \ldots, a_{N/2}$ (i.e. the **half-coefficients**, excluding the fixed $a_0 = 1$). |
| $H(P) = \max_k \|a_k\|$ | **Height** — maximum coefficient magnitude (over the full polynomial). |
| $L(P) = \sum_k \|a_k\|$ | **Length** — sum of coefficient magnitudes (over the full polynomial). |
| $K$ | Count of (individual) roots strictly outside the unit disk ($\|z\| > 1$). |
| $U$ | Count of roots on the unit circle. They come in complex-conjugate pairs $\\{z, \bar z\\}$, so $U$ is even. |
| $Q$ | Count of complex non-unity roots ($\|z\| \neq 1$, $\mathrm{Im}(z) \neq 0$). They come in Salem quadruplets $\\{z, \bar z, 1/z, 1/\bar z\\}$, so $Q$ is a multiple of 4. |
| $R$ | Count of real non-unity roots ($\mathrm{Im}(z) = 0$, $\|z\| \neq 1$). They come in reciprocal pairs $\\{z, 1/z\\}$, so $R$ is even. |

For a reciprocal polynomial the roots satisfy $z \leftrightarrow 1/z$, so the
number of roots inside the unit disk equals the number outside.
This gives two useful identities:

$$2K + U = N, \qquad Q + R = 2K.$$

Storage in `AllKnownAdvanpix` is by half-coefficients $(a_0, a_1, \ldots, a_{N/2})$;
the full polynomial is recovered by reciprocity $a_k = a_{N-k}$.

### Record holders

| Record | Value | Origin | Degree | M(P) | NNZ | H | L | K | U | Q | R |
|---|---:|:---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| Smallest M (Lehmer) | M &asymp; 1.17628 | Known180 | 10 | 1.17628… | 4 | 1 | 9 | 1 | 8 | 0 | 2 |
| Most real non-unity | R = 4 | Known180 | 20 | 1.25363… | 8 | 1 | 17 | 2 | 16 | 0 | 4 |
| Largest height | H = 29 | **New** | 348 | 1.25431… | 168 | 29 | 3657 | 40 | 268 | 80 | 0 |
| Most non-zero coeffs | NNZ = 212 | **New** | 432 | 1.25541… | 212 | 23 | 4173 | 42 | 348 | 84 | 0 |
| Largest length | L = 4173 | **New** | 432 | 1.25541… | 212 | 23 | 4173 | 42 | 348 | 84 | 0 |
| Most roots outside &#124;z&#124;=1 | K = 76 | **New** | 452 | 1.28530… | 151 | 1 | 303 | 76 | 300 | 152 | 0 |
| Most complex non-unity | Q = 152 | **New** | 452 | 1.28530… | 151 | 1 | 303 | 76 | 300 | 152 | 0 |
| Most roots on &#124;z&#124;=1 | U = 396 | **New** | 456 | 1.25491… | 3 | 1 | 7 | 30 | 396 | 60 | 0 |

Six of the eight category records are new discoveries — they live at
degrees 348, 432, 452, and 456, well above the 180-degree boundary of the
prior literature. The two records inside the classical regime (Lehmer's
M-record at degree 10, and the R-record at degree 20) are pre-existing
historical entries included for completeness.

The structural analysis of these polynomials surfaced a concrete
**parametric construction** that reproduces the second-smallest known
Salem polynomial (at $m = 21$) and converges, by the Boyd–Lawton theorem,
to a specific bivariate Mahler measure. The full development — including
PARI/GP-verified decompositions, the parametric scan over odd $m$ up to
1001, and the analytic limit — is in
[`doc/champion-analysis.md`](doc/champion-analysis.md).

### Sparsest extremal polynomial — the Max-U record (New)

A degree-456 polynomial with only **three non-zero half-coefficients** packs
**396 roots onto the unit circle** in conjugate pairs:

$$P(x) = x^{456} + x^{455} - x^{229} - x^{228} - x^{227} + x + 1$$

$$M(P) \approx 1.25491475757884793378\ldots, \quad \deg P = 456, \quad U = 396, \quad K = 30, \quad Q = 60$$

**Exact structural decomposition.** This polynomial is **the cyclotomic
product $(x+1)(x^{455}+1)$ perturbed by a single sparse $\Phi_3$ term:**

$$P(x) = (x+1)(x^{455} + 1) - x^{227} \Phi_3(x), \qquad \Phi_3(x) = x^2 + x + 1.$$

The first factor is a product of cyclotomic polynomials — all 456 of its
roots sit on the unit circle. The single perturbation term $-x^{227}\Phi_3(x)$,
placed at the half-degree, pushes exactly 60 of those roots off the circle
into 15 Salem quadruplets (each at distance $\sim 0.008$ from $|z|=1$),
producing the irreducible polynomial above. **Generalising the construction
to odd $m$,**

$$P_m(x) = (x+1)(x^m + 1) - x^{(m-1)/2} \Phi_3(x),$$

rediscovers, after factoring out cyclotomic content, the second-smallest
known Salem polynomial at $m = 21$
($M \approx 1.18836\ldots$, degree 18). Boyd–Lawton predicts
$\lim_{m\to\infty} M(P_m) = \exp(m(F)) \approx 1.24936\ldots$, where $F$
is the bivariate "lift" of the family. See
[`doc/champion-analysis.md`](doc/champion-analysis.md) for the full
analysis.

### Polynomial with the most real non-unity roots (Known180)

The degree-20 polynomial with the smallest Mahler measure among those with
$R = 4$ real non-unity roots (two reciprocal pairs):

$$P(x) = \sum_{k=17}^{20} x^k - \sum_{k=6}^{14} x^k + \sum_{k=0}^{3} x^k$$

$$M(P) \approx 1.25363556570886317997\ldots, \quad \deg P = 20, \quad K = 2, \quad U = 16, \quad R = 4$$

Equivalently, with all 21 coefficients written out:

$$P(x) = x^{20} + x^{19} + x^{18} + x^{17} - x^{14} - x^{13} - x^{12} - x^{11} - x^{10} - x^9 - x^8 - x^7 - x^6 + x^3 + x^2 + x + 1$$

This polynomial is part of the historical `Known180` list (it predates
the present project), and remains the R-record in the combined database.

### Densest extremal polynomials (New)

The polynomials maximising height ($H$), length ($L$), non-zero count
($\mathrm{NNZ}$), and root count outside the unit disk ($K$) all live near
degrees 350–452 and have hundreds of non-zero coefficients with intricate
combinatorial structure. Click any thumbnail below for the full polynomial
typeset in LaTeX:

<table>
  <tr>
    <td align="center" width="33%">
      <a href="images/champion_nnz212.png">
        <img src="images/champion_nnz212.png" width="100%" alt="NNZ=212 champion polynomial"/>
      </a>
      <br/>
      <sub><b>NNZ = 212 / L = 4173 / H = 23</b><br/>
      degree 432, M &asymp; 1.25541<br/>
      212 non-zero half-coefficients,<br/>
      height 23, length 4173</sub>
    </td>
    <td align="center" width="33%">
      <a href="images/champion_h29.png">
        <img src="images/champion_h29.png" width="100%" alt="H=29 champion polynomial"/>
      </a>
      <br/>
      <sub><b>H = 29</b><br/>
      degree 348, M &asymp; 1.25431<br/>
      168 non-zero half-coefficients,<br/>
      max single coefficient 29</sub>
    </td>
    <td align="center" width="33%">
      <a href="images/champion_k76.png">
        <img src="images/champion_k76.png" width="100%" alt="K=76 champion polynomial"/>
      </a>
      <br/>
      <sub><b>K = 76 / Q = 152</b><br/>
      degree 452, M &asymp; 1.28530<br/>
      striking periodic (1, &minus;1, 0)<br/>
      coefficient signature</sub>
    </td>
  </tr>
</table>

## How it works

### Search space

A monic reciprocal polynomial of even degree $N$ has the form

$$P(x) = x^N + a_1 x^{N-1} + a_2 x^{N-2} + \cdots + a_2 x^2 + a_1 x + 1, \qquad a_k = a_{N-k},$$

so the polynomial is fully determined by the $N/2$ free **half-coefficients**
$a_1, a_2, \ldots, a_{N/2}$ (the leading and constant coefficient $a_0 = a_N = 1$
are fixed). PSMM enumerates all such polynomials with:

- **integer coefficients** drawn from a user-supplied alphabet
  $\mathcal{A} = \\{c_0, c_1, \ldots, c_{b-1}\\}$ (e.g. $\\{-1, 1\\}$, $b = 2$);
- exactly **$k$ non-zero half-coefficients** (sparsity constraint).

#### The bijection that makes everything fast

The two constraints — *which* of the $N/2$ positions are non-zero, and
*what* values they take — decouple cleanly. Pick the **sparsity pattern**
(a subset of size $k$ chosen from $\\{1, 2, \ldots, N/2\\}$), then assign
each chosen position a value from $\mathcal{A}$. The total count is

$$Q(N, k, b) = \binom{N/2}{k} \cdot b^{k}.$$

The enumeration is implemented as a **bijection** between
$\\{0, 1, \ldots, Q(N,k,b) - 1\\}$ and the candidate polynomial set:

- The pattern index runs through all $\binom{N/2}{k}$ subsets via
  `std::next_permutation`.
- The value-assignment index is a **mixed-radix counter** with base $b$ and
  $k$ digits — literally `m_Number[k-1] m_Number[k-2] ... m_Number[0]` in
  the iterator, incremented by add-with-carry.

So the entire search space is, structurally, a single integer counter from
$0$ to $Q - 1$. Three consequences:

1. **Resumability**: a search can be paused and resumed by recording the
   counter value.
2. **Parallelism**: any contiguous interval $[i, j)$ of the counter is an
   independent sub-search — workers can be handed disjoint intervals with no
   coordination beyond merging results at the end. PSMM uses this internally
   (batches of 256 polynomials feed a thread pool, see [`psmm.cpp`](psmm.cpp)),
   and externally for cluster-scale searches (split the counter by degree or
   by nnz across machines, then `-merge` the result files).
3. **Symmetry pruning**: polynomials related by the substitution $x \mapsto -x$
   (which preserves Mahler measure) are detected from the pattern alone and
   skipped, halving the effective work for symmetric alphabets like $\\{-1, 1\\}$.

#### Scale considerations

| Degree $N$ | nnz $k$ | $\mathcal{A}$ | $Q(N,k,b)$ |
|---:|---:|---|---:|
| 100 | 3 | $\{-1,1\}$ | $\binom{50}{3} \cdot 2^3 \approx 1.6 \times 10^5$ |
| 200 | 3 | $\{-1,1\}$ | $\binom{100}{3} \cdot 2^3 \approx 1.3 \times 10^6$ |
| 400 | 3 | $\{-1,1\}$ | $\binom{200}{3} \cdot 2^3 \approx 1.0 \times 10^7$ |
| 200 | 4 | $\{-1,1\}$ | $\binom{100}{4} \cdot 2^4 \approx 6.3 \times 10^7$ |
| 200 | 5 | $\{-1,1\}$ | $\binom{100}{5} \cdot 2^5 \approx 2.4 \times 10^9$ |
| 100 | 5 | $\{-1,0,1\}$ | $\binom{50}{5} \cdot 3^5 \approx 5.1 \times 10^8$ |

Exhaustive search is comfortable at low nnz with a binary alphabet up to
degree several hundred, but grows quickly with both. Most production runs
fix $\mathcal{A} = \{-1, 1\}$ and sweep over $k \in \{1, 2, 3, 4\}$.

### Mahler measure computation

For each candidate polynomial, the Mahler measure is computed by finding
all roots using [MPSolve](https://github.com/robol/MPSolve) (Aberth
method, arbitrary precision) and taking the product of the absolute values
of roots outside the unit circle:

    M(p) = |a_N| * prod_{|z_i| > 1} |z_i|

Two modes are available:

- **Fast estimator** (`USE_FAST_MAHLER_ESTIMATOR`, default ON): a patched
  MPSolve that aborts root-finding early as soon as any root's inclusion
  disk lies entirely outside |z| = threshold. Returns immediately with
  "over threshold" for polynomials that clearly exceed the bound.
- **Full computation**: finds all roots to double precision (64-bit GMP
  limb), computes M(p) exactly.

### Deduplication and verification

Candidates with M(p) &le; threshold are checked against a list of
previously known polynomials (`-known` file) and against polynomials
found earlier in the current session. Matching is by Mahler measure
within a tolerance of 1e-14, restricted to polynomials of equal or
lower degree (a higher-degree polynomial with the same M is a multiple
of a known one).

Surviving candidates are then:

1. **Factored** over Z using [NTL](https://libntl.org/). Reducible
   polynomials are split into irreducible factors; each factor with
   M(p) &le; threshold is checked independently.
2. **Verified** by recomputing the Mahler measure in extended precision
   (72 decimal digits, ~240 bits) and cross-checking against the known
   list at that precision.

Only primitively new, irreducible polynomials with verified Mahler measure
are written to the output.

### Parallelism

The outer search loop is parallelized across `-threads` worker threads.
Each worker computes the Mahler measure for a batch of polynomials
independently (one MPSolve context per call, single-threaded internally).
Deduplication and logging run sequentially on the main thread after each
batch completes.

## Building

### Prerequisites

| Dependency | Purpose | Install (Ubuntu/Debian) |
|---|---|---|
| CMake &ge; 3.20 | Build system | `sudo apt install cmake` |
| GCC &ge; 10 or Clang &ge; 12 | C++17 compiler | `sudo apt install g++` |
| GMP + GMPXX | Arbitrary-precision arithmetic | `sudo apt install libgmp-dev` |
| NTL | Polynomial factorization over Z | `sudo apt install libntl-dev` |
| autoconf, automake, libtool | MPSolve build (autotools) | `sudo apt install autoconf automake libtool` |
| bison, flex | MPSolve parser generator | `sudo apt install bison flex` |
| pkg-config | MPSolve configure | `sudo apt install pkg-config` |
| pthreads | Threading | (included with glibc) |

One-liner for Ubuntu/Debian:

```sh
sudo apt install cmake g++ libgmp-dev libntl-dev autoconf automake libtool bison flex pkg-config
```

### Build

```sh
cmake -B build -S .
cmake --build build -j
```

MPSolve 3.2.2 is downloaded, patched, and built automatically as part of the
CMake configure step (via `ExternalProject`). No system-wide MPSolve
installation is needed.

### Run tests

```sh
cd build && ctest --output-on-failure
```

## Usage

PSMM operates in three modes: **search**, **merge**, and **analyze**.

### Search

Enumerate all reciprocal polynomials of a given degree and report those
with Mahler measure below the threshold.

```sh
./build/psmm \
    -degree=N \
    -coeffs=-1,1 \
    -nnz=1,2,3 \
    -threshold=1.3 \
    -threads=8 \
    -period=3600 \
    -known=AllKnownAdvanpix \
    -addto=AllKnownAdvanpix
```

| Argument | Description |
|---|---|
| `-degree=N` | Polynomial degree (must be even and positive). |
| `-coeffs=c0,c1,...` | Comma-separated list of allowed coefficient values. |
| `-nnz=n1,n2,...` | Comma-separated list of non-zero half-coefficient counts to search. |
| `-threshold=T` | Upper bound on Mahler measure (typically 1.3). |
| `-threads=N` | Number of parallel worker threads (default 1). |
| `-period=S` | Progress report interval in seconds (default 5). |
| `-known=FILE` | File of previously known polynomials to skip (optional). |
| `-addto=FILE` | Append verified results to this file (optional). |

**Example — continue the search from existing results:**

```sh
# Search degree 100, coefficients {-1, 1}, nnz 1-3, using 8 threads.
# Skip polynomials already in AllKnownAdvanpix; append new finds to it.
./build/psmm \
    -degree=100 \
    -coeffs=-1,1 \
    -nnz=1,2,3 \
    -threshold=1.3 \
    -threads=8 \
    -period=3600 \
    -known=AllKnownAdvanpix \
    -addto=AllKnownAdvanpix \
    > 100_psmm.txt
```

The log file (`100_psmm.txt`) contains:

- `***` lines: newly found polynomials (not in `-known`, not seen earlier in this session).
- `+++` lines: polynomials found again in this session (duplicate within the run).
- `---` lines: polynomials already in the `-known` file (with degree of the match in parentheses).
- Periodic progress lines: polynomials per second, estimated time remaining.

**Resuming from a crashed search log:**

If a search crashes mid-run (e.g. out of memory), the intermediate `***`
results can be extracted from the log and re-verified:

```sh
./build/psmm \
    -degree=100 \
    -coeffs=-1,1 \
    -nnz=1,2,3 \
    -threshold=1.3 \
    -threads=8 \
    -fromlog=100_psmm.txt \
    -known=AllKnownAdvanpix \
    -addto=AllKnownAdvanpix
```

### Merge

Combine multiple result files, remove duplicates, and produce a single
sorted output. To fold new search results into the master file, include
`AllKnownAdvanpix` as one of the inputs:

```sh
./build/psmm \
    -merge=AllKnownAdvanpix,results_100.txt,results_200.txt \
    -output=AllKnownAdvanpix
```

All inputs are read first, then deduplicated by extended-precision Mahler
measure and sorted by degree. The `-output` file is rewritten with the
clean union. This is the standard way to consolidate results from parallel
searches at different degrees back into a single master file.

### Analyze

Load a results file and display statistics — minimum Mahler measures,
nearest pair, extremal values of NNZ, H, L, K, U, Q, R, and any
non-primitive polynomials:

```sh
./build/psmm -analyze=AllKnownAdvanpix
```

## Data files

| File | Description |
|---|---|
| `AllKnownAdvanpix` | Extended set of known polynomials with M(p) < 1.3, including all entries from `Known180` plus ~40k new finds at degrees > 180 from this project. 72-digit Mahler measures. |
| `Known180` | Michael Mossinghoff's historical list through degree 180 ([source](http://www.cecm.sfu.ca/~mjm/Lehmer/lists/)). Legacy format with 13-digit precision. |

**File format** (AllKnownAdvanpix):

```
N M NNZ H L K U Q R c_0 c_1 ... c_{N/2}
```

| Field | Meaning |
|---|---|
| N | Polynomial degree |
| M | Mahler measure (72 decimal digits) |
| NNZ | Non-zero half-coefficients (excluding the leading 1) |
| H | Height (max coefficient magnitude) |
| L | Length (sum of coefficient magnitudes) |
| K | Roots outside the unit circle |
| U | Roots on the unit circle (complex conjugate pairs) |
| Q | Complex non-unity roots (quadruplets) |
| R | Real non-unity roots (pairs) |
| c_0 ... c_{N/2} | Half-coefficients (c_0 = 1 always) |

## References

- D. H. Lehmer, *Factorization of certain cyclotomic functions*, Ann. of Math. **34** (1933), 461--479.
- M. J. Mossinghoff, *Polynomials with small Mahler measure*, Math. Comp. **67** (1998), 1697--1715.
- R. Breusch, *On the distribution of the roots of a polynomial with integral coefficients*, Proc. Amer. Math. Soc. **2** (1951), 939--941.
- C. J. Smyth, *On the product of the conjugates outside the unit circle of an algebraic integer*, Bull. London Math. Soc. **3** (1971), 169--175.
- MPSolve: D. A. Bini and G. Fiorentino, *Design, analysis, and implementation of a multiprecision polynomial rootfinder*, Numer. Algorithms **23** (2000), 127--173.

## License & citation

This repository is dual-licensed:

- **Source code** (everything under the repository root except the items
  listed below): [GPL-3.0-or-later](LICENSE). Each source file carries
  a license header indicating this.
- **Polynomial database** (`AllKnownAdvanpix`, `Known180`, files under
  `tests/data/`) and **figures** (`images/*`): [Creative Commons
  Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
  See [LICENSE-DATA](LICENSE-DATA) for the full terms.

If you use any of these resources in published work, please cite:

> Pavel Holoborodko, *PSMM: Polynomials with Small Mahler Measure*,
> Advanpix LLC, 2020–2026. <https://github.com/advanpix/PSMM>

A machine-readable citation entry is provided in [`CITATION.cff`](CITATION.cff)
(GitHub's "Cite this repository" widget reads this file).
