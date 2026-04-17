[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# PSMM — Polynomials with Small Mahler Measure

Exhaustive search for primitive, irreducible, reciprocal integer polynomials
whose [Mahler measure](https://en.wikipedia.org/wiki/Mahler_measure) is below
a given threshold (typically M(p) < 1.3). The smallest known Mahler measure
greater than 1 is Lehmer's number (M &asymp; 1.17628), achieved by a degree-10
polynomial found in 1933. Whether smaller values exist remains one of the
central open problems in number theory.

## How it works

### Search space

A monic reciprocal polynomial of degree N has the form

    p(x) = x^N + a_1 x^{N-1} + ... + a_{N/2} x^{N/2} + ... + a_1 x + 1

so only the N/2 "half-coefficients" a_1, ..., a_{N/2} are free. PSMM
enumerates all such polynomials with:

- **integer coefficients** drawn from a given alphabet (e.g. {-1, 1}),
- exactly **nnz** non-zero half-coefficients (sparsity constraint).

The enumeration is driven by `reciprocal_polynomials_iterator`, which
iterates over all sparsity patterns (via `std::next_permutation`) and all
coefficient combinations (mixed-radix counter). Polynomials whose roots
are identical by symmetry (e.g. p(-x) vs p(x)) are skipped automatically.

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
    -known=AllKnownParallel \
    -addto=AllKnownParallel
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
# Skip polynomials already in AllKnownParallel; append new finds to it.
./build/psmm \
    -degree=100 \
    -coeffs=-1,1 \
    -nnz=1,2,3 \
    -threshold=1.3 \
    -threads=8 \
    -period=3600 \
    -known=AllKnownParallel \
    -addto=AllKnownParallel \
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
    -known=AllKnownParallel \
    -addto=AllKnownParallel
```

### Merge

Combine multiple result files, remove duplicates, and produce a single
sorted output. To fold new search results into the master file, include
`AllKnownParallel` as one of the inputs:

```sh
./build/psmm \
    -merge=AllKnownParallel,results_100.txt,results_200.txt \
    -output=AllKnownParallel
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
./build/psmm -analyze=AllKnownParallel
```

## Data files

| File | Description |
|---|---|
| `AllKnownParallel` | Extended set of known polynomials with M(p) < 1.3, found by Advanpix. 72-digit Mahler measures. ~49k entries. |
| `Known180` | Michael Mossinghoff's historical list through degree 180 ([source](http://www.cecm.sfu.ca/~mjm/Lehmer/lists/)). Legacy format with 13-digit precision. |

**File format** (AllKnownParallel):

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

## License

GPL v3 or later. See [LICENSE](LICENSE).
