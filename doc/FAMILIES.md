# PSMM — parametric families with Boyd-Lawton convergence

This document catalogues the parametric polynomial families identified
in `AllKnownAdvanpix`. Each family is an infinite sequence whose Mahler
measure converges to a limit $L = \exp(m(F))$ given by the Boyd-Lawton
theorem applied to a bivariate companion $F(x, u)$.

For each family below:

- **Univariate form** — the family of polynomials $P_N$ as the
  parameters vary.
- **Bivariate companion** — the polynomial that, by monomial
  substitution $u = x^{N/2-c}$, recovers $P_N$.
- **Boyd-Lawton limit** — $L = \exp(m(F))$, the value the univariate
  Mahler measure approaches as the degree grows.
- **DB membership** — how many DB entries come from this family.

Family naming scheme: `C[a, d, s]` for cyclotomic-perturbation
families (parameters identify the bivariate companion), `S5[\Phi_d]`
for the 5-term reciprocal Salem family.

---

## C[2, 3, −] — 7-term cyclotomic perturbation

The original PSMM family. Historically called "Max-U" because its
members have the highest fraction of roots on the unit circle in the
DB.

**Univariate form** (parameter: odd $m \ge 5$; $N = m + 1$ even):

$$P_N(x) = (x + 1)(x^{N-1} + 1) - x^{(N-2)/2}\,\Phi_3(x).$$

**Bivariate companion** (substitution $u = x^{(N-2)/2}$):

$$F_{2,3}(x, u) = x(x + 1)\,u^2 - \Phi_3(x)\,u + (x + 1).$$

**Boyd-Lawton limit:**

$$L_{(2,3)} = \exp(m(F_{2,3})) = 1.255434037727251\ldots$$

**DB membership.** Every odd-$m$ value in $[5, 1294]$ realised. The
irreducible $P_N$ contributes one entry; reducible $P_N$'s
non-cyclotomic factor contributes one entry at the factor's degree
(typically still close to $N$). Distribution within the family:
roughly 13 of every 15 roots sit on the unit circle.

**Twin.** `C[2, 3, +]` is the same family with the sign of the
$\Phi_3$ perturbation flipped. By Mahler-measure sign-invariance, it
shares the same $L$. Both signs are merged into the DB.

---

## C[2, 12, −]

**Univariate form** (odd $m \ge 5$; $N = m + 1$ even):

$$P_N(x) = (x + 1)(x^{N-1} + 1) - x^{(N-4)/2}\,\Phi_{12}(x).$$

**Bivariate companion:**

$$F_{2,12}(x, u) = x(x + 1)\,u^2 - \Phi_{12}(x)\,u + (x + 1)\,x^{\phi(12) - \phi(2)},$$

with $\phi(12) - \phi(2) = 3$.

**Boyd-Lawton limit:**

$$L_{(2,12)} = 1.30911\ldots$$

**DB membership.** Above 1.3 — the bulk of the family is included in
`AllKnownAdvanpix` opportunistically (494 entries; see
[SCOPE.md](SCOPE.md) for the above-1.3 cataloguing policy). All 16
low-$N$ members with $M < 1.3$ are also in the DB.

---

## S5[Φ_3] — 5-term reciprocal Salem (Φ_3 perturbation)

Covers the sign combinations $(s_1, s_2) \in \{(-, -),\ (+, +)\}$.

**Univariate form** (parameters: even $N \ge 8$, integer $a$ with
$1 \le a < N/2$, signs $s_1, s_2 \in \{-1, +1\}$):

$$P_{N, a, s_1, s_2}(x) = 1 + s_1 x^a + s_2 x^{N/2} + s_1 x^{N-a} + x^N.$$

For $(s_1, s_2) = (-, -)$ explicitly:

$$P(x) = x^N - x^{N-a} - x^{N/2} - x^a + 1.$$

**Bivariate companion** (substitution $u = x^{N/2 - a}$, then
$y = x^a$):

$$G_{\Phi_3}(y, u) = y^2 u^2 - y\,\Phi_3(u) + 1.$$

The torus-measure-preserving substitution makes $m(F_a) =
m(G_{\Phi_3})$ for every $a \ge 1$ — the Boyd-Lawton limit is
**independent of $a$**. The same $L$ also covers the twin $(+, +)$
via $y \to -y$.

**Boyd-Lawton limit:**

$$L_{S5[\Phi_3]} = \exp(m(G_{\Phi_3})) = 1.285714475811306572\ldots$$

---

## S5[Φ_6] — 5-term reciprocal Salem (Φ_6 perturbation)

Covers $(s_1, s_2) \in \{(-, +),\ (+, -)\}$.

**Univariate form** — same as S5[Φ_3] with the other sign pair.

**Bivariate companion:**

$$G_{\Phi_6}(y, u) = y^2 u^2 - y\,\Phi_6(u) + 1,$$

where $\Phi_6(u) = u^2 - u + 1$.

**Boyd-Lawton limit:**

$$L_{S5[\Phi_6]} = 1.285741055747102200\ldots$$

This sits only $2.66 \times 10^{-5}$ above $L_{S5[\Phi_3]}$ — the two
S5 sub-families are indistinguishable on a coarse DB-wide $M$-vs-$N$
plot, jointly appearing as a single band at $M \approx 1.2857$.

---

## Combined S5 contribution to the DB

Across all 4 sign combinations and the high-$N$ factor-extract
backfill, the S5 family accounts for the bulk of the DB's
$M \approx 1.2857$ band, including reducible-case dense factor extracts
at $N > 456$ (catalogued by `tools/factor_5term_reducibles.py`).

---

## Provenance

- The C[2, 3, −] family is the original PSMM finding; first surfaced
  as `N = 456` brute-force entries, then extended analytically.
- The C[2, 12, −] family was identified as the next-smallest $L$
  predicted by the $(a, d)$ cyclotomic-perturbation table; it served
  as the empirical check that $L > 1.3$ families are also catalogued.
- The S5 family was identified on 2026-05-18 from a DB-wide
  $M$-vs-$N$ scatter that revealed a second accumulation band at
  $M \approx 1.286$, distinct from the known $L = 1.2554$ locus.
  Structural analysis of the entries near $M = 1.286$
  (NNZ = 2, H = 1, L = 5) recovered the 5-term form. The bivariate
  companion follows by direct substitution and the Jensen-Lemma
  integral gives the exact $L$ values above.
