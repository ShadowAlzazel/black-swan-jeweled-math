# Weights from idempotents: a check of the Gelfand–Graev weights

ShadowAlzazel shared two results from §14 of a text (source not identified here). This note checks them numerically on GL₂(F_q) and records one correction. Tags as elsewhere: **proven**, **open**, **observed**. The script is [gelfand_graev_check.py](gelfand_graev_check.py) (`python docs/gelfand_graev_check.py`, about a second, standard library plus NumPy).

## The two results

Notation: Γ_N = GL_N(F_q) (the text writes g for q). Labels c are characteristic polynomials of the semisimple parameter s of a Deligne–Lusztig series, and E_c is the central idempotent of the series with label c. For an idempotent e, the weight of a label is wt(e, c) = (e E_c)(1) / e(1).

* **Lemma 14.33.** If e = e* = e² and e(1) > 0, then wt(e, c) ≥ 0 for every label and Σ_c wt(e, c) = 1. The proof writes (e E_c)(1) as a squared norm, ‖λ(e)λ(E_c)δ₁‖², in the left regular representation.
* **Proposition 14.34.** For the Gelfand–Graev idempotent e_U = |U|⁻¹ Σ_u conj(ψ₀(Σ_i u_{i,i+1})) u, and every squarefree label c with irreducible factors of degrees d₁, …, d_j, wt(e_U, c) = sq(c) = 1/Π_i (q^{d_i} − 1). (The last line of the proof shows sq(c) is that reciprocal.)

## What was checked (**observed**, q = 2, 3, 5, 7)

For each q the script builds GL₂(F_q), U, ψ₀(x) = e^{2πix/q}, and the Gelfand–Graev character Φ = Ind_U^G ψ, evaluated from the definition of induction on one element per conjugacy class. It also builds the full classical character table: θ∘det, Steinberg twists, principal series and cuspidal characters. Every item below passed for all four q.

* **The character table:** complete (Σ dim² = |G|) and orthonormal (|Gram − I| ≤ 3 × 10⁻¹⁶).
* **Labels:** each irreducible's central character matches the constant term of its label under the fixed duality.
* **Φ:** Φ(1) = [G:U], the multiplicities are integers, and Φ is recovered exactly as the sum of its constituents.
* **Multiplicity one:** the constituents are exactly the generic (non-one-dimensional) irreducibles, and each label has exactly one constituent.
* **The proof's main identity:** (e_U E_c)(1) equals ‖λ(e_U)λ(E_c)δ₁‖², computed as an explicit convolution on G (deviation ≤ 3 × 10⁻¹⁷). It also equals the Frobenius-reciprocity weight, Σ m_π dim π / [G:U].
* **Lemma 14.33:** every weight is ≥ 0 and the weights sum to exactly 1 (exact rational arithmetic).
* **Proposition 14.34:** on squarefree labels, wt = 1/(q − 1)² for split c and 1/(q² − 1) for irreducible quadratic c, exactly as stated.

Weights for q = 7 (42 labels):

| Label type | Count | Gelfand–Graev constituent | Weight |
| --- | ---: | --- | ---: |
| (X − a)(X − b), a ≠ b | 15 | principal series, degree q + 1 | 1/36 = 1/(q − 1)² |
| irreducible quadratic | 21 | cuspidal, degree q − 1 | 1/48 = 1/(q² − 1) |
| (X − a)² | 6 | Steinberg twist, degree q | 7/288 = q/((q − 1)²(q + 1)) |

## One correction

The proof says each series contains exactly one Gelfand–Graev constituent χ_{s,reg}, with χ_{s,reg}(1) = [Γ_N : C(s)]_{g′}. Stated for every s, that degree formula is **false**. It holds for regular semisimple s, which is exactly the squarefree labels the proposition uses, so **the proposition and its conclusion stand**.

For central s on GL₂, C(s) = G and the formula would give degree 1. The actual constituent is the Steinberg twist, of degree q. Using the stated formula for every label would make the weights sum to q/(q + 1) instead of 1: 2/3, 3/4, 5/6 and 7/8 for q = 2, 3, 5, 7, as the script reports.

The general formula is the standard one, since C(s) is connected for GL_N (**proven**, Deligne–Lusztig theory):

  χ_{s,reg}(1) = [Γ_N : C(s)]_{g′} · |C(s)|_g,

the second factor being the degree of the Steinberg character of C(s). It agrees with the table for every label type of GL₂.

A counting check for larger N needs no character table. For s with characteristic polynomial Π f_i^{m_i}, deg f_i = d_i, C(s) = Π_i GL_{m_i}(F_{q^{d_i}}), and the corrected weight is |C(s)|_g / |C(s)|_{g′}. Summed exactly over every monic polynomial of degree N with nonzero constant term (**observed**):

| q | N = 2 | N = 3 | N = 4 |
| ---: | --- | --- | --- |
| 2 | 1 (uncorrected 2/3) | 1 (2/3) | 1 (28/45) |
| 3 | 1 (3/4) | 1 (3/4) | 1 (117/160) |
| 5 | 1 (5/6) | 1 (5/6) | 1 (775/936) |
| 7 | 1 (7/8) | 1 (7/8) | 1 (2793/3200) |

The corrected weights always sum to 1, as multiplicity one and dim Φ = [G:U] require. For N = 2 and q ≤ 3 they match the character-table weights label by label.

**Every N** (**proven**). Group the labels by their factorizations. With Q = q^d, a factor f^m with deg f = d has corrected weight

  Q^{m(m−1)/2} / Π_{k=1}^{m} (Q^k − 1) = Q^{−m} / (Q^{−1}; Q^{−1})_m.

Euler's identity Σ_m z^m/(y; y)_m = 1/(z; y)_∞ turns each irreducible f ≠ X into the factor Π_{j≥1} 1/(1 − (t q^{−j})^d). Monic polynomials with nonzero constant term have generating function Π_{f≠X} 1/(1 − u^{deg f}) = (1 − u)/(1 − qu). So the generating function of the weight sums is

  Π_{j≥1} (1 − t q^{−j}) / (1 − t q^{1−j}) = 1/(1 − t).

The product telescopes, and every coefficient is 1. The script also expands this series to N = 16 for each q without listing labels, and every sum is exactly 1.

**What the proof's formula loses** (**proven**; the numbers are **observed**). Dropping the Steinberg factor swaps Euler's identity for its other form, Σ_m y^{m(m−1)/2} z^m/(y; y)_m = (−z; y)_∞. The generating function becomes H(t²)/(1 − t), with H(s) = (s/q; q^{−2})_∞ / (s/q²; q^{−2})_∞.

* **The sums repeat in pairs.** N = 2K and N = 2K + 1 give the same value. For q = 2 the sums run 2/3, 2/3, 28/45, 28/45, 248/405, ….
* **They tend to H(1).** Euler's theorem that partitions into distinct parts and into odd parts are equinumerous gives H(1) = 1/((p; p)_∞ (−p; p)_∞²) with p = 1/q. That is 0.6091, 0.7290, 0.8278 and 0.8728 for q = 2, 3, 5 and 7.

So the stated degree formula loses a fixed share of the mass at every size, about 39% for q = 2.

## The bridge to Black Swan

* **Positivity for free** (**proven**). A self-adjoint idempotent turns any partition of unity by central projections into a probability distribution, because each weight is a squared norm. That is the same source of positivity as Black Swan's boost 2(1 + cos a)² and the discrete-log paper's filter Π(X − 1 − jM)². Each is a square, so non-negativity needs no tuning.
* **Experiment 05's energy distribution is such a weight** (**proven**, elementary). On ℓ²({0,1}ⁿ), the uniform unit vector u = 2^{−n/2} Σ_x δ_x is, up to normalization, the image of δ under the averaging idempotent of the group (Z/2)ⁿ. The energy-level projections P_k (multiplication by the indicator of E(x) = k) are commuting self-adjoint idempotents summing to the identity, and ‖P_k u‖² = N_k/2ⁿ. So N_k/2ⁿ ≥ 0 sums to 1 for the same reason as the Gelfand–Graev weights. The projections here live in the multiplication algebra, not the group algebra; the shared structure is "a positive state and a partition of unity by projections".
* **Fourier levels of a landscape** (**proven**, elementary). Group the characters χ_S of (Z/2)ⁿ by level |S| = k. The level projections are again commuting self-adjoint idempotents summing to the identity, so ‖P_k E‖²/‖E‖² is a probability distribution on levels: the same mechanism a third time.
  * **The look-alike lives on level 3.** For an XOR equation in its exact CNF encoding, the number of violated clauses is (1 − (−1)^r χ_abc)/2. So the look-alike's energy has weight on levels 0 and 3 only, while a 3-SAT clause also has weight on levels 1 and 2.
  * **Why the starts failed on XOR.** A conditional expectation reads exactly the coefficients Ê(S) whose sets S lie inside the fixed variables. On XOR every such coefficient apart from the constant is zero until some equation's last variable is fixed. That explains Experiment 05's finding that conditional-expectation starts give no help on XOR.
  * **Literature.** Ricci-Tersenghi, Weigt and Zecchina (Phys. Rev. E 63, 026702, 2001) introduced this model as "hyper-SAT", a perfectly balanced random 3-SAT. Complementing variables to move the planted solution to all zeros (Experiment 05's gauge rewrite) turns it into a diluted ferromagnet with three-spin couplings. Franz, Mézard, Ricci-Tersenghi, Weigt and Zecchina showed that this ferromagnet has a glass transition ("A ferromagnet with a glass transition", Europhys. Lett. 55, 465, 2001): local dynamics from a random start stay in the glass, though the ordered ground state is known. [ROADMAP R15](../ROADMAP.md) measures these level weights.
* **What this could offer** (**open**). Black Swan's clause weights were chosen by hand (golden phases, tile orientations) and only had to be positive. Idempotent-induced weights are canonical probability distributions with that positivity built in. Whether any such canonical weighting helps search is untested, and the controls in [CLAUDE.md](../CLAUDE.md) would apply first.
