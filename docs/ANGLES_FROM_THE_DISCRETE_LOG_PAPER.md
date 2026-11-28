# Angles from the discrete-log paper

In November 2026 we read *A Deterministic Subexponential-time Algorithm for the Discrete Logarithm Problem* (Liu, Zhang and Ren, dated 23 November 2026). This note records what Black Swan can use from it. The experiments are in [experiments/05_counting_jewels](../experiments/05_counting_jewels/README.md), and the companion entries for the Fang & Fox Lemma Office are in that repository's IDEAS.md.

## How to read the paper

* It is a conditional result. Its only unproved premise is a Boolean tensor-network contraction theorem from an October 2026 OpenAI preprint, which it cites as refuting the Exponential Time Hypothesis. A Claude (Anthropic) preprint with a Lean formalization is cited for a 2^o(n) 3-SAT algorithm (**open** for us: these preprints postdate what we could check).
* Its own reference list includes S. Szeider's note "OpenAI's claimed refutation of ETH: exhaustive search below 2^131071", whose title says the claimed algorithm does not beat brute force below that size.
* Even if everything holds, 2^o(n) is not polynomial, so P versus NP is untouched. At 256-bit size the exponent improves by a factor of about √(log₂ log₂ 256) ≈ 1.7 before unknown constants (**proven** arithmetic on the stated bound).
* The paper says all of its lemmas pass Lean, yet its main theorem assumes an outside result. A Lean check of a conditional theorem certifies only the implication.

The value for us is the toolkit, which is classical and checkable on its own.

## Seven angles for Black Swan

1. **Count instead of search.** A formula is a tensor network: one factor per clause, 1 when satisfied and z when violated. Contracting it gives the energy polynomial Z(z) = Σ N_k z^k, where N_k counts assignments with k violated clauses and N_0 counts solutions. The paper computes in the same truncated ring (Z/2^P)[z]/(z^ℓ). Because Z is a polynomial, Cauchy's integral and the argument principle genuinely apply, which Experiment 04's field did not allow.
2. **Counting finds solutions.** Pin variables one at a time and keep a value whose completions still count above zero (self-reducibility, the paper's Proposition 2.6 for factoring).
3. **Deterministic starts by conditional expectations.** The paper derandomizes by fixing bits one at a time while keeping a conditional expectation small. In SAT this is Johnson's 1974 method, which guarantees at most m/8 violated three-literal clauses (**proven**; Håstad 2001 shows 7/8 is optimal unless P = NP). It is the principled replacement for irrational sweeps.
4. **Count-preserving (parsimonious) transformations.** The paper's gadgets preserve the number of witnesses exactly, a stronger invariant than Black Swan's "positive weights keep the solutions". Count-preserving rewrites reshape the landscape without changing the problem, which is the Black Swan question in its sharpest form.
5. **Positivity as a detector, and sums into products.** The paper's filter Π(X − 1 − jM)² ≥ 0 vanishes exactly at the allowed values, so one sum certifies a "for all". Black Swan's boost 2(1 + cos a)² and the Mertens identity 3 + 4 cos a + cos 2a ≥ 0 use the same move. Its binomial generating function (1+z)^(a+b) = (1+z)^a(1+z)^b turns a function of a sum into a product of local factors, the same rule as e^(a+b) = e^a e^b in the conformal jewels (`web/06`).
6. **Boundary compression is renormalization.** "Terminal-preserving contraction" sums out a region and keeps only its boundary wires. That is the move of separator-based dynamic programming, of the tensor renormalization group, and of Feigenbaum's doubling operator. Tilings could only help counting through such separators.
7. **The look-alike test.** The Lemma Office checks every RH-adjacent claim on Davenport–Heilbronn. SAT's look-alike is XOR-SAT: the same sparse random structure and a phase transition, but Gaussian elimination solves it, the same linear algebra as the paper's coset updates. A landscape trick that succeeds on planted 3-SAT but fails on planted 3-XOR-SAT is reading local hints, not structure.

A further remark: Kronecker substitution (packing exponents U^i z^j into W^(i + Kj)) is a lossless cut-and-project. Experiment 04's plane slice was lossy, and the two share one conservation law: an injective packing of 2^n states into one coordinate costs n·log K bits of precision.

## What the experiments found

Measured in [Experiment 05](../experiments/05_counting_jewels/README.md): exact counting reproduces Experiment 03's labels, Johnson starts versus random starts, count-preserving transformations under WalkSAT, the XOR look-alike, contraction widths by family, and whether the hat substitution hierarchy is a good contraction order.
