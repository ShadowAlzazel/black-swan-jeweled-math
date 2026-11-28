# Experiment 07 — Fourier levels

Experiment 05 found that every landscape method fails the XOR look-alike. [docs/IDEMPOTENT_WEIGHTS.md](../../docs/IDEMPOTENT_WEIGHTS.md) noticed why: in Walsh–Fourier terms, the look-alike's energy lives entirely on level 3. This experiment uses that lens to redo four older experiments:

* Experiment 05's search tests, the conditional-expectation starts and the look-alike, along a ladder that moves the hidden solution's hint one level at a time;
* Experiment 05's starts and rewrites, measured formula by formula on Experiment 03's formulas;
* Experiment 05's energy polynomial, done exactly for XOR;
* the clause weights of Experiments 02 and 03, tried on the look-alike.

Claims are tagged **proven**, **open** or **observed**, as elsewhere in the repository.

## The lens

Write χ_S(x) = (−1)^(Σ_{i∈S} x_i). Every energy has a Walsh expansion E(x) = Σ_S Ê(S) χ_S(x), and its level k is the part with |S| = k.

* **Clause coefficients** (**proven**). A clause is violated when every literal is false. The literal x_i is false with indicator (1 + χ_i)/2, and ¬x_i with (1 − χ_i)/2. So a clause's coefficients are ±2^(−|C|) on the subsets of its variables, and a formula's are exact sums over clauses, for any n, with no enumeration. A Walsh–Hadamard transform reproduces them exactly on every family tested (F0).
* **The hint at a solution** (**proven**). At a solution x*, E(x*) = 0, so Ê(∅) + Σ_k L_k = 0 with L_k = Σ_{|S|=k} Ê(S) χ_S(x*). The **hint shares** h_k = L_k / (−Ê(∅)) sum to 1. They say which levels carry the dip of the landscape at the solution. Averaging over all solutions gives the same split for formulas with many solutions.
* **One clause at a hidden assignment A** (**proven**). Let s_i = +1 if a literal is true under A and −1 if not. Then the clause's dip of 1/8 splits as (L₁, L₂, L₃) = (−Σ s_i, Σ s_i s_j, −s₁s₂s₃)/8.
* **Even levels cannot tell A from its complement** (**proven**). χ_S(¬x) = (−1)^|S| χ_S(x), so E(¬A) = Ê(∅)(1 + Σ_k (−1)^(k+1) h_k) exactly. For 3-CNF that is Ê(∅)(1 + h₁ − h₂ + h₃).
* **What local methods read** (**proven**).
  * The majority vote reads level 1.
  * A conditional expectation reads exactly the coefficients Ê(S) whose sets S lie inside the fixed variables, so Johnson's first steps see only level 1.
  * Experiment 04's phasor energy, Σ_C Π (1 − s cos θ_j)/2, is the Walsh expansion with each χ_i replaced by cos θ_i. Its gradient at the equator (all cos θ_i = 0) is therefore minus the level-1 vector, and its Hessian there is the level-2 matrix. On XOR both vanish.

### The hint ladder

Five families of formulas with a hidden solution A. Their expected shares follow from the clause rule above:

| Family | Clauses kept | h₁ | h₂ | h₃ | E(¬A)/Ê(∅) |
| --- | --- | ---: | ---: | ---: | ---: |
| planted | every clause A satisfies | 3/7 | 3/7 | 1/7 | 8/7 |
| golden-hidden | weight q^t for t true literals, q = 1/φ | 0 | 3/φ³ ≈ 0.708 | 2/φ⁴ ≈ 0.292 | 4/φ⁴ ≈ 0.584 |
| deceptive | the same with q = 1/2 | −0.158 | 0.789 | 0.368 | 0.421 |
| two-hidden | clauses that A and ¬A both satisfy | 0 | 1 | 0 | 0 (¬A is a solution) |
| planted 3-XOR | the exact 4-clause encoding of each equation | 0 | 0 | 1 | 2 |

* **The golden ratio hides the solution** (**proven**). With clause weights q^t, the expected level-1 dip is proportional to Σ_t C(3,t) q^t (2t − 3) = 3q(q² + q − 1). It vanishes exactly at q = 1/φ.
  * Jia, Moore and Strain (2005) found this balance point (q ≈ 0.618 for 3-SAT) as the way to stop the literals pointing at the hidden solution. They present evidence that the resulting formulas are exponentially hard for WalkSAT.
  * Black Swan began by hoping golden weights would reveal solutions. In planted 3-SAT, the golden weighting is exactly the one that hides them from the vote.
* **Two families are exact for every formula, not just on average** (**proven**).
  * Each two-hidden clause has one or two true literals, so its split is (±1, −1, ∓1)/8. Hence h₂ = 1 and h₁ = −h₃ in every formula.
  * Each XOR clause has one or three true literals, so L₃ = −w/8 for any clause weight w. Hence h₃ = 1 and h₁ = −h₂ for every positive weighting.

## F2 — The hint ladder (Experiment 05's S1 and S3, redone)

Settings: 20 formulas per cell, α = 4.26 for the 3-SAT families, one equation per variable for XOR, and Experiment 05's WalkSAT with noise 0.5 and a budget of 300,000 flips. Cells show median flips and the number solved; a median above the budget means fewer than half solved.

WalkSAT from a random start:

| Family | n = 50 | n = 100 | n = 200 | n = 400 | n = 800 |
| --- | ---: | ---: | ---: | ---: | ---: |
| planted | 53 (20/20) | 176 (20/20) | 517 (20/20) | 1,930 (20/20) | 4,087 (20/20) |
| two-hidden | 128 (20/20) | 382 (20/20) | 3,296 (20/20) | 16,470 (20/20) | 72,162 (18/20) |
| deceptive | 250 (20/20) | 1,159 (19/20) | 4,955 (17/20) | 29,305 (13/20) | 177,002 (13/20) |
| golden-hidden | 218 (20/20) | 1,905 (20/20) | 3,434 (18/20) | 217,442 (11/20) | >300,000 (8/20) |
| planted 3-XOR | 5,844 (20/20) | >300,000 (9/20) | >300,000 (0/20) | >300,000 (0/20) | >300,000 (0/20) |

WalkSAT from Johnson's conditional-expectation start (index order):

| Family | n = 50 | n = 100 | n = 200 | n = 400 | n = 800 |
| --- | ---: | ---: | ---: | ---: | ---: |
| planted | 16 (20/20) | 50 (20/20) | 357 (20/20) | 868 (20/20) | 1,861 (20/20) |
| two-hidden | 64 (20/20) | 578 (20/20) | 4,143 (20/20) | 21,926 (20/20) | 79,034 (17/20) |
| deceptive | 218 (20/20) | 1,672 (19/20) | 3,656 (17/20) | 37,778 (13/20) | 82,686 (16/20) |
| golden-hidden | 364 (20/20) | 2,470 (20/20) | 16,278 (18/20) | >300,000 (10/20) | >300,000 (7/20) |
| planted 3-XOR | 5,757 (20/20) | >300,000 (7/20) | >300,000 (0/20) | >300,000 (0/20) | >300,000 (0/20) |

Agreement with A of the vote, of Johnson's start, and of the assignment where each WalkSAT run ended. The first four columns are means over all sizes; runs from both starts are pooled.

| Family | measured h₁, h₂, h₃ | vote | Johnson start | solved runs end, mean from n = 50 to 800 | unsolved runs end, middle 90% |
| --- | --- | ---: | ---: | --- | --- |
| planted | 0.432, 0.425, 0.143 | 0.694 | 0.690 | 0.75 → 0.72 | none unsolved |
| two-hidden | 0.003, 1.000, −0.003 | 0.525 * | 0.550 * | two-sided: 97 of 195 below 0.5, only 12 within 0.45–0.55 | 0.29–0.49 (5 runs) |
| deceptive | −0.148, 0.786, 0.362 | 0.435 | 0.415 | 0.52 → 0.39 | 0.33–0.47 |
| golden-hidden | 0.010, 0.705, 0.285 | 0.508 | 0.496 | 0.70 → 0.54 | 0.42–0.59 |
| planted 3-XOR | 0, 0, 1 | 0.500 | 0.504 | 0.94–0.96 (n ≤ 100) | 0.46–0.58 |

\* For two-hidden formulas, the better of the agreements with A and with ¬A.

What this shows (**observed**):

* **The measured shares match the predicted ones** at every size, and the vote follows the sign of h₁: 0.69 for planted, 0.44 for deceptive, and one half where h₁ = 0.
* **Johnson's start helps exactly when h₁ > 0.** It cuts the planted family's work 1.4–3.5× at every size. With h₁ = 0 or negative it gives no consistent help, and on golden-hidden formulas it is worse than a random start at every size.
* **WalkSAT goes where the low levels point.**
  * Planted runs end on A's side.
  * Deceptive runs end on the far side from n = 100 on, where q = 1/2 formulas at this density still have solutions.
  * Two-hidden runs leave the equator in both directions and find solutions on either side, as an even-level hint predicts.
  * Golden-hidden runs drift toward the equator as n grows: the mean agreement of their solutions falls from 0.70 to 0.54.
  * Unsolved XOR runs stay at the equator, while the solutions found at small n sit near A.
* **Difficulty orders by the lowest level that carries a usable hint, at every size.**
  * Planted formulas (level 1) are easiest.
  * Two-hidden formulas come next: level 2 is enough for them, because ¬A is also a solution.
  * Deceptive and golden-hidden formulas are harder. In both, level 1 points nowhere or the wrong way, and only level 3 favours A.
  * XOR, with level 3 alone, is hardest.
* **One prediction failed.** Golden-hidden formulas have a false valley at ¬A, at about 0.58 of the mean energy, but failing runs do not drift there. They stall near the equator, at agreement 0.42–0.59, with 5–31 violated clauses at n = 800.
* **Scaling.**
  * The planted family's work grows 2–4× per doubling of n, and the two-hidden family's 3–9×. Neither ratio grows with n, as it would for exponential growth. This agrees with Jia, Moore and Strain's remark that WalkSAT solves two-hidden formulas (Achlioptas, Jia and Moore 2004) in polynomial time.
  * Golden-hidden work jumps 63× between n = 200 and 400.

These are planted formulas, which leak hints by construction, so none of this bears on the tripwire in [CLAUDE.md](../../CLAUDE.md).

## F1 — Experiment 03's formulas and Experiment 05's rewrites

For each of Experiment 03's 89 satisfiable formulas (n = 12–32), every solution was enumerated (checked against exact counts) and the hint shares were averaged over all of them. WalkSAT ran 10 times from random starts and 10 times from the index-order conditional-expectation start.

| Family | Formulas | h₁ | h₂ | h₃ | Level-1 cosine | Median flips, random | Median flips, Johnson |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| planted | 32 | 0.444 | 0.423 | 0.134 | 0.62 | 17 | 3 |
| geometric | 32 | 0.432 | 0.431 | 0.137 | 0.61 | 19 | 2.5 |
| random | 25 | 0.298 | 0.529 | 0.173 | 0.44 | 55 | 14 |

* **Per-formula correlations** (**observed**).
  * Formulas with more level-1 hint need fewer WalkSAT flips (Spearman ρ = −0.56).
  * The size of Johnson's gain does not track h₁ from formula to formula (ρ = 0.02; ρ = 0.21 for the level-1 cosine). Every family here has a substantial level-1 hint, and Johnson's start nearly solves these small formulas, at 1.4 violated clauses on average.
  * The ladder above is the sharper test: there the gain vanishes when h₁ does.
* **Experiment 05's rewrites, measured the same way** (**observed**; mapped onto each formula's own solutions):

| Rewrite | Variables | h₁ | h₂ | h₃ | Level-1 cosine | Experiment 05 median flips |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| original | 21.8 | 0.399 | 0.456 | 0.146 | 0.57 | 29.5 |
| complemented variables | 21.8 | 0.399 | 0.456 | 0.146 | 0.57 | 21 |
| duplicated clauses | 21.8 | 0.397 | 0.458 | 0.145 | 0.53 | 28.5 |
| implied resolvents added | 21.8 | 0.512 | 0.416 | 0.072 | 0.59 | 17 |
| definitional extension | 113.1 | 0.509 | 0.427 | 0.064 | 0.66 | 392 |
| Tovey splitting | 273.9 | 0.057 | 0.922 | 0.021 | 0.17 | 7,866 (157/178 solved) |

  * **Complementing variables changes nothing** (**proven**). Each coefficient changes sign exactly when χ_S(x*) does.
  * **Tovey splitting moves the hint to level 2** (**proven**, then measured). Around each cycle of implications c₁ → c₂ → … → c₁, the level-1 terms cancel, leaving (d − Σ χ_{c_i} χ_{c_{i+1}})/4, a pure level-2 chain. Each copy keeps only its own clause's vote. The share of the hint on level 1 falls from 0.40 to 0.06 and the cosine from 0.57 to 0.17, in line with the 270-fold slowdown.
  * **Adding resolvents raises h₁ and helped search** (17 flips against 29.5).
  * **The definitional extension is a counterexample.** It raises h₁ and the cosine, through its strongly voted y variables, yet WalkSAT needed 13 times more flips. Level 1 alone does not predict the effect of rewrites that add variables.

## F3 — The XOR energy polynomial (Experiment 05's C4, redone)

For an XOR system with m equations on n variables, let K be the set of dependencies: subsets T of equations whose left sides add to zero. Then, exactly (**proven**):

  Σ_k N_k z^k = 2^(n−m) Σ_{T∈K} (−1)^(rhs(T)) (1 − z)^|T| (1 + z)^(m−|T|).

For one equation, z^u = (1 + z)/2 + (1 − z)/2 · (−1)^r χ_S. Multiplying over the equations and summing over x leaves only the dependencies. This is the high-temperature (loop) expansion of an Ising-type model (van der Waerden 1941) in the form of the MacWilliams identity (1963), with w = (1 − z)/(1 + z) = tanh β. That change of variable is a Möbius map, like the one added to `web/06`; it sends the unit disk to the right half-plane.

* **Corollaries** (**proven**).
  * N₀ = 2^(n − rank) if every dependency has even right side, and 0 otherwise: the same 0 = 1 certificate that Experiment 06's kernel finds.
  * A system of full rank has the energy histogram of m independent coins, 2^(n−m)(1 + z)^m. All its zeros sit at z = −1, so the polynomial cannot tell it from m unrelated equations.
  * In general, z = −1 is a zero of multiplicity m minus the size of the largest dependency, unless the signed count of those largest dependencies cancels.
* **Checks** (**observed**). On 32 systems (n = 20–60), the formula matches Experiment 05's exact contraction coefficient for coefficient. The 12 systems with one right side flipped are unsatisfiable and give N₀ = 0.

| Systems (4 each) | Dependencies | Zeros at z = −1 | Other zeros inside the unit disk |
| --- | --- | --- | --- |
| n = 20, m = 20 | 1–3 | 8–10 of 20 | 4–6 |
| n = 30, m = 30 | 1–3 | 12–28 of 30 | 0–8 |
| n = 40, m = 40 | 1–3 | 18–26 of 40 | 6–10 |
| n = 40, m = 20 and n = 60, m = 30 | none | all | none |

Random 3-SAT spreads its zeros across the disk (Experiment 05). The look-alike piles 40–93% of them onto the single point z = −1.

## F4 — Clause weights on the look-alike (Experiments 02–03, redone)

Black Swan's clause weights act on the Walsh coefficients, Ê_w(S) = Σ_j w_j û_j(S). Weights were tried on 2,000 planted XOR systems (n = m = 200) and 400 golden-hidden formulas (n = 200). Values after ± are standard errors.

| Formulas | Weights | Level 1 in the energy's variance | h₁ | h₂ | h₃ | Level-1 cosine with A | Vote agreement |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| XOR, n = m = 200 (2,000) | equal | 0.00 | 0.000 | 0.000 | 1.000 | 0.0000 | 0.5000 |
| XOR, n = m = 200 (2,000) | golden boost 1 + (1 + cos 2πjφ)² | 0.19 | 0.000 | 0.000 | 1.000 | +0.0003 ± 0.0016 | 0.4991 ± 0.0008 |
| XOR, n = m = 200 (2,000) | random in [0.5, 1.5] | 0.05 | 0.000 | 0.000 | 1.000 | +0.0002 ± 0.0016 | 0.5004 ± 0.0008 |
| XOR, n = m = 200 (2,000) | (1, φ, φ², φ³) within each equation | 0.41 | +0.003 | −0.003 | 1.000 | +0.0029 ± 0.0016 | 0.5011 ± 0.0007 |
| XOR, n = m = 200 (2,000) | oracle: weight 2 when every literal is true under A | 0.30 | +0.600 | −0.600 | 1.000 | +0.8679 ± 0.0002 | 0.9755 ± 0.0002 |
| golden-hidden, n = 200 (400) | equal | 0.43 | +0.001 | +0.708 | 0.291 | +0.0004 ± 0.0025 | 0.5009 ± 0.0013 |
| golden-hidden, n = 200 (400) | golden boost 1 + (1 + cos 2πjφ)² | 0.43 | +0.001 | +0.708 | 0.291 | +0.0009 ± 0.0026 | 0.5000 ± 0.0015 |
| golden-hidden, n = 200 (400) | random in [0.5, 1.5] | 0.43 | −0.001 | +0.709 | 0.292 | −0.0012 ± 0.0025 | 0.5020 ± 0.0015 |
| golden-hidden, n = 200 (400) | oracle: weight 2 when every literal is true under A | 0.44 | +0.204 | +0.457 | 0.339 | +0.2290 ± 0.0030 | 0.5872 ± 0.0015 |

* **No weighting moves the look-alike's dip off level 3** (**proven**): h₃ = 1 and h₁ = −h₂ for all positive weights (see the ladder above).
* **Solution-blind weights create level-1 variation but no level-1 hint** (**proven** in expectation, **observed**).
  * Weights fixed by clause position (golden boosts, random, a (1, φ, φ², φ³) pattern within each equation) give the energy a large level-1 part: up to 41% of its variance.
  * For one equation, averaging over the 8 local values of A gives an inner product of exactly 0 between that part and A, for any weights that depend only on the equation and clause positions. The same holds for golden-hidden clauses, whose expected Σ s_i is 0 whatever their weight.
  * The pattern weights showed a 2.4σ excess in cosine at 40 systems. At 2,000 it is +0.0029 ± 0.0016, with vote agreement 0.5011 ± 0.0007, both within 2σ of zero.
* **The measurement can see a hint when one exists.** An oracle that doubles the weight of clauses with every literal true under A gives cosine 0.87 on XOR, where the vote recovers 97.6% of A, and 0.23 on golden-hidden formulas. On XOR the oracle's level-1 inner product with A is exactly 3/8 per equation, which is 75 for every system. Its level-2 part cancels it, as the proven rule requires.

## What it means for Black Swan

Black Swan reshapes landscapes by reweighting and rewriting. In Fourier terms, that changes the coefficients Ê(S), and local methods consume the lowest levels. The vote and Johnson's start read level 1, the gradient and curvature at the equator read levels 1 and 2, and WalkSAT goes where those levels point.

* **The look-alike gives them nothing to read** (**proven**). Its whole dip is on level 3, no positive reweighting can move it, and the vote, Johnson's start and the equator gradient all see zero. That WalkSAT then fails is **observed**, here and in Experiment 05.
* **The golden ratio's true role here** (**proven**) is the opposite of Black Swan's first hope: it is the clause weighting that hides a planted solution from level 1.
* **Rewrites help when they add level-1 signal and hurt when they push the hint upward** (**observed**). Implied resolvents add signal; Tovey's copy cycles push it to level 2. The definitional extension shows the rule has exceptions.
* **Open:** whether a count-preserving rewrite can move signal downward on unplanted formulas at the threshold. That is the hint version of ROADMAP R1 and R15, and the question that would matter for the tripwire.

## Files

| File | Contents |
| --- | --- |
| `levels.py` | Exact Walsh coefficients of any CNF, level weights, hint shares averaged over solutions, level-1 cosine and vote, Krawtchouk sphere averages, brute-force checks, the hint-ladder generators, the XOR energy polynomial from dependencies, an all-solutions DPLL, and forward maps for Experiment 05's rewrites. |
| `run_levels.py` → `levels_results.json` | Sections F0 (checks), F1, F2, F3 and F4. Rerun one with `python run_levels.py f2`. |

## Run

```sh
python experiments/07_fourier_levels/run_levels.py
```

NumPy is the only dependency, together with Experiment 05's modules. Everything is seeded. F2 (on four cores) and F4 take about two minutes each; the other sections take seconds.
