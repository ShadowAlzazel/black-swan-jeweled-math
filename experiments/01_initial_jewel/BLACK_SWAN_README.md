# Black Swan: SAT jewels, quasicrystal geometry, and zeta weights

This is a reproducible mathematical playground. RH is assumed for the global interpretation of the zeta spectrum. The finite numerical experiments do not require RH, make no claim about P versus NP, and do not establish an advantage for number-theoretic weights.

## Files and reproduction

- `black_swan_experiments.py`: generate formulas, compute weights, compare search methods, find illustrative examples, and validate invariants.
- `black_swan_results.json`: full configuration, formulas, starting states, weights, all 2,304 run records, results, and selected examples.
- `quasicrystal_geometry.py`: generate and validate the octagonal cut-and-project geometry.
- `jewel_geometry.json`: finite lattice candidates, window, accepted patch, and 16 Boolean projections.
- `BLACK_SWAN_README.md`: this guide.

Dependencies: Python 3.12, NumPy, and mpmath. The recorded versions are in the results. From a directory containing the scripts, run:

```sh
python black_swan_experiments.py
python quasicrystal_geometry.py
```

Each script writes its JSON next to itself. Seed: 20261125. No downloaded data or network connection is needed.

## 1. Exact energy and spin transformations

For a Boolean assignment x, let u_j(x) be 1 when clause j is violated, and 0 otherwise. Define

E_w(x) = sum_j w_j u_j(x).

All weights in this experiment satisfy 0.4 <= w_j <= 1.6. Therefore

0.4 E(x) <= E_w(x) <= 1.6 E(x),

where E is the unweighted number of violated clauses. In particular E_w(x)=0 if and only if x satisfies the original formula. No satisfying assignment is added or removed. An unsatisfying Boolean assignment has weighted energy at least 0.4. Search success is checked with exact Boolean clauses, never a floating-point zero test.

This also extends the earlier complex-variable jewel: multiplying each squared clause residual by a strictly positive real weight preserves the complete zero set, including the Boolean-forcing residuals z_i(z_i-1).

For spins s_i=2x_i-1 in {-1,+1}, and literal sign sigma=+1 for a positive literal and -1 for a negative literal,

u_C(s) = (1/8) product_(i in C) (1-sigma_(C,i) s_i).

This is an exact conversion from 3SAT to an energy with interactions of up to three spins. Its ground energy is zero exactly when the formula is satisfiable. It is not generally a pairwise Ising model without auxiliary variables.

An optional exact pairwise conversion is available. For three Boolean failure bits q1,q2,q3 and an auxiliary Boolean y, replace w*q1*q2*q3 by

w*y*q3 + M*(q1*q2 - 2*q1*y - 2*q2*y + 3*y).

For M>w, minimizing over y recovers exactly w*q1*q2*q3. In our weight range, M=2 is sufficient. The nonnegative penalty is zero precisely at y=q1*q2; violating that equality costs at least M and cannot gain more than w in the first term. This auxiliary-variable conversion is described here for further experimentation; it is not used in the benchmark.

## 2. Our trial weighting rules

Let phi=(1+sqrt(5))/2, and let p_j be prime j (starting at j=1).

Golden-ratio rule:

w_j = 1 + 0.6 cos(2*pi*j/phi + theta).

This is a quasiperiodic modulation inspired by the irrational rotations that occur in quasicrystal mathematics. It is not itself a cut-and-project point set.

Zeta rule:

w_j = 1 + (0.6/K) sum_(k=1)^K cos(gamma_k log(p_j) + theta_k), K=10.

Here gamma_k is the imaginary part of zeta zero k, computed using mpmath at 40 decimal digits and then converted to binary64 for the search. The formula is an experimental choice, not an identity or consequence of RH. Under RH every nontrivial zero has real part 1/2; the first ten ordinates can already be computed without assuming the full hypothesis. RH alone implies neither quasicrystallinity of the zeros nor a SAT speedup.

Independent phases are sampled for each formula. Each structured weight vector is also randomly permuted as a control. The shuffled vector preserves the exact multiset, its amplitude distribution, and its mean while destroying the chosen weight-to-clause order. Clause order is shuffled separately for every generated instance. Clause numbering is arbitrary and is not an intrinsic mathematical property of SAT.

## 3. Search methods

Strict greedy examines all n one-bit neighbors and accepts a neighbor only when its fixed weighted energy strictly decreases (numerical comparison tolerance 1e-12). Ties among best candidates use a seeded random choice. At a plateau or local minimum it stops. This intentionally simple comparator has no sideways moves or restarts.

The WalkSAT baseline samples an unsatisfied clause, evaluates the flips of its three variables, and chooses a zero-break flip if one exists. Otherwise it chooses randomly with probability 1/2, or chooses a minimum-break flip with probability 1/2. A break is a previously satisfied clause becoming unsatisfied. Ties are random. This method can cross plateaus and climb uphill.

Each method has a cap of 1,025 full-assignment evaluations per run. Every such evaluation scans every clause, including the initial assignment. Strict greedy consumes n candidate evaluations per decision; WalkSAT consumes three. This is a common evaluation cap, not equal consumed work: greedy often stops early, and unused evaluations are not recycled into restarts. Zeta/weight generation overhead is not charged to this search-only cap. These results are not wall-clock speed measurements.

## 4. Recorded experiment

48 planted satisfiable formulas; 16 variables and 68 distinct clauses per formula; eight shared random starting assignments per formula; 384 runs per method. Signs and distinct variable triples are sampled conditional on a planted satisfying assignment. The planted assignment is used only in generation and is never passed to search.

| Method | Solved / runs | Success | Mean assignment evaluations |
| --- | ---: | ---: | ---: |
| Strict uniform greedy | 87 / 384 | 22.66% | 70.88 |
| Golden-ratio weights | 102 / 384 | 26.56% | 82.67 |
| Shuffled golden weights | 117 / 384 | 30.47% | 84.29 |
| Zeta weights | 111 / 384 | 28.91% | 82.46 |
| Shuffled zeta weights | 119 / 384 | 30.99% | 83.42 |
| WalkSAT | 382 / 384 | 99.48% | 97.45 |

Reweighting sometimes improves a weak strict-greedy search by resolving plateaus, and sometimes creates new traps. The shuffled controls did better than the structured versions in this batch. Ordinary WalkSAT had the highest success rate. This provides no evidence of a special benefit from RH or the golden ratio. No confidence intervals or significance claims are made; starts on the same formula are correlated, and this small planted distribution is not representative of all SAT instances. This is not a comparison with modern complete SAT solvers.

Failure to find a satisfying assignment is not an UNSAT certificate. A complete solver such as exhaustive search, DPLL, or CDCL would be needed to establish unsatisfiability in general.

## 5. Selected four-variable examples

Examples were deliberately selected to illustrate both directions; they are not independent benchmark evidence. An integer state uses bit j for variable x_(j+1). Text labels list x1,x2,x3,x4, so integer 5 is labeled 1010.

- `golden_wins`: uniform stops at state 6 with one violation. Golden weights follow 6 -> 2 -> 0. Unweighted violations are 1 -> 1 -> 0: the first move crosses a plateau.
- `zeta_wins`: uniform stops at state 5. Zeta weights follow 5 -> 1 -> 0, again with violations 1 -> 1 -> 0.
- `golden_loses`: uniform follows 13 -> 9 -> 1 to a solution; golden weights stop immediately at 13 with two violations.
- `zeta_loses`: uniform follows 8 -> 0 -> 1 to a solution; zeta weights follow 8 -> 12 and stop with one violation.

Full formulas, phases, weights, and paths are included in the JSON. All 16 states are checked exhaustively for each selected example.

## 6. The quasicrystal construction

Let r=1/sqrt(2). Our physical projection has columns

P = [(1,0), (r,r), (0,1), (-r,r)],

and our internal projection has columns

Q = [(1,0), (-r,r), (0,-1), (r,r)].

For k in Z^4, use physical coordinate Pk and internal coordinate Qk. Their combined coordinates recover the original lattice point via

k = (1/2)*(P^T(Pk) + Q^T(Qk)).

The acceptance window is the convex hull of Q applied to the 16 corners of [-1/2,1/2]^4, a regular octagon, shifted by (0.071,0.113). Keep k precisely when Qk lies in this window, then display Pk. A finite physical-radius crop produces the supplied patch. The enumeration bounds are checked in the geometry script.

There are 182 candidate sites, 427 candidate unit lattice edges, and 81 accepted sites/144 edges at the initial shift and radius 4.5. Candidate padding permits shifting the window horizontally by up to 0.35 in either direction. This changes accepted sites through phason changes. The generated model set is an octagonal quasicrystal construction; a finite image alone is not a proof of aperiodicity.

The SAT jewel uses the same P to display the 16 points x in {0,1}^4. These 16 projections are distinct and retain their Boolean labels and original 32 one-bit edges. Every assignment is retained regardless of acceptance-window membership. This finite projected cube is not itself an infinite quasicrystal, and its drawing is not used to choose the benchmark moves. Window filtering cannot safely discard arbitrary SAT assignments without an additional correctness argument.

## 7. Promising next experiments

1. Adaptive Breakout: at a trap, increase the weights of currently violated clauses, then resume descent. This is an established algorithm family, unlike our arbitrary fixed weight schedules.
2. Geometric separator search: use a layout to propose a small separator of the variable-interaction graph, verify that it really separates that graph, enumerate separator assignments, and solve disconnected pieces independently. The separator size, not the drawing alone, determines the potential savings.
3. Pairwise Ising/QUBO conversion: implement the auxiliary-variable construction above and verify its minima exhaustively on small instances.
4. Stronger controls: multiple clause permutations per formula, random frequencies matched to the zeta frequency distribution, tuned phases assessed on held-out formulas, and complete SAT baselines. Pre-register larger comparisons before interpreting small numerical differences.

## Primary references

- Selman, Kautz, Cohen (1994), *Noise Strategies for Improving Local Search*: https://cdn.aaai.org/AAAI/1994/AAAI94-051.pdf
- Morris (1993), *The Breakout Method for Escaping from Local Minima*: https://cse.unl.edu/~choueiry/Documents/BreakoutMorris1993.pdf
- Grimm (2015), *Aperiodic crystals and beyond*: https://journals.iucr.org/b/issues/2015/03/00/dq5011/index.html
- Clay Mathematics Institute, *Riemann Hypothesis*: https://www.claymath.org/millennium/riemann-hypothesis/

## Validation

The experiment compares the vectorized evaluator with an independent direct Boolean evaluator, checks positive weight bounds and exact zero sets, checks one-bit paths and weighted descent, verifies every reported success independently, and enforces evaluation caps. The geometry checks inverse reconstruction, window membership, complete enumeration bounds, edge lengths, and injectivity of all Boolean projections. Instance zero and all its method/start outcomes were regenerated exactly. The pairwise penalty identity was separately enumerated over all Boolean triples and auxiliary choices.
