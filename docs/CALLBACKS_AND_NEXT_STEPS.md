# Callbacks and next steps

Checking Experiment 04 turned up several connections to other fields. Each one below is either verified in this repository, with the script that checks it, or stated as a proposal with a prediction that a later experiment can test.

## Verified here

### 1. Optical vortices, not Cauchy zeros (wave physics)

Expanding the clause products shows that F_Φ is a finite sum of plane waves e^{i(K·r + Φ)} in r = (x, y), where each K is a sum of up to three variable wavevectors. Its zeros are phase singularities, the "wave dislocations" of Nye and Berry, carrying charge ±1, as in laser speckle and in random-wave models (Berry and Dennis). Like speckle, both signs occur in nearly equal numbers: 27 vortices and 31 antivortices across 40 formulas. A boundary winding number measures the net charge, by degree theory and Poincaré–Hopf; the argument principle does not apply. Every located zero's charge matched the sign of the Jacobian. Checked by `experiments/04_complex_spectre/audit_exp04.py`, sections 3–4.

### 2. Line arrangements and Cover's counting theorem (geometry, learning theory)

Rounding sgn cos θ_k(z) cuts a 2-D slice of the n-torus along L lines. L lines make at most 1 + L + L(L−1)/2 regions (Steiner); Zaslavsky's theorem handles general arrangements, and Cover's 1965 function-counting theorem is the same count read as perceptron capacity. A d-dimensional slice therefore sees O(L^d) assignments. Seeing all 2^n needs exponentially many lines, hence some exponentially large wavevector: about 1.2 × 10⁴ in the unit disk at n = 35, and 2.5 × 10¹³ at n = 100. Measured: the unit disk shows 128 of the 262,144 assignments of `missed_sat_instance_35.cnf` and none of its 9 solutions. This is the rigorous form of Experiment 04's "phase-frequency tax". Checked in audit section 5.

### 3. Feigenbaum's δ as an eigenvalue (renormalization)

δ is the only expanding eigenvalue of the linearized period-doubling operator at its fixed point g(x) = α g(g(x/α)). With 32 even terms, the eigenvalue equals δ to 45 digits and α = 1/g(1) to 48. The next eigenvalues are 0.15963 = α⁻², −0.12365, −0.05731 and 0.02548 = α⁻⁴. Superstable ratios through period 16,384 give only 8 digits: each period doubling divides the error by about δ but doubles the work, so the ratio method costs exponential work per digit, while the eigenvalue method converges exponentially in the number of terms. "Compose twice and rescale" does converge. The rescaled iterates of z² + c∞ approach g, from an error of 0.034 at level 0 to 1.7 × 10⁻⁹ at level 9. This is the same eigenvector-of-an-operator idea as the cycle eigenvectors in Fleurdelys. Checked by `experiments/04_complex_spectre/feigenbaum_renormalization.py`.

### 4. One ruler function behind three "frontier" directions

Let ν₂(j) be the 2-adic valuation, the ruler function 0, 1, 0, 2, 0, 1, 0, 3, …. Then:

* the kneading sequence at the Feigenbaum accumulation point, RLRRRLRLRLRRRLRR…, is ν₂(j) mod 2 (checked for 1,023 symbols from the period-1,024 superstable orbit);
* the binary reflected Gray code flips bit ν₂(j) at step j;
* Luby's restart sequence is the concatenation, over u = 1, 2, 3, …, of the runs 1, 2, …, 2^{ν₂(u)}. Knuth's "reluctant doubling" generates it with the test u & −u = 2^{ν₂(u)};
* the 2-adic absolute value is |j|₂ = 2^{−ν₂(j)}.

So FURTHER_EXPLORATIONS #3 (2-adic trees) and #8 (kneading restarts) share one object, and a "Feigenbaum kneading restart schedule" is already a relative of Luby's schedule. All three identities are checked in `feigenbaum_renormalization.py`, section 4.

### 5. The sign of α is Euler's half-turn

α < 0, so log α = log|α| + iπ. Under w = e^z, one level of period doubling (x ↦ αx) is the translation z ↦ z + log|α| + iπ: a ring step plus e^{iπ} = −1. Two levels make one complete wrap, α² = |α|²e^{2πi}. The superstable data show the same thing on the real line: the signs of d_k alternate − + − + …. In `web/06_conformal_jewels.html`, the α jewel has this twist built in and is invariant under w ↦ αw. Its step button reports f(z + log α)/f(z) = −2.502908.

### 6. Where SL(2, ℤ) belongs

Experiment 04 folds a point into the modular fundamental domain as a "relabeling" diagnostic, and the solver never uses the result. SL(2, ℤ) acts naturally in the elliptic direction: a change of basis of the period lattice leaves ℘ unchanged, and the lattice shape τ = ω₂/ω₁ lives in the modular fundamental domain. The explorer's ℘ satisfies ℘′² = 4℘³ − g₂℘ − g₃ to about 10⁻⁶. For periods 1 and i, it gives ℘(1/2) = 6.875176, against Γ(1/4)⁴/(8π) = 6.875186. It covers its image twice, as an even function should: the measured multiplicity is 2.00.

## Update: Experiment 05

[Experiment 05](../experiments/05_counting_jewels/README.md) took up the counting angle with methods from a November 2026 discrete-log paper ([angles and caveats](ANGLES_FROM_THE_DISCRETE_LOG_PAPER.md)). It covers:

* exact energy polynomials by contraction;
* their zeros and Cauchy's count of solutions (proposal E below, in [web/07](../web/07_energy_jewels.html));
* conditional-expectation starts;
* count-preserving rewrites;
* the XOR-SAT look-alike test.

## Proposals (not yet run)

### A. A field whose winding number does count zeros

Use complex frequencies, F(z) = Σ_j a_j e^{λ_j z}. Then F is entire, every zero has index +1, and the argument principle holds. By Pólya's theorem, the zeros lie in strips perpendicular to the sides of the convex hull of the conj(λ_j), with density proportional to side length. They track the corner locus of max_j (Re(λ_j z) + log|a_j|). That is the one-variable case of the tropical-amoeba direction (#1). This repairs the mathematics of Stage 3 but not its reach: the rounding still sees only a 2-D slice.

### B. Make the plane search exact

The O(L²) cells of the line arrangement can be enumerated exactly with a topological sweep. Crossing one line flips one variable, so each cell's violated-clause count updates incrementally. This turns "screen the plane" into an exact, polynomial-time search of a 2-D family of assignments. It is a fair test of whether golden, plastic or pinwheel wavevectors choose better slices than random ones. *Prediction:* no better than random 2-planes.

### C. Continuous-time dynamics instead of phasor flow

The one complex-plane ingredient that helped, on small formulas, was the gradient of the clause energy E. The principled relative is Ercsey-Ravasz and Toroczkai's continuous-time SAT dynamics. It has analog variables and clause weights that grow exponentially while a clause stays violated: Breakout's adaptive weights in continuous time, the mechanism Experiment 03 found useful. *Prediction:* WalkSAT-level success on Experiment 03's formulas, with solution times whose transient chaos grows near α ≈ 4.2.

### D. Kneading versus Luby restarts

Following item 4, compare Luby, the ruler sequence 2^{ν₂(j)} and a kneading-derived schedule on heavy-tailed instances at equal total budget. *Prediction:* no measurable difference, because all three are ruler-driven.

### E. An honest SAT jewel

Map the arrangement cells of one formula's slice (B) through e^z and color them by violated clauses. The result is a correct picture of what one slice can see, made as a picture and not as a solver.

## References

* J. F. Nye and M. V. Berry, "Dislocations in wave trains", *Proc. R. Soc. Lond. A* 336 (1974) 165–190.
* M. V. Berry and M. R. Dennis, "Phase singularities in isotropic random waves", *Proc. R. Soc. Lond. A* 456 (2000) 2059–2079.
* T. M. Cover, "Geometrical and statistical properties of systems of linear inequalities with applications in pattern recognition", *IEEE Trans. Electronic Computers* EC-14 (1965) 326–334.
* M. J. Feigenbaum, "Quantitative universality for a class of nonlinear transformations", *J. Stat. Phys.* 19 (1978) 25–52.
* O. E. Lanford III, "A computer-assisted proof of the Feigenbaum conjectures", *Bull. AMS* 6 (1982) 427–434.
* M. Luby, A. Sinclair and D. Zuckerman, "Optimal speedup of Las Vegas algorithms", *Inf. Process. Lett.* 47 (1993) 173–180.
* D. E. Knuth, *The Art of Computer Programming*, Vol. 4B, section 7.2.2.2 (reluctant doubling).
* G. Pólya, "Geometrisches über die Verteilung der Nullstellen gewisser ganzer transzendenter Funktionen", *Sitzungsber. Bayer. Akad. Wiss.* (1920) 285–290.
* M. Ercsey-Ravasz and Z. Toroczkai, "Optimization hardness as transient chaos in an analog approach to constraint satisfaction", *Nature Physics* 7 (2011) 966–970.
* OEIS A096268, the period-doubling sequence.
