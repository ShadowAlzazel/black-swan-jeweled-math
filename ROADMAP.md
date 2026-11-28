# Roadmap

Black Swan asks whether a satisfiability landscape can be reshaped so that its solutions stay fixed while the paths to them become easier to find. Seven experiments in, this page says what we now know, and what to explore next. Every item has a question, a method that reuses code already here, a prediction, and a test that would tell us we were wrong.

Sizes: **S** is one sitting, **M** a few sessions, **L** a project. Tags follow the repository: **proven**, **open**, **observed**. If anything on this list ever looks like evidence that P = NP, the tripwire in [CLAUDE.md](CLAUDE.md) applies first.

## Where we stand

| Experiment | What it established |
| --- | --- |
| 01 Initial jewel | Golden and zeta-zero weights do not beat their own shuffled controls; WalkSAT dominates (**observed**). |
| 02 Tiling jewels | Positive boosts such as 2(1 + cos a)² never move solutions (**proven**); adaptive weights escape traps that fixed weights cannot (**observed**, selected examples). |
| 03 Golden sweep | 58,176 searches show no constant or tiling advantage (**observed**). Fixed clause-position weights can always be trapped by seven clauses (**proven**). |
| 04 Complex Spectre | The two-plane solver's geometry adds nothing measurable. A plane slice sees only polynomially many assignments, so covering all 2^n needs exponentially large frequencies (**proven**). δ is the renormalization eigenvalue, computed to 45 digits (**observed**). |
| 05 Counting jewels | Exact counting by contraction works, at about 2^width; width is ≈ 0.69 n for random formulas at the threshold (**observed**). Conditional-expectation starts cut search 2–9×. Count-preserving rewrites that add variables make search much harder. Every landscape method fails the XOR look-alike (**observed**). |
| 06 Black swan hybrid | Exact XOR recognition inside CNF (Gaussian elimination, checkable 0 = 1 certificates) removes nearly all DPLL work on parity formulas. Every constant and tiling schedule matches a random order and its shuffled control (**observed**). Recovered and verified 2026-11-28. |
| 07 Fourier levels | The look-alike's dip lies wholly on level 3, and no positive reweighting moves it (**proven**). Along a hint ladder (planted, two-hidden, deceptive, golden-hidden, XOR), Johnson's starts help only when level 1 carries a hint, and WalkSAT goes where the low levels point (**observed**). Tovey splitting pushes the hint to level 2. The golden ratio is the weighting that hides a planted solution from level 1 (**proven**). |

The pattern so far: structure helps only when the formula itself has it (geometry, planted hints). Methods that read local hints fail on formulas without hints. Exact counting tells the truth about a landscape at a price set by its width.

## Now: small, high value

### R1. A tripwire harness and a certified hard suite (**S–M**)
* **Question:** does a new method still look good on formulas that leak nothing?
* **Method:** generate unplanted random 3-SAT at α = 4.26 for n = 20 to about 36. Label every formula SAT or UNSAT with exact counts from `counting.py`, slicing as needed. Add quietly planted formulas below the condensation threshold (≈ 3.86 for random 3-SAT; Krzakala and Zdeborová 2009), where planting leaves no trace, and the planted 3-XOR look-alike. Then write one command, `benchmarks/tripwire.py`, that runs any solver across sizes and fits its work against n.
* **Prediction:** every method we have scales exponentially on this suite.
* **How we would know otherwise:** a fitted exponent near zero that survives the five checks in CLAUDE.md.

### R2. Reweighting moves the zeros but never the solutions (**S**)
* **Question:** what does Black Swan's reweighting do to the energy polynomial?
* **Method:** with integer clause weights (duplicated clauses), Z_w(z) = Σ_x z^(Σ w_j u_j(x)). Add a weight control to `web/07_energy_jewels.html`: the zeros move, the flower changes shape, and Z_w(0) = N₀ stays fixed (**proven**, because positive weights keep the zero set).
* **Prediction:** reweightings that help search push the near-solution zeros outward. That is **open**, and worth measuring against WalkSAT flip counts.

### R3. Lee–Yang zeros as n grows (**M**)
* **Question:** do the zeros of Z pinch toward the positive real axis near the threshold as n grows, the way partition-function zeros signal phase transitions (Lee and Yang 1952; Fisher 1965)?
* **Method:** full energy polynomials by contraction at roots of unity plus interpolation, for n up to about 32, across α. Track the smallest angle to the positive axis and the zeros inside the unit disk (Experiment 05 saw 64° → 39° over α = 2 → 6 at n = 20).
* **Prediction:** the angle shrinks with n most sharply near α ≈ 4.2. Positive coefficients keep every zero off that axis (**proven**), so any approach must be a limit.

## Next: search with honest controls

### R4. A derandomization ladder (**M**)
Johnson's method fixes bits using the expected number of violated clauses (first moments). Exact self-reduction fixes bits using exact solution counts. Between them: conditional expectations of pair correlations, then counts on small neighbourhoods.
* **Measure:** solve rate and work at each rung on the R1 suite.
* **Prediction:** each rung costs more and helps less on unplanted formulas. The exact rung always succeeds, at 2^width.

### R5. Marginals and frozen variables (**M**)
* **Method:** exact marginals P(x_i = 1 | satisfiable) from counts. Measure the fraction of frozen variables (fixed in every solution) against α. Near the threshold, solutions split into clusters with frozen cores (Achlioptas and Ricci-Tersenghi 2006).
* **Picture:** color each tile of a hat-patch formula by its marginal.
* **Prediction:** a sharp rise in frozen variables above α ≈ 4 that local search feels as long plateaus.

### R6. Continuous-time SAT dynamics (**M**)
* **Why:** Experiment 04's only useful ingredient was a gradient on the clause energy. Its principled form is Ercsey-Ravasz and Toroczkai's analog solver (Nature Physics 2011): spins plus clause weights that grow exponentially while a clause stays violated. That is Breakout's adaptive weighting in continuous time.
* **Test:** run it on Experiment 03, R1 and the XOR look-alike.
* **Prediction:** analog time looks tame, but the weights grow exponentially. The cost moves into the auxiliary variables rather than vanishing.

### R7. Restart schedules (**S**)
* **Method:** compare Luby's sequence, the ruler sequence 2^ν₂(j), a Feigenbaum-kneading schedule, and fixed restarts, at equal budgets on heavy-tailed instances.
* **Why:** all three structured schedules are driven by the same ruler function (**proven** identities, checked in Experiment 04).
* **Prediction:** no measurable difference among them.

### R8. Which implied clauses help? (**M**)
* **Background:** adding implied resolvents of length ≤ 3 helped WalkSAT slightly in Experiment 05, while rewrites that add variables hurt.
* **Method:** test families of implied clauses (binary resolvents, failed-literal implications, clauses learned by a small CDCL run) as landscape reshapers on R1.
* **Prediction:** modest constant-factor gains, no change in scaling.

### R15. Fourier levels of a landscape (**S**)
* **Background:** group the Walsh characters χ_S by level |S| = k. The level weights ‖P_k E‖²/‖E‖² of the energy form a probability distribution, by the mechanism of the Gelfand–Graev lemma ([docs/IDEMPOTENT_WEIGHTS.md](docs/IDEMPOTENT_WEIGHTS.md)).
  * The XOR look-alike's energy lives on levels 0 and 3 only, so conditional expectations see nothing on it until an equation's last variable is fixed (**proven**). That explains why Johnson's starts gave no help on XOR in Experiment 05.
  * A 3-SAT clause also has weight on levels 1 and 2.
* **Method:** for n ≤ 24, compute these level weights by a fast Walsh–Hadamard transform, along with the agreement between the signs of the level-1 coefficients and the best-matching solution. Run planted and unplanted 3-SAT, the XOR look-alike and Experiment 05's rewrites.
* **Prediction:** that agreement predicts the gain from conditional-expectation starts and WalkSAT's work. It is zero on XOR and small on unplanted formulas at the threshold.
* **How we would know otherwise:** a family with no low-level signal that local search still solves quickly.
* **Literature first:** the planted look-alike is gauge-equivalent to a diluted three-spin ferromagnet with a glass transition (Franz, Mézard, Ricci-Tersenghi, Weigt and Zecchina 2001).
* **Status:** the first pass is [Experiment 07](experiments/07_fourier_levels/README.md), with exact coefficients for any n and hints averaged over every solution. The prediction held for conditional-expectation starts along the ladder, but per formula on Experiment 03's set, h₁ tracks WalkSAT's work (ρ = −0.56) and not Johnson's gain. The definitional extension raises h₁ yet slows WalkSAT. What remains is the unplanted suite of R1: solution-averaged hints for unplanted formulas at the threshold.

## Later: geometry and pictures

### R9. Exact plane search (**M**)
* **Method:** the lines cos θ_k(z) = 0 cut a plane into at most 1 + L + L(L − 1)/2 cells (**proven**). Enumerate them exactly with a sweep, score each cell by its violated clauses, and compare golden, plastic and pinwheel wavevectors with random ones as slice choosers.
* **Picture:** an honest "plane jewel" colored by violated clauses.
* **Prediction:** the constants are no better than random slices.

### R10. A field whose winding number counts zeros (**S–M**)
* **Method:** replace Experiment 04's phasors by complex frequencies, F(z) = Σ a_j e^(λ_j z). Then F is entire and the argument principle holds. Pólya's theorem puts its zeros in strips perpendicular to the sides of the hull of the conj(λ_j), the one-variable case of the tropical-amoeba idea.
* **Prediction:** the mathematics is repaired, but the zeros still say nothing about satisfiability, because the slice bound still applies.

### R11. Tiling-aware contraction orders (**M**)
* **Background:** the hat substitution hierarchy was a poor elimination order (width 73 against 25).
* **Method:** try nested dissection (George 1973) on the tile adjacency graph, and orders that eliminate whole supertiles while keeping their boundaries. Then attempt exact counts on larger patches, such as the 625-triangle pinwheel, by slicing.
* **Prediction:** nested dissection matches min-fill. The hierarchy helps only if its supertiles have short boundaries.

### R12. XOR-SAT and the sheaf obstruction (**S**)
* **Background:** for XOR systems, an unsatisfiable core is a set of equations whose left sides cancel and whose right sides sum to 1. That is a first-cohomology obstruction, and Gaussian elimination finds it in polynomial time (**proven**). FURTHER_EXPLORATIONS #4 asks for the same thing in 3-SAT, where no such shortcut is known.
* **Method:** draw the XOR obstruction on the look-alike instances, and show how the same picture fails for 3-SAT. Experiment 06's `xor_kernel.py` already computes the certificates (rows XOR-ing to 0 = 1); what remains is the picture and the 3-SAT contrast.

## Bridges to the Lemma Office

### R13. Lean receipts for the lemmas used here (**M**, with Kumi-San and Alo)
Candidates:
* Johnson's m/8 bound.
* The positivity filter.
* The Tovey cycle lemma: d implications around a cycle have exactly two solutions.
* The face bound for L lines.
* "Positive weights keep the zero set."
* Lemma 14.33: the weights of a self-adjoint idempotent are non-negative and sum to 1 (short once written as squared norms).

Each is short and has a picture. Check Left to Solve for matching targets first.

### R14. The doubling story (**S**)
* **Content:** Feigenbaum kneading = ν₂(j) mod 2, which the Lemma Office proved for A127752. Add Gray-code flips, Luby runs, the Frobenius identity (1+z)² ≡ 1 + z² (mod 2), and Kumi's orbit v ↦ (v − 2)².
* **Format:** one illustrated page, linked from both sites.

## Already done from earlier lists

* The honest SAT jewel: the energy polynomial and its zeros, Experiment 05 and `web/07`.
* Counting instead of searching: Experiment 05.
* The XOR look-alike: Experiment 05, now a required control (CLAUDE.md), with its algebraic side in Experiment 06's kernel.
* Exact XOR recognition with contradiction certificates: Experiment 06.
* The look-alike's energy polynomial as a signed weight enumerator of its dependencies: Experiment 07.

Still open from [docs/FURTHER_EXPLORATIONS.md](docs/FURTHER_EXPLORATIONS.md): GUE-spaced samplers, 2-adic trees, plastic Rauzy fractals, survey propagation on Spectre factor graphs, and braid holonomy. They are not scheduled; each needs a falsifiable first test before it joins this list.

## Order of attack

1. R1, the tripwire harness and certified hard suite. Every later item is measured against it.
2. R2 and R3: the energy jewels, as both picture and measurement.
3. R4, R5 and R6, the search experiments, each run on R1 and the XOR look-alike.
4. R7, R8 and R15 are cheap and can run alongside.
5. R9–R14 as interest leads, with R13 coordinated with the Lemma Office.
