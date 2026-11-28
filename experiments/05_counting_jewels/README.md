# Experiment 05 — Counting jewels

Experiment 05 tries the methods of a November 2026 discrete-log paper on Black Swan's question: can a formula be reshaped so that solutions become easier to find without the problem changing? The paper's headline is conditional and unverified; its toolkit is classical. See [the angles and caveats](../../docs/ANGLES_FROM_THE_DISCRETE_LOG_PAPER.md). Five ideas were tested:

1. count solutions exactly by tensor contraction instead of searching;
2. derandomize starting points with conditional expectations;
3. rewrite formulas without changing their solutions, or their solution count;
4. check every trick on XOR-SAT, the look-alike that Gaussian elimination solves at once;
5. use tilings as contraction orders.

Claims are tagged **proven**, **open** or **observed**, as elsewhere in the repository.

## Summary

* **Exact counting works and is checked** (**observed**). Contraction reproduces brute force exactly, including the whole energy polynomial of `missed_sat_instance_35.cnf`: N₀..N₂₂ = 9, 94, 538, …, 4. Exact counts reproduce all 96 of Experiment 03's DPLL labels, with 48 cross-checked by brute force. On XOR-SAT the counts equal 2^(n − rank), as linear algebra demands. Counting alone also finds solutions (self-reduction): one witness in 0.02 s, the lexicographically smallest, which is the witness Experiment 03 supplied.
* **The price of counting is the width** (**observed**). Contraction costs about 2^width.
  * Random and planted 3-SAT at α = 4.26 need width ≈ 0.69·n at every size from 20 to 200, so counting beats brute force (2^n) but stays exponential. Planting is invisible to it.
  * Planted 3-XOR-SAT with one equation per variable needs ≈ 0.25·n.
  * Formulas on a square grid grow like ≈ 2.2·√n, and formulas on the 169-hat patch reach 25 (1.9·√n).
* **Tiles help counting through geometry, not through their hierarchy** (**observed**). On formulas living on the 169-hat patch, the mean widths of the elimination orders are:

  | Order | Mean width |
  | --- | ---: |
  | hat substitution hierarchy | 72.8 |
  | hierarchy, reversed | 43.4 |
  | left-to-right sweep | 30.0 |
  | min-degree | 27.0 |
  | min-fill | 25.2 |
  | random | 96.0 |

  The substitution addresses do not trace small boundaries.
* **Deterministic starts by conditional expectations cut the work 2–9×** (**observed**). Johnson's method starts Experiment 03's formulas with 1.4 violated clauses on average, against 11.5 for a random start (it guarantees at most m/8, **proven**), and solves 26 of the 89 with no search at all. WalkSAT's median falls from 28 flips to 3; on 100 planted n = 35 formulas from 48 to 8–13, and on 20 planted n = 200 formulas from 535 to 174–284. A golden-ratio variable order is no better than index, degree or random orders.
* **Rewrites that keep the solutions are harmless or mildly helpful; rewrites that add variables hurt badly** (**observed**). Median WalkSAT flips on Experiment 03's formulas, all solved unless noted:

  | Rewrite | Median flips |
  | --- | ---: |
  | original | 29.5 |
  | complemented variables | 21 |
  | duplicated clauses | 28.5 |
  | implied resolvents added | 17 |
  | definitional extension (count-preserving) | 392 |
  | Tovey splitting (count-preserving) | 7,866, solving only 157 of 178 |

  Tovey splitting also blows the contraction width up from 13 to 56 (min-degree), or 101 if each variable's copies are eliminated together.
* **The look-alike test separates the methods** (**observed**). On planted 3-XOR-SAT, WalkSAT's median grows from 109 flips at n = 20 to beyond 300,000 at n = 100, where it solves only 3 of 8. The same-size planted 3-SAT stays below 200 flips. Gaussian elimination solves all 40 XOR formulas in under 6 ms per batch of 8. Conditional-expectation starts give no help on XOR, whose clauses carry no local bias. Landscape methods read local hints; XOR has none. Experiment 06's kernel recognizes these exact CNF encodings as XOR equations and solves every one, 8/8 at n = 400 in about 12 ms per formula.
* **The energy polynomial is an honest complex-plane jewel** (**observed**, n ≤ 20). Z(z) = Σ N_k z^k is a genuine polynomial, so Cauchy's integral returns N₀ exactly on any circle, and the argument principle counts its zeros inside the unit disk.
  * Across clause density, the median number of zeros inside |z| < 1 rises steadily: 4, 6, 7, 10, 10.5, 11, 11.5 and 14 for α = 2 to 6 at n = 20.
  * The nearest zero moves toward the origin as solutions thin out, from 0.43 at α = 2 to 0.09 at α = 4.25.
  * The zeros also lean toward the positive real axis: the median smallest angle falls from 64° at α = 2 to 39° at α = 6 (n = 20). Coefficients are positive, so no zero ever lies on that axis (**proven**).
  * Unsatisfiable formulas have their zero exactly at the origin, because Z(0) = N₀ = 0. The satisfiability transition is a zero arriving at the origin.
  * The `web/07_energy_jewels.html` lab draws these zeros and the flower Z(re^{iθ}) for 91 formulas. Across all of them, at radii 0.5, 1 and 2, the flower's winding number equals the number of zeros inside, and its average equals N₀.
  * The nearest zero sits on the negative real axis for 43 of the 80 satisfiable formulas. It always lies farther out than the first-order guess N₀/N₁ (median ratio 1.8).
* **Exact count beyond brute force** (**observed**). A planted formula living on the 169-hat patch (169 variables, 710 clauses) has exactly 33,639,161,723,328 solutions (≈ 2^44.9), and 2,713,059,408,490,640 assignments that violate exactly one clause.
  * The count took 125 s: an elimination order of width 25, four pinned variables giving 16 slices of width 21, and six 31-bit primes.
  * A second count with a different order and slicing (7 pinned variables, 128 slices of width 20) gives the same N₀.
  * Brute force would enumerate 2^169 ≈ 7.5 × 10^50 assignments.
  * The same contraction cannot touch a random 169-variable formula at the threshold, whose width is about 117.

## Files

| File | Contents |
| --- | --- |
| `formulas.py` | Generators (random and planted 3-SAT, planted 3-XOR-SAT with its exact CNF, formulas on hat or grid points), Gaussian elimination over GF(2), and five rewrites: gauge, duplication, implied resolvents, Tovey splitting, definitional extension. Each rewrite returns a map that lifts solutions back. |
| `counting.py` | Brute-force energy counts (n ≤ 26); elimination orders and widths; exact bucket elimination over truncated polynomials modulo 31-bit primes with Chinese remaindering; slicing by pinned variables; count-guided self-reduction. |
| `search.py` | Incremental WalkSAT (Experiment 03's rule, about 200,000 flips per second in Python) from any starting assignment; Johnson's conditional-expectation starts under several variable orders. |
| `run_counting.py` → `counting_results.json` | C1 validation, C2 widths, C3 tile orders, C4 zeros across density, C5 one exact count beyond brute force. Rerun single sections with `python run_counting.py c4`. |
| `run_search.py` → `search_results.json` | S1 conditional-expectation starts, S2 rewrites, S3 the XOR look-alike. About 50 seconds. |
| `make_lab_data.py` | Energy polynomials and zeros embedded in `web/07_energy_jewels.html`. |

## Run

```sh
python experiments/05_counting_jewels/run_search.py
python experiments/05_counting_jewels/run_counting.py
python experiments/05_counting_jewels/make_lab_data.py
```

NumPy is the only dependency. Everything is seeded, and every claimed witness is rechecked against the formula, including witnesses lifted back through rewrites.

## What it means for the Black Swan question

Counting shows exactly what a landscape is: how many solutions and near-solutions there are. That makes it the right instrument for judging reweightings, but its cost is set by the width, not by any constant or tiling. Random formulas at the threshold have width proportional to n, so geometry can only help formulas that are themselves geometric.

For search, the reshaping that helped is the most ordinary: start from a conditional-expectation assignment and add implied short clauses. Rewrites that preserve the count but introduce new variables make the same problem look much harder to local search. The XOR control shows that local guidance, not structure, is what all of these methods consume. That is the precise sense in which Black Swan's "jewel" transformations cannot beat the exponential barrier: they reshape the hints, and some formulas have no hints to reshape. (**Open**: whether any count-preserving rewrite helps local search on unplanted formulas at the threshold.)
