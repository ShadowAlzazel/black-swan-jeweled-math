# Golden-ratio jewel sweep: findings and reproduction

No candidate establishing P=NP emerged. The experiments found useful adaptive local-search behavior, and a rigorous obstruction to a narrower family of fixed-weight searches.

## What ran

- 202 configurations, covering 11 constant labels, three tilings, and six schedules, plus four ordinary search controls.
- 96 formulas with 12,18,24,32 variables; three starting assignments per formula: **58,176 main search runs**.
- Formula families: planted, unplanted random, and planted local/geometric.
- Training/testing split fixed before any heuristic ran.
- Independent exact classification: 89 SAT and seven UNSAT. Testing contains 43 confirmed SAT formulas, giving 129 SAT starts per method.
- A locked relabeling stress test added 1,548 run records: 387 exact identity replays and 1,161 fresh clause/variable/complement permutations.
- A separate exhaustive toy experiment ran **136,680 strict-descent searches**, testing 67 weight templates on all eight starts of all 255 satisfiable subsets of the eight distinct full-width three-variable clauses. The remaining subset is UNSAT. This is a precisely specified small family, not all larger 3SAT.

## Main held-out results

| Method | Solved / 129 confirmed-SAT starts |
| --- | ---: |
| Strict greedy | 14 |
| Greedy with restarts | 97 |
| Plain Breakout | 116 |
| Training-selected pinwheel / 13/8 / normalized 3-4-1 | 120 |
| WalkSAT | 129 |

The candidate was selected using training results only:133/138 training SAT starts. Its nine held-out failures remain failures in the report. Looking at testing to select a different configuration would be tuning on the test set; no such replacement was made.

On three fresh relabelings per held-out SAT formula, with three starts each:

| Frozen method | Solved / 387 |
| --- | ---: |
| Selected pinwheel / 13/8 candidate | 354 |
| Plain Breakout | 352 |
| WalkSAT | 384 |

These results do not establish a special golden-ratio advantage, or a statistically reliable advantage over ordinary Breakout. Runs sharing formulas are correlated. The complete exact solver used no more than40 nodes on any of these modest generated instances; this is not a difficult production-SAT benchmark or evidence about asymptotic scaling. Heuristic failures never certify UNSAT.

## Constants and transformations

Constants: phi,phi²,1/phi,sqrt(2),sqrt(3),pi,e,pi²/6,the plastic constant,13/8,and a seeded random real. Golden powers and reciprocals are related rotation schedules, not independent mathematical discoveries. The main engine uses distinct tie-breaking random seeds for configuration IDs, so small empirical differences between equivalent underlying rotations must not be interpreted as new mathematics.

Tilings: actual169-hat and625-triangle pinwheel substitution patches, and a periodic square-grid control. The constant and tile orientations enter cosine phases. Adaptive schedules change positive penalties at a trap. The diffusion schedule also averages weights over the actual tile-adjacency graph. Mapping clause j to tile j is an arbitrary experimental attachment, not a theorem about SAT locality.

Every adaptive boost is normalized to average one over the currently violated clauses. This controls the raw boost magnitude, addressing the confound in the earlier small demonstration. Diffusion uses a nonnegative row-stochastic matrix, so positive weights stay positive. SAT solutions remain exactly the zero-energy assignments.

All methods use a cap of2049 full-assignment evaluations,1000 flips,and100 weight updates. Methods consume different candidates per decision. Trigonometry, geometry preparation, and weight updates are outside this evaluation proxy, so these are success-under-cap comparisons, not total-runtime speed measurements. See ENGINE_README.md for the full equations and protocol.

## A small obstruction with a proof

Consider any seven fixed positive weights w1,...,w7 assigned to clause positions independently of clause contents. For each three-bit state a other than111, form the clause C_a that is violated only at a: use x_i if a_i=0 and not-x_i if a_i=1.

The seven clauses have unique satisfying assignment111. Place C_000 at an index with minimum weight. At000 the weighted energy is the minimum positive weight. Each one-bit neighbor violates exactly one other clause and therefore has energy at least that minimum. Strict one-bit descent cannot move from000, even though111 solves the formula.

This applies to every fixed content-independent constant/tiling weight template when arbitrary clause ordering is allowed. It is not an obstruction to adaptive updates, sideways moves, nonlocal moves, or content-dependent weighting. Ordinary Breakout escapes the displayed toy examples by changing the weights.

exhaustive_small.json includes the explicit seven clauses and energies for every tested template. Its representative phi/hat example escapes through000 ->100 ->110 ->111 after two boosts (labels list x1,x2,x3).

## Concrete missed SAT instance

missed_sat_instance_35.cnf is an18-variable76-clause SAT instance that the frozen candidate missed under its cap. missed_sat_instance_35.witness supplies a complete independently verified satisfying assignment. This is a budget-limited failure of the tested heuristic, not a lower bound on the work necessary to solve the instance.

## Audit

audit_checks.py imports no engine evaluator. The audit independently:

- classified all96 formulas;
- checked every main-run final state and35,957 SAT witnesses;
- checked all1548 stress-run terminal states and inverse-mapped witnesses;
- verified172 relabeling maps and5676 assignment equivalences;
- reproduced all387 identity replays exactly;
- checked positive zero-set preservation, normalized positive boosts, stochastic diffusion, golden identities, and the fixed-weight obstruction.

No mismatch was found. Detailed results are in audit_report.json. The small experiment and the analytical obstruction are separately recorded so that neither is confused with a proof about all NP problems.

## Reproduce

Python's standard library and a C++17 compiler called g++ suffice for the main engine:

```sh
python run_sweep.py
python permutation_stress.py
python audit_checks.py
python exhaustive_small.py
```

The cached tiling_templates.json makes the main engine self-contained. The archive also contains the original geometry JSON files under tiling_research/ because exhaustive_small.py uses those orientation snapshots. The main CSV and compressed stress-test CSV are included. Reproduction regenerates the compiled engine as needed.

Source attribution and the upstream hat BSD license are included. See ENGINE_README.md for detailed data formats, phases, stochastic decisions, and reproducibility limits.
