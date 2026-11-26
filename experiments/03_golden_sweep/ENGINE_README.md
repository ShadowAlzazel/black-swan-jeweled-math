# Broad jewel experiment: reproducible engine and protocol

This directory implements 202 configurations on 96 formulas, three starting assignments each: 58,176 main searches. A locked follow-up adds 1,548 searches under formula relabelings, for 59,724 total search runs. These are finite heuristic experiments, not an argument for P=NP.

## Run

Requires Python 3.9+ and a C++17 compiler named `g++`; no Python packages or network connection are needed.

```sh
python run_sweep.py
python permutation_stress.py
```

`run_sweep.py` regenerates formulas and configuration definitions deterministically, computes exact DPLL labels, compiles `sweep_engine.cpp`, runs the search, checks every terminal state independently, and writes `summary.json`, `raw_runs.csv`, and `raw_runs.csv.gz`. `permutation_stress.py` freezes the training-selected noncontrol, transforms only held-out known-SAT formulas, and writes the stress summary and compressed raw records. The source tiling templates are retained in `tiling_templates.json`; their geometries came from the previously constructed actual Hat and pinwheel substitutions. Hat source attribution and BSD license must accompany redistribution of that source material.

## Dataset fixed before running heuristics

- Variables: 12,18,24,32.
- Clauses: round(4.2*n), every clause has three distinct variables.
- Families: planted random, unplanted random, and planted local/geometric (variables on a square grid; clauses use a variable and nearby variables).
- Eight formulas per size/family: 96 total.
- Even replicate numbers form training; odd replicate numbers form held-out testing. The partition is stratified and fixed before search.
- Three shared random starting assignments per formula.
- Exact independent DPLL labeling found 89 SAT and seven UNSAT formulas. Training has 46 SAT/two UNSAT; testing has 43 SAT/five UNSAT. Heuristics never receive the DPLL witness or a planted answer.
- Report SAT-success fractions using confirmed-SAT runs only. A failed heuristic search says nothing about UNSAT.

## Configurations

Four controls: strict uniform greedy, uniform greedy with random restarts, plain adaptive Breakout, and WalkSAT.

Eleven constants: phi,phi^2,1/phi,sqrt(2),sqrt(3),pi,e,pi^2/6,plastic constant,13/8,and one fixed-seed random real. These labels are not eleven independent rotation families: some golden-ratio variants have identical or reversed fractional increments.

Three tiling templates: 169 Hat tiles,625 pinwheel triangles,and a 13x13 periodic square-grid control. A clause uses the same-index tile; only the first m tiles participate. This attachment is arbitrary and does not assume spatial proximity equals logical interaction. Actual positive-length shared polygon boundaries determine the Hat/pinwheel adjacency. Square-grid control neighbors share square edges.

Six schedules for each constant/tiling pair:

1. Static positive cosine weights, amplitude .25.
2. Static positive cosine weights, amplitude .75.
3. Adaptive positive cosine boosts, amplitude .25.
4. Adaptive positive cosine boosts, amplitude .75.
5. Adaptive positive 3-4-1 boosts.
6. Adaptive 3-4-1 boosts plus positive neighbor diffusion.

This is 4+11*3*6=202 configurations.

## Equations

Let u_j(x) be the exact Boolean clause-violation indicator. Energy is E_w(x)=sum_j w_j*u_j(x). Write

a_j(r)=theta_j+2*pi*(j+1)/c+phase0-2*pi*r/c,

where j starts at zero, theta_j is its tile orientation, c is the tested constant, r counts weight updates, and phase0 is seeded per formula/start and shared across configurations. Static weights use r=0 and w_j=1+A*cos(a_j(0)).

At a nonsolution trap, adaptive cosine uses g_j=1+A*cos(a_j(r)). The 3-4-1 rule uses g_j=.25+2*(1+cos(a_j(r)))^2. Only violated clauses are boosted. Divide their g_j by their mean before adding, so the average increment over violated clauses is exactly one. This prevents large raw 3-4-1 boost sizes from masquerading as a geometry effect.

For diffusion, first replace w by .8*w+.2*P*w, where P averages each retained tile's own weight and its retained neighbors' weights. P is nonnegative and row stochastic, even on boundary or isolated sites. Then apply the normalized boost. All weights remain positive. Thus the exact SAT zero set is preserved at every step. Energy need not decrease across weight updates; it strictly decreases on accepted greedy flips for fixed weights.

## Search and cost cap

Every run permits at most 2,049 full-assignment evaluations, 1,000 flips, and 100 weight updates. Every assignment evaluation scans all m clauses using bit masks. The unweighted cost proxy is therefore m times the assignment-evaluation count; 3*m times it is an equivalent literal-inspection count, not a wall-clock or instruction count. Constant calculation, trigonometric work, geometry preparation, and weight/diffusion updates are not charged to this proxy. These are not claims about total runtime speedups.

Greedy variants examine all n one-bit neighbors. Strict greedy stops at a plateau. The restart control spends remaining budget on fresh uniform starts. Adaptive rules change weights at a trap. WalkSAT selects an unsatisfied clause and examines its three flips: zero-break choice if available, otherwise p=.5 random flip versus minimum-break flip. Ties use seeded random choices. Methods consume different numbers of candidate evaluations per decision.

SAT stopping uses the exact Boolean clause test. The floating tolerance 1e-10 is used only to compare weighted move scores. All final states are rechecked by a separate Python Boolean evaluator. This same cap is not an equal wall-clock budget and does not reproduce production SAT-solver implementations.

## Selection and results

Choose the best noncontrol by training SAT success fraction; break ties by smaller training mean assignment-evaluation count, then smaller fixed config ID. Testing does not participate in this selection.

The selected configuration is ID176, normalized3-4-1 with constant13/8 and pinwheel orientations:

- Training:133/138 confirmed-SAT starts solved.
- Held-out:120/129 solved.
- Held-out controls: strictgreedy14/129, restarts97/129, Breakout116/129, WalkSAT129/129.

The held-out-best noncontrol is different (pi/square3-4-1,123/129), but it was not selected because that would use testing for tuning. No special golden-ratio advantage emerges from the finite sweep. These small correlated samples do not establish statistical superiority or asymptotic complexity.

## Locked permutation stress

The frozen winner, Breakout, and WalkSAT are rerun on the43 held-out SAT formulas. Variant0 is identity; variants1-3 independently permute clause order, bijectively relabel variables, and complement selected variables (Boolean gauge flips). The original three starts and a SAT witness transform exactly. Original instance IDs preserve RNG/phase seeds. Identity records reproduce all original fields exactly.

For the three fresh variants combined, 387 runs per method:

- Frozen pinwheel13/8 candidate:354/387.
- Plain Breakout:352/387.
- WalkSAT:384/387.

These transformations preserve satisfiability exactly while exposing arbitrary labeling/geometry-attachment sensitivity. No parameters are retuned after this test. Formula metadata and all transformations are retained for independent audit.

## Golden rotation relationships

The increments1/phi and phi differ by an integer, so c=phi and c=1/phi produce the same j/c rotation modulo one. Also1/phi^2=1-1/phi, so the phi^2 increment is the reversed golden increment. With tile-dependent theta_j, reversing just the frequency need not make the full cosine field identical; the spatial phase must reverse as well for that equivalence. Count these as related representations rather than independent mathematical discoveries.
