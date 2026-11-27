# Experiment 04 — Complex Spectre (`HOLOMORPHIC SPECTRE`)

A prototype that lifts 3-SAT into two linked complex planes, ℂ_in → ℂ_out, and searches with:

* an **Euler-phasor field** F_Φ(z) = Σ_c e^{2πic/ρ} Π_j (1 − s_{c,j} e^{iθ_{v(c,j)}(z)})/2, where θ(z) is a real-linear function of z and ρ is the plastic number;
* **winding numbers** around 12–16 Spectre-labelled outlines;
* a **Möbius lens** and **pinwheel** arctan(1/2) rotations, with an **SL(2, ℤ)** folding diagnostic;
* a **logistic map** just above the Feigenbaum accumulation point μ∞ as a deterministic noise source;
* **Fiedler-vector canonicalization** of the variable order.

The delivered solver, benchmark script and recorded results are kept here. This directory also adds an independent audit, controls and switches for the audit, and a separate check of the Feigenbaum numbers. Read the audit summary before the recorded tables.

## Audit summary

* **Reproduction.** Every recorded solve rate and step count reproduces exactly. Only wall-clock times differ, because they depend on the machine.
* **The missed instance has 18 variables, not 35.** `missed_sat_instance_35.cnf` is Experiment 03's instance id 35 (random family, held-out split): 18 variables, 76 clauses, α = 4.22. Experiment B's n = 35 does not mirror it. HOLOMORPHIC SPECTRE solves it: the CLI command below finishes in 35 steps, and 100 seeds give 100/100 with a median of 47 steps. WalkSAT does the same: 100/100, median 44. The formula has exactly 9 satisfying assignments, in two one-flip clusters of 6 and 3.
* **F_Φ is not holomorphic.** Each phase θ_k(z) = k_re x + k_im y + offset depends on both z and z̄. Measured |∂F/∂z̄| / |∂F/∂z| has a median of 0.83–1.09, and 45–56% of the plane reverses orientation. The boundary integral is therefore a topological degree: it counts vortices minus antivortices, not zeros, and Cauchy's argument principle does not apply. Every located zero's charge matched the sign of the Jacobian, and the square's dense boundary winding equalled N₊ − N₋.
* **Zeros of F_Φ are not satisfying assignments.** Across 40 n = 35 formulas, the solver's disk holds 27 vortices and 31 antivortices; 7 formulas have none. No zero rounds to a solution. Zeros round to slightly worse than typical assignments: 19.3 violated clauses at zeros, 18.6 on average over the disk, and 10.1 at the best point. |F_Φ| does not track violated clauses (Spearman 0.04 and −0.02); the nonnegative clause energy E(θ(z)) does (0.81 and 0.72). The missed instance's field has no zeros at all in [−0.9, 0.9]², so every tile's winding is 0 there.
* **One z-plane shows very few assignments.** Rounding sgn cos θ_k(z) cuts the plane along L lines into at most 1 + L + L(L−1)/2 cells. For the missed instance, the unit disk shows 128 of 262,144 assignments; that bound is 211. None of the 9 solutions is among them, and the best point violates 2 clauses. The solver leaves the plane with its first flip (θ_k += π), so Stages 3–4 only choose six starting points.
* **Ablations.** Every ablation keeps the step budget and runs on more formulas than the recorded benchmark.
  * *Random starts* (no tiles, Möbius lens or winding screen) solve as many formulas, run 3–4× faster, and use fewer steps on the 100 n = 35 formulas (sign test p = 0.005).
  * *A seeded PRNG in place of the logistic map* makes no detectable difference anywhere.
  * *The true Tile(1,1) contour* makes no detectable difference either.
  * *The continuous gradient step* helps on Experiment 03's formulas (p < 0.001) but not detectably at n = 35.
* **Controls.** The delivered `golden_sweep` is a weak control (58/100 at n = 35). Against WalkSAT with Experiment 03's rule, step counts are statistically indistinguishable (p = 0.23 at n = 35; p = 0.62 on Experiment 03's formulas), and WalkSAT runs 1.5–3.4× faster. The "rescue" in Experiment B is the familiar effect of adding noise to a deterministic local search.
* **Permutations.** The recorded 5.3× improvement is the most favorable of 20 base formulas. Across all 20, canonicalization lowers the mean on 13 and raises it on 7 (sign test p = 0.26). The remaining variance comes from clause order: canonicalizing clauses as well makes all ten relabelings run identically on 20/20 formulas. The headline formula then takes 79 steps, not 23.3.
* **Phase-frequency tax.** Most of Experiment C's growth comes from the imposed √n frequency scaling. At n = 50 with the solver's fixed frequency, aliasing occurs in 6.7% of trials, against 26.7% with √n scaling. A 12-sample winding number never disagreed with a 1,500-sample one in 360 trials. The rigorous form of the tax is the line-arrangement bound: to show all 2^n assignments, the unit disk needs about 2^{(n+1)/2} lines and some wavevector |k| ≥ (π/2)(L/n − 1). That is about 1.2 × 10⁴ at n = 35.
* **Geometry.** The delivered "Spectre" contour takes 13 unit steps and closes implicitly with an edge of length √7 at 10.9°. It is a simple 14-gon, not the Spectre. The code now also provides the true Tile(1,1) (area 3 + 3√3; its hat sibling has area 8√3 as expected). The "mesh" is 12–16 separate outlines on a golden-angle spiral. It covers about 45% of the disk, with small overlaps, and is not a tiling.
* **Constants.** All eleven constants agree with mpmath to double precision, including μ∞ and c∞ = μ∞/2 − μ∞²/4.
* **Feigenbaum.** The quoted ratios 4.668548581 → 4.669060661 → 4.669171555 are exactly δ_k for the superstable cycles through periods 128, 256 and 512. The quoted −2.502906472 is the signed ratio d₈/d₉ (periods 256/512). The renormalization fixed point gives δ as the only eigenvalue outside the unit circle, correct to 45 digits, and α = 1/g(1) to 48 digits.

None of this bears on P versus NP beyond the observation that a 2-D slice cannot see 2^n assignments without exponentially large frequencies.

## Files

| File | Status | Contents |
| --- | --- | --- |
| `constants_and_tilings.py` | delivered + audit additions | Constants, DIMACS loader, canonicalization, Möbius / Poincaré / SL(2, ℤ) helpers, contours, F_Φ, winding integrator, Feigenbaum driver. |
| `holomorphic_spectre_solver.py` | delivered + audit additions | The six-stage solver, the `uniform_random_walk` and `golden_sweep` controls, and a CLI. |
| `run_experiments.py` | delivered, unchanged | Experiments A–D; rewrites `spectre_benchmark_results.json`. |
| `spectre_benchmark_results.json` | delivered, unchanged | Results recorded on the delivering machine. |
| `audit_exp04.py` | added | The independent audit (sections 0–9 below); writes `audit_results.json`. |
| `audit_results.json` | added | Output of the full audit run reported here. |
| `feigenbaum_renormalization.py` | added | Superstable cycles, renormalization eigenvalue, kneading / ruler checks; writes `feigenbaum_results.json`. |
| `feigenbaum_results.json` | added | Output of that script. |

Changes to delivered code. With default arguments the recorded benchmark still reproduces exactly; the audit checks this on every run.

* Witnesses are rechecked with the pure-Python evaluator before a run reports success.
* Clauses with more than three literals are rejected; they used to be silently truncated to three.
* `spectre_prototile_vertices` keeps its behavior but now documents what it is. Added: `tile_ab_vertices`, `spectre_tile11_vertices`, `canonicalize_cnf` and `PRNGPhasonDriver`.
* Opt-in switches on `solve_holomorphic_spectre`: `init="random"`, `use_gradient=False`, `driver="prng"`, `polygon="tile11"`, `canonical_clauses=True`. New diagnostics: `explicit_phase_flips` and `rounded_bit_changes`.
* A standard `solve_walksat` control, the same rule as Experiment 03's WalkSAT, and a `--solver` option on the CLI.

## Run

Python 3.11+ with NumPy; the audit and the Feigenbaum check also need mpmath (`requirements.txt`). The recorded benchmark reproduced exactly under Python 3.13 with NumPy 2.5.3 and with the pinned NumPy 2.3.5; the audit run reported here used NumPy 2.5.3 and mpmath 1.3.0.

```sh
python experiments/04_complex_spectre/run_experiments.py      # recorded benchmark A-D (seconds)
python experiments/04_complex_spectre/holomorphic_spectre_solver.py experiments/03_golden_sweep/missed_sat_instance_35.cnf --steps 4000 --tiles 16
python experiments/04_complex_spectre/holomorphic_spectre_solver.py experiments/03_golden_sweep/missed_sat_instance_35.cnf --solver walksat
python experiments/04_complex_spectre/audit_exp04.py          # full audit, about 3 minutes; --quick for a smaller run
python experiments/04_complex_spectre/feigenbaum_renormalization.py   # about 5 seconds
```

`run_experiments.py` overwrites the recorded JSON with this machine's timings; `git checkout` restores the delivered file. The solver CLI prints the result and then an independent witness check against the file:

```text
"solved": true, "steps": 35, "n_vars": 18, "m_clauses": 76, ... "winding_hits": 0
Independent witness check against experiments/03_golden_sweep/missed_sat_instance_35.cnf: PASS (76 clauses)
```

## Recorded benchmark (as delivered)

The tables below are the delivered `spectre_benchmark_results.json`. Each one is followed by what the audit found.

### Experiment A — phase transition (n = 28, planted 3-SAT, 10 trials per α, budget 900 steps)

| Density α | Uniform walk (solve / steps) | Golden sweep (solve / steps) | Holomorphic Spectre (solve / steps / bursts) |
| :--- | :--- | :--- | :--- |
| 2.50 | 100% / 13.9 | 100% / 10.5 | 100% / 5.5 (0.0 bursts) |
| 3.50 | 100% / 28.4 | 90% / 218.2 | 100% / 21.5 (1.3 bursts) |
| 4.00 | 90% / 144.1 | 60% / 404.8 | 100% / 101.4 (11.3 bursts) |
| 4.26 | 100% / 59.5 | 70% / 340.1 | 90% / 152.4 (17.6 bursts) |

*Audit:* reproduced exactly. With 10 trials per cell, these solve rates have 95% intervals about ±20 points wide.

### Experiment B — critical stress test (n = 35, α = 4.26, 12 planted trials, budget 1,600)

| Solver | Solve rate | Mean steps | Median steps | Mean time (ms) |
| :--- | :--- | :--- | :--- | :--- |
| Uniform random walk | 100.0% (12/12) | 320.4 | 239.5 | 28.75 |
| Golden sweep | 58.3% (7/12) | 891.3 | 739.0 | 60.33 |
| Holomorphic Spectre | 100.0% (12/12) | 287.9 | 185.0 | 74.58 |

*Delivered takeaway:* Möbius zoom, pinwheel kicks and Feigenbaum bursts rescue deterministic sweeps from critical-threshold traps, slightly beating WalkSAT in step count.

*Audit:* reproduced exactly, apart from wall time on this machine (11.3 / 27.1 / 33.1 ms). On 100 such formulas the step counts are not distinguishable from WalkSAT's. Random starts beat the tile-based starts (see the ablation table). The "uniform random walk" is a WalkSAT-like walk that scores flips by the total violated count and restarts every budget/8 steps; it is not WalkSAT.

### Experiment C — phase-frequency tax (fixed contour, 12 samples per edge, frequency scaled by √n/2.5)

| n | Mean ∫∂T \|d arg F_Φ\| (rad) | Mean max step jump (rad) | Aliasing rate (jump > π/2) |
| :--- | :--- | :--- | :--- |
| 6 | 2.72 | 0.093 | 0.0% |
| 10 | 3.37 | 0.181 | 0.0% |
| 16 | 5.62 | 0.289 | 0.0% |
| 24 | 6.60 | 0.486 | 0.0% |
| 35 | 9.81 | 0.770 | 16.7% |
| 50 | 11.79 | 1.032 | 16.7% |

*Delivered takeaway:* exponential work re-emerges as boundary sampling frequency, Ω(2^{n/2}) under bijective packing.

*Audit:* reproduced exactly. Each 16.7% is one trial out of six. See "Experiment C revisited" below for 30 trials, a fixed-frequency comparison, and the counting bound that makes the Ω(2^{n/2}) statement rigorous for the rounding map.

### Experiment D — permutation audit (one n = 30, α = 4.1 formula, 10 relabelings, budget 1,100)

| Mode | Solve rate | Mean steps | Step std |
| :--- | :--- | :--- | :--- |
| Without canonicalization | 100% | 123.7 | 143.6 |
| With spectral canonicalization | 100% | 23.3 | 31.4 |

*Audit:* reproduced exactly; it is the most favorable of 20 formulas (see below).

## Audit details

All runs use the delivered code paths with the switches listed above. Steps are recorded as delivered: one search step can flip several bits during a burst. "Fewer / more" counts the runs in which a variant used fewer or more steps than the delivered solver on the same formula and seed, with failures counted as budget + 1; p is an exact two-sided sign test.

### Ablations and controls

Experiment 03's 89 confirmed-SAT formulas (n = 12–32), three seeds each, budget 1,000 steps:

| Solver | Solved / 267 | Wilson 95% | Median steps | Mean ms | Fewer / more vs delivered | p |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| Holomorphic Spectre (delivered) | 262 | 0.957–0.992 | 20 | 21.4 | — | — |
| random starts (no tiles, lens, winding) | 264 | 0.967–0.996 | 20 | 5.3 | 121 / 132 | 0.53 |
| no gradient steps | 262 | 0.957–0.992 | 39 | 21.6 | 75 / 168 | <0.001 |
| PRNG instead of logistic map | 263 | 0.962–0.994 | 18 | 21.9 | 112 / 131 | 0.25 |
| true Tile(1,1) contour | 261 | 0.952–0.990 | 17 | 21.3 | 81 / 70 | 0.42 |
| no canonicalization | 265 | 0.973–0.998 | 18 | 20.8 | 126 / 134 | 0.66 |
| golden sweep (delivered control) | 193 | 0.666–0.773 | 144 | 10.7 | 55 / 206 | <0.001 |
| uniform random walk (delivered control) | 264 | 0.967–0.996 | 28 | 3.3 | 97 / 162 | <0.001 |
| WalkSAT p = 0.5 (Experiment 03 rule) | 264 | 0.967–0.996 | 21 | 6.3 | 125 / 134 | 0.62 |

On the 43 held-out formulas (129 runs), the solvers solved respectively 125, 126, 127, 129, 125, 128, 99, 128 and 127 runs. For context only: Experiment 03's own held-out table (a different engine, with a 1,000-flip and 2,049-evaluation cap) has WalkSAT 129/129, the frozen pinwheel candidate 120/129 and Breakout 116/129.

Experiment 04's family: 100 planted n = 35, α = 4.26 formulas (the first 12 are Experiment B's), budget 1,600:

| Solver | Solved / 100 | Wilson 95% | Median steps | Mean ms | Fewer / more vs delivered | p |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| Holomorphic Spectre (delivered) | 97 | 0.915–0.990 | 130.5 | 33.6 | — | — |
| random starts (no tiles, lens, winding) | 98 | 0.930–0.994 | 100.5 | 11.2 | 64 / 35 | 0.005 |
| no gradient steps | 95 | 0.888–0.978 | 159.5 | 30.6 | 43 / 55 | 0.27 |
| PRNG instead of logistic map | 98 | 0.930–0.994 | 121 | 33.7 | 52 / 46 | 0.61 |
| true Tile(1,1) contour | 99 | 0.946–0.998 | 130 | 37.8 | 45 / 42 | 0.83 |
| no canonicalization | 97 | 0.915–0.990 | 108 | 35.0 | 61 / 38 | 0.027 |
| golden sweep (delivered control) | 58 | 0.482–0.672 | 832.5 | 31.1 | 20 / 78 | <0.001 |
| uniform random walk (delivered control) | 99 | 0.946–0.998 | 192 | 11.6 | 47 / 51 | 0.76 |
| WalkSAT p = 0.5 (Experiment 03 rule) | 97 | 0.915–0.990 | 111 | 23.0 | 56 / 43 | 0.23 |

When the winding screen finds a nonzero tile (0.82 tiles per run here), it ranks that tile first. That tile starts at a zero of F_Φ, and zeros round to worse-than-typical assignments. This is a plausible reason random starts do better.

### Permutation audit across 20 formulas

The protocol is Experiment D's (n = 30, α = 4.1, planted, ten relabelings, budget 1,100) on base seeds 31415 + 17b.

| | Pooled mean | Pooled median | Pooled std |
| --- | ---: | ---: | ---: |
| without canonicalization | 122.5 | 49 | 183.7 |
| variable canonicalization (delivered) | 89.2 | 29 | 120.0 |
| variables and clauses canonical (`canonical_clauses=True`) | 97.8 | 75.5 | — |

Variable canonicalization lowers the per-formula mean on 13 formulas and raises it on 7 (p = 0.26). It lowers the per-formula std on 13 and raises it on 7 (p = 0.26). The canonical variable order itself was exactly invariant in every test. Step counts still varied, because clause weights e^{2πic/ρ} and violated-clause picks depend on clause positions. With clauses canonical too, all ten relabelings of every formula give identical runs (std 0 on 20/20). This is the invariance the delivered note aimed for, and it removes variance without improving the average.

### The missed instance

100 seeds, budget 4,000 steps, 16 tiles:

| Solver | Solved | Median steps | Max steps |
| --- | ---: | ---: | ---: |
| Holomorphic Spectre (delivered) | 100 | 47 | 995 |
| random starts | 100 | 57.5 | 738 |
| no gradient steps | 100 | 84.5 | 869 |
| PRNG driver | 100 | 38 | 1,370 |
| true Tile(1,1) contour | 100 | 55 | 781 |
| no canonicalization | 100 | 87.5 | 990 |
| golden sweep | 46 | 4,000 | 4,000 |
| uniform random walk | 100 | 126 | 1,544 |
| WalkSAT | 100 | 44 | 679 |

The 9 solutions, listed x1…x18, are:

* cluster A: `000000010000000000`, `000100010000000000`, `000010010000000000`, `000110010000000000`, `000100010000000100`, `000110010000000100`
* cluster B: `001000010100000110`, `001100010100000110`, `001110010100000110`

Experiment 03's supplied witness is the first one, and the CLI run above finds the last one. Across the 100 seeds the solvers found all nine.

### What one z-plane can see

The table counts distinct rounded assignments in the unit disk on a 1,400 × 1,400 grid, which is a lower bound. The upper bound is the arrangement bound from the line count.

| n | freq 1 (solver): found / bound | freq √n/2.5 (Exp C): found / bound | 2^n |
| ---: | ---: | ---: | ---: |
| 6 | 16 / 22 | 15 / 22 | 64 |
| 10 | 38 / 56 | 64 / 121 | 1,024 |
| 16 | 98 / 154 | 246 / 466 | 65,536 |
| 24 | 249 / 466 | 967 / 1,892 | 1.7 × 10⁷ |
| 35 | 642 / 1,226 | 3,467 / 7,022 | 3.4 × 10¹⁰ |
| 50 | 1,517 / 3,004 | 11,545 / 24,091 | 1.1 × 10¹⁵ |
| 100 | 9,967 / 20,707 | 128,665 / 345,697 | 1.3 × 10³⁰ |

To show all 2^n assignments, the arrangement needs L lines with 1 + L + L(L−1)/2 ≥ 2^n. Since n families share those lines, some |k| must reach (π/2)(L/n − 1): about 62 for n = 18, 1.2 × 10⁴ for n = 35, 1.5 × 10⁶ for n = 50 and 2.5 × 10¹³ for n = 100. This is the rigorous version of the delivered "phase-frequency tax": exponential frequency is forced by counting cells, before any question of contour sampling arises.

### Zeros of F_Φ (vortex census)

* `missed_sat_instance_35`, field as in the solver: 0 vortices and 0 antivortices in [−0.9, 0.9]². |F| against violated clauses has Spearman 0.04; E(θ(z)) has 0.81.
* Experiment B formula 0: 1 vortex and 1 antivortex. Both round to 21 violated clauses, while the region's best point has 12. The single tile with nonzero winding (W = −1) contains the antivortex.
* 40 formulas of the n = 35 family: 27 vortices, 31 antivortices, 7 formulas without zeros, no satisfying zero. Zeros average 19.3 violated clauses, against 18.6 over the disk.

### Experiment C revisited (30 trials per n)

| n | √n/2.5 frequency: max jump, aliasing (95% CI) | fixed frequency 1: max jump, aliasing (95% CI) |
| ---: | --- | --- |
| 6 | 0.19 rad, 3.3% (0.6–16.7%) | 0.16 rad, 0% (0–11.4%) |
| 16 | 0.33 rad, 3.3% (0.6–16.7%) | 0.15 rad, 0% (0–11.4%) |
| 24 | 0.44 rad, 3.3% (0.6–16.7%) | 0.14 rad, 0% (0–11.4%) |
| 35 | 0.81 rad, 10.0% (3.5–25.6%) | 0.33 rad, 6.7% (1.8–21.3%) |
| 50 | 1.21 rad, 26.7% (14.2–44.4%) | 0.41 rad, 6.7% (1.8–21.3%) |

In all 360 trials, the 12-sample winding number agreed with a 1,500-sample one.

### Feigenbaum check (`feigenbaum_renormalization.py`, `feigenbaum_results.json`)

| Quantity | Value | Error vs literature |
| --- | --- | --- |
| δ_k, superstable cycles through period 16,384 (k = 14) | 4.66920159553749 | 1.4 × 10⁻⁸ |
| signed α_k = d₁₃/d₁₄ | −2.50290787444281 | 6.5 × 10⁻¹⁰ |
| δ from the renormalization eigenvalue (32 terms in x²) | 4.669201609102990671853204… | 2.1 × 10⁻⁴⁵ |
| α = 1/g(1) at the fixed point | −2.502907875095892822283903… | 9.7 × 10⁻⁴⁹ |

* **Convergence.** Each period doubling divides the error in δ_k by about δ while doubling the work, so the ratio method costs exponential work per digit. The eigenvalue method converges exponentially in the number of terms.
* **The rest of the spectrum.** After δ the eigenvalues are 0.15963, −0.12365, −0.05731, 0.02548, −0.01015, …; 0.15963 = α⁻² and 0.02548 = α⁻⁴.
* **"Compose twice and rescale" works.** The rescaled iterates h_k(x) = f^{2^k}(s_k x)/s_k of z² + c∞ converge to the fixed point g: max |h_k − g| on [−1, 1] falls from 0.034 at k = 0 to 1.7 × 10⁻⁹ at k = 9. The two-level error ratios 0.018, 0.021, 0.022, 0.023 approach α⁻⁴ ≈ 0.0255. The scale ratios s_k/s_{k+1} converge to α.
* **The kneading sequence.** At the accumulation point it is RLRRRLRLRLRRRLRR…; for the first 1,023 symbols it equals the parity of the ruler function ν₂(j), the 2-adic valuation. The same function gives the Gray-code flip positions and Luby's restart runs (1, 2, …, 2^{ν₂(u)}). All three identities are checked in the script.
* **The sign of α is Euler's half-turn.** log α = log|α| + iπ = 0.917453205977598 + iπ.

## What the experiment shows

The useful parts of the solver are ordinary: focused flips on violated clauses, noise, and on small formulas a continuous relaxation of the clause energy E. The complex-plane stages choose six starting points. Those starting points are no better than random ones and cost most of the running time. The two-plane picture is still worth keeping as a picture, and as a clean demonstration of why a 2-D slice cannot see 2^n states. The interactive version of that picture is `web/06_conformal_jewels.html`.
