# Black Swan Jeweled Math

A mathematical playground for SAT landscapes, golden-ratio rotations, quasicrystals, aperiodic tilings, and solution-preserving transformations.

The central question is visual and computational: **can we reshape an energy landscape so its satisfying points stay fixed while useful paths become easier to find?** This repository records the experiments, including failures, controls, and independent checks. It does not claim a proof of P=NP or RH.

## Start with the pictures

Open [web/index.html](web/index.html) in a browser. Everything runs locally, including the controls; no server, account, or API key is required.

The same journey is also a web page: [`index.html`](index.html) at the repository root walks through all four experiments with pictures of each lab, every claim tagged proven, open, or observed. With GitHub Pages turned on for this repository (Settings → Pages → Deploy from a branch → `main`, `/ (root)`), it is served at https://shadowalzazel.github.io/black-swan-jeweled-math/ and the labs at `/web/`. `.nojekyll` keeps Pages from running Jekyll. `node make_site_images.js` redraws the page's screenshots and `images/social_card.png` (needs Playwright).

| Explorer | What to try |
| --- | --- |
| [3SAT jewel](web/01_3sat_jewel.html) | Add clauses to a projected Boolean cube and compare solutions with a Hamiltonian walk. |
| [Quasicrystal lab](web/02_quasicrystal_lab.html) | Compare weight schedules and move a cut-and-project window. |
| [Tiling lab](web/03_tiling_lab.html) | Refine actual hat and pinwheel patches, inspect adaptive penalties, and explore coprime lattice points. |
| [Sweep results](web/04_sweep_results.html) | Compare constants and tilings on held-out formulas; inspect a fixed-weight trap and its adaptive escape. |
| [Complex Julia-SAT lab](web/05_complex_julia_sat.html) | The delivered two-plane view of Experiment 04: phasor field, winding probe, Möbius zoom. |
| [Conformal jewels](web/06_conformal_jewels.html) | Carry grids, hats, and pinwheels through eᶻ, z², tanh(z/2), or ℘; follow a tile between planes and apply Euler's identity. |
| [Energy jewels](web/07_energy_jewels.html) | The energy polynomial of a formula: its zeros, the flower it draws on the unit circle, and Cauchy's count of the solutions. |

## Repository layout

```text
experiments/
  01_initial_jewel/    First fixed-weight comparison and quasicrystal construction
  02_tiling_jewels/    Hat/pinwheel geometry and adaptive orientation weights
  03_golden_sweep/     Broad sweep, raw data, exact labels, audit, and toy obstruction
  04_complex_spectre/  Two-plane phasor solver, recorded benchmark, independent audit, Feigenbaum check
  05_counting_jewels/  Exact counting by contraction, deterministic starts, count-preserving rewrites, XOR look-alike
index.html            GitHub Pages front page (the journey)
images/               Lab screenshots and the social card used by index.html
web/                  Standalone interactive explorers
docs/                 Mathematical ideas, findings, and GitHub instructions
provenance.json       Source archive hashes and consolidation checks
requirements.txt      Python dependencies (NumPy, mpmath)
```

Each experiment directory preserves its original scripts, data, and notes. The five stages are kept separately so the development of the ideas remains reproducible. Later audits add files or opt-in switches beside the delivered code and say so; provenance.json records every change. Generated binaries and temporary files are ignored by Git.

## Reproduce the experiments

Use Python 3.11 or newer. The initial experiment and Experiment 04 need NumPy; their checks also use mpmath:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python experiments/01_initial_jewel/black_swan_experiments.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. The tiling/adaptive script uses Python's standard library:

```sh
python experiments/02_tiling_jewels/tiling_adaptive.py
```

The full sweep uses a C++17 compiler named `g++`, plus Python's standard library:

```sh
python experiments/03_golden_sweep/run_sweep.py
python experiments/03_golden_sweep/permutation_stress.py
python experiments/03_golden_sweep/audit_checks.py
python experiments/03_golden_sweep/exhaustive_small.py
```

The complex-plane experiment, its audit, and the Feigenbaum check:

```sh
python experiments/04_complex_spectre/run_experiments.py
python experiments/04_complex_spectre/holomorphic_spectre_solver.py experiments/03_golden_sweep/missed_sat_instance_35.cnf --steps 4000 --tiles 16
python experiments/04_complex_spectre/audit_exp04.py
python experiments/04_complex_spectre/feigenbaum_renormalization.py
python experiments/05_counting_jewels/run_counting.py
python experiments/05_counting_jewels/run_search.py
```

These commands regenerate results next to their scripts. To redraw the hat substitution from the retained authors' source, use Node.js:

```sh
node experiments/02_tiling_jewels/tiling_research/generate_hat_patch.js
python experiments/02_tiling_jewels/tiling_research/validate_hat_patches.py
python experiments/02_tiling_jewels/tiling_research/generate_pinwheel.py
```

## What the broad sweep found

The largest experiment contains 202 configurations, 96 formulas with 12–32 variables, and three starts each: **58,176 main searches**. A locked follow-up adds 1,548 relabeling records, including 387 identity replays. An exhaustive three-variable experiment adds 136,680 toy searches.

The training-selected experimental candidate used **13/8, pinwheel orientations, and normalized adaptive 3–4–1 weights**.

| Method | Held-out SAT starts | Fresh relabelings |
| --- | ---: | ---: |
| Selected pinwheel candidate | 120/129 | 354/387 |
| Plain Breakout | 116/129 | 352/387 |
| WalkSAT | 129/129 | 384/387 |

Every reported solution was checked independently. Exact labeling found 89 SAT and seven UNSAT formulas. These are modest generated cases; the results do not establish a special golden-ratio advantage, a runtime improvement, or a worst-case complexity bound.

The accompanying mathematical obstruction is narrower and exact: any fixed positive clause-position weights can be paired with seven three-variable clauses so strict one-bit descent stalls despite a unique solution. Adaptive weighting and escape moves are outside that obstruction.

Read [the full sweep report](experiments/03_golden_sweep/START_HERE.md) and [the experimental protocol](experiments/03_golden_sweep/ENGINE_README.md) before interpreting the numerical rankings.

## What the complex-plane experiment found

Experiment 04 lifts each formula to a complex phasor field F_Φ on a plane, screens tiles by winding numbers, zooms with a Möbius lens, and adds a Feigenbaum logistic driver. Its recorded benchmark reproduces exactly. The independent audit found:

* **Wrong instance size.** `missed_sat_instance_35.cnf` has 18 variables; the 35 is its instance id. The new solver and plain WalkSAT both solve it in all 100 seeds. It has exactly 9 solutions.
* **Not holomorphic.** F_Φ depends on z̄, so a winding number counts vortices minus antivortices; Cauchy's argument principle does not apply. Its zeros are not solutions: none of 58 zeros found on 40 formulas rounded to one.
* **A plane is a thin slice.** The disk shows 128 of the 262,144 assignments of the missed instance, and none of its solutions. In general, L lines make at most 1 + L + L(L−1)/2 cells, so seeing all 2^n assignments needs exponentially large frequencies.
* **The geometry stages do not help.** Random starts are as good or better, and 3–4× faster. A plain PRNG works as well as the Feigenbaum driver. Step counts match WalkSAT's. The reported 5.3× permutation gain is the best of 20 formulas and not significant overall.
* **The Feigenbaum numbers check out.** The quoted δ ratios are exactly those through period 512. δ, computed as the renormalization eigenvalue, is correct to 45 digits.

See [the Experiment 04 report](experiments/04_complex_spectre/README.md), [the delivered notes with audit marks](docs/COMPLEX_SPECTRE_NOTES.md), and [callbacks and next steps](docs/CALLBACKS_AND_NEXT_STEPS.md).

## Counting instead of searching

Experiment 05 takes its methods from a November 2026 discrete-log paper ([angles and caveats](docs/ANGLES_FROM_THE_DISCRETE_LOG_PAPER.md)): count solutions exactly by tensor contraction, derandomize with conditional expectations, rewrite formulas without changing their solution count, and test every trick on XOR-SAT, the look-alike that algebra solves at once. Read [the Experiment 05 report](experiments/05_counting_jewels/README.md) for the results.

## Related journey

[Fleurdelys](https://shadowalzazel.github.io/fleurdelys-math/) is the related exploration of zeta flowers, Newton dynamics, prime spectra, and cycle eigenvectors. Black Swan Jeweled Math keeps the SAT and tiling experiments in their own repository.

## Attribution and rights

The hat substitution source is by David Smith, Joseph Samuel Myers, Craig S. Kaplan, and Chaim Goodman-Strauss and collaborators, and retains its BSD-3-Clause license. Copies of that license accompany the relevant experiment data and source. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Project-wide license: not yet selected. Upstream third-party licenses remain in force.
