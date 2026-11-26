# Black Swan Jeweled Math

A mathematical playground for SAT landscapes, golden-ratio rotations, quasicrystals, aperiodic tilings, and solution-preserving transformations.

The central question is visual and computational: **can we reshape an energy landscape so its satisfying points stay fixed while useful paths become easier to find?** This repository records the experiments, including failures, controls, and independent checks. It does not claim a proof of P=NP or RH.

## Start with the pictures

Open [web/index.html](web/index.html) in a browser. Everything runs locally, including the controls; no server, account, or API key is required.

| Explorer | What to try |
| --- | --- |
| [3SAT jewel](web/01_3sat_jewel.html) | Add clauses to a projected Boolean cube and compare solutions with a Hamiltonian walk. |
| [Quasicrystal lab](web/02_quasicrystal_lab.html) | Compare weight schedules and move a cut-and-project window. |
| [Tiling lab](web/03_tiling_lab.html) | Refine actual hat and pinwheel patches, inspect adaptive penalties, and explore coprime lattice points. |
| [Sweep results](web/04_sweep_results.html) | Compare constants and tilings on held-out formulas; inspect a fixed-weight trap and its adaptive escape. |

## Repository layout

```text
experiments/
  01_initial_jewel/    First fixed-weight comparison and quasicrystal construction
  02_tiling_jewels/    Hat/pinwheel geometry and adaptive orientation weights
  03_golden_sweep/     Broad sweep, raw data, exact labels, audit, and toy obstruction
web/                  Standalone interactive explorers
docs/                 Mathematical ideas, findings, and GitHub instructions
provenance.json       Source archive hashes and consolidation checks
requirements.txt      Python dependencies for the initial experiment
```

Each experiment directory preserves its original scripts, data, and notes. The three stages are kept separately so the development of the ideas remains reproducible. Generated binaries and temporary files are ignored by Git.

## Reproduce the experiments

Use Python 3.11 or newer. Only the initial experiment needs NumPy and mpmath:

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

## Related journey

[Fleurdelys](https://shadowalzazel.github.io/fleurdelys-math/) is the related exploration of zeta flowers, Newton dynamics, prime spectra, and cycle eigenvectors. Black Swan Jeweled Math keeps the SAT and tiling experiments in their own repository.

## Attribution and rights

The hat substitution source is by David Smith, Joseph Samuel Myers, Craig S. Kaplan, and Chaim Goodman-Strauss and collaborators, and retains its BSD-3-Clause license. Copies of that license accompany the relevant experiment data and source. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Project-wide license: not yet selected. Upstream third-party licenses remain in force.
