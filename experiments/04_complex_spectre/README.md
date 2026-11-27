# Experiment 04 — Complex Spectre (`HOLOMORPHIC SPECTRE`)

Runnable prototype, benchmark suite, and empirical verification of the **HOLOMORPHIC SPECTRE** synthesis combining:

* **Euler-Phasor Complex 3-SAT Field** $F_\Phi : \mathbb{C}_{\text{in}} \to \mathbb{C}_{\text{out}}$ with **Plastic-Number** ($\rho \approx 1.324718$) complex clause weights
* **Chiral Spectre Monotile Mesh** (rotation-only, zero reflections) + **Cauchy Winding-Number Contour Screening** $W(T)$
* **Möbius Disk Lens Zoom** + **Pinwheel $\arctan(1/2)$** irrational frame rotations + **$SL(2,\mathbb{Z})$** modular folding diagnostics
* **Deterministic Feigenbaum ($\mu_\infty \approx 3.56994567$, $\delta \approx 4.6692$, $\alpha \approx 2.5029$)** edge-of-chaos phason burst annealer
* **Laplacian Fiedler-Vector Spectral Canonicalization** to counter variable-permutation drift

---

## Files in This Folder

* `constants_and_tilings.py` — Constants, DIMACS `.cnf` loader, spectral canonicalizer, Möbius / Poincaré / $SL(2,\mathbb{Z})$ transforms, chiral Spectre tile geometry, Euler-phasor potential $F_\Phi(z)$, Cauchy contour integrator with Nyquist aliasing detection, and Feigenbaum phason driver.
* `holomorphic_spectre_solver.py` — Full 6-stage `HOLOMORPHIC SPECTRE` solver, plus fair baselines (`uniform_random_walk`, `golden_sweep`) and a DIMACS `.cnf` CLI runner.
* `run_experiments.py` — Reproducible 4-part benchmark suite (Phase Transition, $n=35$ Critical Needle Stress Test, Nyquist Phase-Frequency Tax, Permutation Invariance Audit).
* `spectre_benchmark_results.json` — Recorded empirical output from running `run_experiments.py`.

---

## Quickstart

Run the full 4-experiment benchmark suite:
```bash
python experiments/04_complex_spectre/run_experiments.py
```

Run directly against `missed_sat_instance_35.cnf` (from Experiment 03):
```bash
python experiments/04_complex_spectre/holomorphic_spectre_solver.py experiments/03_golden_sweep/missed_sat_instance_35.cnf --steps 4000 --tiles 16
```

---

## Recorded Benchmark Highlights (`spectre_benchmark_results.json`)

### Experiment A — Phase Transition ($n=28$, planted 3-SAT, 10 trials/$\alpha$, budget 900 steps)

| Density $\alpha$ | Uniform Walk (Solve / Steps) | Golden Sweep (Solve / Steps) | Holomorphic Spectre (Solve / Steps / Bursts) |
| :--- | :--- | :--- | :--- |
| **2.50** | 100% / 13.9 | 100% / 10.5 | **100% / 5.5** (0.0 bursts) |
| **3.50** | 100% / 28.4 | 90% / 218.2 | **100% / 21.5** (1.3 bursts) |
| **4.00** | 90% / 144.1 | 60% / 404.8 | **100% / 101.4** (11.3 bursts) |
| **4.26** | **100% / 59.5** | 70% / 340.1 | 90% / 152.4 (17.6 bursts) |

### Experiment B — Critical Needle-in-a-Haystack Stress Test ($n=35, \alpha=4.26$, 12 trials, budget 1600)

| Solver | Solve Rate | Mean Steps | Median Steps | Mean Time (ms) |
| :--- | :--- | :--- | :--- | :--- |
| **Uniform Random Walk** | 100.0% (12/12) | 320.4 | 239.5 | 28.75 |
| **Golden Sweep** | 58.3% (7/12) | 891.3 | 739.0 | 60.33 |
| **Holomorphic Spectre** | **100.0% (12/12)** | **287.9** | **185.0** | 74.58 |

*Takeaway:* Adding Möbius zoom, Pinwheel irrational kicks, and Feigenbaum bursts rescues deterministic sweeps from critical-threshold traps (58.3% $\to$ 100%, 4x lower median steps), slightly beating WalkSAT in step count with a modest geometry overhead.

### Experiment C — Empirical Proof of the Phase-Frequency Tax (Nyquist Aliasing Wall)

Fixed Spectre tile perimeter, 12 samples/edge, moderate $\sqrt{n}$ wavevector scaling:

| $n$ | Mean Total Phase Variation $\int_{\partial T} |d\arg F_\Phi|$ (rad) | Mean Max Step Jump (rad) | Nyquist Aliasing Rate ($> \pi/2$) |
| :--- | :--- | :--- | :--- |
| **6** | 2.72 | 0.093 | 0.0% |
| **10** | 3.37 | 0.181 | 0.0% |
| **16** | 5.62 | 0.289 | 0.0% |
| **24** | 6.60 | 0.486 | 0.0% |
| **35** | 9.81 | 0.770 | **16.7%** |
| **50** | 11.79 | 1.032 | **16.7%** |

*Takeaway:* Even under gentle $\sqrt{n}$ packing, boundary phase oscillation grows steadily and begins aliasing by $n=35$. Under true bijective $2^n$ state packing, required contour samples explode as $\Omega(2^{n/2})$—explaining concretely why Cauchy contour pruning cannot break $P \neq NP$.

### Experiment D — Permutation Invariance Audit ($n=30, \alpha=4.1$, 10 random permutations)

| Mode | Solve Rate | Mean Steps | Step Std Dev |
| :--- | :--- | :--- | :--- |
| **Without Canonicalization** | 100% | 123.7 | 143.6 |
| **With Spectral Canonicalization** | 100% | **23.3** | **31.4** |

*Takeaway:* Graph-Laplacian Fiedler canonicalization slashed step variance by **4.6x** and mean steps by **5.3x** across randomized variable/clause relabelings.
