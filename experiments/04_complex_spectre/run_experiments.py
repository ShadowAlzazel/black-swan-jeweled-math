"""
run_experiments.py
------------------
Reproducible benchmark & empirical verification suite for Experiment 04:
  A. Phase-Transition Benchmark (alpha = 2.5, 3.5, 4.0, 4.26)
  B. Planted Needle-in-a-Haystack Stress Test (simulating missed_sat_instance_35.cnf)
  C. Empirical Measurement of the Phase-Frequency Tax (Nyquist Aliasing Wall)
  D. Permutation Invariance Audit (callback to 03_golden_sweep/permutation_stress.py)

Writes full results to spectre_benchmark_results.json.
"""

from __future__ import annotations
import json
import math
import os
import sys
import numpy as np

# Ensure local imports work regardless of working directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from constants_and_tilings import (
    EulerPhasorField,
    PINWHEEL_ANGLE,
    cauchy_winding_number,
    generate_random_3sat,
    permute_cnf,
    spectre_prototile_vertices,
)
from holomorphic_spectre_solver import (
    solve_golden_sweep,
    solve_holomorphic_spectre,
    solve_uniform_random,
)


def experiment_a_phase_transition():
    """Compare Uniform Walk, Golden Sweep, and Holomorphic Spectre across alpha."""
    alphas = [2.5, 3.5, 4.0, 4.26]
    n_vars = 28
    trials = 10
    budget = 900
    summary = {}

    for alpha in alphas:
        stats = {
            "uniform_random_walk": {"solved": 0, "steps": [], "best_unsat": []},
            "golden_sweep": {"solved": 0, "steps": [], "best_unsat": []},
            "holomorphic_spectre": {"solved": 0, "steps": [], "best_unsat": [], "bursts": []},
        }
        for t in range(trials):
            cnf = generate_random_3sat(n_vars=n_vars, alpha=alpha, seed=1000 + int(alpha * 100) + t, planted=True)
            r_u = solve_uniform_random(cnf, max_steps=budget, seed=t)
            r_g = solve_golden_sweep(cnf, max_steps=budget, seed=t)
            r_h = solve_holomorphic_spectre(cnf, max_steps=budget, n_spectre_tiles=10, seed=t)

            for r in (r_u, r_g, r_h):
                bucket = stats[r.solver_name]
                if r.solved:
                    bucket["solved"] += 1
                bucket["steps"].append(r.steps)
                bucket["best_unsat"].append(r.best_unsat)
                if r.solver_name == "holomorphic_spectre" and r.diagnostics:
                    bucket["bursts"].append(r.diagnostics["feigenbaum_bursts"])

        summary[f"alpha_{alpha:.2f}"] = {
            name: {
                "solve_rate": round(b["solved"] / trials, 3),
                "solved_count": b["solved"],
                "trials": trials,
                "mean_steps": round(float(np.mean(b["steps"])), 1),
                "mean_best_unsat": round(float(np.mean(b["best_unsat"])), 2),
                **({"mean_feigenbaum_bursts": round(float(np.mean(b["bursts"])), 1)} if "bursts" in b else {}),
            }
            for name, b in stats.items()
        }
    return summary


def experiment_b_needle_stress():
    """
    Harder planted 3-SAT instances right at the critical threshold (n=35, alpha=4.26),
    mirroring the exact variable count & profile of missed_sat_instance_35.cnf.
    """
    n_vars = 35
    alpha = 4.26
    trials = 12
    budget = 1600
    out = {
        "n_vars": n_vars,
        "alpha": alpha,
        "trials": trials,
        "budget": budget,
        "solvers": {},
    }
    buckets = {
        "uniform_random_walk": [],
        "golden_sweep": [],
        "holomorphic_spectre": [],
    }
    for t in range(trials):
        cnf = generate_random_3sat(n_vars=n_vars, alpha=alpha, seed=7700 + t, planted=True)
        for fn in (solve_uniform_random, solve_golden_sweep, solve_holomorphic_spectre):
            res = fn(cnf, max_steps=budget, seed=t)
            buckets[res.solver_name].append(res)

    for name, res_list in buckets.items():
        solved = sum(1 for r in res_list if r.solved)
        out["solvers"][name] = {
            "solve_rate": round(solved / trials, 3),
            "solved_count": solved,
            "mean_steps": round(float(np.mean([r.steps for r in res_list])), 1),
            "median_steps": round(float(np.median([r.steps for r in res_list])), 1),
            "mean_best_unsat": round(float(np.mean([r.best_unsat for r in res_list])), 2),
            "mean_elapsed_ms": round(float(np.mean([r.elapsed_ms for r in res_list])), 2),
        }
    return out


def experiment_c_nyquist_phase_tax():
    """
    Empirically measure the Phase-Frequency Tax (Nyquist Aliasing Wall):
    as n grows (and high-frequency quasicrystal wavevectors pack more states into C_in),
    track total boundary phase variation int |d arg F_Phi| and aliasing rate around
    a fixed-size Spectre tile at fixed sample count.
    """
    ns = [6, 10, 16, 24, 35, 50]
    poly = spectre_prototile_vertices(scale=0.35, angle=PINWHEEL_ANGLE, center=0.15 + 0.1j)
    samples_per_edge = 12
    rows = []
    for n in ns:
        variations = []
        max_jumps = []
        aliased_count = 0
        trials = 6
        for t in range(trials):
            cnf = generate_random_3sat(n_vars=n, alpha=4.2, seed=900 + n * 10 + t, planted=True)
            # Scale frequency moderately with sqrt(n) to reflect packing density in 2D
            field = EulerPhasorField(cnf, freq_scale=math.sqrt(n) / 2.5, canonicalize=True)
            rep = cauchy_winding_number(field, poly, samples_per_edge=samples_per_edge)
            variations.append(rep.total_abs_phase_variation)
            max_jumps.append(rep.max_step_phase_jump)
            if rep.nyquist_aliased:
                aliased_count += 1
        rows.append({
            "n_vars": n,
            "samples_per_edge": samples_per_edge,
            "mean_total_phase_variation_rad": round(float(np.mean(variations)), 2),
            "mean_max_step_jump_rad": round(float(np.mean(max_jumps)), 3),
            "nyquist_aliasing_rate": round(aliased_count / trials, 3),
        })
    return rows


def experiment_d_permutation_invariance():
    """
    Callback to 03_golden_sweep/permutation_stress.py:
    Compare solution stability across 10 random variable/clause permutations of the same
    hard n=30, alpha=4.1 planted CNF with vs. without spectral canonicalization.
    """
    base_cnf = generate_random_3sat(n_vars=30, alpha=4.1, seed=31415, planted=True)
    n_perms = 10
    budget = 1100

    raw_steps, canon_steps = [], []
    raw_solved, canon_solved = 0, 0

    for p in range(n_perms):
        p_cnf, _ = permute_cnf(base_cnf, seed=500 + p)
        r_raw = solve_holomorphic_spectre(p_cnf, max_steps=budget, canonicalize=False, seed=0)
        r_can = solve_holomorphic_spectre(p_cnf, max_steps=budget, canonicalize=True, seed=0)
        raw_steps.append(r_raw.steps)
        canon_steps.append(r_can.steps)
        if r_raw.solved:
            raw_solved += 1
        if r_can.solved:
            canon_solved += 1

    return {
        "n_vars": base_cnf.n_vars,
        "alpha": round(base_cnf.alpha, 3),
        "n_permutations": n_perms,
        "without_canonicalization": {
            "solve_rate": round(raw_solved / n_perms, 3),
            "mean_steps": round(float(np.mean(raw_steps)), 1),
            "std_steps": round(float(np.std(raw_steps)), 1),
        },
        "with_spectral_canonicalization": {
            "solve_rate": round(canon_solved / n_perms, 3),
            "mean_steps": round(float(np.mean(canon_steps)), 1),
            "std_steps": round(float(np.std(canon_steps)), 1),
        },
    }


def main():
    print("Running Experiment A: Phase-Transition Benchmark...")
    exp_a = experiment_a_phase_transition()
    print("Running Experiment B: Planted Needle-in-a-Haystack Stress Test (n=35, alpha=4.26)...")
    exp_b = experiment_b_needle_stress()
    print("Running Experiment C: Phase-Frequency Tax (Nyquist Aliasing Wall)...")
    exp_c = experiment_c_nyquist_phase_tax()
    print("Running Experiment D: Permutation Invariance Audit...")
    exp_d = experiment_d_permutation_invariance()

    results = {
        "experiment": "04_complex_spectre",
        "algorithm": "HOLOMORPHIC SPECTRE",
        "experiment_a_phase_transition": exp_a,
        "experiment_b_needle_stress_n35": exp_b,
        "experiment_c_nyquist_phase_tax": exp_c,
        "experiment_d_permutation_invariance": exp_d,
    }
    out_path = os.path.join(SCRIPT_DIR, "spectre_benchmark_results.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    print(f"Saved results to {out_path}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
