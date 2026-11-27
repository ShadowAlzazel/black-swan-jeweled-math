"""
holomorphic_spectre_solver.py
-----------------------------
Prototype implementation of the HOLOMORPHIC SPECTRE 3-SAT solver (Stages 1-6)
plus comparative baselines (Uniform Random Sweep, Golden-Ratio Sweep) and CLI
support for running directly on DIMACS .cnf files (including
../03_golden_sweep/missed_sat_instance_35.cnf).
"""

from __future__ import annotations
import argparse
import cmath
import json
import math
import time
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple
import numpy as np

from constants_and_tilings import (
    CNFInstance,
    EulerPhasorField,
    FeigenbaumPhasonDriver,
    GOLDEN_ANGLE,
    PINWHEEL_ANGLE,
    PHI,
    PLASTIC_RHO,
    cauchy_winding_number,
    chiral_spectre_mesh_centers,
    disk_to_upper_half_plane,
    load_dimacs_cnf,
    mobius_disk_automorphism,
    reduce_sl2z_fundamental_domain,
    spectre_prototile_vertices,
)


@dataclass
class SolveResult:
    solver_name: str
    instance_name: str
    n_vars: int
    m_clauses: int
    alpha: float
    solved: bool
    best_unsat: int
    steps: int
    elapsed_ms: float
    witness: Optional[List[bool]] = None
    diagnostics: Optional[Dict] = None


def _fast_clause_arrays(cnf: CNFInstance):
    m = cnf.m_clauses
    c_vars = np.zeros((m, 3), dtype=np.int32)
    c_pos = np.zeros((m, 3), dtype=bool)
    for ci, cl in enumerate(cnf.clauses):
        padded = list(cl) + [cl[-1]] * max(0, 3 - len(cl))
        for j in range(3):
            c_vars[ci, j] = abs(padded[j]) - 1
            c_pos[ci, j] = padded[j] > 0
    return c_vars, c_pos


def _eval_unsat_fast(assign: np.ndarray, c_vars: np.ndarray, c_pos: np.ndarray) -> np.ndarray:
    """Return boolean mask of UNSATISFIED clauses."""
    lit_ok = assign[c_vars] == c_pos
    return ~np.any(lit_ok, axis=1)


# ---------------------------------------------------------------------------
# Baseline 1: Uniform Random WalkSAT-lite / Sweep
# ---------------------------------------------------------------------------
def solve_uniform_random(
    cnf: CNFInstance,
    max_steps: int = 800,
    restarts: int = 8,
    noise: float = 0.35,
    seed: int = 0,
) -> SolveResult:
    t0 = time.perf_counter()
    rng = np.random.default_rng(seed)
    c_vars, c_pos = _fast_clause_arrays(cnf)
    best_unsat = cnf.m_clauses
    best_assign: Optional[List[bool]] = None
    steps_used = 0
    steps_per_restart = max(1, max_steps // max(1, restarts))

    for _ in range(restarts):
        assign = rng.integers(0, 2, size=cnf.n_vars).astype(bool)
        for _ in range(steps_per_restart):
            steps_used += 1
            unsat_mask = _eval_unsat_fast(assign, c_vars, c_pos)
            u_count = int(unsat_mask.sum())
            if u_count < best_unsat:
                best_unsat = u_count
                best_assign = assign.tolist()
            if u_count == 0:
                return SolveResult(
                    solver_name="uniform_random_walk",
                    instance_name=cnf.name,
                    n_vars=cnf.n_vars,
                    m_clauses=cnf.m_clauses,
                    alpha=cnf.alpha,
                    solved=True,
                    best_unsat=0,
                    steps=steps_used,
                    elapsed_ms=(time.perf_counter() - t0) * 1000.0,
                    witness=assign.tolist(),
                )
            unsat_idx = np.flatnonzero(unsat_mask)
            chosen_cl = int(rng.choice(unsat_idx))
            cand_vars = c_vars[chosen_cl]
            if rng.random() < noise:
                flip_v = int(rng.choice(cand_vars))
            else:
                # Greedy break-count choice
                best_v = int(cand_vars[0])
                best_score = cnf.m_clauses + 1
                for v in cand_vars:
                    assign[v] = not assign[v]
                    sc = int(_eval_unsat_fast(assign, c_vars, c_pos).sum())
                    assign[v] = not assign[v]
                    if sc < best_score:
                        best_score = sc
                        best_v = int(v)
                flip_v = best_v
            assign[flip_v] = not assign[flip_v]

    return SolveResult(
        solver_name="uniform_random_walk",
        instance_name=cnf.name,
        n_vars=cnf.n_vars,
        m_clauses=cnf.m_clauses,
        alpha=cnf.alpha,
        solved=False,
        best_unsat=best_unsat,
        steps=steps_used,
        elapsed_ms=(time.perf_counter() - t0) * 1000.0,
        witness=best_assign if best_unsat == 0 else None,
    )


# ---------------------------------------------------------------------------
# Baseline 2: Deterministic Golden-Ratio Quasicrystal Sweep
# ---------------------------------------------------------------------------
def solve_golden_sweep(
    cnf: CNFInstance,
    max_steps: int = 800,
    seed: int = 0,
) -> SolveResult:
    """
    Deterministic golden-ratio phase sweep across the assignment hypercube
    with simple greedy descent polish (mirrors fixed-stride sweep baselines).
    """
    t0 = time.perf_counter()
    c_vars, c_pos = _fast_clause_arrays(cnf)
    best_unsat = cnf.m_clauses
    best_assign: Optional[List[bool]] = None
    steps_used = 0

    # Irrational golden offsets per variable
    var_freqs = np.array([(i + 1) / PHI for i in range(cnf.n_vars)], dtype=float)
    n_restarts = 8
    steps_per = max(1, max_steps // n_restarts)

    for r in range(n_restarts):
        # Golden-sequence seed point
        t_param = (r + 1 + seed * 0.137) * PHI
        assign = (np.sin(2.0 * math.pi * (t_param * var_freqs)) >= 0.0).astype(bool)

        for s in range(steps_per):
            steps_used += 1
            unsat_mask = _eval_unsat_fast(assign, c_vars, c_pos)
            u_count = int(unsat_mask.sum())
            if u_count < best_unsat:
                best_unsat = u_count
                best_assign = assign.tolist()
            if u_count == 0:
                return SolveResult(
                    solver_name="golden_sweep",
                    instance_name=cnf.name,
                    n_vars=cnf.n_vars,
                    m_clauses=cnf.m_clauses,
                    alpha=cnf.alpha,
                    solved=True,
                    best_unsat=0,
                    steps=steps_used,
                    elapsed_ms=(time.perf_counter() - t0) * 1000.0,
                    witness=assign.tolist(),
                )
            unsat_idx = np.flatnonzero(unsat_mask)
            # Pick unsatisfied clause via golden-ratio stride
            pick_idx = int(((s + 1) * PHI * len(unsat_idx)) % len(unsat_idx))
            cand_vars = c_vars[unsat_idx[pick_idx]]

            # Deterministic golden-threshold kick vs greedy pick
            frac = ((steps_used * PHI) % 1.0)
            if frac < 0.25:
                flip_v = int(cand_vars[int((steps_used * PHI * 3) % 3)])
            else:
                best_v = int(cand_vars[0])
                best_score = cnf.m_clauses + 1
                for v in cand_vars:
                    assign[v] = not assign[v]
                    sc = int(_eval_unsat_fast(assign, c_vars, c_pos).sum())
                    assign[v] = not assign[v]
                    if sc < best_score:
                        best_score = sc
                        best_v = int(v)
                flip_v = best_v
            assign[flip_v] = not assign[flip_v]

    return SolveResult(
        solver_name="golden_sweep",
        instance_name=cnf.name,
        n_vars=cnf.n_vars,
        m_clauses=cnf.m_clauses,
        alpha=cnf.alpha,
        solved=False,
        best_unsat=best_unsat,
        steps=steps_used,
        elapsed_ms=(time.perf_counter() - t0) * 1000.0,
        witness=best_assign if best_unsat == 0 else None,
    )


# ---------------------------------------------------------------------------
# HOLOMORPHIC SPECTRE Solver (All 6 Stages Combined)
# ---------------------------------------------------------------------------
def solve_holomorphic_spectre(
    cnf: CNFInstance,
    max_steps: int = 800,
    n_spectre_tiles: int = 12,
    canonicalize: bool = True,
    seed: int = 0,
) -> SolveResult:
    """
    Stages:
      1. Spectral canonicalization + Poincare disk / SL(2, Z) modular folding diagnostic
      2. Euler-Phasor complex field F_Phi(z) with Plastic-number clause weights
      3. Chiral Spectre mesh + Cauchy winding-number tile screening
      4. Mobius lens zoom onto highest-winding / lowest-|F| tile + Pinwheel frame kicks
      5. Continuous Euler-phasor descent + Feigenbaum edge-of-chaos phason bursts
      6. Discrete rounding + focused conflict-guided repair & exact verification
    """
    t0 = time.perf_counter()
    c_vars, c_pos = _fast_clause_arrays(cnf)

    # Stage 1 & 2: Canonical Euler-Phasor Field
    field = EulerPhasorField(cnf, freq_scale=1.0, canonicalize=canonicalize, seed=seed)

    # Stage 3: Chiral Spectre tiling screen + Cauchy winding probe
    mesh = chiral_spectre_mesh_centers(n_spectre_tiles, radius=0.88, twist=PINWHEEL_ANGLE)
    tile_scores = []
    aliased_tiles = 0
    total_winding_hits = 0

    for center_z, tile_r, orient in mesh:
        poly = spectre_prototile_vertices(scale=tile_r, angle=orient, center=center_z)
        report = cauchy_winding_number(field, poly, samples_per_edge=6)
        if report.nyquist_aliased:
            aliased_tiles += 1
        if abs(report.rounded_winding) > 0:
            total_winding_hits += 1
        # Rank tiles favoring nonzero winding and small minimum |F_Phi|
        priority = (0 if abs(report.rounded_winding) > 0 else 1, report.min_abs_F)
        tile_scores.append((priority, center_z, tile_r, orient, report))

    tile_scores.sort(key=lambda x: x[0])

    # Stage 4: Mobius zoom onto top-ranked candidate tiles + SL(2, Z) diagnostic
    top_k = min(6, len(tile_scores))
    feig = FeigenbaumPhasonDriver(r0=(1.0 / PHI + 0.07 * seed) % 0.9 + 0.05, level_k=1)

    best_unsat = cnf.m_clauses
    best_assign: Optional[List[bool]] = None
    steps_used = 0
    feigenbaum_bursts = 0
    sl2z_folds = 0
    mobius_zooms = 0

    budget_per_tile = max(1, max_steps // max(1, top_k))

    for tile_idx in range(top_k):
        _, center_z, tile_r, orient, rep = tile_scores[tile_idx]
        z_focus = rep.best_z
        if abs(z_focus) >= 0.95:
            z_focus = z_focus / abs(z_focus) * 0.95

        # Mobius lens re-centering + Pinwheel irrational frame rotation
        mobius_zooms += 1
        pinwheel_rot = (tile_idx + 1) * PINWHEEL_ANGLE
        z_warped = mobius_disk_automorphism(center_z, z_focus, rot=pinwheel_rot)

        # Fold via Cayley transform into SL(2, Z) fundamental domain for diagnostics
        tau = disk_to_upper_half_plane(z_warped)
        _, fold_steps = reduce_sl2z_fundamental_domain(tau)
        sl2z_folds += fold_steps

        # Initialize continuous phasors from zoomed complex coordinate
        theta = field.phases_at(z_focus + 0.15 * z_warped)
        stagnant_streak = 0
        local_best_u = cnf.m_clauses

        # Stage 5 & 6: Coupled continuous phasor flow + Feigenbaum phason bursts + discrete polish
        for local_s in range(budget_per_tile):
            steps_used += 1

            # Continuous Euler-phasor energy gradient step (every 3rd step, or on stagnation)
            if local_s % 3 == 0:
                _, grad = field.eval_energy_and_grad_theta(theta)
                # Plastic-ratio damped step size
                lr = 0.85 / (1.0 + 0.02 * local_s / PLASTIC_RHO)
                theta = (theta - lr * grad) % (2.0 * math.pi)

            # Snap phasors to discrete Boolean assignment
            assign_list = field.assignment_from_phases(theta)
            assign = np.array(assign_list, dtype=bool)
            unsat_mask = _eval_unsat_fast(assign, c_vars, c_pos)
            u_count = int(unsat_mask.sum())

            if u_count < best_unsat:
                best_unsat = u_count
                best_assign = assign_list
            if u_count < local_best_u:
                local_best_u = u_count
                stagnant_streak = 0
                feig.adapt_level(stagnant=False)
            else:
                stagnant_streak += 1

            if u_count == 0:
                return SolveResult(
                    solver_name="holomorphic_spectre",
                    instance_name=cnf.name,
                    n_vars=cnf.n_vars,
                    m_clauses=cnf.m_clauses,
                    alpha=cnf.alpha,
                    solved=True,
                    best_unsat=0,
                    steps=steps_used,
                    elapsed_ms=(time.perf_counter() - t0) * 1000.0,
                    witness=assign_list,
                    diagnostics={
                        "aliased_tiles": aliased_tiles,
                        "winding_hits": total_winding_hits,
                        "feigenbaum_bursts": feigenbaum_bursts,
                        "mobius_zooms": mobius_zooms,
                        "sl2z_folds": sl2z_folds,
                        "canonicalized": canonicalize,
                    },
                )

            # Conflict-directed phasor flip + Feigenbaum edge-of-chaos phason kick
            r_feig, z_kick, is_burst = feig.next_kick()
            unsat_idx = np.flatnonzero(unsat_mask)

            if stagnant_streak >= 6 and is_burst:
                # Long-range deterministic Feigenbaum phason flip across multiple violated clauses!
                feigenbaum_bursts += 1
                feig.adapt_level(stagnant=True)
                n_kicks = min(len(unsat_idx), max(2, int(round(1.0 + abs(z_kick) * 2.5))))
                for k_i in range(n_kicks):
                    cl_i = unsat_idx[int(((local_s + k_i + 1) * PLASTIC_RHO * len(unsat_idx)) % len(unsat_idx))]
                    v_orig = int(c_vars[cl_i, int((local_s + k_i) % 3)])
                    rank = field.inv_order[v_orig]
                    # Phase rotation by pi + irrational Pinwheel offset (breaks reflection lock)
                    theta[rank] = (theta[rank] + math.pi + 0.25 * z_kick.real) % (2.0 * math.pi)
                stagnant_streak = 0
            else:
                # Targeted conflict-clause variable flip with Pinwheel-modulated tie breaking
                cl_pick = unsat_idx[int((r_feig * len(unsat_idx)) % len(unsat_idx))]
                cand_vars = c_vars[cl_pick]
                if r_feig < 0.28:
                    # Irrational Pinwheel exploration pick
                    v_orig = int(cand_vars[int(((local_s + 1) * PINWHEEL_ANGLE * 10.0) % 3)])
                else:
                    best_v = int(cand_vars[0])
                    best_sc = cnf.m_clauses + 1
                    for v in cand_vars:
                        assign[v] = not assign[v]
                        sc = int(_eval_unsat_fast(assign, c_vars, c_pos).sum())
                        assign[v] = not assign[v]
                        if sc < best_sc:
                            best_sc = sc
                            best_v = int(v)
                    v_orig = best_v
                rank = field.inv_order[v_orig]
                # Flip phasor across pi boundary (True <-> False) with subtle continuous drift
                theta[rank] = (theta[rank] + math.pi + 0.08 * math.sin(local_s * PINWHEEL_ANGLE)) % (2.0 * math.pi)

            if steps_used >= max_steps:
                break
        if steps_used >= max_steps:
            break

    return SolveResult(
        solver_name="holomorphic_spectre",
        instance_name=cnf.name,
        n_vars=cnf.n_vars,
        m_clauses=cnf.m_clauses,
        alpha=cnf.alpha,
        solved=False,
        best_unsat=best_unsat,
        steps=steps_used,
        elapsed_ms=(time.perf_counter() - t0) * 1000.0,
        witness=best_assign if best_unsat == 0 else None,
        diagnostics={
            "aliased_tiles": aliased_tiles,
            "winding_hits": total_winding_hits,
            "feigenbaum_bursts": feigenbaum_bursts,
            "mobius_zooms": mobius_zooms,
            "sl2z_folds": sl2z_folds,
            "canonicalized": canonicalize,
        },
    )


def main():
    parser = argparse.ArgumentParser(description="Run HOLOMORPHIC SPECTRE on a DIMACS .cnf file.")
    parser.add_argument("cnf_path", help="Path to DIMACS .cnf file (e.g. ../03_golden_sweep/missed_sat_instance_35.cnf)")
    parser.add_argument("--steps", type=int, default=4000, help="Max step budget")
    parser.add_argument("--tiles", type=int, default=16, help="Number of chiral Spectre tiles")
    parser.add_argument("--seed", type=int, default=0, help="Deterministic seed offset")
    args = parser.parse_args()

    cnf = load_dimacs_cnf(args.cnf_path)
    res = solve_holomorphic_spectre(cnf, max_steps=args.steps, n_spectre_tiles=args.tiles, seed=args.seed)
    print(json.dumps(asdict(res), indent=2))


if __name__ == "__main__":
    main()
