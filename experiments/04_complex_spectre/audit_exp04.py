"""
audit_exp04.py
--------------
Independent audit of Experiment 04 (HOLOMORPHIC SPECTRE).

    python experiments/04_complex_spectre/audit_exp04.py          # full audit
    python experiments/04_complex_spectre/audit_exp04.py --quick  # fewer trials

Writes audit_results.json next to this script. Requires numpy and mpmath.

Sections
  0. Reproduction of the recorded benchmark (every non-timing field)
  1. Constants, recomputed with mpmath
  2. Geometry: the delivered contour, the Tile(a, b) family, and the "mesh"
  3. Holomorphy: Cauchy-Riemann residual of F_Phi
  4. Zeros of F_Phi: vortex/antivortex census, degree identity, SAT relevance
  5. Plane reachability: how many assignments one z-plane slice can show
  6. Ablations and controls (Experiment 03 formulas; Experiment 04 n=35 family)
  7. Permutation audit across many base formulas
  8. missed_sat_instance_35.cnf: exact solution set and multi-seed runs
  9. Experiment C revisited: is the aliasing driven by n or by the imposed frequency?
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import os
import platform
import sys
import time
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

import constants_and_tilings as ct  # noqa: E402
import holomorphic_spectre_solver as hs  # noqa: E402
import run_experiments as rx  # noqa: E402

EXP03_DIR = os.path.join(REPO_DIR, "experiments", "03_golden_sweep")
MISSED_CNF = os.path.join(EXP03_DIR, "missed_sat_instance_35.cnf")
RECORDED = os.path.join(SCRIPT_DIR, "spectre_benchmark_results.json")
OUT_PATH = os.path.join(SCRIPT_DIR, "audit_results.json")


# ---------------------------------------------------------------------------
# Small statistics helpers (standard library + numpy only)
# ---------------------------------------------------------------------------
def wilson(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (round(max(0.0, mid - half), 3), round(min(1.0, mid + half), 3))


def sign_test_p(wins: int, losses: int) -> float:
    """Exact two-sided binomial sign test, ties excluded."""
    n = wins + losses
    if n == 0:
        return 1.0
    k = min(wins, losses)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / 2 ** n
    return round(min(1.0, 2 * tail), 4)


def rankdata(a: np.ndarray) -> np.ndarray:
    """Average ranks (ties share their mean rank)."""
    a = np.asarray(a)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(len(a), dtype=float)
    _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
    sums = np.bincount(inv, weights=ranks)
    return (sums / counts)[inv]


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean()
    rb -= rb.mean()
    return float(np.dot(ra, rb) / math.sqrt(np.dot(ra, ra) * np.dot(rb, rb)))


# ---------------------------------------------------------------------------
# Independent CNF helpers (no solver code)
# ---------------------------------------------------------------------------
def read_dimacs(path: str) -> Tuple[int, List[List[int]]]:
    n, clauses, cur = 0, [], []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if not parts or parts[0] == "c":
                continue
            if parts[0] == "p":
                n = int(parts[2])
                continue
            for tok in parts:
                lit = int(tok)
                if lit == 0:
                    clauses.append(cur)
                    cur = []
                else:
                    cur.append(lit)
    return n, clauses


def all_solutions(n: int, clauses: Sequence[Sequence[int]]) -> np.ndarray:
    """Exhaustive satisfying assignments (bit i of the code = variable i+1). n <= 24."""
    codes = np.arange(1 << n, dtype=np.uint32)
    sat = np.ones(1 << n, dtype=bool)
    for cl in clauses:
        ok = np.zeros(1 << n, dtype=bool)
        for lit in cl:
            bit = (codes >> np.uint32(abs(lit) - 1)) & np.uint32(1)
            ok |= (bit == 1) if lit > 0 else (bit == 0)
        sat &= ok
    return np.flatnonzero(sat)


def code_to_bits(code: int, n: int) -> List[bool]:
    return [bool((int(code) >> i) & 1) for i in range(n)]


def exp03_instances() -> List[dict]:
    with open(os.path.join(EXP03_DIR, "dataset.json"), "r", encoding="utf-8") as fh:
        return json.load(fh)["instances"]


def instance_to_cnf(inst: dict) -> ct.CNFInstance:
    return ct.CNFInstance(
        n_vars=int(inst["n"]),
        clauses=[tuple(int(x) for x in cl) for cl in inst["clauses"]],
        name=f"exp03_id{inst['id']}_{inst['family']}_n{inst['n']}",
    )


# ---------------------------------------------------------------------------
# Vectorized field evaluation on point sets (chunked to bound memory)
# ---------------------------------------------------------------------------
def field_values(field: ct.EulerPhasorField, zs: np.ndarray, chunk: int = 60000) -> np.ndarray:
    zs = np.asarray(zs, dtype=complex).ravel()
    out = np.empty(zs.shape, dtype=complex)
    for s in range(0, zs.size, chunk):
        z = zs[s:s + chunk]
        th = field.k_real[:, None] * z.real[None, :] + field.k_imag[:, None] * z.imag[None, :] + field.phase_offset[:, None]
        ph = np.exp(1j * th)
        acc = np.zeros(z.size, dtype=complex)
        for c in range(field.m):
            prod = np.ones(z.size, dtype=complex)
            for j in range(3):
                prod *= 0.5 * (1.0 - field.c_signs[c, j] * ph[field.c_vars[c, j]])
            acc += field.clause_weights[c] * prod
        out[s:s + chunk] = acc
    return out


def rank_codes(field: ct.EulerPhasorField, zs: np.ndarray, chunk: int = 60000) -> np.ndarray:
    """Rounded assignments sgn(cos theta(z)) packed as uint64 words, canonical rank order."""
    zs = np.asarray(zs, dtype=complex).ravel()
    words = (field.n + 63) // 64
    out = np.zeros((zs.size, words), dtype=np.uint64)
    for s in range(0, zs.size, chunk):
        z = zs[s:s + chunk]
        th = field.k_real[:, None] * z.real[None, :] + field.k_imag[:, None] * z.imag[None, :] + field.phase_offset[:, None]
        bits = (np.cos(th) >= 0.0)
        for r in range(field.n):
            w, b = divmod(r, 64)
            out[s:s + chunk, w] |= bits[r].astype(np.uint64) << np.uint64(b)
    return out


def unsat_counts_from_rank_bits(field: ct.EulerPhasorField, zs: np.ndarray, chunk: int = 60000) -> np.ndarray:
    zs = np.asarray(zs, dtype=complex).ravel()
    out = np.empty(zs.size, dtype=np.int32)
    for s in range(0, zs.size, chunk):
        z = zs[s:s + chunk]
        th = field.k_real[:, None] * z.real[None, :] + field.k_imag[:, None] * z.imag[None, :] + field.phase_offset[:, None]
        truth = np.cos(th) >= 0.0  # rank order
        lit_true = np.where(field.c_signs[:, :, None] > 0, truth[field.c_vars], ~truth[field.c_vars])
        out[s:s + chunk] = np.count_nonzero(~np.any(lit_true, axis=1), axis=0)
    return out


def boundary_winding(field: ct.EulerPhasorField, polygon: Sequence[complex], samples_per_edge: int) -> Tuple[float, float]:
    pts = []
    nv = len(polygon)
    for i in range(nv):
        a, b = polygon[i], polygon[(i + 1) % nv]
        t = np.arange(samples_per_edge) / samples_per_edge
        pts.append(a + t * (b - a))
    pts = np.concatenate(pts + [np.array([polygon[0]])])
    vals = field_values(field, pts)
    d = np.angle(vals[1:] / vals[:-1])
    return float(d.sum() / (2 * math.pi)), float(np.abs(d).max())


def points_in_polygon(px: np.ndarray, py: np.ndarray, polygon: Sequence[complex]) -> np.ndarray:
    inside = np.zeros(px.shape, dtype=bool)
    nv = len(polygon)
    for i in range(nv):
        a, b = polygon[i], polygon[(i + 1) % nv]
        xa, ya, xb, yb = a.real, a.imag, b.real, b.imag
        cond = (ya > py) != (yb > py)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = xa + (py - ya) * (xb - xa) / (yb - ya)
        inside ^= cond & (px < xint)
    return inside


# ---------------------------------------------------------------------------
# 0. Reproduction
# ---------------------------------------------------------------------------
def audit_reproduction() -> dict:
    with open(RECORDED, "r", encoding="utf-8") as fh:
        rec = json.load(fh)
    new = {
        "experiment_a_phase_transition": rx.experiment_a_phase_transition(),
        "experiment_b_needle_stress_n35": rx.experiment_b_needle_stress(),
        "experiment_c_nyquist_phase_tax": rx.experiment_c_nyquist_phase_tax(),
        "experiment_d_permutation_invariance": rx.experiment_d_permutation_invariance(),
    }
    diffs, timing = [], []

    def walk(x, y, path):
        if isinstance(x, dict):
            for k in sorted(set(x) | set(y)):
                walk(x.get(k), y.get(k), f"{path}/{k}")
        elif isinstance(x, list):
            if len(x) != len(y):
                diffs.append((path, "length", len(x), len(y)))
            for i, (u, v) in enumerate(zip(x, y)):
                walk(u, v, f"{path}[{i}]")
        elif x != y:
            (timing if path.endswith("mean_elapsed_ms") else diffs).append((path, x, y))

    for key, value in new.items():
        walk(rec[key], value, key)
    return {
        "non_timing_mismatches": diffs,
        "exact_match_excluding_wall_clock": not diffs,
        "wall_clock_fields_recorded_vs_this_machine": timing,
    }


# ---------------------------------------------------------------------------
# 1. Constants
# ---------------------------------------------------------------------------
def audit_constants() -> dict:
    import mpmath as mp

    mp.mp.dps = 40
    phi = (1 + mp.sqrt(5)) / 2
    rho = mp.findroot(lambda x: x ** 3 - x - 1, 1.3)
    roots = mp.polyroots([1, 0, -1, -1], maxsteps=200, extraprec=200)
    lam = [r for r in roots if mp.im(r) > 0][0]
    # Literature values (OEIS A006890, A006891, A098587), 27 significant digits.
    delta_ref = mp.mpf("4.66920160910299067185320382")
    alpha_ref = mp.mpf("2.50290787509589282228390287")
    mu_ref = mp.mpf("3.56994567187094490184200515")
    checks = {
        "PHI": (ct.PHI, phi),
        "SILVER_RATIO": (ct.SILVER_RATIO, 1 + mp.sqrt(2)),
        "PLASTIC_RHO (real root of x^3-x-1)": (ct.PLASTIC_RHO, rho),
        "PLASTIC_CONJ real part": (ct.PLASTIC_CONJ.real, mp.re(lam)),
        "PLASTIC_CONJ imag part": (ct.PLASTIC_CONJ.imag, mp.im(lam)),
        "FEIGENBAUM_DELTA": (ct.FEIGENBAUM_DELTA, delta_ref),
        "FEIGENBAUM_ALPHA (|alpha|)": (ct.FEIGENBAUM_ALPHA, alpha_ref),
        "LOGISTIC_MU_INF": (ct.LOGISTIC_MU_INF, mu_ref),
        "MANDELBROT_C_INF = mu/2 - mu^2/4": (ct.MANDELBROT_C_INF, mu_ref / 2 - mu_ref ** 2 / 4),
        "PINWHEEL_ANGLE = atan(1/2)": (ct.PINWHEEL_ANGLE, mp.atan(mp.mpf(1) / 2)),
        "GOLDEN_ANGLE = 2 pi / phi^2": (ct.GOLDEN_ANGLE, 2 * mp.pi / phi ** 2),
    }
    rows = {}
    for name, (module_value, ref) in checks.items():
        err = abs(mp.mpf(module_value) - ref)
        rows[name] = {
            "module": module_value,
            "reference": mp.nstr(ref, 25),
            "abs_error": float(err),
            "ok_to_double_precision": bool(err < 1e-14 * max(1, abs(ref))),
        }
    identities = {
        "rho * |lambda|^2 = 1 (product of the three roots)": float(abs(rho * abs(lam) ** 2 - 1)),
        "|lambda| = rho^(-1/2)": float(abs(abs(lam) - rho ** mp.mpf(-0.5))),
        "phi^2 - phi - 1": float(abs(phi ** 2 - phi - 1)),
        "golden angle in degrees": float(360 / phi ** 2),
        "pinwheel angle in degrees": float(mp.degrees(mp.atan(mp.mpf(1) / 2))),
    }
    return {"constants": rows, "identities": identities}


# ---------------------------------------------------------------------------
# 2. Geometry
# ---------------------------------------------------------------------------
def polygon_area(p: Sequence[complex]) -> float:
    return 0.5 * sum((p[i].conjugate() * p[(i + 1) % len(p)]).imag for i in range(len(p)))


def segments_cross(a, b, c, d) -> bool:
    def orient(p, q, r):
        return ((q - p).conjugate() * (r - p)).imag
    return orient(a, b, c) * orient(a, b, d) < 0 and orient(c, d, a) * orient(c, d, b) < 0


def self_intersections(p: Sequence[complex]) -> int:
    n, hits = len(p), 0
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            if segments_cross(p[i], p[(i + 1) % n], p[j], p[(j + 1) % n]):
                hits += 1
    return hits


def audit_geometry() -> dict:
    headings = [0, 60, 150, 90, 180, 120, 210, 150, 240, 300, 210, 270, 330, 270]
    raw = [0j]
    for h in headings[:-1]:
        raw.append(raw[-1] + cmath.exp(1j * math.radians(h)))
    closing = raw[0] - raw[-1]
    delivered = {
        "listed_unit_steps_used": len(headings) - 1,
        "implicit_closing_edge_length": round(abs(closing), 6),
        "implicit_closing_edge_heading_deg": round(math.degrees(cmath.phase(closing)) % 360, 3),
        "intended_closing_edge": "length 1 at 270 deg (the unused last heading)",
        "closing_edge_equals_sqrt7": abs(abs(closing) - math.sqrt(7)) < 1e-12,
        "area_unit_edges": round(abs(polygon_area(raw)), 6),
        "self_intersections": self_intersections(raw),
    }
    family = {}
    for label, (a, b) in {
        "Tile(1, sqrt3) hat": (1.0, math.sqrt(3)),
        "Tile(1, 1) spectre base": (1.0, 1.0),
        "Tile(sqrt3, 1) turtle": (math.sqrt(3), 1.0),
    }.items():
        p = ct.tile_ab_vertices(a, b)
        last_kind, last_heading = ct.TILE_AB_EDGES[-1]
        last_len = a if last_kind == "a" else b
        gap = abs(p[-1] + last_len * cmath.exp(1j * math.radians(last_heading)) - p[0])
        family[label] = {
            "closure_gap": gap,
            "area": round(abs(polygon_area(p)), 9),
            "self_intersections": self_intersections(p),
        }
    family["hat area equals 8*sqrt(3) (8 kites)"] = abs(family["Tile(1, sqrt3) hat"]["area"] - 8 * math.sqrt(3)) < 1e-9
    family["Tile(1,1) area equals 3 + 3*sqrt(3)"] = abs(family["Tile(1, 1) spectre base"]["area"] - (3 + 3 * math.sqrt(3))) < 1e-9

    # Is the "chiral Spectre mesh" a tiling? Rasterize the tiles the solver screens.
    g = 700
    xs = np.linspace(-1, 1, g)
    px, py = np.meshgrid(xs, xs)
    disk = px ** 2 + py ** 2 <= 0.88 ** 2
    mesh_report = {}
    for n_tiles in (12, 16):
        for label, fn in (("delivered 14-gon", ct.spectre_prototile_vertices), ("Tile(1,1)", ct.spectre_tile11_vertices)):
            cover = np.zeros(px.shape, dtype=np.int32)
            for cz, tr, orient in ct.chiral_spectre_mesh_centers(n_tiles, radius=0.88, twist=ct.PINWHEEL_ANGLE):
                cover += points_in_polygon(px, py, fn(scale=tr, angle=orient, center=cz))
            covered = cover[disk] >= 1
            mesh_report[f"{n_tiles} tiles, {label}"] = {
                "fraction_of_disk_covered": round(float(covered.mean()), 3),
                "fraction_of_disk_uncovered": round(float(1 - covered.mean()), 3),
                "fraction_of_covered_area_overlapped": round(float((cover[disk] >= 2).sum() / max(1, covered.sum())), 3),
            }
    return {"delivered_contour": delivered, "tile_ab_family": family, "solver_mesh_coverage": mesh_report}


# ---------------------------------------------------------------------------
# 3. Holomorphy
# ---------------------------------------------------------------------------
def audit_holomorphy(fields: Dict[str, ct.EulerPhasorField]) -> dict:
    out = {}
    rng = np.random.default_rng(11)
    h = 1e-6
    for name, field in fields.items():
        zs = rng.uniform(-0.85, 0.85, 400) + 1j * rng.uniform(-0.85, 0.85, 400)
        fx = (field_values(field, zs + h) - field_values(field, zs - h)) / (2 * h)
        fy = (field_values(field, zs + 1j * h) - field_values(field, zs - 1j * h)) / (2 * h)
        dz = 0.5 * (fx - 1j * fy)
        dzbar = 0.5 * (fx + 1j * fy)
        ratio = np.abs(dzbar) / np.abs(dz)
        out[name] = {
            "median_|dF/dzbar|/|dF/dz|": round(float(np.median(ratio)), 3),
            "fraction_orientation_reversing (Jacobian < 0)": round(float(np.mean(ratio > 1)), 3),
            "holomorphic_would_give": "ratio 0 everywhere and Jacobian >= 0",
        }
    return out


# ---------------------------------------------------------------------------
# 4. Zeros: vortex census and degree identity
# ---------------------------------------------------------------------------
def zero_census(field: ct.EulerPhasorField, half_width: float, g: int) -> dict:
    xs = np.linspace(-half_width, half_width, g)
    X, Y = np.meshgrid(xs, xs)
    Z = X + 1j * Y
    F = field_values(field, Z).reshape(Z.shape)
    d_right = np.angle(F[:, 1:] / F[:, :-1])
    d_up = np.angle(F[1:, :] / F[:-1, :])
    wind = (d_right[:-1, :] + d_up[:, 1:] - d_right[1:, :] - d_up[:, :-1]) / (2 * math.pi)
    wind = np.rint(wind).astype(int)
    centers = (Z[:-1, :-1] + Z[1:, 1:]) / 2
    return {
        "grid": g,
        "max_edge_phase_step_rad": float(max(np.abs(d_right).max(), np.abs(d_up).max())),
        "wind": wind,
        "centers": centers,
        "n_plus": int((wind > 0).sum()),
        "n_minus": int((wind < 0).sum()),
        "higher_charge_cells": int((np.abs(wind) > 1).sum()),
    }


def audit_zeros(field: ct.EulerPhasorField, label: str, n_tiles: int = 16) -> dict:
    hw = 0.9
    coarse = zero_census(field, hw, 361)
    fine = zero_census(field, hw, 721)
    square = [complex(-hw, -hw), complex(hw, -hw), complex(hw, hw), complex(-hw, hw)]
    w_square, _ = boundary_winding(field, square, 20000)

    # Jacobian sign at each located zero must match its plaquette charge.
    idx = np.argwhere(fine["wind"] != 0)
    zc = fine["centers"][fine["wind"] != 0]
    charges = fine["wind"][fine["wind"] != 0]
    h = 1e-6
    fx = (field_values(field, zc + h) - field_values(field, zc - h)) / (2 * h)
    fy = (field_values(field, zc + 1j * h) - field_values(field, zc - 1j * h)) / (2 * h)
    jac = 0.25 * (np.abs(fx - 1j * fy) ** 2 - np.abs(fx + 1j * fy) ** 2)
    agree = float(np.mean(np.sign(jac) == np.sign(charges))) if len(charges) else 1.0

    # Rounded assignments at the zeros versus everywhere.
    unsat_at_zeros = unsat_counts_from_rank_bits(field, zc)
    sample = fine["centers"].ravel()[:: 7]
    unsat_all = unsat_counts_from_rank_bits(field, sample)
    absF = np.abs(field_values(field, sample))
    energy = []
    for z in sample[:: 5]:
        e, _ = field.eval_energy_and_grad_theta(field.phases_at(complex(z)))
        energy.append(e)

    # The solver's own tiles: dense winding versus the zeros they contain.
    tiles = []
    cx, cy = fine["centers"].real, fine["centers"].imag
    for k, (cz, tr, orient) in enumerate(ct.chiral_spectre_mesh_centers(n_tiles, radius=0.88, twist=ct.PINWHEEL_ANGLE)):
        poly = ct.spectre_prototile_vertices(scale=tr, angle=orient, center=cz)
        w6, jump6 = boundary_winding(field, poly, 6)
        w_dense, _ = boundary_winding(field, poly, 400)
        inside = points_in_polygon(cx, cy, poly)
        plus = int(((fine["wind"] > 0) & inside).sum())
        minus = int(((fine["wind"] < 0) & inside).sum())
        tiles.append({
            "tile": k,
            "W_6_samples_per_edge": round(w6, 3),
            "W_dense": round(w_dense, 3),
            "vortices_inside": plus,
            "antivortices_inside": minus,
        })
    tiles_zero_w_with_zeros = sum(1 for t in tiles if round(t["W_dense"]) == 0 and t["vortices_inside"] + t["antivortices_inside"] > 0)
    return {
        "instance": label,
        "region": f"square [-{hw}, {hw}]^2",
        "vortices(+1)": fine["n_plus"],
        "antivortices(-1)": fine["n_minus"],
        "higher_charge_cells": fine["higher_charge_cells"],
        "coarse_grid_counts(+,-)": [coarse["n_plus"], coarse["n_minus"]],
        "max_edge_phase_step_rad_fine_grid": round(fine["max_edge_phase_step_rad"], 3),
        "dense_boundary_winding_of_square": round(w_square, 4),
        "degree_identity_holds (W = N+ - N-)": int(round(w_square)) == fine["n_plus"] - fine["n_minus"],
        "jacobian_sign_matches_charge_fraction": round(agree, 4),
        "unsat_at_zeros": {
            "min": int(unsat_at_zeros.min()) if len(unsat_at_zeros) else None,
            "mean": round(float(unsat_at_zeros.mean()), 2) if len(unsat_at_zeros) else None,
            "satisfying_zeros": int((unsat_at_zeros == 0).sum()),
        },
        "unsat_over_region": {
            "min": int(unsat_all.min()),
            "mean": round(float(unsat_all.mean()), 2),
        },
        "spearman(|F|, unsat count)": round(spearman(absF, unsat_all), 3),
        "spearman(E(theta(z)), unsat count)": round(spearman(np.array(energy), unsat_all[:: 5]), 3),
        "solver_tiles": tiles,
        "tiles_with_W_0_that_contain_zeros": tiles_zero_w_with_zeros,
        "tiles_with_nonzero_W_at_6_samples": sum(1 for t in tiles if round(t["W_6_samples_per_edge"]) != 0),
    }


def audit_zeros_across_formulas(quick: bool) -> dict:
    """Do the zeros of F (what the winding screen looks for) sit at good assignments?"""
    count = 12 if quick else 40
    g = 241
    xs = np.linspace(-0.88, 0.88, g)
    X, Y = np.meshgrid(xs, xs)
    rows, at_zero, region_mean, region_min = [], [], [], []
    for t in range(count):
        cnf = ct.generate_random_3sat(n_vars=35, alpha=4.26, seed=7700 + t, planted=True)
        field = ct.EulerPhasorField(cnf, freq_scale=1.0, canonicalize=True, seed=t)
        census = zero_census(field, 0.88, g)
        inside = np.abs(census["centers"]) < 0.88
        mask = (census["wind"] != 0) & inside
        zc = census["centers"][mask]
        everywhere = unsat_counts_from_rank_bits(field, (X + 1j * Y).ravel()[np.abs((X + 1j * Y).ravel()) < 0.88])
        u0 = unsat_counts_from_rank_bits(field, zc) if len(zc) else np.array([], dtype=int)
        rows.append({
            "formula": t,
            "vortices": int(((census["wind"] > 0) & inside).sum()),
            "antivortices": int(((census["wind"] < 0) & inside).sum()),
        })
        at_zero += [int(u) for u in u0]
        region_mean.append(float(everywhere.mean()))
        region_min.append(int(everywhere.min()))
    v = sum(r["vortices"] for r in rows)
    a = sum(r["antivortices"] for r in rows)
    return {
        "formulas": count,
        "field_setting": "solver setting: freq_scale=1, canonicalized, jitter seed = run seed",
        "region": "disk |z| < 0.88 (where the solver's tiles live)",
        "total_vortices": v,
        "total_antivortices": a,
        "formulas_without_any_zero": sum(1 for r in rows if r["vortices"] + r["antivortices"] == 0),
        "mean_unsat_at_zeros": round(float(np.mean(at_zero)), 2) if at_zero else None,
        "mean_unsat_over_region": round(float(np.mean(region_mean)), 2),
        "mean_best_unsat_over_region": round(float(np.mean(region_min)), 2),
        "satisfying_zeros": sum(1 for u in at_zero if u == 0),
        "per_formula": rows,
    }


# ---------------------------------------------------------------------------
# 5. Plane reachability
# ---------------------------------------------------------------------------
def line_bound(field: ct.EulerPhasorField, radius: float) -> Tuple[int, int]:
    """Lines cos(theta_k) = 0 meeting the disk, and the arrangement face bound 1 + L + C(L, 2)."""
    total = 0
    kmag = np.hypot(field.k_real, field.k_imag)
    for c, k in zip(field.phase_offset, kmag):
        lo, hi = c - k * radius, c + k * radius
        total += max(0, math.floor((hi - math.pi / 2) / math.pi) - math.ceil((lo - math.pi / 2) / math.pi) + 1)
    return total, 1 + total + total * (total - 1) // 2


def disk_points(g: int, radius: float = 1.0) -> np.ndarray:
    xs = np.linspace(-radius, radius, g)
    X, Y = np.meshgrid(xs, xs)
    Z = (X + 1j * Y).ravel()
    return Z[np.abs(Z) <= radius]


def audit_reachability(quick: bool) -> dict:
    g = 900 if quick else 1400
    zs = disk_points(g)
    n_file, clauses = read_dimacs(MISSED_CNF)
    cnf = ct.load_dimacs_cnf(MISSED_CNF)
    field = ct.EulerPhasorField(cnf, freq_scale=1.0, canonicalize=True)
    codes = rank_codes(field, zs)
    distinct = np.unique(codes, axis=0)
    sols = all_solutions(n_file, clauses)
    sol_rank_codes = set()
    for s in sols:
        bits = code_to_bits(int(s), n_file)
        sol_rank_codes.add(sum(int(bits[orig]) << r for r, orig in enumerate(field.order)))
    reachable = set(int(c[0]) for c in distinct)
    L, faces = line_bound(field, 1.0)
    unsat = unsat_counts_from_rank_bits(field, zs)
    missed = {
        "n": n_file,
        "assignments_total": 2 ** n_file,
        "distinct_assignments_in_unit_disk": int(len(distinct)),
        "line_arrangement_upper_bound": faces,
        "lines_cos_theta_zero_in_disk": L,
        "satisfying_assignments_total": int(len(sols)),
        "satisfying_assignments_reachable_in_disk": len(reachable & sol_rank_codes),
        "min_unsat_in_disk": int(unsat.min()),
    }
    scaling = []
    for n in (6, 10, 16, 24, 35, 50, 100):
        inst = ct.generate_random_3sat(n_vars=n, alpha=4.2, seed=4242 + n, planted=True)
        row = {"n": n, "log2_assignments": n}
        for tag, fs in (("freq_scale_1", 1.0), ("freq_scale_sqrt_n_over_2.5", math.sqrt(n) / 2.5)):
            f = ct.EulerPhasorField(inst, freq_scale=fs, canonicalize=True)
            d = len(np.unique(rank_codes(f, zs), axis=0))
            L, faces = line_bound(f, 1.0)
            row[tag] = {"distinct_in_disk": int(d), "upper_bound": int(faces), "lines": int(L)}
        scaling.append(row)
    needed = {}
    for n in (18, 35, 50, 100):
        # smallest L with 1 + L + L(L-1)/2 >= 2^n, i.e. L(L+1) >= 2^(n+1) - 2
        target = 2 ** (n + 1) - 2
        L = max(1, math.isqrt(target) - 2)
        while L * (L + 1) < target:
            L += 1
        # n line families share L lines, so some |k_j| must reach (pi / 2R)(L/n - 1)
        needed[str(n)] = {
            "lines_needed_for_2^n_faces": L,
            "max_wavevector_needed_unit_disk": float(math.pi / 2 * (L / n - 1)),
        }
    return {
        "missed_instance_unit_disk": missed,
        "scaling_planted_alpha_4.2": scaling,
        "minimum_frequency_for_full_coverage": needed,
        "note": "Faces of an arrangement of L lines number at most 1 + L + C(L,2), so a 2-D slice of the n-torus shows O(L^2) assignments.",
    }


# ---------------------------------------------------------------------------
# 6. Ablations and controls
# ---------------------------------------------------------------------------
SolverFn = Callable[[ct.CNFInstance, int, int], hs.SolveResult]


def solver_table(n_tiles: int) -> Dict[str, SolverFn]:
    H = hs.solve_holomorphic_spectre
    return {
        "holomorphic_spectre (delivered)": lambda c, b, s: H(c, max_steps=b, n_spectre_tiles=n_tiles, seed=s),
        "  ablation: random starts, no tiles/Mobius": lambda c, b, s: H(c, max_steps=b, n_spectre_tiles=n_tiles, seed=s, init="random"),
        "  ablation: no gradient steps": lambda c, b, s: H(c, max_steps=b, n_spectre_tiles=n_tiles, seed=s, use_gradient=False),
        "  ablation: PRNG instead of logistic map": lambda c, b, s: H(c, max_steps=b, n_spectre_tiles=n_tiles, seed=s, driver="prng"),
        "  variant: true Tile(1,1) contour": lambda c, b, s: H(c, max_steps=b, n_spectre_tiles=n_tiles, seed=s, polygon="tile11"),
        "  variant: no canonicalization": lambda c, b, s: H(c, max_steps=b, n_spectre_tiles=n_tiles, seed=s, canonicalize=False),
        "golden_sweep (delivered control)": lambda c, b, s: hs.solve_golden_sweep(c, max_steps=b, seed=s),
        "uniform_random_walk (delivered control)": lambda c, b, s: hs.solve_uniform_random(c, max_steps=b, seed=s),
        "walksat p=0.5 (Exp 03 control rule)": lambda c, b, s: hs.solve_walksat(c, max_steps=b, seed=s),
    }


def run_raw(jobs: List[Tuple[str, ct.CNFInstance, int]], budget: int, n_tiles: int) -> Dict[str, List[hs.SolveResult]]:
    table = solver_table(n_tiles)
    raw: Dict[str, List[hs.SolveResult]] = {name: [] for name in table}
    for name, fn in table.items():
        for _, cnf, seed in jobs:
            raw[name].append(fn(cnf, budget, seed))
    return raw


def summarize(raw: Dict[str, List[hs.SolveResult]], budget: int, keep: Sequence[bool] = None) -> Dict[str, dict]:
    if keep is not None:
        raw = {name: [r for r, k in zip(res, keep) if k] for name, res in raw.items()}
    ref = raw["holomorphic_spectre (delivered)"]
    summary = {}
    for name, results in raw.items():
        solved = [r.solved for r in results]
        k = int(sum(solved))
        censored = [r.steps if r.solved else budget + 1 for r in results]
        row = {
            "solved": k,
            "runs": len(results),
            "solve_rate": round(k / len(results), 3),
            "wilson_95": wilson(k, len(results)),
            "median_steps_failures_as_budget+1": float(np.median(censored)),
            "mean_steps_when_solved": round(float(np.mean([r.steps for r in results if r.solved])), 1) if k else None,
            "mean_ms": round(float(np.mean([r.elapsed_ms for r in results])), 2),
        }
        if name != "holomorphic_spectre (delivered)":
            only_this = sum(1 for a, b in zip(results, ref) if a.solved and not b.solved)
            only_ref = sum(1 for a, b in zip(results, ref) if b.solved and not a.solved)
            fewer = sum(1 for a, b in zip(censored, [r.steps if r.solved else budget + 1 for r in ref]) if a < b)
            more = sum(1 for a, b in zip(censored, [r.steps if r.solved else budget + 1 for r in ref]) if a > b)
            row["vs_delivered"] = {
                "solved_only_by_this": only_this,
                "solved_only_by_delivered": only_ref,
                "runs_with_fewer_steps": fewer,
                "runs_with_more_steps": more,
                "sign_test_p_steps": sign_test_p(fewer, more),
            }
        diags = [r.diagnostics for r in results if r.diagnostics]
        if diags and "explicit_phase_flips" in diags[0]:
            flips = sum(d["explicit_phase_flips"] for d in diags)
            changes = sum(d["rounded_bit_changes"] for d in diags)
            row["explicit_phase_flips_total"] = flips
            row["rounded_bit_changes_total"] = changes
            row["winding_hits_per_run"] = round(float(np.mean([d["winding_hits"] for d in diags])), 3)
        summary[name] = row
    return summary


def audit_ablations(quick: bool) -> dict:
    instances = [i for i in exp03_instances() if i["sat_label"] == "SAT"]
    seeds = (0, 1) if quick else (0, 1, 2)
    jobs_all = [(i["split"], instance_to_cnf(i), s) for i in instances for s in seeds]
    t0 = time.perf_counter()
    raw03 = run_raw(jobs_all, budget=1000, n_tiles=12)
    exp03_all = summarize(raw03, 1000)
    exp03_test = summarize(raw03, 1000, [j[0] == "test" for j in jobs_all])
    by_family = {
        fam: summarize(raw03, 1000, [j[1].name.split("_")[2] == fam for j in jobs_all])
        for fam in ("planted", "random", "geometric")
    }
    t1 = time.perf_counter()
    count = 40 if quick else 100
    jobs_b = [("exp04", ct.generate_random_3sat(n_vars=35, alpha=4.26, seed=7700 + t, planted=True), t) for t in range(count)]
    exp04_family = summarize(run_raw(jobs_b, budget=1600, n_tiles=12), 1600)
    t2 = time.perf_counter()
    return {
        "exp03_sat_formulas_all": {
            "formulas": len(instances), "seeds": list(seeds), "budget_steps": 1000, "results": exp03_all,
            "context": "Experiment 03 (held-out, 1000-flip cap): WalkSAT 129/129, frozen pinwheel 120/129, Breakout 116/129",
        },
        "exp03_sat_formulas_heldout_test": {"budget_steps": 1000, "results": exp03_test},
        "exp03_by_family": by_family,
        "exp04_planted_n35_alpha4.26": {
            "instances": count, "budget_steps": 1600, "results": exp04_family,
            "note": "Instances 0-11 are exactly Experiment B's twelve formulas.",
        },
        "seconds": {"exp03": round(t1 - t0, 1), "exp04_family": round(t2 - t1, 1)},
    }


# ---------------------------------------------------------------------------
# 7. Permutation audit
# ---------------------------------------------------------------------------
def canonical_invariance(cnf: ct.CNFInstance, perms: int) -> float:
    """Fraction of random variable relabelings for which the canonical order is exactly invariant."""
    base = ct.spectral_canonical_order(cnf)
    same = 0
    for p in range(perms):
        pcnf, perm = ct.permute_cnf(cnf, seed=900 + p)
        order = ct.spectral_canonical_order(pcnf)
        # perm maps new index -> old index (see permute_cnf); express order in old labels
        mapped = [perm[v] for v in order]
        same += int(mapped == base)
    return same / perms


def audit_permutations(quick: bool) -> dict:
    bases = 8 if quick else 20
    rows = []
    pooled_raw, pooled_can, pooled_full = [], [], []
    for b in range(bases):
        base = ct.generate_random_3sat(n_vars=30, alpha=4.1, seed=31415 + 17 * b, planted=True)
        raw, can, full = [], [], []
        for p in range(10):
            pcnf, _ = ct.permute_cnf(base, seed=500 + p)
            raw.append(hs.solve_holomorphic_spectre(pcnf, max_steps=1100, canonicalize=False, seed=0).steps)
            can.append(hs.solve_holomorphic_spectre(pcnf, max_steps=1100, canonicalize=True, seed=0).steps)
            full.append(hs.solve_holomorphic_spectre(pcnf, max_steps=1100, canonicalize=True, seed=0,
                                                     canonical_clauses=True).steps)
        pooled_raw += raw
        pooled_can += can
        pooled_full += full
        rows.append({
            "base_seed": 31415 + 17 * b,
            "mean_steps_raw": round(float(np.mean(raw)), 1),
            "mean_steps_canonical": round(float(np.mean(can)), 1),
            "mean_steps_canonical_vars_and_clauses": round(float(np.mean(full)), 1),
            "std_raw": round(float(np.std(raw)), 1),
            "std_canonical": round(float(np.std(can)), 1),
            "std_canonical_vars_and_clauses": round(float(np.std(full)), 1),
            "distinct_canonical_step_counts": len(set(can)),
            "distinct_step_counts_vars_and_clauses": len(set(full)),
            "canonical_order_invariance_rate": canonical_invariance(base, 10),
        })
    lower_mean = sum(1 for r in rows if r["mean_steps_canonical"] < r["mean_steps_raw"])
    higher_mean = sum(1 for r in rows if r["mean_steps_canonical"] > r["mean_steps_raw"])
    lower_std = sum(1 for r in rows if r["std_canonical"] < r["std_raw"])
    higher_std = sum(1 for r in rows if r["std_canonical"] > r["std_raw"])
    return {
        "base_formulas": bases,
        "permutations_per_formula": 10,
        "per_formula": rows,
        "pooled": {
            "mean_raw": round(float(np.mean(pooled_raw)), 1),
            "mean_canonical": round(float(np.mean(pooled_can)), 1),
            "median_raw": float(np.median(pooled_raw)),
            "median_canonical": float(np.median(pooled_can)),
            "std_raw": round(float(np.std(pooled_raw)), 1),
            "std_canonical": round(float(np.std(pooled_can)), 1),
            "mean_canonical_vars_and_clauses": round(float(np.mean(pooled_full)), 1),
            "median_canonical_vars_and_clauses": float(np.median(pooled_full)),
        },
        "formulas_fully_invariant_with_clause_canonicalization": sum(
            1 for r in rows if r["distinct_step_counts_vars_and_clauses"] == 1),
        "formulas_where_canonical_mean_is_lower": lower_mean,
        "formulas_where_canonical_mean_is_higher": higher_mean,
        "sign_test_p_mean": sign_test_p(lower_mean, higher_mean),
        "formulas_where_canonical_std_is_lower": lower_std,
        "formulas_where_canonical_std_is_higher": higher_std,
        "sign_test_p_std": sign_test_p(lower_std, higher_std),
        "mean_canonical_order_invariance_rate": round(float(np.mean([r["canonical_order_invariance_rate"] for r in rows])), 3),
    }


# ---------------------------------------------------------------------------
# 8. The missed instance
# ---------------------------------------------------------------------------
def audit_missed(quick: bool) -> dict:
    n, clauses = read_dimacs(MISSED_CNF)
    sols = all_solutions(n, clauses)
    bits = np.array([[(int(s) >> i) & 1 for i in range(n)] for s in sols])
    ham = (bits[:, None, :] != bits[None, :, :]).sum(-1)
    # single-linkage clusters at Hamming distance 1 (one-flip connectivity)
    parent = list(range(len(sols)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(sols)):
        for j in range(i + 1, len(sols)):
            if ham[i, j] == 1:
                parent[find(i)] = find(j)
    clusters: Dict[int, List[str]] = {}
    for i, s in enumerate(sols):
        clusters.setdefault(find(i), []).append("".join(str(b) for b in bits[i]))
    match = None
    for inst in exp03_instances():
        if int(inst["n"]) == n and [list(c) for c in inst["clauses"]] == clauses:
            match = {k: inst[k] for k in ("id", "family", "split", "n", "m", "sat_label")}
    cnf = ct.load_dimacs_cnf(MISSED_CNF)
    seeds = list(range(30 if quick else 100))
    raw = run_raw([("missed", cnf, s) for s in seeds], budget=4000, n_tiles=16)
    runs = {}
    for name, res in raw.items():
        k = sum(r.solved for r in res)
        runs[name] = {
            "solved": k,
            "runs": len(res),
            "median_steps": float(np.median([r.steps for r in res])),
            "max_steps": int(max(r.steps for r in res)),
            "distinct_witnesses_found": len({tuple(r.witness) for r in res if r.solved}),
        }
    return {
        "file": os.path.relpath(MISSED_CNF, REPO_DIR),
        "n_vars": n,
        "m_clauses": len(clauses),
        "alpha": round(len(clauses) / n, 4),
        "matches_exp03_dataset_instance": match,
        "exact_satisfying_assignments": int(len(sols)),
        "solutions_x1_to_xn": ["".join(str(b) for b in row) for row in bits],
        "one_flip_clusters": [sorted(v) for v in clusters.values()],
        "pairwise_hamming_min_max": [int(ham[np.triu_indices(len(sols), 1)].min()), int(ham.max())],
        "runs_budget_4000_tiles_16": runs,
        "note": "The file name's 35 is the Experiment 03 instance id; the formula has 18 variables, so Experiment B (n=35) does not mirror its size.",
    }


# ---------------------------------------------------------------------------
# 9. Experiment C revisited
# ---------------------------------------------------------------------------
def audit_nyquist(quick: bool) -> dict:
    trials = 12 if quick else 30
    poly = ct.spectre_prototile_vertices(scale=0.35, angle=ct.PINWHEEL_ANGLE, center=0.15 + 0.1j)
    rows = []
    for n in (6, 10, 16, 24, 35, 50):
        row = {"n": n}
        for tag, fs in (("delivered freq sqrt(n)/2.5", math.sqrt(n) / 2.5), ("fixed freq 1.0", 1.0)):
            jumps, variation, aliased, wind_err = [], [], 0, 0
            for t in range(trials):
                cnf = ct.generate_random_3sat(n_vars=n, alpha=4.2, seed=900 + n * 10 + t, planted=True)
                field = ct.EulerPhasorField(cnf, freq_scale=fs, canonicalize=True)
                rep = ct.cauchy_winding_number(field, poly, samples_per_edge=12)
                dense, _ = boundary_winding(field, poly, 1500)
                jumps.append(rep.max_step_phase_jump)
                variation.append(rep.total_abs_phase_variation)
                aliased += int(rep.nyquist_aliased)
                wind_err += int(rep.rounded_winding != int(round(dense)))
            row[tag] = {
                "mean_max_step_jump_rad": round(float(np.mean(jumps)), 3),
                "mean_total_phase_variation_rad": round(float(np.mean(variation)), 2),
                "aliasing_rate": round(aliased / trials, 3),
                "aliasing_wilson_95": wilson(aliased, trials),
                "winding_wrong_at_12_samples": wind_err,
            }
        rows.append(row)
    return {"trials_per_n": trials, "rows": rows}


# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="smaller trial counts")
    args = parser.parse_args()
    t_start = time.perf_counter()
    results = {
        "audit": "Experiment 04 independent audit",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "quick_mode": args.quick,
    }
    steps = [
        ("0_reproduction", audit_reproduction),
        ("1_constants", audit_constants),
        ("2_geometry", audit_geometry),
    ]
    for key, fn in steps:
        print(f"[audit] {key} ...", flush=True)
        results[key] = fn()

    missed = ct.load_dimacs_cnf(MISSED_CNF)
    planted35 = ct.generate_random_3sat(n_vars=35, alpha=4.26, seed=7700, planted=True)
    fields = {
        "missed_sat_instance_35 (n=18)": ct.EulerPhasorField(missed, freq_scale=1.0, canonicalize=True),
        "Experiment B formula 0 (n=35)": ct.EulerPhasorField(planted35, freq_scale=1.0, canonicalize=True),
    }
    print("[audit] 3_holomorphy ...", flush=True)
    results["3_holomorphy"] = audit_holomorphy(fields)
    print("[audit] 4_zeros ...", flush=True)
    results["4_zeros"] = {name: audit_zeros(f, name) for name, f in fields.items()}
    results["4_zeros"]["across_exp04_formulas"] = audit_zeros_across_formulas(args.quick)
    print("[audit] 5_reachability ...", flush=True)
    results["5_reachability"] = audit_reachability(args.quick)
    print("[audit] 6_ablations ...", flush=True)
    results["6_ablations"] = audit_ablations(args.quick)
    print("[audit] 7_permutations ...", flush=True)
    results["7_permutations"] = audit_permutations(args.quick)
    print("[audit] 8_missed_instance ...", flush=True)
    results["8_missed_instance"] = audit_missed(args.quick)
    print("[audit] 9_nyquist_revisited ...", flush=True)
    results["9_nyquist_revisited"] = audit_nyquist(args.quick)
    results["seconds_total"] = round(time.perf_counter() - t_start, 1)

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    print(f"[audit] wrote {OUT_PATH} in {results['seconds_total']} s")


if __name__ == "__main__":
    main()
