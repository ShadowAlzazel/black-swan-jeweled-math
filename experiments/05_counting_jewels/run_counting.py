"""
run_counting.py -- Experiment 05 counting experiments. Writes counting_results.json after each section.

C1  Validation: contraction against brute force; exact counts against Experiment 03's DPLL labels.
C2  Contraction width by formula family and size (the cost of counting is about 2^width).
C3  Tiles as contraction orders: the hat substitution hierarchy versus simple sweeps and heuristics.
C4  Zeros of the energy polynomial across clause density (brute force, n <= 20).
C5  One exact count far beyond brute force: a 169-variable formula living on the hat patch.

    python experiments/05_counting_jewels/run_counting.py
"""

from __future__ import annotations

import json
import math
import os
import random
import statistics
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import counting as C  # noqa: E402
import formulas as F  # noqa: E402

OUT = os.path.join(HERE, "counting_results.json")
results: dict = {}


def save():
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)


def c1_validation():
    n, cl = F.read_dimacs(os.path.join(F.EXP03, "missed_sat_instance_35.cnf"))
    t = time.time()
    brute = C.energy_counts_bruteforce(n, cl)
    tb = time.time() - t
    t = time.time()
    con = C.energy_polynomial(n, cl)
    tc = time.time() - t
    sliced = C.energy_polynomial(n, cl, max_width=8)
    t = time.time()
    sol = C.solve_by_counting(n, cl)
    ts = time.time() - t
    labels, brute_checked = 0, 0
    t = time.time()
    for inst in F.exp03_instances():
        c = [tuple(x) for x in inst["clauses"]]
        N0 = C.count_solutions(inst["n"], c)
        labels += int(("SAT" if N0 > 0 else "UNSAT") == inst["sat_label"])
        if inst["n"] <= 18:
            assert C.energy_counts_bruteforce(inst["n"], c)[0] == N0
            brute_checked += 1
    results["C1_validation"] = {
        "missed_instance_N_k": {str(k): v for k, v in enumerate(brute) if v},
        "contraction_equals_brute_force": con["N"] == brute,
        "contraction_width": con["width"],
        "sliced_contraction_equals_brute_force (width cap 8)": sliced["N"] == brute,
        "slices_pinned_variables": sliced["pins"],
        "seconds_brute_vs_contraction": [round(tb, 2), round(tc, 2)],
        "count_guided_witness": "".join("1" if b else "0" for b in sol),
        "count_guided_witness_satisfies": F.satisfies(cl, sol),
        "count_guided_seconds": round(ts, 3),
        "exp03_labels_matched_by_exact_count": f"{labels}/96",
        "exp03_counts_cross_checked_by_brute_force": brute_checked,
        "exp03_seconds": round(time.time() - t, 1),
    }
    save()
    print("[C1] done", flush=True)


def c2_widths():
    rows = []

    def width_of(n, cl, heavy=True):
        adj = C.primal_graph(n, cl)
        w = C.induced_width(adj, C.greedy_order(adj, "min_degree", 0))
        if heavy and n <= 170:
            w = min(w, C.induced_width(adj, C.greedy_order(adj, "min_fill", 0)))
        return w

    for n in (20, 40, 60, 80, 100, 150, 200):
        for fam in ("random 3-SAT, alpha 4.26", "planted 3-SAT, alpha 4.26", "planted 3-XOR-SAT, 1 equation per variable"):
            ws = []
            for s in range(3):
                if fam.startswith("random"):
                    cl, _ = F.random_3sat(n, round(4.26 * n), seed=40 + s)
                elif fam.startswith("planted 3-SAT"):
                    cl, _ = F.random_3sat(n, round(4.26 * n), seed=60 + s, planted=True)
                else:
                    _, cl, _ = F.planted_3xor(n, n, seed=80 + s)
                ws.append(width_of(n, cl, heavy=n <= 100))
            rows.append({"family": fam, "n": n, "mean_width": round(statistics.mean(ws), 1),
                         "width_over_n": round(statistics.mean(ws) / n, 3)})
    for n in (25, 49, 100, 169, 256, 400):
        ws = []
        for s in range(3):
            cl, _ = F.tile_formula(F.grid_points(n), 4.2, seed=20 + s)
            ws.append(width_of(n, cl, heavy=n <= 169))
        rows.append({"family": "square-grid geometric, alpha 4.2", "n": n, "mean_width": round(statistics.mean(ws), 1),
                     "width_over_sqrt_n": round(statistics.mean(ws) / math.sqrt(n), 2)})
    for level, label in ((1, "hat-patch geometric (25 tiles), alpha 4.2"), (2, "hat-patch geometric (169 tiles), alpha 4.2")):
        pts, _ = F.hat_points(level)
        ws = [width_of(len(pts), F.tile_formula(pts, 4.2, seed=30 + s)[0]) for s in range(3)]
        rows.append({"family": label, "n": len(pts), "mean_width": round(statistics.mean(ws), 1),
                     "width_over_sqrt_n": round(statistics.mean(ws) / math.sqrt(len(pts)), 2)})
    results["C2_widths"] = rows
    save()
    print("[C2] done", flush=True)


def c3_tile_orders():
    pts, addrs = F.hat_points(2)
    n = len(pts)
    by = {}
    for s in range(5):
        cl, _ = F.tile_formula(pts, 4.2, seed=s)
        adj = C.primal_graph(n, cl)
        hier = sorted(range(1, n + 1), key=lambda v: (addrs[v - 1], v))
        orders = {
            "hat substitution hierarchy (address order)": hier,
            "hierarchy, reversed": hier[::-1],
            "left-to-right sweep": sorted(range(1, n + 1), key=lambda v: (pts[v - 1][0], pts[v - 1][1])),
            "min-degree heuristic": C.greedy_order(adj, "min_degree", 0),
            "min-fill heuristic": C.greedy_order(adj, "min_fill", 0),
            "random order": random.Random(s).sample(range(1, n + 1), n),
        }
        for name, o in orders.items():
            by.setdefault(name, []).append(C.induced_width(adj, o))
    results["C3_hat_orders"] = {"formulas": 5, "n": n, "mean_width": {k: round(statistics.mean(v), 1) for k, v in by.items()},
                                "all_widths": by}
    save()
    print("[C3] done", flush=True)


def c4_zeros():
    rows = []
    for n in (12, 16, 20):
        for alpha in (2.0, 3.0, 3.5, 4.0, 4.25, 4.5, 5.0, 6.0):
            stats = {"sat": 0, "N0": [], "ratio": [], "near": [], "near_negative_real": 0, "inside": [],
                     "min_angle": [], "degree": []}
            for s in range(12):
                cl, _ = F.random_3sat(n, round(alpha * n), seed=1000 * n + int(100 * alpha) + s)
                N = C.energy_counts_bruteforce(n, cl)
                d = max(k for k, v in enumerate(N) if v)
                zs = np.roots(np.array(N[:d + 1][::-1], dtype=float)) if d > 0 else np.array([])
                stats["degree"].append(d)
                stats["inside"].append(int(np.sum(np.abs(zs) < 1)))
                # A zero at the origin (exactly the unsatisfiable formulas) has no angle; skip it.
                pos = [abs(math.atan2(z.imag, z.real)) for z in zs if abs(z) > 1e-9]
                stats["min_angle"].append(min(pos) if pos else math.pi)
                if N[0] > 0:
                    stats["sat"] += 1
                    near = zs[np.argmin(np.abs(zs))]
                    stats["N0"].append(N[0])
                    stats["ratio"].append(N[0] / N[1] if N[1] else float("inf"))
                    stats["near"].append(abs(near))
                    stats["near_negative_real"] += int(abs(near.imag) < 1e-7 and near.real < 0)
            rows.append({
                "n": n, "alpha": alpha, "formulas": 12, "sat": stats["sat"],
                "unsat_zero_at_origin": 12 - stats["sat"],
                "median_solutions": statistics.median(stats["N0"]) if stats["N0"] else 0,
                "median_N0_over_N1": round(statistics.median(stats["ratio"]), 4) if stats["ratio"] else None,
                "median_nearest_zero_modulus": round(statistics.median(stats["near"]), 4) if stats["near"] else None,
                "nearest_zero_on_negative_real_axis": stats["near_negative_real"],
                "median_zeros_inside_unit_disk": statistics.median(stats["inside"]),
                "median_degree": statistics.median(stats["degree"]),
                "median_smallest_angle_to_positive_axis_deg": round(math.degrees(statistics.median(stats["min_angle"])), 1),
            })
        print(f"[C4] n={n} done", flush=True)
        results["C4_energy_zeros"] = rows
        save()


def c5_big_count():
    pts, _ = F.hat_points(2)
    n = len(pts)
    cl, hidden = F.tile_formula(pts, 4.2, seed=1)
    t = time.time()
    order, width = C.best_order(n, cl, tries=1)
    r = C.energy_polynomial(n, cl, L=2, order=order, max_width=21)
    results["C5_hat169_exact"] = {
        "n": n, "m": len(cl), "planted": True, "order_width": width, "pinned_variables": r["pins"],
        "sliced_width": r["width"], "primes": r["primes"],
        "N0_solutions": str(r["N"][0]), "N1_one_violated": str(r["N"][1]),
        "log2_N0": round(math.log2(r["N"][0]), 3) if r["N"][0] else None,
        "planted_assignment_is_counted": F.satisfies(cl, hidden),
        "seconds": round(time.time() - t, 1),
    }
    save()
    print("[C5] done", flush=True)


if __name__ == "__main__":
    t0 = time.time()
    steps = {"c1": c1_validation, "c2": c2_widths, "c3": c3_tile_orders, "c4": c4_zeros, "c5": c5_big_count}
    chosen = sys.argv[1:] or list(steps)
    if os.path.exists(OUT) and sys.argv[1:]:
        with open(OUT, "r", encoding="utf-8") as fh:
            results.update(json.load(fh))      # rerun some sections, keep the others
    for name in chosen:
        steps[name]()
    results["seconds"] = round(time.time() - t0, 1)
    save()
    print(f"[counting] wrote {OUT} in {results['seconds']} s")
