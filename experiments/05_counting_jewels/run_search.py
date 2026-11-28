"""
run_search.py -- Experiment 05 search experiments. Writes search_results.json after each section.

S1  Deterministic starts by conditional expectations (Johnson) versus random starts, under WalkSAT.
S2  Transformations that keep the solutions (gauge, duplicated clauses, implied resolvents) or the
    solution count (Tovey splitting, definitional extension), under WalkSAT.
S3  The XOR-SAT look-alike: planted 3-XOR-SAT, which Gaussian elimination solves at once.

    python experiments/05_counting_jewels/run_search.py
"""

from __future__ import annotations

import json
import os
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import counting as C  # noqa: E402
import formulas as F  # noqa: E402
import search as S  # noqa: E402

OUT = os.path.join(HERE, "search_results.json")
results: dict = {}


def save():
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)


def summary(runs, budget):
    """runs: list of (solved, flips)."""
    k = sum(1 for s, _ in runs if s)
    censored = [f if s else budget + 1 for s, f in runs]
    return {"solved": k, "runs": len(runs), "median_flips": statistics.median(censored),
            "mean_flips_solved": round(statistics.mean([f for s, f in runs if s]), 1) if k else None}


def families():
    exp03 = [(f"exp03_{i['id']}", i["n"], [tuple(c) for c in i["clauses"]])
             for i in F.exp03_instances() if i["sat_label"] == "SAT"]
    n35 = [(f"planted35_{t}", 35, F.random_3sat(35, 149, seed=7700 + t, planted=True)[0]) for t in range(100)]
    n200 = [(f"planted200_{t}", 200, F.random_3sat(200, 840, seed=9900 + t, planted=True)[0]) for t in range(20)]
    return {"exp03_sat (89, n 12-32)": exp03, "planted n=35 m=149 (100)": n35, "planted n=200 m=840 (20)": n200}


def s1_starts():
    budget, seeds = 100_000, (0, 1, 2)
    out = {}
    for fam, forms in families().items():
        rows = {}
        for name, n, cl in forms:
            ords = S.orders(n, cl)
            starts = {"random start": None}
            for oname in ("index", "degree", "golden", "random"):
                starts[f"conditional expectations, {oname} order"] = S.conditional_expectation_start(n, cl, ords[oname])
            for sname, st in starts.items():
                r = rows.setdefault(sname, {"runs": [], "start_violated": [], "solved_by_start": 0})
                for s in seeds:
                    start = st if st is not None else S.random_start(n, s)
                    v0 = F.violated(cl, start)
                    ok, flips, a = S.walksat(n, cl, start, budget, seed=s)
                    assert not ok or F.satisfies(cl, a)
                    r["runs"].append((ok, flips))
                    r["start_violated"].append(v0)
                if st is not None and F.violated(cl, st) == 0:
                    r["solved_by_start"] += 1
        out[fam] = {k: {**summary(v["runs"], budget),
                        "mean_start_violated": round(statistics.mean(v["start_violated"]), 2),
                        "formulas_solved_by_start_alone": v["solved_by_start"] if k != "random start" else None}
                    for k, v in rows.items()}
        out[fam]["budget_flips"] = budget
        print(f"[S1] {fam} done", flush=True)
    results["S1_conditional_expectation_starts"] = out
    save()


def s2_transformations():
    budget = 50_000
    fams = families()
    sets = {"exp03_sat (89, n 12-32), seeds 0-1": (fams["exp03_sat (89, n 12-32)"], (0, 1)),
            "planted n=35 m=149 (100), seed 0": (fams["planted n=35 m=149 (100)"], (0,))}
    transforms = {
        "original": lambda n, cl, s: (n, list(cl), lambda a: a),
        "gauge (complement half the variables)": lambda n, cl, s: F.gauge(n, cl, 100 + s),
        "duplicate each clause 1-3 times": lambda n, cl, s: F.duplicate(n, cl, 100 + s),
        "add implied resolvents (length <= 3)": lambda n, cl, s: F.resolvents(n, cl),
        "Tovey split (count-preserving)": lambda n, cl, s: F.tovey_split(n, cl),
        "definitional extension (count-preserving)": lambda n, cl, s: F.definitional(n, cl),
    }
    out = {}
    for label, (forms, seeds) in sets.items():
        rows = {}
        sizes = {}
        for name, n, cl in forms:
            for tname, tr in transforms.items():
                for s in seeds:
                    n2, c2, lift = tr(n, cl, s)
                    ok, flips, a = S.walksat(n2, c2, S.random_start(n2, s), budget, seed=s)
                    if ok:
                        assert F.satisfies(cl, lift(a)), (name, tname)
                    rows.setdefault(tname, []).append((ok, flips))
                    sz = sizes.setdefault(tname, [0, 0, 0])
                    sz[0] += n2
                    sz[1] += len(c2)
                    sz[2] += 1
        out[label] = {t: {**summary(r, budget), "mean_variables": round(sizes[t][0] / sizes[t][2], 1),
                          "mean_clauses": round(sizes[t][1] / sizes[t][2], 1)} for t, r in rows.items()}
        out[label]["budget_flips"] = budget
        print(f"[S2] {label} done", flush=True)
    results["S2_transformations"] = out
    save()


def gf2_rank(n, eqs):
    rows = []
    for vs, _ in eqs:
        r = 0
        for v in vs:
            r ^= 1 << v
        rows.append(r)
    rank = 0
    for bit in range(n):
        piv = next((i for i in range(rank, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] >> bit & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def s3_xor():
    budget = 300_000
    rows = []
    for n in (20, 40, 60, 80, 100):
        xr, xj, sr, gauss, t_gauss = [], [], [], 0, 0.0
        counted = []
        for s in range(8):
            eqs, cl, hidden = F.planted_3xor(n, n, 500 + s)
            t0 = time.perf_counter()
            g = F.gauss_xor(n, eqs)
            t_gauss += time.perf_counter() - t0
            gauss += int(g is not None and F.satisfies(cl, g))
            xr.append(S.walksat(n, cl, S.random_start(n, s), budget, seed=s)[:2])
            js = S.conditional_expectation_start(n, cl, list(range(1, n + 1)))
            xj.append(S.walksat(n, cl, js, budget, seed=s)[:2])
            sat_cl, _ = F.random_3sat(n, 4 * n, seed=800 + s, planted=True)
            sr.append(S.walksat(n, sat_cl, S.random_start(n, s), budget, seed=s)[:2])
            if n <= 40:
                counted.append(C.count_solutions(n, cl, max_width=20) == 2 ** (n - gf2_rank(n, eqs)))
        rows.append({
            "n": n, "xor_equations": n, "cnf_clauses": 4 * n,
            "xor_walksat_random_start": summary(xr, budget),
            "xor_walksat_conditional_expectation_start": summary(xj, budget),
            "xor_gaussian_elimination_solved": gauss, "xor_gaussian_ms_total": round(1000 * t_gauss, 2),
            "planted_3sat_same_size_walksat": summary(sr, budget),
            "exact_count_equals_2^(n-rank)": (all(counted) if counted else None),
        })
        print(f"[S3] n={n} done", flush=True)
        results["S3_xor_lookalike"] = {"budget_flips": budget, "formulas_per_row": 8, "rows": rows}
        save()


if __name__ == "__main__":
    t0 = time.time()
    s1_starts()
    s3_xor()
    s2_transformations()
    results["seconds"] = round(time.time() - t0, 1)
    save()
    print(f"[search] wrote {OUT} in {results['seconds']} s")
