"""
run_levels.py -- Experiment 07: older experiments redone through the Fourier levels of the energy.

    python run_levels.py              # every section
    python run_levels.py f2 f3        # chosen sections

  F0  checks: closed-form coefficients against a Walsh-Hadamard transform, sphere averages and the XOR
      energy polynomial against brute force, the solution enumerator, the rewrites' forward maps
  F1  Experiment 05's conditional-expectation starts, redone per formula on Experiment 03's 89
      satisfiable formulas, and the levels of Experiment 05's rewrites
  F2  the hint ladder: planted, golden-hidden, deceptive, two-hidden and XOR formulas (Experiment 05's
      S1 and S3 redone with the hint moved level by level)
  F3  Experiment 05's energy polynomial redone for XOR: an exact signed weight enumerator
  F4  Experiments 02-03's clause weights redone on the look-alike and on golden-hidden formulas

Results go to levels_results.json, saved after each section.
"""

from __future__ import annotations

import json
import math
import multiprocessing as mp
import random
import statistics
import sys
import time

import numpy as np

import levels as lv
import counting
from formulas import (exp03_instances, planted_3xor, random_3sat, gauge, duplicate, resolvents,
                      tovey_split, definitional, satisfies, violated)
from search import walksat, conditional_expectation_start, orders, random_start

OUT = lv.HERE / "levels_results.json"
CAP = 300_000                          # WalkSAT flip budget, as in Experiment 05's S3


def save(section: str, data) -> None:
    results = json.loads(OUT.read_text()) if OUT.exists() else {}
    results[section] = data
    OUT.write_text(json.dumps(results, indent=2) + "\n")


def r(x, k=4):
    return None if x is None else round(float(x), k)


def spearman(a, b) -> float:
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    return float(ra @ rb / math.sqrt((ra @ ra) * (rb @ rb)))


def bits_of(mask: int, n: int):
    return [bool((mask >> i) & 1) for i in range(n)]


# ---------------------------------------------------------------------------
# F0: checks
# ---------------------------------------------------------------------------
def f0():
    t0 = time.time()
    out = {}
    cases = []
    n = 12
    cases.append(("random 3-SAT", n, random_3sat(n, 51, 3)[0]))
    for name, kind, q in (("planted", "q", 1.0), ("golden-hidden", "q", lv.Q_GOLDEN),
                          ("deceptive q = 0.5", "q", 0.5), ("two-hidden", "two", 1.0)):
        cases.append((name, n, lv.hidden_3sat(n, 51, 3, kind, q)[0]))
    eqs, xcl, _ = planted_3xor(n, n, 3)
    cases.append(("planted 3-XOR", n, xcl))
    inst = [d for d in exp03_instances() if d["n"] == 12 and d["sat_label"] == "SAT"][0]
    ic = [tuple(c) for c in inst["clauses"]]
    cases.append(("gauge", n, gauge(n, ic, 0)[1]))
    cases.append(("duplicate", n, duplicate(n, ic, 0)[1]))
    nr, rc, _ = resolvents(n, ic)
    cases.append(("resolvents", nr, rc))
    tiny = random_3sat(5, 5, 1)[0]              # 15 occurrences, so the split formula has 15 copies
    nt, tc, _ = tovey_split(5, tiny)
    cases.append(("Tovey split (n = 5, m = 5)", nt, tc))
    nd, dc, _ = definitional(5, tiny)
    cases.append(("definitional (n = 5, m = 5)", nd, dc))
    worst = {}
    for name, nn, cl in cases:
        w = lv.walsh_transform(lv.energy_vector(nn, cl))
        coef = lv.coefficients(cl)
        dense = np.zeros(2 ** nn)
        for S, c in coef.items():
            dense[sum(1 << i for i in S)] += c
        worst[name] = float(np.max(np.abs(dense - w)))
    out["coefficients_vs_walsh_transform_max_error"] = worst

    # sphere averages around a planted solution, against enumeration
    cl, A = lv.hidden_3sat(12, 51, 5, "q", lv.Q_GOLDEN)
    E = lv.energy_vector(12, cl)
    a_mask = sum(1 << i for i in range(12) if A[i])
    dist = np.array([bin(x ^ a_mask).count("1") for x in range(2 ** 12)])
    coef = lv.coefficients(cl)
    out["sphere_average_max_error"] = max(abs(float(E[dist == d].mean()) - lv.sphere_average(coef, A, d)) for d in range(13))

    # XOR energy polynomial against enumeration, satisfiable and not
    errs, kinds = [], []
    for seed in range(6):
        eqs, xcl, hidden = planted_3xor(13, 13, seed)
        for flip in (False, True):
            e2 = list(eqs)
            if flip:
                basis = lv.xor_kernel(e2)
                if not basis:
                    continue
                j = (basis[0][0] & -basis[0][0]).bit_length() - 1
                vs, t = e2[j]
                e2[j] = (vs, 1 - t)
            clauses = []
            for vs, t in e2:
                for bitsv in __import__("itertools").product((0, 1), repeat=3):
                    if sum(bitsv) % 2 != t:
                        clauses.append(tuple(v + 1 if b == 0 else -(v + 1) for v, b in zip(vs, bitsv)))
            Ev = lv.energy_vector(13, clauses).astype(int)
            brute = np.bincount(Ev, minlength=len(e2) + 1).tolist()
            got = lv.xor_energy_polynomial(13, e2)
            errs.append(int(np.max(np.abs(np.array(brute[:len(got["N"])]) - np.array(got["N"])))))
            kinds.append("unsat" if flip else "sat")
    out["xor_polynomial_vs_enumeration"] = {"cases": len(errs), "unsat_cases": kinds.count("unsat"), "max_error": max(errs)}

    # solution enumerator and rewrites' forward maps
    ok_count, ok_forward = 0, 0
    small = [d for d in exp03_instances() if d["n"] == 12 and d["sat_label"] == "SAT"][:6]
    for d in small:
        cl = [tuple(c) for c in d["clauses"]]
        sols = lv.all_solutions(12, cl)
        ok_count += len(sols) == int((lv.energy_vector(12, cl) == 0).sum())
        x = bits_of(sols[0], 12)
        nt, tc, lift_t = tovey_split(12, cl)
        nd, dc, lift_d = definitional(12, cl)
        xt, xd = lv.tovey_forward(12, cl, x), lv.definitional_forward(12, cl, x)
        ok_forward += (len(xt) == nt and satisfies(tc, xt) and lift_t(xt) == x
                       and len(xd) == nd and satisfies(dc, xd) and lift_d(xd) == x)
    big = [d for d in exp03_instances() if d["n"] == 32 and d["sat_label"] == "SAT"][:2]
    count_matches = []
    for d in big:
        cl = [tuple(c) for c in d["clauses"]]
        count_matches.append(len(lv.all_solutions(32, cl)) == counting.energy_polynomial(32, cl, L=1)["N"][0])
    out["solution_enumerator"] = {"n12_matches_enumeration": f"{ok_count}/{len(small)}",
                                  "n32_matches_exact_count": f"{sum(count_matches)}/{len(big)}"}
    out["rewrite_forward_maps_ok"] = f"{ok_forward}/{len(small)}"
    out["seconds"] = round(time.time() - t0, 1)
    save("F0_checks", out)
    print(json.dumps(out, indent=1))


# ---------------------------------------------------------------------------
# F1: Experiment 03's formulas, Johnson starts per formula, and Experiment 05's rewrites
# ---------------------------------------------------------------------------
def f1_instance(d):
    n = d["n"]
    cl = [tuple(c) for c in d["clauses"]]
    sols = lv.all_solutions(n, cl)
    X = lv.chi_matrix(sols, n)
    coef = lv.coefficients(cl)
    shares = lv.hint_shares(coef, X)
    c1 = lv.level1_cosine(coef, X)
    start = conditional_expectation_start(n, cl, orders(n, cl)["index"])
    flips_r, flips_j = [], []
    for s in range(10):
        ok, f, _ = walksat(n, cl, random_start(n, 7000 + 97 * d["id"] + s), 100_000, seed=s)
        flips_r.append(f if ok else None)
        ok, f, _ = walksat(n, cl, start, 100_000, seed=s)
        flips_j.append(f if ok else None)
    rec = {"id": d["id"], "family": d["family"], "n": n, "m": len(cl), "solutions": len(sols),
           "hint_shares": {k: r(v) for k, v in shares.items()}, "level1_cosine": r(c1),
           "johnson_start_violated": violated(cl, start),
           "median_flips_random": statistics.median([f if f is not None else 10 ** 9 for f in flips_r]),
           "median_flips_johnson": statistics.median([f if f is not None else 10 ** 9 for f in flips_j]),
           "unsolved_runs": sum(f is None for f in flips_r + flips_j)}
    # rewrites: levels at every solution, mapped forward
    rew = {}
    xs = [bits_of(s, n) for s in sols]
    for name, (nn, cc, fwd) in {
        "original": (n, cl, lambda x: x),
        "gauge (complement half the variables)": (lambda g: (g[0], g[1], g[2]))(gauge(n, cl, 0)),
        "duplicate each clause 1-3 times": (n, duplicate(n, cl, 0)[1], lambda x: x),
        "add implied resolvents (length <= 3)": (n, resolvents(n, cl)[1], lambda x: x),
        "definitional extension (count-preserving)": (definitional(n, cl)[0], definitional(n, cl)[1], lambda x: lv.definitional_forward(n, cl, x)),
        "Tovey split (count-preserving)": (tovey_split(n, cl)[0], tovey_split(n, cl)[1], lambda x: lv.tovey_forward(n, cl, x)),
    }.items():
        Xn = lv.chi_matrix([fwd(x) for x in xs], nn)
        cf = lv.coefficients(cc)
        sh = lv.hint_shares(cf, Xn)
        rew[name] = {"variables": nn, "clauses": len(cc), "hint_shares": {k: r(v) for k, v in sh.items()},
                     "level1_cosine": r(lv.level1_cosine(cf, Xn))}
    rec["rewrites"] = rew
    return rec


def f1():
    t0 = time.time()
    inst = [d for d in exp03_instances() if d["sat_label"] == "SAT"]
    with mp.Pool(4) as pool:
        recs = pool.map(f1_instance, inst, chunksize=2)
    gain = [math.log((x["median_flips_random"] + 1) / (x["median_flips_johnson"] + 1)) for x in recs]
    h1 = [x["hint_shares"].get(1, 0.0) or 0.0 for x in recs]
    h2 = [x["hint_shares"].get(2, 0.0) or 0.0 for x in recs]
    c1 = [x["level1_cosine"] for x in recs]
    summary = {
        "formulas": len(recs),
        "spearman_h1_vs_log_johnson_gain": r(spearman(h1, gain)),
        "spearman_level1_cosine_vs_log_johnson_gain": r(spearman(c1, gain)),
        "spearman_h2_vs_log_johnson_gain": r(spearman(h2, gain)),
        "spearman_h1_vs_log_random_start_flips": r(spearman(h1, [math.log(x["median_flips_random"] + 1) for x in recs])),
        "by_family": {},
    }
    for fam in sorted({x["family"] for x in recs}):
        sub = [x for x, g in zip(recs, gain) if x["family"] == fam]
        gs = [g for x, g in zip(recs, gain) if x["family"] == fam]
        summary["by_family"][fam] = {
            "formulas": len(sub),
            "mean_h1": r(statistics.mean(x["hint_shares"].get(1, 0.0) for x in sub)),
            "mean_h2": r(statistics.mean(x["hint_shares"].get(2, 0.0) for x in sub)),
            "mean_h3": r(statistics.mean(x["hint_shares"].get(3, 0.0) for x in sub)),
            "mean_level1_cosine": r(statistics.mean(x["level1_cosine"] for x in sub)),
            "median_flips_random": statistics.median(x["median_flips_random"] for x in sub),
            "median_flips_johnson": statistics.median(x["median_flips_johnson"] for x in sub),
            "median_log_gain": r(statistics.median(gs)),
        }
    stored = json.loads((lv.HERE.parent / "05_counting_jewels" / "search_results.json").read_text())
    s2 = stored["S2_transformations"]["exp03_sat (89, n 12-32), seeds 0-1"]
    table = {}
    for name in recs[0]["rewrites"]:
        rows = [x["rewrites"][name] for x in recs]
        key = name if name in s2 else None
        table[name] = {
            "mean_variables": r(statistics.mean(x["variables"] for x in rows), 1),
            "mean_h1": r(statistics.mean(x["hint_shares"].get(1, 0.0) for x in rows)),
            "mean_h2": r(statistics.mean(x["hint_shares"].get(2, 0.0) for x in rows)),
            "mean_h3": r(statistics.mean(x["hint_shares"].get(3, 0.0) for x in rows)),
            "mean_level1_cosine": r(statistics.mean(x["level1_cosine"] for x in rows)),
            "experiment05_median_flips": s2[key]["median_flips"] if key else None,
            "experiment05_solved": f'{s2[key]["solved"]}/{s2[key]["runs"]}' if key else None,
        }
    summary["rewrites_vs_experiment05_S2"] = table
    summary["seconds"] = round(time.time() - t0, 1)
    save("F1_experiment03_and_rewrites", {"summary": summary, "formulas": recs})
    print(json.dumps(summary, indent=1))


# ---------------------------------------------------------------------------
# F2: the hint ladder
# ---------------------------------------------------------------------------
FAMILIES = [("planted (q = 1)", "q", 1.0), ("golden-hidden (q = 1/phi)", "q", lv.Q_GOLDEN),
            ("deceptive (q = 0.5)", "q", 0.5), ("two-hidden (A and its complement)", "two", 1.0),
            ("planted 3-XOR", "xor", None)]
SIZES = (50, 100, 200, 400, 800)
PER_CELL = 20


def make_family(fam, n, seed):
    name, kind, q = fam
    if kind == "xor":
        eqs, cl, A = planted_3xor(n, n, 10_000 + seed)
        return cl, A
    return lv.hidden_3sat(n, round(4.26 * n), 10_000 + seed, kind, q)


def f2_task(args):
    fi, n, seed = args
    fam = FAMILIES[fi]
    cl, A = make_family(fam, n, seed)
    notA = [not b for b in A]
    coef = lv.coefficients(cl)
    shares = lv.hint_shares(coef, lv.chi_matrix([A], n))
    vote = lv.majority_vote(coef, n)
    targets = [A, notA] if fam[1] == "two" else [A]
    start = conditional_expectation_start(n, cl, orders(n, cl)["index"])
    ok_r, f_r, end_r = walksat(n, cl, random_start(n, 50_000 + seed), CAP, seed=seed)
    ok_j, f_j, end_j = walksat(n, cl, start, CAP, seed=seed)
    prof = None
    if seed == 0 and n == 200:
        half = lv.sphere_average(coef, A, n // 2)
        prof = [r(lv.sphere_average(coef, A, d) / half) for d in range(0, n + 1, n // 20)]
    return {"family": fam[0], "n": n, "seed": seed, "clauses": len(cl),
            "hint_shares": {k: r(v) for k, v in shares.items()},
            "majority_vote_agreement": r(max(lv.agreement(vote, t) for t in targets)),
            "johnson_start_violated": violated(cl, start),
            "johnson_start_agreement": r(max(lv.agreement(start, t) for t in targets)),
            "walksat_random": f_r if ok_r else None, "walksat_johnson": f_j if ok_j else None,
            "end_agreement_with_A_random": r(lv.agreement(end_r, A)), "end_agreement_with_A_johnson": r(lv.agreement(end_j, A)),
            "end_violated_random": violated(cl, end_r), "end_violated_johnson": violated(cl, end_j),
            "sphere_profile": prof}


def cell_median(vals):
    """Median with unsolved runs as +infinity; None when fewer than half are solved."""
    vs = sorted(v if v is not None else math.inf for v in vals)
    med = statistics.median(vs)
    return None if math.isinf(med) else med


def f2():
    t0 = time.time()
    tasks = [(fi, n, s) for fi in range(len(FAMILIES)) for n in SIZES for s in range(PER_CELL)]
    with mp.Pool(4) as pool:
        recs = pool.map(f2_task, tasks, chunksize=1)
    rows = []
    for fam in FAMILIES:
        for n in SIZES:
            sub = [x for x in recs if x["family"] == fam[0] and x["n"] == n]
            rows.append({
                "family": fam[0], "n": n, "formulas": len(sub), "clauses": sub[0]["clauses"],
                "mean_hint_shares": {k: r(statistics.mean(x["hint_shares"].get(k, 0.0) for x in sub)) for k in (1, 2, 3)},
                "mean_majority_vote_agreement": r(statistics.mean(x["majority_vote_agreement"] for x in sub)),
                "mean_johnson_start_violated": r(statistics.mean(x["johnson_start_violated"] for x in sub), 2),
                "mean_johnson_start_agreement": r(statistics.mean(x["johnson_start_agreement"] for x in sub)),
                "walksat_random_solved": sum(x["walksat_random"] is not None for x in sub),
                "walksat_random_median_flips": cell_median([x["walksat_random"] for x in sub]),
                "walksat_johnson_solved": sum(x["walksat_johnson"] is not None for x in sub),
                "walksat_johnson_median_flips": cell_median([x["walksat_johnson"] for x in sub]),
                "unsolved_runs_end_agreement_with_A": r(statistics.mean(
                    [x[f"end_agreement_with_A_{k}"] for x in sub for k in ("random", "johnson") if x[f"walksat_{k}"] is None]))
                    if any(x[f"walksat_{k}"] is None for x in sub for k in ("random", "johnson")) else None,
                "unsolved_runs_mean_end_violated": r(statistics.mean(
                    [x[f"end_violated_{k}"] for x in sub for k in ("random", "johnson") if x[f"walksat_{k}"] is None]), 2)
                    if any(x[f"walksat_{k}"] is None for x in sub for k in ("random", "johnson")) else None,
                "solved_runs_end_agreement_with_A": r(statistics.mean(
                    [x[f"end_agreement_with_A_{k}"] for x in sub for k in ("random", "johnson") if x[f"walksat_{k}"] is not None]))
                    if any(x[f"walksat_{k}"] is not None for x in sub for k in ("random", "johnson")) else None,
            })
    expected = {fam[0]: ({1: 0.0, 2: 0.0, 3: 1.0} if fam[1] == "xor" else
                         {k: r(v) for k, v in lv.expected_hint_shares(fam[1], fam[2]).items()}) for fam in FAMILIES}
    profiles = {x["family"]: x["sphere_profile"] for x in recs if x["sphere_profile"]}
    for x in recs:
        x.pop("sphere_profile", None)
    data = {"budget_flips": CAP, "alpha_3sat": 4.26, "xor_equations_per_variable": 1,
            "expected_hint_shares": expected, "rows": rows,
            "sphere_profiles_n200": {"distance_fraction": [i / 20 for i in range(21)], "relative_mean_energy": profiles},
            "seconds": round(time.time() - t0, 1), "formulas": recs}
    save("F2_hint_ladder", data)
    for row in rows:
        print(row["family"][:26].ljust(26), row["n"], row["mean_hint_shares"], "vote", row["mean_majority_vote_agreement"],
              "J-start", row["mean_johnson_start_agreement"], "WS rnd", row["walksat_random_solved"], row["walksat_random_median_flips"],
              "WS J", row["walksat_johnson_solved"], row["walksat_johnson_median_flips"],
              "stuck at agreement", row["unsolved_runs_end_agreement_with_A"], "solved at", row["solved_runs_end_agreement_with_A"])
    print("seconds", data["seconds"])


# ---------------------------------------------------------------------------
# F3: the XOR energy polynomial
# ---------------------------------------------------------------------------
def xor_cnf(eqs):
    import itertools
    out = []
    for vs, t in eqs:
        for bitsv in itertools.product((0, 1), repeat=3):
            if sum(bitsv) % 2 != t:
                out.append(tuple(v + 1 if b == 0 else -(v + 1) for v, b in zip(vs, bitsv)))
    return out


def f3():
    t0 = time.time()
    rows = []
    for n, m in ((20, 20), (30, 30), (40, 40), (40, 20), (60, 30)):
        for seed in range(4):
            eqs, cl, _ = planted_3xor(n, m, seed)
            for variant in ("planted", "one right side flipped"):
                e2 = list(eqs)
                if variant != "planted":
                    basis = lv.xor_kernel(e2)
                    if not basis:
                        continue
                    j = (basis[0][0] & -basis[0][0]).bit_length() - 1
                    vs, t = e2[j]
                    e2[j] = (vs, 1 - t)
                got = lv.xor_energy_polynomial(n, e2)
                exact = counting.energy_polynomial(n, xor_cnf(e2))
                N_exact = list(exact["N"]) + [0] * (len(got["N"]) - len(exact["N"]))
                mult = lv.zero_multiplicity_at_minus_one(got["N"])
                quotient = lv.divide_by_one_plus_z(got["N"], mult)
                top = max(i for i, c in enumerate(quotient) if c)
                rest = np.roots(list(reversed(quotient[:top + 1]))) if top > 0 else np.array([])
                roots = list(rest) + [-1.0] * mult
                rows.append({"n": n, "m": m, "seed": seed, "variant": variant, "rank": got["rank"],
                             "kernel_dim": got["kernel_dim"], "satisfiable": got["satisfiable"],
                             "N0": got["N"][0], "matches_contraction": N_exact[:len(got["N"])] == got["N"],
                             "contraction_width": exact["width"],
                             "zeros_at_minus_one": mult, "degree": len(roots),
                             "zeros_inside_unit_disk": int(sum(abs(z) < 1 - 1e-9 for z in rest)),
                             "zeros_at_origin": next(i for i, c in enumerate(got["N"]) if c),
                             "binomial": got["kernel_dim"] == 0})
    data = {"rows": rows,
            "all_match": all(x["matches_contraction"] for x in rows),
            "mean_share_of_zeros_at_minus_one": {f"n={n}, m={m}": r(statistics.mean(x["zeros_at_minus_one"] / x["degree"] for x in rows if x["n"] == n and x["m"] == m and x["variant"] == "planted"))
                                                 for n, m in ((20, 20), (30, 30), (40, 40), (40, 20), (60, 30))},
            "seconds": round(time.time() - t0, 1)}
    save("F3_xor_energy_polynomial", data)
    for x in rows:
        print(x)
    print("all match:", data["all_match"], data["mean_share_of_zeros_at_minus_one"], "seconds", data["seconds"])


# ---------------------------------------------------------------------------
# F4: clause weights on the look-alike and on golden-hidden formulas
# ---------------------------------------------------------------------------
def weight_schemes(cl, A, kind_of_formula, seed):
    rng = random.Random(seed)
    m = len(cl)
    t = [sum((l > 0) == A[abs(l) - 1] for l in c) for c in cl]
    schemes = {
        "equal weights": [1.0] * m,
        "golden boost 1 + (1 + cos 2 pi j phi)^2": [1 + (1 + math.cos(2 * math.pi * j * lv.PHI)) ** 2 for j in range(m)],
        "random in [0.5, 1.5]": [rng.uniform(0.5, 1.5) for _ in range(m)],
    }
    if kind_of_formula == "xor":
        schemes["pattern (1, phi, phi^2, phi^3) within each equation"] = [lv.PHI ** (j % 4) for j in range(m)]
    schemes["oracle control: weight 2 when every literal is true under A"] = [2.0 if tt == 3 else 1.0 for tt in t]
    return schemes


def f4():
    t0 = time.time()
    out = {}
    for label, kind in (("planted 3-XOR, n = 200, m = 200 (2,000 systems)", "xor"), ("golden-hidden 3-SAT, n = 200, m = 852 (400 formulas)", "golden")):
        acc = {}
        for seed in range(2000 if kind == "xor" else 400):
            if kind == "xor":
                eqs, cl, A = planted_3xor(200, 200, 20_000 + seed)
            else:
                cl, A = lv.hidden_3sat(200, 852, 20_000 + seed, "q", lv.Q_GOLDEN)
            X = lv.chi_matrix([A], 200)
            for name, w in weight_schemes(cl, A, kind, seed).items():
                coef = lv.coefficients(cl, w)
                sh = lv.hint_shares(coef, X)
                lw = lv.level_weights(coef)
                tot = sum(lw.values())
                acc.setdefault(name, []).append({
                    "h1": sh.get(1, 0.0), "h2": sh.get(2, 0.0), "h3": sh.get(3, 0.0),
                    "inner": float(-(lv.level1_vector(coef, 200) @ X[0])),
                    "cos": lv.level1_cosine(coef, X),
                    "vote": lv.agreement(lv.majority_vote(coef, 200), A),
                    "spectrum_level1_share": lw.get(1, 0.0) / tot})
        out[label] = {name: {**{k: r(statistics.mean(x[k] for x in v)) for k in v[0]},
                             **{f"{k}_standard_error": r(statistics.stdev(x[k] for x in v) / math.sqrt(len(v)), 5)
                                for k in ("inner", "cos", "vote")},
                             "formulas": len(v)}
                      for name, v in acc.items()}
    out["seconds"] = round(time.time() - t0, 1)
    save("F4_clause_weights", out)
    print(json.dumps(out, indent=1))


SECTIONS = {"f0": f0, "f1": f1, "f2": f2, "f3": f3, "f4": f4}

if __name__ == "__main__":
    chosen = [a.lower() for a in sys.argv[1:]] or list(SECTIONS)
    for name in chosen:
        print(f"== {name} ==", flush=True)
        SECTIONS[name]()
