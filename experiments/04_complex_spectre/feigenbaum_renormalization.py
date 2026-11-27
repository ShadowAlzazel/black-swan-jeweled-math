"""
feigenbaum_renormalization.py
-----------------------------
Independent check of the Feigenbaum numbers quoted with Experiment 04, and the
"eigenvalue" route to delta.

    python experiments/04_complex_spectre/feigenbaum_renormalization.py

Writes feigenbaum_results.json next to this script. Requires mpmath.

1. Superstable parameters c_k of f_c(z) = z^2 + c (0 has period 2^k), through
   period 2^14, by Newton's method on f_c^(2^k)(0) = 0.
   delta_k = (c_{k-1} - c_{k-2}) / (c_k - c_{k-1})     -> 4.669201609...
   alpha_k = d_k / d_{k+1},  d_k = f_{c_k}^(2^(k-1))(0)  -> -2.502907875...
2. The Feigenbaum-Cvitanovic fixed point g(x) = alpha g(g(x / alpha)),
   g(0) = 1, alpha = 1 / g(1), as an even polynomial solved by collocation.
   delta is the only eigenvalue of the linearized doubling operator outside
   the unit circle.
3. "Compose twice and rescale": the rescaled iterates of z^2 + c_inf converge to g.
4. The kneading sequence at the accumulation point is the parity of the ruler
   function nu_2(j), the same function that schedules Gray-code flips and
   Luby restarts.
5. The sign of alpha as an Euler half-turn: alpha = |alpha| e^{i pi}.
"""

from __future__ import annotations

import json
import os
import time

import mpmath as mp
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.join(SCRIPT_DIR, "feigenbaum_results.json")

# Values quoted in the Experiment 04 conversation (through period 512).
QUOTED_DELTA_RATIOS = ["4.668548581", "4.669060661", "4.669171555"]
QUOTED_SIGNED_ALPHA = "-2.502906472"
# Literature (OEIS A006890, A006891), for comparison only. Kept as strings and
# converted at the working precision (an mpf made at import time would carry only
# mpmath's default 15 digits).
DELTA_REF_STR = "4.669201609102990671853203820466201617258185577475768632745651"
ALPHA_REF_STR = "-2.502907875095892822283902873218215786381271376727149977336192"


def DELTA_REF():
    return mp.mpf(DELTA_REF_STR)


def ALPHA_REF():
    return mp.mpf(ALPHA_REF_STR)


def s(x, digits=20) -> str:
    return mp.nstr(x, digits)


# ---------------------------------------------------------------------------
# 1. Superstable cycles
# ---------------------------------------------------------------------------
def critical_orbit(c, steps):
    z, dz = mp.mpf(0), mp.mpf(0)
    for _ in range(steps):
        dz = 2 * z * dz + 1
        z = z * z + c
    return z, dz


def superstable_parameters(k_max: int):
    cs = [mp.mpf(0), mp.mpf(-1)]
    newton_steps = []
    for k in range(2, k_max + 1):
        ratio = (cs[-2] - cs[-3]) / (cs[-1] - cs[-2]) if k >= 3 else mp.mpf("4.669")
        x = cs[-1] + (cs[-1] - cs[-2]) / ratio
        for it in range(200):
            g, dg = critical_orbit(x, 2 ** k)
            step = g / dg
            x -= step
            if abs(step) < mp.mpf(10) ** (-(mp.mp.dps - 12)):
                break
        newton_steps.append(it + 1)
        # check the period is exactly 2^k (not a divisor)
        z = mp.mpf(0)
        for j in range(1, 2 ** k):
            z = z * z + x
            if abs(z) < mp.mpf(10) ** (-(mp.mp.dps // 2)):
                raise RuntimeError(f"Newton converged to a lower period at k={k}")
        cs.append(x)
    return cs, newton_steps


def section_superstable(k_max: int = 14) -> dict:
    cs, newton_steps = superstable_parameters(k_max)
    deltas = []
    for k in range(2, len(cs)):
        deltas.append((k, (cs[k - 1] - cs[k - 2]) / (cs[k] - cs[k - 1])))
    d = {}
    for k in range(1, len(cs)):
        z = mp.mpf(0)
        for _ in range(2 ** (k - 1)):
            z = z * z + cs[k]
        d[k] = z
    alphas = [(k, d[k] / d[k + 1]) for k in range(1, len(cs) - 1)]
    # Aitken extrapolation of the accumulation point
    a, b, c = cs[-3], cs[-2], cs[-1]
    c_inf = c - (c - b) ** 2 / ((c - b) - (b - a))

    def matches(value, quoted: str) -> bool:
        return mp.nstr(value, len(quoted.replace("-", "").replace(".", ""))) == quoted

    through_512 = [(k, v) for k, v in deltas if k <= 9]
    quoted_match = {}
    for q in QUOTED_DELTA_RATIOS:
        hit = [k for k, v in through_512 if matches(v, q)]
        quoted_match[q] = {"matches_delta_k_with_k": hit, "period_of_last_parameter": [2 ** k for k in hit]}
    alpha_hit = [k for k, v in alphas if k + 1 <= 9 and matches(v, QUOTED_SIGNED_ALPHA)]
    return {
        "c_k": {str(2 ** k): s(v, 30) for k, v in enumerate(cs)},
        "newton_iterations": newton_steps,
        "delta_k": {f"k={k} (periods {2 ** (k - 2)},{2 ** (k - 1)},{2 ** k})": s(v, 15) for k, v in deltas},
        "signed_alpha_k": {f"k={k} (d_{k}/d_{k + 1}, periods {2 ** k}/{2 ** (k + 1)})": s(v, 15) for k, v in alphas},
        "quoted_delta_ratios_reproduced": quoted_match,
        "quoted_signed_alpha_reproduced_at_k": alpha_hit,
        "best_delta_estimate": s(deltas[-1][1], 15),
        "best_delta_error_vs_literature": s(abs(deltas[-1][1] - DELTA_REF()), 3),
        "best_alpha_estimate": s(alphas[-1][1], 15),
        "best_alpha_error_vs_literature": s(abs(alphas[-1][1] - ALPHA_REF()), 3),
        "delta_error_ratios": [s((deltas[i - 1][1] - DELTA_REF()) / (deltas[i][1] - DELTA_REF()), 6)
                               for i in range(1, len(deltas))],
        "c_inf_aitken": s(c_inf, 25),
        "sign_pattern_of_d_k": "".join("+" if d[k] > 0 else "-" for k in sorted(d)),
        "_cs": cs,
        "_c_inf": c_inf,
    }


# ---------------------------------------------------------------------------
# 2. Renormalization fixed point and its linearization
# ---------------------------------------------------------------------------
def make_g(coeffs):
    def g(x):
        x2 = x * x
        total, power = mp.mpf(1), mp.mpf(1)
        for a in coeffs:
            power *= x2
            total += a * power
        return total
    return g


def doubling(coeffs, xs):
    """T[g](x) = alpha g(g(x / alpha)), alpha = 1 / g(1), at the collocation points."""
    g = make_g(coeffs)
    alpha = 1 / g(mp.mpf(1))
    return [alpha * g(g(x / alpha)) for x in xs], alpha


def section_fixed_point(n_terms: int) -> dict:
    xs = [mp.sqrt((1 + mp.cos(mp.pi * (2 * j - 1) / (2 * n_terms))) / 2) for j in range(1, n_terms + 1)]
    V = mp.matrix(n_terms, n_terms)
    for j, x in enumerate(xs):
        for i in range(n_terms):
            V[j, i] = x ** (2 * (i + 1))
    guess = [-1.52763, 0.104815, 0.0267057, -0.00352741, 0.0000816, 0.0000254, -0.0000027]
    coeffs = [mp.mpf(guess[i]) if i < len(guess) else mp.mpf(0) for i in range(n_terms)]
    h = mp.mpf(10) ** (-(mp.mp.dps // 2))

    def residual(a):
        tv, _ = doubling(a, xs)
        g = make_g(a)
        return mp.matrix([tv[j] - g(xs[j]) for j in range(n_terms)])

    for it in range(40):
        r0 = residual(coeffs)
        J = mp.matrix(n_terms, n_terms)
        for i in range(n_terms):
            bumped = list(coeffs)
            bumped[i] += h
            ri = residual(bumped)
            for j in range(n_terms):
                J[j, i] = (ri[j] - r0[j]) / h
        step = mp.lu_solve(J, -r0)
        coeffs = [coeffs[i] + step[i] for i in range(n_terms)]
        if mp.norm(step) < mp.mpf(10) ** (-(mp.mp.dps // 2 - 5)):
            break
    newton_iterations = it + 1

    # Linearized doubling operator in coefficient space: a' = V^{-1}(T[g](x) - 1)
    t0, alpha = doubling(coeffs, xs)
    DT = mp.matrix(n_terms, n_terms)
    Vinv = mp.inverse(V)
    for i in range(n_terms):
        bumped = list(coeffs)
        bumped[i] += h
        ti, _ = doubling(bumped, xs)
        col = Vinv * mp.matrix([(ti[j] - t0[j]) / h for j in range(n_terms)])
        for j in range(n_terms):
            DT[j, i] = col[j]
    spectrum = np.linalg.eigvals(np.array(DT.tolist(), dtype=float))
    spectrum = sorted(spectrum, key=lambda z: -abs(z))
    # Refine the dominant eigenvalue with power iteration in full precision.
    v = mp.matrix([1] * n_terms)
    lam = mp.mpf(0)
    for _ in range(400):
        w = DT * v
        lam_new = w[0] / v[0]
        v = w / mp.norm(w)
        if abs(lam_new - lam) < mp.mpf(10) ** (-(mp.mp.dps - 20)):
            lam = lam_new
            break
        lam = lam_new
    return {
        "terms_in_x^2": n_terms,
        "newton_iterations": newton_iterations,
        "alpha = 1/g(1)": s(alpha, 25),
        "alpha_error_vs_literature": s(abs(alpha - ALPHA_REF()), 3),
        "delta = leading eigenvalue of DT": s(lam, 25),
        "delta_error_vs_literature": s(abs(lam - DELTA_REF()), 3),
        "eigenvalues_outside_unit_circle": [complex(z) for z in spectrum if abs(z) > 1],
        "next_eigenvalues": [complex(round(z.real, 8), round(z.imag, 8)) for z in spectrum[1:6]],
        "next_eigenvalues_matching_alpha^-2j": {
            f"alpha^-{2 * j}": next((round(z.real, 8) for z in spectrum[1:8]
                                     if abs(z - float(alpha ** (-2 * j))) < 1e-6), None)
            for j in (1, 2)
        },
        "g_first_coefficients": [s(c, 15) for c in coeffs[:8]],
        "_coeffs": coeffs,
        "_alpha": alpha,
    }


# ---------------------------------------------------------------------------
# 3. Compose twice, rescale, compare with the previous level
# ---------------------------------------------------------------------------
def section_rescaled_iterates(c_inf, g_coeffs, levels: int = 9) -> dict:
    g = make_g(g_coeffs)
    xs = [mp.mpf(i) / 40 for i in range(-40, 41)]

    def iterate(z, n):
        for _ in range(n):
            z = z * z + c_inf
        return z

    rows = []
    prev = None
    for k in range(levels + 1):
        p = 2 ** k
        scale = iterate(mp.mpf(0), p)
        hk = [iterate(scale * x, p) / scale for x in xs]
        err_g = max(abs(hk[i] - g(xs[i])) for i in range(len(xs)))
        err_prev = max(abs(hk[i] - prev[i]) for i in range(len(xs))) if prev else None
        rows.append({
            "level_k": k,
            "iterate": f"f^{p}",
            "rescale_by": s(scale, 12),
            "max|h_k - g| on [-1,1]": s(err_g, 4),
            "max|h_k - h_(k-1)|": s(err_prev, 4) if err_prev is not None else None,
        })
        prev = hk
    errs = [mp.mpf(r["max|h_k - g| on [-1,1]"]) for r in rows]
    return {
        "map": "f(z) = z^2 + c_inf, h_k(x) = f^(2^k)(s_k x) / s_k with s_k = f^(2^k)(0)",
        "levels": rows,
        "successive_error_ratios": [s(errs[i] / errs[i - 1], 5) for i in range(2, len(errs))],
        "scale_ratios_s_k/s_(k+1)": [s(mp.mpf(rows[i]["rescale_by"]) / mp.mpf(rows[i + 1]["rescale_by"]), 10)
                                      for i in range(len(rows) - 1)],
    }


# ---------------------------------------------------------------------------
# 4. Kneading sequence, ruler function, Gray code, Luby
# ---------------------------------------------------------------------------
def nu2(j: int) -> int:
    return (j & -j).bit_length() - 1


def section_ruler(cs) -> dict:
    k = 10
    c = cs[k]
    z = mp.mpf(0)
    symbols = []
    for _ in range(1, 2 ** k):
        z = z * z + c
        # z^2 + c with c < 0 is conjugate to the logistic map with orientation reversed:
        # z < 0 <-> x > 1/2 (R)
        symbols.append("R" if z < 0 else "L")
    ruler = ["R" if nu2(j) % 2 == 0 else "L" for j in range(1, 2 ** k)]
    gray_flips = [((j ^ (j >> 1)) ^ ((j - 1) ^ ((j - 1) >> 1))).bit_length() - 1 for j in range(1, 2 ** k)]

    def luby(i: int) -> int:
        kk = 1
        while (1 << kk) - 1 < i:
            kk += 1
        if i == (1 << kk) - 1:
            return 1 << (kk - 1)
        return luby(i - (1 << (kk - 1)) + 1)

    reluctant, u, v = [], 1, 1
    while len(reluctant) < 2 ** k:
        reluctant.append(v)
        if (u & -u) == v:
            u, v = u + 1, 1
        else:
            v = 2 * v
    runs, u = [], 1
    while len(runs) < 2 ** k:
        runs.extend(2 ** e for e in range(nu2(u) + 1))
        u += 1
    n_check = 2 ** k
    return {
        "superstable_parameter_used": f"period {2 ** k}",
        "kneading_first_32": "".join(symbols[:32]),
        "kneading_equals_ruler_parity_for_first": f"{2 ** k - 1} symbols" if symbols == ruler else "MISMATCH",
        "gray_code_flip_position_equals_nu2(j)": gray_flips == [nu2(j) for j in range(1, 2 ** k)],
        "luby_recursive_equals_knuth_reluctant_doubling": [luby(i) for i in range(1, n_check + 1)] == reluctant[:n_check],
        "luby_equals_runs_1_2_..._2^nu2(u)": reluctant[:n_check] == runs[:n_check],
        "luby_first_31": reluctant[:31],
        "two_adic_absolute_value": "|j|_2 = 2^(-nu2(j)): the same ruler function defines the 2-adic ultrametric",
    }


# ---------------------------------------------------------------------------
# 5. The sign of alpha as an Euler half-turn
# ---------------------------------------------------------------------------
def section_half_turn(alpha) -> dict:
    log_alpha = mp.log(mp.mpc(alpha))  # principal branch: log|alpha| + i pi
    return {
        "log(alpha)": f"{s(mp.re(log_alpha), 15)} + {s(mp.im(log_alpha), 15)} i",
        "imag_part_equals_pi": bool(abs(mp.im(log_alpha) - mp.pi) < mp.mpf(10) ** -30),
        "meaning": "In w = e^z coordinates, x -> alpha x is the translation z -> z + log|alpha| + i pi: "
                   "a radial step of log|alpha| plus Euler's half-turn e^{i pi} = -1.",
        "two_levels": f"alpha^2 = {s(alpha ** 2, 15)} = |alpha|^2 e^(2 pi i): two period doublings make one full wrap",
        "log|alpha|": s(mp.log(abs(alpha)), 15),
        "log(phi)": s(mp.log((1 + mp.sqrt(5)) / 2), 15),
        "log(delta)": s(mp.log(DELTA_REF()), 15),
    }


def strip_private(obj):
    if isinstance(obj, dict):
        return {k: strip_private(v) for k, v in obj.items() if not k.startswith("_")}
    if isinstance(obj, list):
        return [strip_private(v) for v in obj]
    if isinstance(obj, complex):
        return [round(obj.real, 10), round(obj.imag, 10)]
    return obj


def main() -> None:
    t0 = time.perf_counter()
    mp.mp.dps = 80
    print("[feigenbaum] superstable cycles through period 2^14 ...", flush=True)
    sup = section_superstable(14)
    print("[feigenbaum] renormalization fixed point ...", flush=True)
    mp.mp.dps = 90
    fixed = {}
    for n in (8, 16, 24, 32):
        fixed[str(n)] = section_fixed_point(n)
        print(f"  {n:2d} terms: delta = {fixed[str(n)]['delta = leading eigenvalue of DT']}", flush=True)
    best = fixed["32"]
    print("[feigenbaum] rescaled iterates ...", flush=True)
    mp.mp.dps = 80
    rescaled = section_rescaled_iterates(sup["_c_inf"], best["_coeffs"])
    print("[feigenbaum] kneading / ruler / Gray / Luby ...", flush=True)
    ruler = section_ruler(sup["_cs"])
    half = section_half_turn(best["_alpha"])
    results = {
        "1_superstable_cycles": sup,
        "2_fixed_point_and_eigenvalue": fixed,
        "3_compose_twice_and_rescale": rescaled,
        "4_ruler_function": ruler,
        "5_alpha_sign_as_half_turn": half,
        "seconds": round(time.perf_counter() - t0, 1),
    }
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(strip_private(results), fh, indent=2)
    out = strip_private(results)
    print(json.dumps({k: out["1_superstable_cycles"][k] for k in (
        "quoted_delta_ratios_reproduced", "quoted_signed_alpha_reproduced_at_k",
        "best_delta_estimate", "best_alpha_estimate", "c_inf_aitken")}, indent=2))
    print(f"[feigenbaum] wrote {OUT_PATH} in {results['seconds']} s")


if __name__ == "__main__":
    main()
