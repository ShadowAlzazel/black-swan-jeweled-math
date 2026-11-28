"""
make_lab_data.py -- data for web/07_energy_jewels.html (Experiment 05).

For each formula it computes the whole energy polynomial
    Z(z) = sum_x z^E(x) = sum_k N_k z^k      (E = number of violated clauses)
exactly, finds its complex zeros, and writes everything into the page between the
markers /*LAB_DATA_START*/ ... /*LAB_DATA_END*/ so the page stays self-contained.

Exact counts for n <= 26 use a split brute force: write an assignment as (high half,
low half). A clause is violated exactly when all its high-half literals and all its
low-half literals are false, so the violation count of every assignment is the matrix
product A^T B of two 0/1 tables (clauses x half-assignments). One matrix product and
one histogram give every N_k; the result is cross-checked against counting.py.

    python experiments/05_counting_jewels/make_lab_data.py
"""

from __future__ import annotations

import json
import os
import re
import sys
from typing import List, Sequence, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import counting as C  # noqa: E402
import formulas as F  # noqa: E402

REPO = os.path.dirname(os.path.dirname(HERE))
PAGE = os.path.join(REPO, "web", "07_energy_jewels.html")
Clause = Tuple[int, ...]


def energy_counts_split(n: int, clauses: Sequence[Clause]) -> List[int]:
    """Exact N_0..N_m for n <= 26 via one matrix product over the two halves of the variables."""
    m = len(clauses)
    lo_bits = n // 2
    hi_bits = n - lo_bits
    lo = np.arange(1 << lo_bits, dtype=np.int64)
    hi = np.arange(1 << hi_bits, dtype=np.int64)
    A = np.ones((m, hi.size), dtype=np.float32)   # all high-half literals of clause c are false
    B = np.ones((m, lo.size), dtype=np.float32)   # all low-half literals of clause c are false
    for ci, c in enumerate(clauses):
        for l in c:
            v = abs(l) - 1
            if v < lo_bits:
                bit = (lo >> v) & 1
                B[ci] *= (bit == 0) if l > 0 else (bit == 1)
            else:
                bit = (hi >> (v - lo_bits)) & 1
                A[ci] *= (bit == 0) if l > 0 else (bit == 1)
    V = A.T @ B                                     # violated-clause count of (hi, lo)
    counts = np.bincount(np.rint(V).astype(np.int64).ravel(), minlength=m + 1)
    return [int(x) for x in counts]


def trimmed(N: Sequence[int]) -> List[int]:
    last = max((k for k, c in enumerate(N) if c), default=0)
    return list(N[:last + 1])


def zeros_of(N: Sequence[int]) -> List[List[float]]:
    coeffs = np.array(trimmed(N)[::-1], dtype=float)       # highest degree first
    if coeffs.size <= 1:
        return []
    roots = np.roots(coeffs)
    return [[round(float(z.real), 5), round(float(z.imag), 5)] for z in sorted(roots, key=lambda z: (abs(z), np.angle(z)))]


def entry(label: str, group: str, n: int, clauses: Sequence[Clause], **meta) -> dict:
    N = energy_counts_split(n, clauses)
    assert sum(N) == 2 ** n
    Nt = trimmed(N)
    return {"label": label, "group": group, "n": n, "m": len(clauses), "N": Nt, "zeros": zeros_of(Nt), **meta}


def main() -> None:
    # Cross-check the split brute force against counting.py.
    n0, cl0 = F.read_dimacs(os.path.join(F.EXP03, "missed_sat_instance_35.cnf"))
    assert energy_counts_split(n0, cl0) == C.energy_counts_bruteforce(n0, cl0)
    small = [i for i in F.exp03_instances() if i["n"] == 12][:4]
    for inst in small:
        cl = [tuple(c) for c in inst["clauses"]]
        assert energy_counts_split(inst["n"], cl) == C.energy_counts_bruteforce(inst["n"], cl)

    data = []
    data.append(entry("missed_sat_instance_35.cnf (Exp 03 id 35, random, n = 18)", "Experiment 03 missed instance",
                      n0, cl0, family="random", sat_label="SAT", source="experiments/03_golden_sweep/missed_sat_instance_35.cnf"))
    for inst in F.exp03_instances():
        if inst["n"] > 24:
            continue
        cl = [tuple(c) for c in inst["clauses"]]
        e = entry(f"Exp 03 id {inst['id']} · {inst['family']} · n = {inst['n']} · {inst['split']} · {inst['sat_label']}",
                  "Experiment 03 formulas (n ≤ 24)", inst["n"], cl,
                  id=inst["id"], family=inst["family"], split=inst["split"], sat_label=inst["sat_label"])
        assert (e["N"][0] > 0) == (inst["sat_label"] == "SAT")
        data.append(e)
    for alpha in (2, 3, 4, 4.25, 5, 6):
        for seed in range(3):
            m = round(alpha * 16)
            cl, _ = F.random_3sat(16, m, seed=1600 + int(alpha * 100) + seed, planted=False)
            e = entry(f"random 3-SAT · n = 16 · α = {alpha} · seed {seed}", "Random 3-SAT, n = 16, across density",
                      16, cl, alpha=alpha, seed=seed)
            e["sat_label"] = "SAT" if e["N"][0] > 0 else "UNSAT"
            data.append(e)

    payload = json.dumps(data, separators=(",", ":"))
    with open(PAGE, "r", encoding="utf-8") as fh:
        page = fh.read()
    page, k = re.subn(r"/\*LAB_DATA_START\*/.*?/\*LAB_DATA_END\*/",
                      lambda _: "/*LAB_DATA_START*/" + payload + "/*LAB_DATA_END*/", page, flags=re.S)
    if k != 1:
        raise SystemExit("data markers not found in " + PAGE)
    with open(PAGE, "w", encoding="utf-8") as fh:
        fh.write(page)

    # Summary of what the data show
    def inside(e, r=1.0):
        return sum(1 for x, y in e["zeros"] if x * x + y * y < r * r)
    print(f"wrote {len(data)} formulas, {len(payload) / 1024:.0f} KB of data into {os.path.relpath(PAGE, REPO)}")
    near = [min(e["zeros"], key=lambda z: z[0] ** 2 + z[1] ** 2) for e in data if e["zeros"] and e["N"][0] > 0]
    on_neg = sum(1 for x, y in near if x < 0 and abs(y) < 1e-6)
    print(f"nearest zero on the negative real axis for {on_neg}/{len(near)} satisfiable formulas")
    ratios = [abs((x * x + y * y) ** 0.5 / (e["N"][0] / e["N"][1])) for e, (x, y) in
              zip([e for e in data if e["zeros"] and e["N"][0] > 0], near) if len(e["N"]) > 1]
    print(f"|nearest zero| / (N0/N1): median {np.median(ratios):.3f}, range {min(ratios):.3f}-{max(ratios):.3f}")
    for alpha in (2, 3, 4, 4.25, 5, 6):
        es = [e for e in data if e.get("alpha") == alpha]
        print(f"alpha {alpha}: SAT {sum(e['N'][0] > 0 for e in es)}/3, zeros inside |z|<1: {[inside(e) for e in es]}, "
              f"degree {[len(e['N']) - 1 for e in es]}, nearest |z| {[round(min((x*x+y*y)**0.5 for x, y in e['zeros']), 3) for e in es]}")
    unsat = [e for e in data if e["N"][0] == 0]
    print("UNSAT formulas:", len(unsat), "all have a zero at the origin:",
          all(any(abs(x) < 1e-9 and abs(y) < 1e-9 for x, y in e["zeros"]) for e in unsat))


if __name__ == "__main__":
    main()
