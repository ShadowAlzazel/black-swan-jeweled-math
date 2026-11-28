"""
levels.py -- Walsh-Fourier levels of SAT energies.

Conventions. An assignment is x in {0,1}^n and chi_S(x) = (-1)^(sum of x_i over i in S). A clause is
violated when every literal is false. The positive literal x_i is false iff x_i = 0, with indicator
(1 + chi_i)/2; the negative literal is false iff x_i = 1, with indicator (1 - chi_i)/2. A clause's
violation indicator is therefore prod_{l in C} (1 + sigma_l chi_l)/2, with sigma = +1 for a positive
literal and -1 for a negative one, and its Walsh coefficient on a subset S of the clause's variables
is 2^-|C| prod_{l in S} sigma_l. Summing over clauses gives the energy's coefficients exactly, for
any n, without enumerating assignments.

At a solution x*, E(x*) = 0 splits as Ehat(empty) + sum_k L_k(x*) = 0, with
L_k(x*) = sum_{|S| = k} Ehat(S) chi_S(x*). The hint shares h_k = L_k / (-Ehat(empty)) sum to 1:
they say which levels carry the energy's dip at the solution. Averaging over all solutions gives
the same split for a formula with many solutions.
"""

from __future__ import annotations

import itertools
import math
import random
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "05_counting_jewels"))

Clause = Tuple[int, ...]
PHI = (1 + 5 ** 0.5) / 2
Q_GOLDEN = 1 / PHI                      # (sqrt 5 - 1)/2: the q-hidden weight with no level-1 hint


# ---------------------------------------------------------------------------
# Exact Walsh coefficients
# ---------------------------------------------------------------------------
def normalize(clause: Sequence[int]) -> Optional[Clause]:
    """Remove repeated literals; None for a tautology (never violated)."""
    lits = set(clause)
    if any(-l in lits for l in lits):
        return None
    return tuple(sorted(lits, key=abs))


def coefficients(clauses: Sequence[Sequence[int]], weights: Optional[Sequence[float]] = None) -> Dict[Tuple[int, ...], float]:
    """Walsh coefficients of E_w(x) = sum_j w_j [clause j violated], keyed by sorted 0-based variable tuples."""
    coef: Dict[Tuple[int, ...], float] = defaultdict(float)
    for j, raw in enumerate(clauses):
        c = normalize(raw)
        if c is None:
            continue
        w = 1.0 if weights is None else float(weights[j])
        k = len(c)
        vs = [abs(l) - 1 for l in c]
        sg = [1 if l > 0 else -1 for l in c]
        for r in range(k + 1):
            for idx in itertools.combinations(range(k), r):
                coef[tuple(vs[i] for i in idx)] += w * math.prod(sg[i] for i in idx) / 2 ** k
    return dict(coef)


def level_weights(coef: Dict[Tuple[int, ...], float]) -> Dict[int, float]:
    """Sum of squared coefficients per level k >= 1 (Parseval: they add up to the variance of E)."""
    out: Dict[int, float] = defaultdict(float)
    for S, c in coef.items():
        if S:
            out[len(S)] += c * c
    return dict(out)


def chi_matrix(assignments: Sequence, n: int) -> np.ndarray:
    """Rows of chi_i(x) = +1 (x_i = 0) or -1 (x_i = 1). Accepts bool lists or int bitmasks (bit i = x_i)."""
    rows = []
    for a in assignments:
        if isinstance(a, (int, np.integer)):
            bits = [(int(a) >> i) & 1 for i in range(n)]
        else:
            bits = [1 if b else 0 for b in a]
        rows.append([1 - 2 * b for b in bits])
    return np.array(rows, dtype=np.int8).reshape(len(rows), n)


def hint_levels(coef: Dict[Tuple[int, ...], float], X: np.ndarray) -> Tuple[float, Dict[int, float]]:
    """(Ehat(empty), {k: L_k}) with L_k averaged over the rows of X (chi values of solutions)."""
    mean = coef.get((), 0.0)
    L: Dict[int, float] = defaultdict(float)
    Xf = X.astype(np.float64)
    for S, c in coef.items():
        if not S or c == 0.0:
            continue
        prod = Xf[:, S[0]].copy()
        for i in S[1:]:
            prod *= Xf[:, i]
        L[len(S)] += c * float(prod.mean())
    return mean, dict(L)


def hint_shares(coef, X) -> Dict[int, float]:
    mean, L = hint_levels(coef, X)
    total = sum(L.values())
    assert abs(mean + total) < 1e-6 * max(1.0, mean), "rows of X must be solutions"
    return {k: v / total for k, v in sorted(L.items())} if total else {}   # total = -Ehat(empty)


def level1_vector(coef, n: int) -> np.ndarray:
    return np.array([coef.get((i,), 0.0) for i in range(n)])


def level1_cosine(coef, X: np.ndarray) -> float:
    """Cosine between the level-1 descent direction (-Ehat({i})) and the mean chi vector of the solutions."""
    n = X.shape[1]
    v = -level1_vector(coef, n)
    mvec = X.astype(np.float64).mean(axis=0)
    nv, nm = np.linalg.norm(v), np.linalg.norm(mvec)
    return float(v @ mvec / (nv * nm)) if nv > 0 and nm > 0 else 0.0


def majority_vote(coef, n: int) -> List[Optional[bool]]:
    """x_i = True where the level-1 coefficient is positive (more positive literals), None on ties."""
    v = level1_vector(coef, n)
    return [True if c > 0 else (False if c < 0 else None) for c in v]


def agreement(vote: Sequence[Optional[bool]], x: Sequence[bool]) -> float:
    """Fraction of coordinates where vote matches x; ties count one half."""
    return sum(0.5 if a is None else float(a == b) for a, b in zip(vote, x)) / len(x)


# ---------------------------------------------------------------------------
# Sphere averages around a solution
# ---------------------------------------------------------------------------
def krawtchouk(n: int, k: int, d: int) -> int:
    return sum((-1) ** j * math.comb(d, j) * math.comb(n - d, k - j) for j in range(0, min(d, k) + 1))


def sphere_average(coef, x: Sequence[bool], d: int) -> float:
    """Exact mean energy over all assignments at Hamming distance d from x."""
    n = len(x)
    mean, L = hint_levels(coef, chi_matrix([x], n))
    return mean + sum(Lk * krawtchouk(n, k, d) / math.comb(n, k) for k, Lk in L.items())


# ---------------------------------------------------------------------------
# Brute force (small n) for checks
# ---------------------------------------------------------------------------
def energy_vector(n: int, clauses: Sequence[Sequence[int]], weights=None) -> np.ndarray:
    idx = np.arange(2 ** n, dtype=np.int64)
    E = np.zeros(2 ** n)
    for j, c in enumerate(clauses):
        viol = np.ones(2 ** n, dtype=bool)
        for l in c:
            bit = (idx >> (abs(l) - 1)) & 1
            viol &= (bit == 0) if l > 0 else (bit == 1)
        E += viol * (1.0 if weights is None else weights[j])
    return E


def walsh_transform(f: np.ndarray) -> np.ndarray:
    a = f.astype(np.float64).copy()
    h = 1
    while h < len(a):
        a = a.reshape(-1, 2 * h)
        left, right = a[:, :h].copy(), a[:, h:].copy()
        a[:, :h], a[:, h:] = left + right, left - right
        a = a.reshape(-1)
        h *= 2
    return a / len(a)


# ---------------------------------------------------------------------------
# Hidden-solution generators: the hint ladder
# ---------------------------------------------------------------------------
def hidden_3sat(n: int, m: int, seed: int, kind: str, q: float = 1.0) -> Tuple[List[Clause], List[bool]]:
    """
    Random 3-SAT with a hidden assignment A, distinct clauses on three distinct variables.
      kind 'q':   a clause with t literals true under A (t >= 1) is kept with probability q^(t-1), so
                  P(clause) is proportional to q^t (Jia, Moore and Strain 2005). q = 1 is naive planting;
                  q = 1/phi cancels the level-1 hint; q < 1/phi points it away from A.
      kind 'two': clauses satisfied by both A and its complement, i.e. t in {1, 2} (Achlioptas, Jia
                  and Moore 2004).
    """
    rng = random.Random(seed)
    A = [rng.random() < 0.5 for _ in range(n)]
    seen, clauses = set(), []
    while len(clauses) < m:
        vs = sorted(rng.sample(range(1, n + 1), 3))
        c = tuple(v if rng.random() < 0.5 else -v for v in vs)
        if c in seen:
            continue
        t = sum((l > 0) == A[abs(l) - 1] for l in c)
        if t == 0 or (kind == "two" and t == 3):
            continue
        if kind == "q" and q < 1.0 and rng.random() >= q ** (t - 1):
            continue
        seen.add(c)
        clauses.append(c)
    return clauses, A


def expected_hint_shares(kind: str, q: float = 1.0) -> Dict[int, float]:
    """Per-clause expected hint shares at A for the generators above (exact, from the 7 or 6 true/false patterns)."""
    pats = [p for p in itertools.product((0, 1), repeat=3) if sum(p) >= 1]
    if kind == "two":
        pats = [p for p in pats if sum(p) <= 2]
        wts = [1.0] * len(pats)
    else:
        wts = [q ** sum(p) for p in pats]
    Z = sum(wts)
    s = lambda p: [2 * b - 1 for b in p]
    e1 = sum(w * sum(s(p)) for p, w in zip(pats, wts)) / Z
    e2 = sum(w * sum(a * b for a, b in itertools.combinations(s(p), 2)) for p, w in zip(pats, wts)) / Z
    e3 = sum(w * math.prod(s(p)) for p, w in zip(pats, wts)) / Z
    # u_C(A) = (1/8)(1 - sum s + sum ss - sss) = 0, so the deficit 1/8 splits as e1, -e2, e3 (over 8)
    return {1: e1, 2: -e2, 3: e3}


# ---------------------------------------------------------------------------
# The XOR energy polynomial as a signed weight enumerator
# ---------------------------------------------------------------------------
def xor_kernel(eqs) -> List[Tuple[int, int]]:
    """Basis of the dependencies among the equations' left sides: (bitmask of equations, parity of their right sides)."""
    pivots: Dict[int, Tuple[int, int, int]] = {}
    basis = []
    for j, (vs, r) in enumerate(eqs):
        vec = 0
        for v in vs:
            vec ^= 1 << v
        tag, rhs = 1 << j, int(r) & 1
        while vec:
            p = vec.bit_length() - 1
            if p not in pivots:
                pivots[p] = (vec, tag, rhs)
                break
            pv, pt, pr = pivots[p]
            vec, tag, rhs = vec ^ pv, tag ^ pt, rhs ^ pr
        if vec == 0:
            basis.append((tag, rhs))
    return basis


def poly_mul(a: List[int], b: List[int]) -> List[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def xor_energy_polynomial(n: int, eqs, max_dim: int = 24) -> Dict:
    """
    N_k = #{x : exactly k equations violated}, from
        sum_k N_k z^k = 2^(n-m) sum_{T in K} (-1)^(rhs(T)) (1 - z)^|T| (1 + z)^(m - |T|),
    K the dependencies among the left sides. Proof: z^u = (1+z)/2 + (1-z)/2 (-1)^r chi_S for a single
    equation, multiply over equations and sum over x; only T with cancelling left sides survive.
    """
    m = len(eqs)
    basis = xor_kernel(eqs)
    dim = len(basis)
    if dim > max_dim:
        raise ValueError(f"kernel dimension {dim} too large to enumerate")
    signed = [0] * (m + 1)                 # signed count of dependencies by size
    for mask in range(2 ** dim):
        tag, rhs = 0, 0
        for b in range(dim):
            if mask >> b & 1:
                tag ^= basis[b][0]
                rhs ^= basis[b][1]
        signed[bin(tag).count("1")] += -1 if rhs else 1
    total = [0] * (m + 1)
    for w, a in enumerate(signed):
        if a:
            term = [1]
            for _ in range(w):
                term = poly_mul(term, [1, -1])
            for _ in range(m - w):
                term = poly_mul(term, [1, 1])
            for k, c in enumerate(term):
                total[k] += a * c
    scale = Fraction(2) ** (n - m)
    N = [scale * c for c in total]
    assert all(x.denominator == 1 for x in N), "counts must be integers"
    N = [int(x) for x in N]
    return {"N": N, "kernel_dim": dim, "rank": m - dim, "signed_kernel_sizes": signed,
            "satisfiable": all(r == 0 for _, r in basis)}


def zero_multiplicity_at_minus_one(N: List[int]) -> int:
    """How many times (1 + z) divides sum N_k z^k (exact integer division)."""
    p = list(N)
    count = 0
    while len(p) > 1 and sum(c * (-1) ** k for k, c in enumerate(p)) == 0:
        d = len(p) - 1                  # p = (z + 1) q:  p_k = q_(k-1) + q_k
        q = [0] * d
        q[d - 1] = p[d]
        for k in range(d - 1, 0, -1):
            q[k - 1] = p[k] - q[k]
        p = q
        count += 1
    return count


def divide_by_one_plus_z(N: List[int], times: int) -> List[int]:
    """Exact quotient of sum N_k z^k by (1 + z)^times (the division must be exact)."""
    p = list(N)
    for _ in range(times):
        d = len(p) - 1
        q = [0] * d
        q[d - 1] = p[d]
        for k in range(d - 1, 0, -1):
            q[k - 1] = p[k] - q[k]
        assert q[0] == p[0], "not divisible"
        p = q
    return p


# ---------------------------------------------------------------------------
# All solutions of a small formula (for averaging hints over every solution)
# ---------------------------------------------------------------------------
def all_solutions(n: int, clauses: Sequence[Sequence[int]], limit: int = 500_000) -> List[int]:
    """Every satisfying assignment as an int bitmask (bit v-1 = variable v), by DPLL with unit propagation."""
    clauses = [tuple(c) for c in clauses]
    deg = [0] * (n + 1)
    for c in clauses:
        for l in c:
            deg[abs(l)] += 1
    order = sorted(range(1, n + 1), key=lambda v: (-deg[v], v))
    val: List[Optional[bool]] = [None] * (n + 1)
    sols: List[int] = []

    def status():
        unit = None
        for c in clauses:
            nfree, free, sat = 0, 0, False
            for l in c:
                x = val[abs(l)]
                if x is None:
                    nfree += 1
                    free = l
                elif x == (l > 0):
                    sat = True
                    break
            if sat:
                continue
            if nfree == 0:
                return "conflict", 0
            if nfree == 1 and unit is None:
                unit = free
        return ("unit", unit) if unit is not None else ("ok", 0)

    def rec():
        trail = []
        while True:
            st, l = status()
            if st == "conflict":
                for v in trail:
                    val[v] = None
                return
            if st == "unit":
                val[abs(l)] = l > 0
                trail.append(abs(l))
                continue
            break
        free = next((v for v in order if val[v] is None), None)
        if free is None:
            sols.append(sum(1 << (v - 1) for v in range(1, n + 1) if val[v]))
        else:
            for b in (False, True):
                if len(sols) >= limit:
                    break
                val[free] = b
                rec()
                val[free] = None
        for v in trail:
            val[v] = None

    rec()
    if len(sols) >= limit:
        raise RuntimeError("solution limit reached")
    return sols


# ---------------------------------------------------------------------------
# Forward maps of Experiment 05's rewrites (old solution -> new solution)
# ---------------------------------------------------------------------------
def tovey_forward(n: int, clauses: Sequence[Sequence[int]], x: Sequence[bool]) -> List[bool]:
    """Copies numbered as in formulas.tovey_split: one per literal occurrence in clause order, then one free copy per unused variable."""
    out = []
    used = set()
    for c in clauses:
        for l in c:
            out.append(x[abs(l) - 1])
            used.add(abs(l))
    for v in range(1, n + 1):
        if v not in used:
            out.append(x[v - 1])
    return out


def definitional_forward(n: int, clauses: Sequence[Sequence[int]], x: Sequence[bool]) -> List[bool]:
    """y = (a or b) for each clause of length >= 3, in clause order, as in formulas.definitional."""
    lit = lambda l: x[abs(l) - 1] == (l > 0)
    return list(x) + [lit(c[0]) or lit(c[1]) for c in clauses if len(c) >= 3]
