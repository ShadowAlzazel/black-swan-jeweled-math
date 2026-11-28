"""
formulas.py -- formula families and solution-preserving transformations (Experiment 05).

Clauses are tuples of nonzero integers in DIMACS style: literal +v means x_v,
-v means not x_v, variables are 1..n. Every generator is seeded and pure Python.
"""

from __future__ import annotations

import itertools
import json
import math
import os
import random
from typing import Dict, List, Optional, Sequence, Tuple

Clause = Tuple[int, ...]
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
EXP03 = os.path.join(REPO, "experiments", "03_golden_sweep")
HAT_PATCHES = os.path.join(REPO, "experiments", "02_tiling_jewels", "tiling_research", "hat_patches.json")
PHI = (1 + math.sqrt(5)) / 2


# ---------------------------------------------------------------------------
# Loading and checking
# ---------------------------------------------------------------------------
def read_dimacs(path: str) -> Tuple[int, List[Clause]]:
    n, clauses, cur = 0, [], []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if not parts or parts[0] in ("c", "%"):
                continue
            if parts[0] == "p":
                n = int(parts[2])
                continue
            for tok in parts:
                lit = int(tok)
                if lit == 0:
                    clauses.append(tuple(cur))
                    cur = []
                else:
                    cur.append(lit)
    return n, clauses


def exp03_instances() -> List[dict]:
    with open(os.path.join(EXP03, "dataset.json"), "r", encoding="utf-8") as fh:
        return json.load(fh)["instances"]


def satisfies(clauses: Sequence[Clause], assignment: Sequence[bool]) -> bool:
    return all(any((l > 0) == assignment[abs(l) - 1] for l in c) for c in clauses)


def violated(clauses: Sequence[Clause], assignment: Sequence[bool]) -> int:
    return sum(not any((l > 0) == assignment[abs(l) - 1] for l in c) for c in clauses)


# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------
def random_3sat(n: int, m: int, seed: int, planted: bool = False) -> Tuple[List[Clause], Optional[List[bool]]]:
    """Distinct clauses on three distinct variables; planted formulas reject clauses the hidden assignment violates."""
    rng = random.Random(seed)
    hidden = [rng.random() < 0.5 for _ in range(n)] if planted else None
    seen, clauses = set(), []
    while len(clauses) < m:
        vs = sorted(rng.sample(range(1, n + 1), 3))
        c = tuple(v if rng.random() < 0.5 else -v for v in vs)
        if c in seen or (hidden and not any((l > 0) == hidden[abs(l) - 1] for l in c)):
            continue
        seen.add(c)
        clauses.append(c)
    return clauses, hidden


def planted_3xor(n: int, m: int, seed: int) -> Tuple[List[Tuple[Tuple[int, int, int], int]], List[Clause], List[bool]]:
    """Planted x_a + x_b + x_c = t (mod 2) equations and their exact CNF: four clauses forbid the four wrong parities."""
    rng = random.Random(seed)
    hidden = [rng.random() < 0.5 for _ in range(n)]
    eqs, clauses = [], []
    for _ in range(m):
        vs = tuple(sorted(rng.sample(range(n), 3)))
        t = sum(hidden[v] for v in vs) % 2
        eqs.append((vs, t))
        for bits in itertools.product((0, 1), repeat=3):
            if sum(bits) % 2 != t:
                clauses.append(tuple(v + 1 if b == 0 else -(v + 1) for v, b in zip(vs, bits)))
    return eqs, clauses, hidden


def gauss_xor(n: int, eqs) -> Optional[List[bool]]:
    """Gaussian elimination over GF(2); returns a solution (free variables 0) or None."""
    pivots: Dict[int, Tuple[int, int]] = {}
    for vs, t in eqs:
        row = 0
        for v in vs:
            row ^= 1 << v
        for p, (prow, pt) in pivots.items():
            if row >> p & 1:
                row ^= prow
                t ^= pt
        if row == 0:
            if t:
                return None
            continue
        p = row.bit_length() - 1
        for q in list(pivots):
            if pivots[q][0] >> p & 1:
                pivots[q] = (pivots[q][0] ^ row, pivots[q][1] ^ t)
        pivots[p] = (row, t)
    x = [False] * n
    for p, (row, t) in pivots.items():
        x[p] = bool(t)
    return x


def hat_patch(level: int = 2) -> dict:
    with open(HAT_PATCHES, "r", encoding="utf-8") as fh:
        return json.load(fh)["patches"][level]


def tile_formula(points: Sequence[Tuple[float, float]], alpha: float, seed: int,
                 planted: bool = True, neighbors: int = 6) -> Tuple[List[Clause], Optional[List[bool]]]:
    """
    Experiment 03's geometric family on any point set: each clause is a variable plus two of its
    `neighbors` nearest points, with random signs. Variable j lives at points[j].
    """
    n = len(points)
    rng = random.Random(seed)
    near = [sorted((k for k in range(n) if k != j),
                   key=lambda k: ((points[k][0] - points[j][0]) ** 2 + (points[k][1] - points[j][1]) ** 2, k))[:neighbors]
            for j in range(n)]
    hidden = [rng.random() < 0.5 for _ in range(n)] if planted else None
    m = round(alpha * n)
    seen, clauses = set(), []
    while len(clauses) < m:
        j = rng.randrange(n)
        vs = sorted([j] + rng.sample(near[j], 2))
        c = tuple((v + 1) if rng.random() < 0.5 else -(v + 1) for v in vs)
        if c in seen or (hidden and not any((l > 0) == hidden[abs(l) - 1] for l in c)):
            continue
        seen.add(c)
        clauses.append(c)
    return clauses, hidden


def hat_points(level: int = 2) -> Tuple[List[Tuple[float, float]], List[Tuple[int, ...]]]:
    """Centroids and substitution addresses of the hat tiles."""
    patch = hat_patch(level)
    pts, addrs = [], []
    for t in patch["tiles"]:
        poly = t["polygon"]
        pts.append((sum(p[0] for p in poly) / len(poly), sum(p[1] for p in poly) / len(poly)))
        addrs.append(tuple(t["address"]))
    return pts, addrs


def grid_points(n: int) -> List[Tuple[float, float]]:
    side = math.ceil(math.sqrt(n))
    return [(j % side, j // side) for j in range(n)]


# ---------------------------------------------------------------------------
# Transformations. Each returns (n', clauses', lift) where lift maps a solution of the
# new formula to a solution of the old one. "Same solutions" transformations keep n.
# ---------------------------------------------------------------------------
def gauge(n: int, clauses: Sequence[Clause], seed: int):
    """Complement a random half of the variables: an isomorphic copy (control)."""
    rng = random.Random(seed)
    flip = {v for v in range(1, n + 1) if rng.random() < 0.5}
    new = [tuple(-l if abs(l) in flip else l for l in c) for c in clauses]
    return n, new, (lambda a: [(not a[v - 1]) if v in flip else a[v - 1] for v in range(1, n + 1)])


def duplicate(n: int, clauses: Sequence[Clause], seed: int, copies=(1, 2, 3)):
    """Repeat each clause 1-3 times: integer clause weights. Identical solution set."""
    rng = random.Random(seed)
    new = [c for c in clauses for _ in range(rng.choice(copies))]
    return n, new, (lambda a: list(a))


def resolvents(n: int, clauses: Sequence[Clause], max_len: int = 3):
    """Add every non-tautological resolvent of length <= max_len that is not already present: implied clauses, identical solution set."""
    have = {tuple(sorted(set(c))) for c in clauses}
    by_lit: Dict[int, List[Clause]] = {}
    for c in clauses:
        for l in set(c):
            by_lit.setdefault(l, []).append(c)
    added = []
    for l, pos in by_lit.items():
        if l < 0:
            continue
        for a in pos:
            for b in by_lit.get(-l, []):
                r = set(a) - {l} | set(b) - {-l}
                if len(r) > max_len or any(-x in r for x in r):
                    continue
                key = tuple(sorted(r))
                if key and key not in have:
                    have.add(key)
                    added.append(key)
    return n, list(clauses) + added, (lambda a: list(a))


def tovey_split(n: int, clauses: Sequence[Clause]):
    """
    Give every occurrence of a variable its own copy and tie the copies of each variable with a
    cycle of implications c1 -> c2 -> ... -> c1. Copies must be equal, so the number of solutions
    is unchanged; each copy occurs at most three times (Tovey's normal form).
    """
    occ: Dict[int, List[int]] = {}
    new_clauses: List[Clause] = []
    nxt = 0
    for c in clauses:
        out = []
        for l in c:
            nxt += 1
            occ.setdefault(abs(l), []).append(nxt)
            out.append(nxt if l > 0 else -nxt)
        new_clauses.append(tuple(out))
    first = {}
    for v in range(1, n + 1):
        copies = occ.get(v)
        if not copies:
            nxt += 1
            copies = [nxt]            # a variable that never occurs keeps one free copy
        first[v] = copies[0]
        if len(copies) > 1:
            for i, cv in enumerate(copies):
                new_clauses.append((-cv, copies[(i + 1) % len(copies)]))
    return nxt, new_clauses, (lambda a: [a[first[v] - 1] for v in range(1, n + 1)])


def definitional(n: int, clauses: Sequence[Clause]):
    """
    For each clause (a or b or c) add y <-> (a or b), i.e. (not a or y), (not b or y), (a or b or not y),
    and replace the clause by (y or c). Each y is determined, so the number of solutions is unchanged.
    """
    new: List[Clause] = []
    nxt = n
    for c in clauses:
        if len(c) < 3:
            new.append(tuple(c))
            continue
        a, b, rest = c[0], c[1], c[2:]
        nxt += 1
        y = nxt
        new += [(-a, y), (-b, y), (a, b, -y), (y,) + tuple(rest)]
    return nxt, new, (lambda s: list(s[:n]))
