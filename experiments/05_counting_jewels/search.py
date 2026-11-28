"""
search.py -- incremental WalkSAT with chosen starts, and deterministic starts by
conditional expectations (Experiment 05).

WalkSAT follows the rule of Experiment 03's control: pick a uniformly random violated
clause; flip a zero-break variable if one exists, otherwise with probability `noise`
flip a random variable of the clause, else a minimum-break variable (random ties).
Clause truth counts are maintained incrementally, so a flip costs O(occurrences).
"""

from __future__ import annotations

import math
import random
from typing import Callable, Dict, List, Optional, Sequence, Tuple

Clause = Tuple[int, ...]


def walksat(n: int, clauses: Sequence[Clause], start: Sequence[bool], max_flips: int,
            noise: float = 0.5, seed: int = 0) -> Tuple[bool, int, List[bool]]:
    """Returns (solved, flips used, final assignment)."""
    rng = random.Random(seed)
    a = list(start)
    occ_pos: List[List[int]] = [[] for _ in range(n + 1)]
    occ_neg: List[List[int]] = [[] for _ in range(n + 1)]
    for ci, c in enumerate(clauses):
        for l in c:
            (occ_pos if l > 0 else occ_neg)[abs(l)].append(ci)
    true_count = [sum((l > 0) == a[abs(l) - 1] for l in c) for c in clauses]
    unsat: List[int] = []
    where = [-1] * len(clauses)
    for ci, t in enumerate(true_count):
        if t == 0:
            where[ci] = len(unsat)
            unsat.append(ci)

    def drop(ci):
        i = where[ci]
        last = unsat.pop()
        if last != ci:
            unsat[i] = last
            where[last] = i
        where[ci] = -1

    def flip(v):
        now_false = occ_pos[v] if a[v - 1] else occ_neg[v]
        now_true = occ_neg[v] if a[v - 1] else occ_pos[v]
        a[v - 1] = not a[v - 1]
        for ci in now_false:
            true_count[ci] -= 1
            if true_count[ci] == 0:
                where[ci] = len(unsat)
                unsat.append(ci)
        for ci in now_true:
            true_count[ci] += 1
            if true_count[ci] == 1:
                drop(ci)

    def breaks(v):
        cur = occ_pos[v] if a[v - 1] else occ_neg[v]
        return sum(1 for ci in cur if true_count[ci] == 1)

    for flips in range(max_flips + 1):
        if not unsat:
            return True, flips, a
        if flips == max_flips:
            break
        c = clauses[unsat[rng.randrange(len(unsat))]]
        cand = sorted({abs(l) for l in c})
        bs = [breaks(v) for v in cand]
        zero = [v for v, b in zip(cand, bs) if b == 0]
        if zero:
            v = rng.choice(zero)
        elif rng.random() < noise:
            v = rng.choice(cand)
        else:
            low = min(bs)
            v = rng.choice([v for v, b in zip(cand, bs) if b == low])
        flip(v)
    return False, max_flips, a


def random_start(n: int, seed: int) -> List[bool]:
    rng = random.Random(10_000 + seed)
    return [rng.random() < 0.5 for _ in range(n)]


def conditional_expectation_start(n: int, clauses: Sequence[Clause], order: Sequence[int]) -> List[bool]:
    """
    Johnson's method of conditional expectations: give variables their values in `order`,
    each time choosing the value that does not increase the expected number of violated
    clauses when the remaining variables are uniformly random. The result violates at most
    as many clauses as a random assignment does on average (m/8 for three-literal clauses).
    Ties go to True.
    """
    occ: List[List[Tuple[int, bool]]] = [[] for _ in range(n + 1)]
    for ci, c in enumerate(clauses):
        for l in c:
            occ[abs(l)].append((ci, l > 0))
    unassigned = [len({abs(l) for l in c}) for c in clauses]
    satisfied = [False] * len(clauses)
    value: Dict[int, bool] = {}
    for v in order:
        delta = {True: 0.0, False: 0.0}
        seen = set()
        for ci, positive in occ[v]:
            if satisfied[ci] or ci in seen:
                continue
            seen.add(ci)
            lits = [l for l in clauses[ci] if abs(l) == v]
            before = 0.5 ** unassigned[ci]
            for choice in (True, False):
                if any((l > 0) == choice for l in lits):
                    delta[choice] += -before
                else:
                    delta[choice] += before          # 2^-(u-1) - 2^-u
        choice = delta[True] <= delta[False]
        value[v] = choice
        for ci in {ci for ci, _ in occ[v]}:
            if satisfied[ci]:
                continue
            if any(abs(l) == v and (l > 0) == choice for l in clauses[ci]):
                satisfied[ci] = True
            else:
                unassigned[ci] -= 1
    return [value.get(v, True) for v in range(1, n + 1)]


def orders(n: int, clauses: Sequence[Clause], seed: int = 0) -> Dict[str, List[int]]:
    """Variable orders for conditional expectations."""
    deg = [0] * (n + 1)
    for c in clauses:
        for l in c:
            deg[abs(l)] += 1
    rng = random.Random(seed)
    shuffled = list(range(1, n + 1))
    rng.shuffle(shuffled)
    return {
        "index": list(range(1, n + 1)),
        "degree": sorted(range(1, n + 1), key=lambda v: (-deg[v], v)),
        "golden": sorted(range(1, n + 1), key=lambda v: ((v * (math.sqrt(5) - 1) / 2) % 1.0, v)),
        "random": shuffled,
    }
