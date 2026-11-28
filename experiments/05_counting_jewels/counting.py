"""
counting.py -- exact counting by tensor contraction (Experiment 05).

For a CNF formula with m clauses, the energy polynomial is

    Z(z) = sum over assignments x of z^E(x) = sum_k N_k z^k,

where E(x) is the number of violated clauses and N_k counts the assignments with
exactly k violated clauses. N_0 is the number of solutions, and Z(1) = 2^n.

Z is the contraction of a Boolean tensor network with one factor per clause: 1 when
the clause is satisfied and z when it is violated. Summing variables out one at a time
(bucket elimination) costs about 2^(width) per step, where width is the induced width
of the elimination order. Arithmetic is exact: tensors hold polynomials in z truncated
after z^(L-1), reduced modulo 31-bit primes, and coefficients are rebuilt with the
Chinese remainder theorem. This is the "count, do not search" angle of the 2026
discrete-log paper (its truncated ring (Z/2^P)[z]/(z^l)) carried out with classical,
checkable tools only.
"""

from __future__ import annotations

import itertools
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

Clause = Tuple[int, ...]


def _is_prime(q: int) -> bool:
    if q < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if q % p == 0:
            return q == p
    d, s = q - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):           # deterministic for q < 3.4e14
        x = pow(a, d, q)
        if x in (1, q - 1):
            continue
        for _ in range(s - 1):
            x = x * x % q
            if x == q - 1:
                break
        else:
            return False
    return True


PRIMES: List[int] = []
_q = (1 << 31) - 1
while len(PRIMES) < 12:
    if _is_prime(_q):
        PRIMES.append(_q)
    _q -= 2


# ---------------------------------------------------------------------------
# Brute force (reference, n <= 26)
# ---------------------------------------------------------------------------
def energy_counts_bruteforce(n: int, clauses: Sequence[Clause]) -> List[int]:
    """N_k for k = 0..m by enumerating all 2^n assignments (vectorized)."""
    if n > 26:
        raise ValueError("brute force is limited to n <= 26")
    m = len(clauses)
    total = np.zeros(m + 1, dtype=np.int64)
    chunk_bits = min(n, 22)
    low = np.arange(1 << chunk_bits, dtype=np.uint32)
    for hi in range(1 << (n - chunk_bits)):
        idx = low | np.uint32(hi << chunk_bits)
        viol = np.zeros(idx.size, dtype=np.int16)
        for c in clauses:
            ok = np.zeros(idx.size, dtype=bool)
            for l in c:
                bit = (idx >> np.uint32(abs(l) - 1)) & np.uint32(1)
                ok |= (bit == 1) if l > 0 else (bit == 0)
            viol += ~ok
        total += np.bincount(viol, minlength=m + 1)
    return [int(x) for x in total]


# ---------------------------------------------------------------------------
# Graphs, elimination orders, widths
# ---------------------------------------------------------------------------
def primal_graph(n: int, clauses: Sequence[Clause]) -> Dict[int, set]:
    adj = {v: set() for v in range(1, n + 1)}
    for c in clauses:
        vs = {abs(l) for l in c}
        for a in vs:
            adj[a] |= vs - {a}
    return adj


def induced_width(adj: Dict[int, set], order: Sequence[int]) -> int:
    live = {v: set(s) for v, s in adj.items()}
    width = 0
    for v in order:
        nb = live.pop(v)
        width = max(width, len(nb))
        for a in nb:
            live[a].discard(v)
            live[a] |= nb - {a}
    return width


def greedy_order(adj: Dict[int, set], rule: str = "min_fill", seed: int = 0) -> List[int]:
    """Min-degree or min-fill elimination order (ties broken by a seeded random key)."""
    rng = np.random.default_rng(seed)
    tie = {v: float(rng.random()) for v in adj}
    live = {v: set(s) for v, s in adj.items()}
    order = []
    while live:
        if rule == "min_degree":
            v = min(live, key=lambda u: (len(live[u]), tie[u]))
        else:
            def fill(u):
                nb = list(live[u])
                return sum(1 for i in range(len(nb)) for j in range(i + 1, len(nb)) if nb[j] not in live[nb[i]])
            v = min(live, key=lambda u: (fill(u), len(live[u]), tie[u]))
        nb = live.pop(v)
        for a in nb:
            live[a].discard(v)
            live[a] |= nb - {a}
        order.append(v)
    return order


def best_order(n: int, clauses: Sequence[Clause], tries: Optional[int] = None) -> Tuple[List[int], int]:
    adj = primal_graph(n, clauses)
    tries = tries if tries is not None else (3 if n <= 80 else 1)
    best = None
    for rule in ("min_fill", "min_degree"):
        for s in range(tries):
            o = greedy_order(adj, rule, seed=s)
            w = induced_width(adj, o)
            if best is None or w < best[1]:
                best = (o, w)
    return best


def choose_pins(adj: Dict[int, set], order: Sequence[int], max_width: int) -> Tuple[List[int], int]:
    """
    Slicing: simulate the elimination order, and while some step's clique is wider than
    max_width, pin the variable that occurs most often in those oversized cliques.
    The order itself is kept; pinned variables are simply removed from it.
    """
    pins: List[int] = []
    while True:
        live = {v: set(s) - set(pins) for v, s in adj.items() if v not in pins}
        width, freq = 0, {}
        for v in order:
            if v in pins:
                continue
            nb = live.pop(v)
            width = max(width, len(nb))
            if len(nb) > max_width:
                for u in nb | {v}:
                    freq[u] = freq.get(u, 0) + 1
            for a in nb:
                live[a].discard(v)
                live[a] |= nb - {a}
        if width <= max_width:
            return pins, width
        pins.append(max(freq, key=lambda u: (freq[u], -u)))


# ---------------------------------------------------------------------------
# Contraction over (Z/p)[z]/(z^L)
# ---------------------------------------------------------------------------
def _poly_mul(a: np.ndarray, b: np.ndarray, p: int, L: int) -> np.ndarray:
    shape = np.broadcast_shapes(a.shape, b.shape)
    out = np.zeros(shape, dtype=np.int64)
    for i in range(L):
        ai = a[..., i:i + 1]
        if not ai.any():
            continue
        out[..., i:] = (out[..., i:] + (ai * b[..., :L - i]) % p) % p
    return out


def _clause_factor(c: Clause, L: int) -> Tuple[Tuple[int, ...], np.ndarray]:
    vs = tuple(sorted({abs(l) for l in c}))
    t = np.zeros((2,) * len(vs) + (L,), dtype=np.int64)
    for bits in itertools.product((0, 1), repeat=len(vs)):
        val = dict(zip(vs, bits))
        sat = any((l > 0) == bool(val[abs(l)]) for l in c)
        if sat:
            t[bits + (0,)] = 1
        elif L > 1:
            t[bits + (1,)] = 1
    return vs, t


def _eliminate(n_vars: Sequence[int], factors, order: Sequence[int], p: int, L: int) -> np.ndarray:
    pos = {v: i for i, v in enumerate(order)}
    buckets: Dict[int, list] = {}
    const = np.zeros(L, dtype=np.int64)
    const[0] = 1
    for vs, arr in factors:
        if not vs:
            const = _poly_mul(const, arr % p, p, L)
            continue
        buckets.setdefault(min(vs, key=pos.__getitem__), []).append((vs, arr % p))
    for v in order:
        fs = buckets.pop(v, [])
        if not fs:
            const = const * 2 % p
            continue
        union = sorted(set().union(*[set(f[0]) for f in fs]))
        prod = None
        for vs, arr in fs:
            a = arr.reshape([2 if u in vs else 1 for u in union] + [L])
            prod = a if prod is None else _poly_mul(prod, a, p, L)
        summed = prod.sum(axis=union.index(v)) % p
        rest = tuple(u for u in union if u != v)
        if rest:
            buckets.setdefault(min(rest, key=pos.__getitem__), []).append((rest, summed))
        else:
            const = _poly_mul(const, summed, p, L)
    return const


def _crt(residues: Sequence[int], primes: Sequence[int]) -> int:
    x, mod = 0, 1
    for r, p in zip(residues, primes):
        t = ((r - x) * pow(mod, -1, p)) % p
        x += mod * t
        mod *= p
    return x


def condition(clauses: Sequence[Clause], pins: Dict[int, bool]) -> Tuple[List[Clause], int]:
    """Apply a partial assignment. Returns the shortened clauses and the number of clauses it already violates."""
    out, dead = [], 0
    for c in clauses:
        if any(abs(l) in pins and pins[abs(l)] == (l > 0) for l in c):
            continue
        rest = tuple(l for l in c if abs(l) not in pins)
        if rest:
            out.append(rest)
        else:
            dead += 1
    return out, dead


def energy_polynomial(n: int, clauses: Sequence[Clause], L: Optional[int] = None,
                      order: Optional[Sequence[int]] = None, max_width: int = 20,
                      pins: Optional[Sequence[int]] = None) -> Dict:
    """
    Exact N_0..N_{L-1} (L defaults to m + 1, the whole polynomial) by contraction.
    If the elimination width exceeds max_width, high-degree variables are pinned and the
    2^k conditioned formulas are contracted separately (slicing).
    """
    m = len(clauses)
    L = m + 1 if L is None else L
    primes = PRIMES[:max(1, math.ceil((n + 2) / 30))]
    adj = primal_graph(n, clauses)
    if order is None:
        order, width = best_order(n, clauses)
    else:
        width = induced_width(adj, order)
    if pins is None:
        pins, width = choose_pins(adj, order, max_width) if width > max_width else ([], width)
    pins = list(pins)
    free_vars = [v for v in order if v not in pins]
    total = [0] * L
    for bits in itertools.product((False, True), repeat=len(pins)):
        sub_clauses, dead = condition(clauses, dict(zip(pins, bits)))
        if dead >= L:
            continue
        factors = [_clause_factor(c, L) for c in sub_clauses]
        residues = []
        for p in primes:
            r = _eliminate(free_vars, factors, free_vars, p, L)
            residues.append(r)
        for k in range(L - dead):
            total[k + dead] += _crt([int(r[k]) for r in residues], primes)
    return {"N": total, "width": width, "pins": len(pins), "primes": len(primes)}


def contract(variables: Sequence[int], clauses: Sequence[Clause], L: int, order: Sequence[int]) -> List[int]:
    """Exact N_0..N_{L-1} summing over `variables` only (every clause must use only these), in a given order."""
    primes = PRIMES[:max(1, math.ceil((len(variables) + 2) / 30))]
    factors = [_clause_factor(c, L) for c in clauses]
    order = [v for v in order if v in set(variables)]
    residues = [_eliminate(order, factors, order, p, L) for p in primes]
    return [_crt([int(r[k]) for r in residues], primes) for k in range(L)]


def count_solutions(n: int, clauses: Sequence[Clause], **kw) -> int:
    return energy_polynomial(n, clauses, L=1, **kw)["N"][0]


def solve_by_counting(n: int, clauses: Sequence[Clause], order: Optional[Sequence[int]] = None) -> Optional[List[bool]]:
    """
    Self-reduction: fix x_1, x_2, ... in turn, keeping a value whose completions still count > 0.
    One count decides satisfiability; n more counts turn it into a witness, with no search and
    no randomness. Every count reuses one elimination order restricted to the free variables.
    """
    if order is None:
        order, _ = best_order(n, clauses)
    if contract(list(range(1, n + 1)), clauses, 1, order)[0] == 0:
        return None
    pins: Dict[int, bool] = {}
    for v in range(1, n + 1):
        trial = {**pins, v: False}
        sub, dead = condition(clauses, trial)
        free = [u for u in range(1, n + 1) if u not in trial]
        ok = dead == 0 and contract(free, sub, 1, order)[0] > 0
        pins[v] = not ok
    return [pins[v] for v in range(1, n + 1)]
