#!/usr/bin/env python3
"""Reproducible small SAT landscape experiment; no search reads planted answers.

Every charged assignment evaluation checks every clause.  Gold/zeta weights
are deliberately exploratory; shuffled controls preserve their multisets.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import random
import statistics
from pathlib import Path

import mpmath as mp
import numpy as np

SEED = 20261125
N = 16
M = 68
INSTANCES = 48
STARTS = 8
BUDGET = 1025
K = 10
AMPLITUDE = 0.6
PHI = (1 + math.sqrt(5)) / 2
METHODS = ('uniform', 'golden', 'golden_shuffled', 'zeta', 'zeta_shuffled', 'walksat')


def seed_for(*parts):
    text = ':'.join(map(str, (SEED,) + parts)).encode('ascii')
    return int.from_bytes(hashlib.blake2b(text, digest_size=8).digest(), 'big')


def label(state, n):
    return ''.join(str((state >> j) & 1) for j in range(n))


def primes(count):
    answer = []
    candidate = 2
    while len(answer) < count:
        if all(candidate % p for p in answer if p * p <= candidate):
            answer.append(candidate)
        candidate += 1
    return answer


def planted_formula(n, m, rng):
    """Uniform triples and uniform signs conditional on one planted solution."""
    plant = rng.randrange(1 << n)
    clauses = []
    seen = set()
    while len(clauses) < m:
        variables = sorted(rng.sample(range(n), 3))
        clause = tuple((j + 1) * rng.choice((-1, 1)) for j in variables)
        if clause in seen or not any(bool((plant >> (abs(lit) - 1)) & 1) == (lit > 0) for lit in clause):
            continue
        seen.add(clause)
        clauses.append(clause)
    rng.shuffle(clauses)  # The weight-to-clause correspondence varies by instance.
    return clauses, plant


class Evaluator:
    def __init__(self, n, clauses):
        self.n = n
        self.clauses = clauses
        self.variables = np.array([[abs(lit) - 1 for lit in clause] for clause in clauses], dtype=np.uint32)
        self.positive = np.array([[lit > 0 for lit in clause] for clause in clauses], dtype=bool)
        self.evaluations = 0

    def evaluate(self, states):
        values = np.asarray(states, dtype=np.uint32)
        self.evaluations += len(values)
        bits = ((values[:, None, None] >> self.variables[None, :, :]) & 1).astype(bool)
        return ~np.any(bits == self.positive[None, :, :], axis=2)


def direct_violations(state, clauses):
    return [not any(bool((state >> (abs(lit) - 1)) & 1) == (lit > 0) for lit in clause) for clause in clauses]


def make_weights(m, ordinates, rng):
    phase = rng.uniform(0, 2 * math.pi)
    phases = np.array([rng.uniform(0, 2 * math.pi) for _ in ordinates])
    j = np.arange(1, m + 1)
    golden = 1 + AMPLITUDE * np.cos(2 * math.pi * j / PHI + phase)
    logarithms = np.log(np.asarray(primes(m), dtype=float))
    zeta = 1 + AMPLITUDE * np.cos(logarithms[:, None] * np.asarray(ordinates)[None, :] + phases[None, :]).mean(axis=1)
    golden_shuffled = list(golden)
    zeta_shuffled = list(zeta)
    rng.shuffle(golden_shuffled)
    rng.shuffle(zeta_shuffled)
    weights = {'uniform': np.ones(m), 'golden': golden, 'golden_shuffled': np.array(golden_shuffled),
               'zeta': zeta, 'zeta_shuffled': np.array(zeta_shuffled), 'walksat': np.ones(m)}
    for weight in weights.values():
        assert np.all(weight >= 0.4 - 1e-12) and np.all(weight <= 1.6 + 1e-12)
    return weights, {'golden_phase': phase, 'zeta_phases': phases.tolist()}


def search(n, clauses, weights, start, rng, budget=BUDGET, walksat=False, keep_path=False):
    ev = Evaluator(n, clauses)
    state = start
    violation = ev.evaluate([state])[0]
    flips = 0
    path = [state] if keep_path else None
    reason = 'budget'
    while True:
        if not np.any(violation):
            reason = 'sat'
            break
        if walksat:
            clause_index = rng.choice(np.flatnonzero(violation).tolist())
            variables = [abs(lit) - 1 for lit in clauses[clause_index]]
        else:
            variables = list(range(n))
        if ev.evaluations + len(variables) > budget:
            break
        neighbors = [state ^ (1 << j) for j in variables]
        candidate_violations = ev.evaluate(neighbors)
        if walksat:
            breaks = np.sum(candidate_violations & ~violation[None, :], axis=1)
            free = np.flatnonzero(breaks == 0)
            if len(free):
                chosen = rng.choice(free.tolist())
            elif rng.random() < 0.5:
                chosen = rng.randrange(len(variables))
            else:
                chosen = rng.choice(np.flatnonzero(breaks == breaks.min()).tolist())
        else:
            score = float(violation @ weights)
            candidate_scores = candidate_violations @ weights
            best = float(candidate_scores.min())
            if best >= score - 1e-12:
                reason = 'trap'
                break
            chosen = rng.choice(np.flatnonzero(np.abs(candidate_scores - best) <= 1e-12).tolist())
        state = neighbors[chosen]
        violation = candidate_violations[chosen]
        flips += 1
        if keep_path:
            path.append(state)
    result = {'success': reason == 'sat', 'reason': reason, 'flips': flips,
              'assignment_evaluations': ev.evaluations, 'clause_evaluations': ev.evaluations * len(clauses),
              'final_violations': int(violation.sum()), 'final_state': state}
    assert ev.evaluations <= budget
    if keep_path:
        assert all((a ^ b).bit_count() == 1 for a, b in zip(path, path[1:]))
        result['path'] = path
        result['path_labels'] = [label(s, n) for s in path]
        result['path_violations'] = [sum(direct_violations(s, clauses)) for s in path]
        result['path_weighted_energy'] = [sum(w for w, bad in zip(weights, direct_violations(s, clauses)) if bad) for s in path]
    assert result['success'] == (not any(direct_violations(state, clauses)))
    return result


def check_instance(n, clauses, weights, sample_states):
    ev = Evaluator(n, clauses)
    violations = ev.evaluate(sample_states)
    for index, state in enumerate(sample_states):
        assert violations[index].tolist() == direct_violations(state, clauses)
    unweighted = violations.sum(axis=1)
    for weight in weights.values():
        energy = violations @ weight
        assert np.array_equal(energy == 0, unweighted == 0)
        assert np.all(energy >= 0.4 * unweighted - 1e-10)
        assert np.all(energy <= 1.6 * unweighted + 1e-10)


def summarize(runs):
    answer = {}
    for method in METHODS:
        records = [r for r in runs if r['method'] == method]
        successes = [r for r in records if r['success']]
        answer[method] = {
            'runs': len(records), 'successes': len(successes), 'success_rate': len(successes) / len(records),
            'traps': sum(r['reason'] == 'trap' for r in records),
            'budget_exhaustions': sum(r['reason'] == 'budget' for r in records),
            'mean_assignment_evaluations_all_runs': statistics.mean(r['assignment_evaluations'] for r in records),
            'mean_clause_evaluations_all_runs': statistics.mean(r['clause_evaluations'] for r in records),
            'median_assignment_evaluations_successes': statistics.median(r['assignment_evaluations'] for r in successes) if successes else None,
            'mean_flips_all_runs': statistics.mean(r['flips'] for r in records),
        }
    return answer


def find_demonstrations(ordinates):
    examples = {}
    for trial in range(3000):
        rng = random.Random(seed_for('demo', trial))
        clauses, plant = planted_formula(4, 14 + trial % 5, rng)
        weights, phases = make_weights(len(clauses), ordinates, rng)
        check_instance(4, clauses, weights, list(range(16)))
        for start in range(16):
            outcomes = {}
            for method in ('uniform', 'golden', 'zeta'):
                outcomes[method] = search(4, clauses, weights[method], start,
                    random.Random(seed_for('demo', trial, start, method)), budget=1025, keep_path=True)
            for method in ('golden', 'zeta'):
                if not outcomes['uniform']['success'] and outcomes[method]['success']:
                    key = method + '_wins'
                elif outcomes['uniform']['success'] and not outcomes[method]['success']:
                    key = method + '_loses'
                else:
                    continue
                if key not in examples:
                    examples[key] = {'trial': trial, 'n': 4, 'm': len(clauses), 'clauses': clauses,
                        'start': start, 'start_label': label(start, 4), 'weights': {k: v.tolist() for k, v in weights.items()},
                        'phases': phases, 'outcomes': outcomes,
                        'satisfying_states': [s for s in range(16) if not any(direct_violations(s, clauses))]}
        if len(examples) == 4:
            break
    assert 'golden_wins' in examples or 'zeta_wins' in examples
    return examples


def main():
    mp.mp.dps = 40
    ordinates = [float(mp.im(mp.zetazero(k))) for k in range(1, K + 1)]
    instances = []
    runs = []
    for index in range(INSTANCES):
        clauses, plant = planted_formula(N, M, random.Random(seed_for('instance', index)))
        weights, phases = make_weights(M, ordinates, random.Random(seed_for('weights', index)))
        starts = [random.Random(seed_for('start', index, j)).randrange(1 << N) for j in range(STARTS)]
        if index < 4:
            check_instance(N, clauses, weights,
                [random.Random(seed_for('check', index, j)).randrange(1 << N) for j in range(64)])
        for start_index, start in enumerate(starts):
            for method in METHODS:
                result = search(N, clauses, weights[method], start,
                    random.Random(seed_for('search', index, start_index, method)), walksat=method == 'walksat')
                result.update(instance=index, start_index=start_index, start=start, method=method)
                runs.append(result)
        instances.append({'index': index, 'clauses': clauses, 'planted_state_for_generation_only': plant,
            'starts': starts, 'weights': {k: v.tolist() for k, v in weights.items()}, 'phases': phases})
    examples = find_demonstrations(ordinates)
    output = {
        'configuration': {'seed': SEED, 'variables': N, 'clauses': M, 'instances': INSTANCES,
            'starts_per_instance': STARTS, 'budget_assignment_evaluations': BUDGET, 'zero_count': K,
            'weight_amplitude': AMPLITUDE, 'golden_ratio': PHI,
            'python_version': platform.python_version(), 'numpy_version': np.__version__, 'mpmath_version': mp.__version__},
        'methodology': {
            'generation': 'Distinct random 3-variable clauses with signs sampled uniformly conditional on one planted solution; clauses shuffled before assigning weights. Planted solution never passed to search.',
            'budget': 'Each run has an upper cap of 1025 full-assignment evaluations. Every evaluated assignment checks all 68 clauses. Strict greedy stops early at a trap; unused budget is not recycled into restarts.',
            'uniform': 'Strict best-improvement greedy over all n one-bit neighbors; random ties.',
            'golden': 'w_j=1+0.6*cos(2*pi*(j+1)/phi+phase), j zero-based.',
            'zeta': 'w_j=1+(0.6/K)*sum_k cos(gamma_k*log(p_j)+phase_k). gamma are first K zeta zero ordinates; p_j are consecutive primes.',
            'shuffled_controls': 'An independent random permutation of exactly the same weight multiset as its structured method.',
            'walksat': 'Choose a random unsatisfied clause, evaluate all 3 variable flips, choose a zero-break flip if available; otherwise with probability 0.5 choose random, else minimum break count. Random ties. Counts full clause scans for every candidate.',
            'ties_and_stopping': 'Independent seeded Python RNG for each instance/start/method. Strict improvement tolerance 1e-12; exact Boolean clause satisfaction determines success.',
            'interpretation': 'Small planted satisfiable distribution, no runtime timing or asymptotic inference. The baseline can make many more flips because it evaluates only 3 candidate neighbors at a time. RH is an optional global interpretation; computing the finite ordinates needs no RH assumption.',
            'demo_selection': 'Four-variable examples were deliberately selected to exhibit both success and failure of fixed reweighting; they are illustrative, not benchmark evidence.'},
        'zeta_ordinates': ordinates, 'summary': summarize(runs), 'instances': instances,
        'runs': runs, 'demonstrations': examples,
        'tests': ['Vectorized evaluator equals independent direct Boolean evaluator on 256 sampled 16-variable assignments.',
            'Every four-variable demonstration examined exhaustively preserves exact SAT zero sets for all six weights.',
            'Positive bounds 0.4*E <= E_w <= 1.6*E checked numerically.',
            'Every retained path consists only of one-bit moves; every reported success independently checks all clauses.',
            'All runs respect the common assignment-evaluation cap.']}
    destination = Path(__file__).with_name('black_swan_results.json')
    destination.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'summary': output['summary'], 'demonstrations': list(examples), 'output': str(destination)}, indent=2))


if __name__ == '__main__':
    main()
