#!/usr/bin/env python3
"""Independent Boolean, parity, and status audit for the hybrid experiment.

The reference evaluator and reference solvers below do not import the engine.
Optional black-box module tests use their outputs only as claims to verify.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import math
import random
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def satisfied(state, clauses):
    return all(any(bool((state >> (abs(lit)-1)) & 1) == (lit > 0)
                   for lit in clause) for clause in clauses)


def independent_parity(n, clauses):
    groups = defaultdict(set)
    normalized = []
    for clause in clauses:
        if len(clause) != 3 or len({abs(v) for v in clause}) != 3:
            normalized.append(None)
            continue
        variables = tuple(sorted(abs(v)-1 for v in clause))
        signs = {abs(v)-1: v < 0 for v in clause}
        forbidden = sum(int(signs[v]) << j for j, v in enumerate(variables))
        groups[variables].add(forbidden)
        normalized.append((variables, forbidden))
    equations, accounted = [], defaultdict(set)
    for variables, patterns in groups.items():
        for parity in (0, 1):
            forbidden_set = {s for s in range(8) if s.bit_count() % 2 == parity}
            if forbidden_set <= patterns:
                equations.append((sum(1 << v for v in variables), 1-parity))
                accounted[variables].update(forbidden_set)
    pure = all(item is not None and item[1] in accounted[item[0]] for item in normalized)
    pivots = {}
    for mask, rhs in equations:
        while mask:
            pivot = mask.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = (mask, rhs)
                break
            other_mask, other_rhs = pivots[pivot]
            mask ^= other_mask
            rhs ^= other_rhs
        if not mask and rhs:
            return {'status': 'UNSAT', 'rank': len(pivots), 'pure': pure, 'equations': equations}
    witness = 0
    for pivot, (mask, rhs) in sorted(pivots.items()):
        bit = rhs ^ ((mask & ~(1 << pivot) & witness).bit_count() % 2)
        witness |= bit << pivot
    if pure:
        assert satisfied(witness, clauses)
    return {'status': 'SAT' if pure else 'NOT_PURE_XOR', 'rank': len(pivots),
            'nullity': n-len(pivots), 'pure': pure, 'witness': witness,
            'equations': equations}


class SearchCap(Exception):
    pass


def independent_dpll(n, clauses, cap=500000):
    """Simple independent reduced-clause solver with exact status semantics."""
    original = clauses
    normalized = []
    for clause in clauses:
        literals = set(clause)
        if any(-lit in literals for lit in literals):
            continue
        normalized.append(tuple(sorted(literals)))
    nodes = 0

    def apply(work, lit):
        return tuple(tuple(x for x in c if x != -lit) for c in work if lit not in c)

    def visit(work, assigned, truth):
        nonlocal nodes
        nodes += 1
        if nodes > cap:
            raise SearchCap
        while True:
            if any(not c for c in work):
                return None
            unit = next((c[0] for c in work if len(c) == 1), None)
            if unit is None:
                break
            bit = 1 << (abs(unit)-1)
            if assigned & bit and bool(truth & bit) != (unit > 0):
                return None
            assigned |= bit
            truth = truth | bit if unit > 0 else truth & ~bit
            work = apply(work, unit)
        if not work:
            return truth
        scores = defaultdict(int)
        shortest = min(len(c) for c in work)
        for c in work:
            if len(c) == shortest:
                for lit in c:
                    scores[abs(lit)] += 1
        variable = min(scores, key=lambda v: (-scores[v], v))
        bit = 1 << (variable-1)
        for lit in (variable, -variable):
            answer = visit(apply(work, lit), assigned | bit,
                           truth | bit if lit > 0 else truth & ~bit)
            if answer is not None:
                return answer
        return None

    try:
        witness = visit(tuple(normalized), 0, 0)
    except SearchCap:
        return {'status': 'UNKNOWN', 'nodes': nodes, 'cap': cap}
    if witness is not None:
        assert satisfied(witness, original)
    return {'status': 'SAT' if witness is not None else 'UNSAT',
            'witness': witness, 'nodes': nodes, 'cap': cap}


def exact_label(n, clauses):
    if n <= 14:
        witnesses = [s for s in range(1 << n) if satisfied(s, clauses)]
        return {'status': 'SAT' if witnesses else 'UNSAT',
                'witness': witnesses[0] if witnesses else None,
                'solution_count': len(witnesses), 'method': 'exhaustive'}
    parity = independent_parity(n, clauses)
    if parity['status'] in ('SAT', 'UNSAT'):
        return dict(parity, method='independent_GF2')
    return dict(independent_dpll(n, clauses), method='independent_DPLL')


def clause_forbidding(state):
    return [-(j+1) if state >> j & 1 else j+1 for j in range(3)]


def foundational_checks():
    rng = random.Random(20261126)
    count = 0
    for n in range(1, 8):
        for _ in range(8):
            clauses = []
            for _ in range(2*n):
                variables = rng.sample(range(1, n+1), min(n, 3))
                clauses.append([v*rng.choice((-1, 1)) for v in variables])
            expected = next((s for s in range(1 << n) if satisfied(s, clauses)), None)
            result = independent_dpll(n, clauses)
            assert (result['status'] == 'SAT') == (expected is not None)
            permutation = list(range(n)); rng.shuffle(permutation)
            gauge = [rng.randrange(2) for _ in range(n)]
            transformed = [[(permutation[abs(lit)-1]+1)*(1 if lit > 0 else -1)
                            * (-1 if gauge[abs(lit)-1] else 1) for lit in c] for c in clauses]
            for state in range(1 << n):
                mapped = sum((((state >> j) & 1) ^ gauge[j]) << permutation[j] for j in range(n))
                assert satisfied(state, clauses) == satisfied(mapped, transformed)
                count += 1
    for _ in range(32):
        weights = [0.01+rng.random()*10 for _ in range(7)]
        lightest = min(range(7), key=weights.__getitem__)
        forbidden = list(range(1, 7)); forbidden.insert(lightest, 0)
        energy = [0.0]*8
        for state, weight in zip(forbidden, weights):
            energy[state] = weight
        clauses = [clause_forbidding(s) for s in forbidden]
        assert [s for s in range(8) if satisfied(s, clauses)] == [7]
        assert all(energy[0] <= energy[1 << bit] for bit in range(3))
    z = 0.37+0.63j
    import cmath
    assert z*z == (-z)*(-z)
    assert abs(cmath.exp(z)-cmath.exp(z+2j*math.pi)) < 1e-12
    return {'independent_DPLL_vs_exhaustive': 56,
            'variable_gauge_assignment_equivalences': count,
            'fixed_weight_traps_checked': 32,
            'complex_many_to_one_examples_checked': 2}


def load_module(name):
    spec = importlib.util.spec_from_file_location('audit_target_'+name, ROOT/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def xor_blackbox_checks():
    if not (ROOT/'xor_kernel.py').exists():
        return {'available': False}
    module = load_module('xor_kernel')
    checks = certificates = equations_checked = affine_sets = 0

    def verify(n, clauses, result):
        nonlocal certificates, equations_checked, affine_sets
        equations = result['equations']
        for equation in equations:
            variables = [j for j in range(n) if equation['mask'] >> j & 1]
            sources = [clauses[j] for j in equation['clauses']]
            for local in range(1 << len(variables)):
                state = sum(((local >> k) & 1) << variable for k, variable in enumerate(variables))
                assert satisfied(state, sources) == ((state & equation['mask']).bit_count() % 2 == equation['rhs'])
            equations_checked += 1
        basis = result['nullspace_basis']
        for vector in basis:
            assert vector and all((vector & eq['mask']).bit_count() % 2 == 0 for eq in equations)
        certificate = result.get('certificate')
        if certificate and certificate['kind'] == 'xor_rows':
            combination = certificate['combination']
            assert 0 < combination < (1 << len(equations))
            mask = rhs = 0
            for j, equation in enumerate(equations):
                if combination >> j & 1:
                    mask ^= equation['mask']; rhs ^= equation['rhs']
            assert mask == 0 and rhs == 1
            certificates += 1
        elif certificate:
            assert certificate['kind'] == 'empty_clause'
            assert not clauses[certificate['clause_index']]
            certificates += 1
        if result['pure'] and result['status'] == 'SAT':
            generated = {result['assignment']}
            for vector in basis:
                generated |= {state ^ vector for state in tuple(generated)}
            expected = {state for state in range(1 << n) if satisfied(state, clauses)}
            assert generated == expected
            assert len(basis) == result['nullity'] == n-result['rank']
            assert len(generated) == 1 << len(basis)
            affine_sets += 1
    for mask in range(256):
        clauses = [clause_forbidding(s) for s in range(8) if mask >> s & 1]
        expected_sat = any(satisfied(s, clauses) for s in range(8))
        result = module.solve_pure_xor(3, clauses)
        verify(3, clauses, result['metadata'])
        if result['status'] == 'SAT':
            assert expected_sat and satisfied(result['witness'], clauses), (mask, result)
        elif result['status'] == 'UNSAT':
            assert not expected_sat, (mask, result)
        else:
            assert result['status'] == 'NOT_PURE_XOR', result
        checks += 1
    rng = random.Random(82317)
    random_checks = 0
    for n in range(1, 9):
        for _ in range(24):
            clauses = []
            for _ in range(rng.randrange(n+3)):
                variables = sorted(rng.sample(range(n), rng.randrange(1, min(n, 3)+1)))
                rhs = rng.randrange(2)
                for forbidden in range(1 << len(variables)):
                    if forbidden.bit_count() % 2 != rhs:
                        clauses.append([-(v+1) if forbidden >> j & 1 else v+1 for j, v in enumerate(variables)])
            if clauses and rng.randrange(2):
                clauses += [clauses[0].copy()]
            for mixed in (False, True):
                trial = [c.copy() for c in clauses]
                if mixed:
                    trial.append([rng.choice((-1, 1))*(v+1) for v in rng.sample(range(n), min(n, 3))])
                result = module.analyze_cnf(n, trial)
                verify(n, trial, result)
                truth = [s for s in range(1 << n) if satisfied(s, trial)]
                if result['status'] == 'SAT':
                    assert result['assignment'] in truth
                elif result['status'] == 'UNSAT':
                    assert not truth
                random_checks += 1
    return {'available': True, 'all_256_three_variable_clause_subsets_checked': checks,
            'random_parity_or_mixed_formulas': random_checks,
            'source_CNF_equation_equivalences': equations_checked,
            'contradiction_certificates': certificates, 'complete_affine_sets': affine_sets}


def engine_blackbox_checks():
    if not (ROOT/'run_hybrid.py').exists():
        return {'available': False}
    engine = load_module('run_hybrid')
    fixtures = [(1, [[1, 1], [-1]], 'UNSAT'),
                (1, [[1, -1], [-1]], 'SAT'),
                (2, [[1, 1, 2], [-1], [-2]], 'UNSAT'),
                (1, [[]], 'UNSAT'), (0, [], 'SAT'),
                (3, [[1, 2, 3], [-1, -2, -3]], 'SAT')]
    for n, clauses, expected in fixtures:
        result = engine.dpll(n, clauses)
        assert result['status'] == expected, (n, clauses, result)
        if expected == 'SAT':
            assert satisfied(result['witness'], clauses)
    rng = random.Random(992613)
    random_checks = 0
    for n in range(1, 8):
        for _ in range(16):
            clauses = [[rng.choice((-1, 1))*rng.randrange(1, n+1)
                        for _ in range(rng.randrange(4))]
                       for _ in range(rng.randrange(2*n+3))]
            truth = next((s for s in range(1 << n) if satisfied(s, clauses)), None)
            result = engine.dpll(n, clauses)
            assert (result['status'] == 'SAT') == (truth is not None), (clauses, result)
            if truth is not None:
                assert satisfied(result['witness'], clauses)
            random_checks += 1
    unsat = [[1, 2], [-1, 2], [1, -2], [-1, -2]]
    assert engine.dpll(2, unsat, node_cap=1)['status'] == 'UNKNOWN'
    assert engine.dpll(2, unsat, node_cap=100)['status'] == 'UNSAT'
    reductions = 0
    for width in range(4, 8):
        for _ in range(8):
            clause = [(j+1)*rng.choice((-1, 1)) for j in range(width)]
            augmented_n, encoded, metadata = engine.to_three_cnf(width, [clause])
            assert augmented_n == width+width-3
            assert all(1 <= len(c) <= 3 for c in encoded)
            assert metadata[0]['raw_clause'] == clause
            for state in range(1 << width):
                exists = any(satisfied(state | (aux << width), encoded)
                             for aux in range(1 << (augmented_n-width)))
                assert exists == satisfied(state, [clause])
                reductions += 1
    multi = [[1, 2, 3, 4], [-1, -2, 4, 5], [1, -5]]
    augmented_n, encoded, metadata = engine.to_three_cnf(5, multi)
    all_aux = [v for group in metadata for v in group['aux_variables']]
    assert len(set(all_aux)) == len(all_aux)
    for state in range(32):
        exists = any(satisfied(state | (aux << 5), encoded)
                     for aux in range(1 << (augmented_n-5)))
        assert exists == satisfied(state, multi)
        reductions += 1
    return {'available': True, 'general_CNF_edge_fixtures': len(fixtures),
            'random_small_formulas_including_duplicates_tautologies_empty_clauses': random_checks,
            'budget_status_checks': 2, 'existential_3CNF_reduction_assignments': reductions}


def verify_dataset_kernel(n, clauses, claim):
    equations = claim['equations']
    for equation in equations:
        variables = [j for j in range(n) if equation['mask'] >> j & 1]
        sources = [clauses[j] for j in equation['clauses']]
        assert all(abs(lit)-1 in variables for c in sources for lit in c)
        for local in range(1 << len(variables)):
            state = sum(((local >> j) & 1) << v for j, v in enumerate(variables))
            assert satisfied(state, sources) == ((state & equation['mask']).bit_count() % 2 == equation['rhs'])
    cert = claim.get('certificate')
    if cert and cert['kind'] == 'xor_rows':
        mask = rhs = 0
        assert 0 < cert['combination'] < 1 << len(equations)
        for j, equation in enumerate(equations):
            if cert['combination'] >> j & 1:
                mask ^= equation['mask']; rhs ^= equation['rhs']
        assert mask == 0 and rhs == 1
    elif cert:
        assert cert['kind'] == 'empty_clause' and not clauses[cert['clause_index']]
    basis = claim['nullspace_basis']
    reduced = {}
    for vector in basis:
        assert 0 < vector < 1 << n
        assert all((vector & eq['mask']).bit_count() % 2 == 0 for eq in equations)
        residual = vector
        while residual:
            pivot = residual.bit_length()-1
            if pivot not in reduced:
                reduced[pivot] = residual
                break
            residual ^= reduced[pivot]
        assert residual, 'Dependent reported nullspace basis'
    assert len(basis) == claim['nullity'] == n-claim['rank']
    return bool(cert)


def main():
    started = time.monotonic()
    report = {'foundational_checks': foundational_checks(), 'xor_blackbox': xor_blackbox_checks(),
              'engine_blackbox': engine_blackbox_checks(),
              'mismatches': [], 'interpretation': 'Independent finite correctness audit; no polynomial worst-case conclusion.'}
    dataset_path = ROOT/'dataset.json'
    if dataset_path.exists():
        dataset = json.loads(dataset_path.read_text())
        labels = {}
        kernel_module = load_module('xor_kernel') if (ROOT/'xor_kernel.py').exists() else None
        kernel_claims, graph_checks, affine_members = {}, [], 0
        dataset_certificates = 0
        for instance in dataset['instances']:
            result = exact_label(instance['n'], instance['clauses'])
            labels[str(instance['id'])] = result
            expected = instance.get('expected_status')
            if expected and result['status'] != 'UNKNOWN':
                assert result['status'] == expected, (instance['id'], expected, result)
            witness = instance.get('expected_witness')
            if witness is not None:
                assert satisfied(witness, instance['clauses'])
            if kernel_module:
                kernel = kernel_module.analyze_cnf(instance['n'], instance['clauses'])
                dataset_certificates += verify_dataset_kernel(instance['n'], instance['clauses'], kernel)
                kernel_claims[instance['id']] = kernel
                if instance['family'] in ('xor_prism', 'xor_cubic'):
                    construction = instance['construction']
                    vertices, edges = construction['vertices'], construction['edges']
                    assert kernel['pure'] and kernel['rank'] == vertices-1
                    assert kernel['nullity'] == len(edges)-vertices+1
                    assert kernel['status'] == ('UNSAT' if construction['charge_sum_mod2'] else 'SAT')
                    graph_checks.append({'instance': instance['id'], 'vertices': vertices,
                                         'edges': len(edges), 'rank': kernel['rank'],
                                         'cycle_nullity': kernel['nullity'], 'status': kernel['status']})
                    if kernel['status'] == 'SAT':
                        generated = {kernel['assignment']}
                        for vector in kernel['nullspace_basis']:
                            generated |= {s ^ vector for s in tuple(generated)}
                        assert len(generated) == 1 << kernel['nullity']
                        assert all(satisfied(s, instance['clauses']) for s in generated)
                        affine_members += len(generated)
        report['independent_labels'] = labels
        report['instances_checked'] = len(labels)
        report['dataset_XOR_certificates_checked'] = dataset_certificates
        report['graph_cycle_kernel_checks'] = graph_checks
        report['structured_affine_family_members_checked'] = affine_members
        results_path = ROOT/'results.json'
        if results_path.exists():
            results = json.loads(results_path.read_text())
            by_id = {str(i['id']): i for i in dataset['instances']}
            checked = sat_count = unsat_count = unknown_count = 0
            for run in results['runs']:
                instance = by_id[str(run['instance'])]
                exact = labels[str(run['instance'])]
                status = run['status']
                if status == 'SAT':
                    assert run.get('witness') is not None
                    assert satisfied(run['witness'], instance['clauses']), run
                    assert exact['status'] != 'UNSAT', run
                    sat_count += 1
                elif status == 'UNSAT':
                    assert exact['status'] == 'UNSAT', (run, exact)
                    unsat_count += 1
                else:
                    assert status == 'UNKNOWN', run
                    unknown_count += 1
                checked += 1
                assert run['dpll_nodes'] <= dataset['node_cap']
                assert run['frontend_evaluations'] <= dataset['frontend_eval_cap']
                if run.get('kernel'):
                    independent_claim = kernel_claims[run['instance']]
                    for key in ('status', 'pure', 'rank', 'nullity', 'contradiction_combination'):
                        assert run['kernel'][key] == independent_claim.get(key)
            report['runs_checked'] = checked
            report['SAT_witnesses_checked'] = sat_count
            report['UNSAT_answers_independently_checked'] = unsat_count
            report['UNKNOWN_answers'] = unknown_count
    report['elapsed_seconds'] = round(time.monotonic()-started, 4)
    destination = ROOT/'audit_results.json'
    destination.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key != 'independent_labels'}, indent=2))


if __name__ == '__main__':
    main()
