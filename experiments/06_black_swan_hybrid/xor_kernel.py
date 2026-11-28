#!/usr/bin/env python3
"""Recognize small CNF XOR groups and solve their exact GF(2) kernel.

This module does not relax OR constraints or merge Boolean states. A mixed
formula is retained unchanged unless a recognized subset is contradictory.
Only groups of width at most three are enumerated by the detector.
"""
from __future__ import annotations

import json
import random
from collections import defaultdict


def parity(value: int) -> int:
    return value.bit_count() & 1


def satisfies(assignment: int, clauses) -> bool:
    """Direct CNF evaluation, independent of the linear elimination code."""
    return all(any(bool(assignment & (1 << (abs(lit) - 1))) == (lit > 0)
                   for lit in clause) for clause in clauses)


def solve_gf2(n: int, equations) -> dict:
    """Solve mask*x=rhs over GF(2), retaining a row-combination certificate.

    An equation is a dict with integer ``mask`` and binary ``rhs``. A returned
    contradiction_combination selects original equations whose XOR is 0=1.
    A nullspace vector is an integer bit mask over the n variables.
    """
    basis = {}
    contradiction = None
    for index, equation in enumerate(equations):
        mask, rhs, combination = int(equation['mask']), int(equation['rhs']), 1 << index
        if mask < 0 or mask >= (1 << n) or rhs not in (0, 1):
            raise ValueError('Invalid GF(2) equation')
        while mask:
            pivot = mask.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = (mask, rhs, combination)
                break
            old_mask, old_rhs, old_combination = basis[pivot]
            mask ^= old_mask
            rhs ^= old_rhs
            combination ^= old_combination
        if not mask and rhs and contradiction is None:
            contradiction = combination

    ordered = sorted(basis)
    particular = 0
    for pivot in ordered:
        mask, rhs, _ = basis[pivot]
        if parity(mask & particular) ^ rhs:
            particular |= 1 << pivot

    nullspace = []
    for free in range(n):
        if free in basis:
            continue
        vector = 1 << free
        for pivot in ordered:
            mask, _, _ = basis[pivot]
            if parity(mask & vector):
                vector |= 1 << pivot
        nullspace.append(vector)

    return {
        'consistent': contradiction is None,
        'rank': len(basis),
        'nullity': n - len(basis),
        'particular_solution': particular if contradiction is None else None,
        'nullspace_basis': nullspace,
        'contradiction_combination': contradiction,
        'pivot_rows': [{'pivot': pivot, 'mask': basis[pivot][0],
                        'rhs': basis[pivot][1], 'combination': basis[pivot][2]}
                       for pivot in ordered],
    }


def analyze_cnf(n: int, clauses) -> dict:
    """Return SAT/UNSAT/UNKNOWN using only exactly recognized XOR structure.

    Variable j is integer bit j; a signed literal +/-k refers to bit k-1.
    Pure, consistent XOR CNFs return a full SAT assignment. A consistent mixed
    CNF returns UNKNOWN, even if its linear subset has an assignment. An XOR
    contradiction proves UNSAT for a mixed CNF as well. Empty clauses are
    detected directly. Wider clauses and partial groups remain uncovered.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError('n must be a nonnegative integer')
    original = [tuple(clause) for clause in clauses]
    groups = defaultdict(lambda: defaultdict(list))
    ignored = []
    active = set()
    for index, clause in enumerate(original):
        if any(not isinstance(lit, int) or isinstance(lit, bool)
               or lit == 0 or abs(lit) > n for lit in clause):
            raise ValueError('A literal is outside the declared variable range')
        values = set(clause)
        if any(-lit in values for lit in values):
            ignored.append(index)
            continue
        if not values:
            return {'status': 'UNSAT', 'pure': False, 'assignment': None,
                    'equations': [], 'rank': 0, 'nullity': n,
                    'particular_solution': None, 'nullspace_basis': [],
                    'contradiction_combination': None,
                    'certificate': {'kind': 'empty_clause', 'clause_index': index},
                    'covered_clauses': [], 'uncovered_clauses': [index],
                    'ignored_tautologies': ignored,
                    'metadata': {'detected_equations': 0, 'empty_clause': True}}
        active.add(index)
        variables = tuple(sorted(abs(lit) - 1 for lit in values))
        if len(variables) > 3:
            continue
        forbidden = sum((1 << j) for j, variable in enumerate(variables)
                        if -(variable + 1) in values)
        groups[variables][forbidden].append(index)

    covered = set()
    equations = []
    for variables, patterns in sorted(groups.items()):
        width = len(variables)
        for forbidden_parity in (0, 1):
            required = [value for value in range(1 << width)
                        if parity(value) == forbidden_parity]
            if not all(value in patterns for value in required):
                continue
            source = sorted(index for value in required for index in patterns[value])
            covered.update(source)
            equations.append({'mask': sum(1 << variable for variable in variables),
                              'rhs': 1 ^ forbidden_parity,
                              'variables': list(variables),
                              'clauses': source,
                              'forbidden_patterns': required})

    linear = solve_gf2(n, equations)
    pure = active <= covered
    if not linear['consistent']:
        status, assignment = 'UNSAT', None
        certificate = {'kind': 'xor_rows',
                       'combination': linear['contradiction_combination']}
    elif pure:
        status, assignment = 'SAT', linear['particular_solution']
        assert satisfies(assignment, original)
        certificate = None
    else:
        status, assignment, certificate = 'UNKNOWN', None, None
    return {
        'status': status, 'pure': pure, 'assignment': assignment,
        'equations': equations, **linear, 'certificate': certificate,
        'covered_clauses': sorted(covered),
        'uncovered_clauses': sorted(active - covered),
        'ignored_tautologies': ignored,
        'metadata': {'detected_equations': len(equations),
                     'covered_clause_count': len(covered),
                     'uncovered_clause_count': len(active - covered),
                     'rank': linear['rank'], 'nullity': linear['nullity']},
    }


def solve_pure_xor(n: int, clauses) -> dict:
    """Small compatibility wrapper; mixed consistent formulas stay untouched."""
    result = analyze_cnf(n, clauses)
    return {'status': result['status'] if result['status'] != 'UNKNOWN' else 'NOT_PURE_XOR',
            'pure': result['pure'], 'witness': result['assignment'], 'metadata': result}


def xor_clauses(variables, rhs: int):
    """Encode an XOR equation as clauses; variables are zero-based indices."""
    variables = tuple(sorted(variables))
    return [[-(variable + 1) if value & (1 << j) else variable + 1
             for j, variable in enumerate(variables)]
            for value in range(1 << len(variables)) if parity(value) != rhs]


def verify_result(n: int, clauses, result: dict) -> None:
    """Check the returned witness, nullspace, and contradiction certificate."""
    equations = result['equations']
    if result['status'] == 'SAT':
        assert satisfies(result['assignment'], clauses)
    for vector in result['nullspace_basis']:
        assert vector and all(parity(eq['mask'] & vector) == 0 for eq in equations)
    if result.get('consistent'):
        assert all(parity(eq['mask'] & result['particular_solution']) == eq['rhs']
                   for eq in equations)
    cert = result['certificate']
    if cert and cert['kind'] == 'xor_rows':
        mask = rhs = 0
        for j, eq in enumerate(equations):
            if cert['combination'] & (1 << j):
                mask ^= eq['mask']
                rhs ^= eq['rhs']
        assert mask == 0 and rhs == 1
    elif cert:
        assert not clauses[cert['clause_index']]


def self_test() -> dict:
    rng = random.Random(20261126)
    checks = {'three_variable_clause_subsets': 0, 'random_parity_systems': 0,
              'mixed_formulas': 0, 'edge_cases': 0, 'large_variable_systems': 0}
    clauses8 = [[-(j + 1) if state & (1 << j) else j + 1 for j in range(3)]
                for state in range(8)]
    for subset in range(256):
        clauses = [clause for j, clause in enumerate(clauses8) if subset & (1 << j)]
        result = analyze_cnf(3, clauses)
        verify_result(3, clauses, result)
        satisfying = [state for state in range(8) if satisfies(state, clauses)]
        if result['status'] == 'UNSAT':
            assert not satisfying
        if result['status'] == 'SAT':
            assert result['assignment'] in satisfying
        for eq in result['equations']:
            assert all(parity(eq['mask'] & state) == eq['rhs'] for state in satisfying)
        checks['three_variable_clause_subsets'] += 1
    for n in range(1, 9):
        for _ in range(100):
            original_equations = []
            clauses = []
            for _ in range(rng.randrange(n + 5)):
                variables = rng.sample(range(n), rng.randrange(1, min(3, n) + 1))
                rhs = rng.randrange(2)
                original_equations.append((sum(1 << j for j in variables), rhs))
                clauses.extend(xor_clauses(variables, rhs))
            rng.shuffle(clauses)
            if clauses and rng.randrange(2):
                clauses.append(clauses[0].copy())
            result = analyze_cnf(n, clauses)
            verify_result(n, clauses, result)
            satisfying = [state for state in range(1 << n)
                          if all(parity(mask & state) == rhs
                                 for mask, rhs in original_equations)]
            assert result['pure'] and (result['status'] == 'SAT') == bool(satisfying)
            assert all(satisfies(state, clauses) == (state in satisfying)
                       for state in range(1 << n))
            if satisfying:
                basis = result['nullspace_basis']
                generated = {result['assignment']}
                for vector in basis:
                    generated |= {state ^ vector for state in list(generated)}
                assert generated == set(satisfying)
            checks['random_parity_systems'] += 1
            if n >= 3:
                mixed = clauses + [[(v + 1) * rng.choice((-1, 1))
                                    for v in rng.sample(range(n), 3)]]
                result = analyze_cnf(n, mixed)
                verify_result(n, mixed, result)
                ground = [state for state in range(1 << n) if satisfies(state, mixed)]
                if result['status'] == 'SAT':
                    assert result['assignment'] in ground
                elif result['status'] == 'UNSAT':
                    assert not ground
                checks['mixed_formulas'] += 1
    for n, clauses in [(0, []), (0, [[]]), (3, [[1, -1, 2]]),
                       (3, [[1, 1, 1]]), (3, [[1, 2, 3]]),
                       (4, [[1, 2, 3, 4]])]:
        result = analyze_cnf(n, clauses)
        verify_result(n, clauses, result)
        checks['edge_cases'] += 1
    assert analyze_cnf(3, [[1, 2, 3]])['status'] == 'UNKNOWN'
    assert analyze_cnf(4, [[1, 2, 3, 4]])['status'] == 'UNKNOWN'
    large_clauses = [[j + 1] for j in range(80)]
    large_result = analyze_cnf(80, large_clauses)
    verify_result(80, large_clauses, large_result)
    assert large_result['rank'] == 80 and large_result['assignment'] == (1 << 80) - 1
    checks['large_variable_systems'] += 1
    large_clauses.append([-80])
    large_result = analyze_cnf(80, large_clauses)
    verify_result(80, large_clauses, large_result)
    assert large_result['status'] == 'UNSAT'
    checks['large_variable_systems'] += 1
    checks['status'] = 'passed'
    return checks


if __name__ == '__main__':
    print(json.dumps(self_test(), indent=2))
