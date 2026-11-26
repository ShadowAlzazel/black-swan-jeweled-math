#!/usr/bin/env python3
"""Exhaust the 256 subsets of the eight full-width three-variable clauses.

This is a restricted, explicitly finite family, not all 3SAT formulas at larger n.
Also construct a strict-local-search trap for every fixed positive weight vector.
"""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'tiling_research'
if not SOURCE.exists():
    SOURCE = ROOT.parent / 'tiling_research'
PHI = (1 + math.sqrt(5)) / 2
plastic = 1.3
for _ in range(12):
    plastic -= (plastic ** 3 - plastic - 1) / (3 * plastic ** 2 - 1)
CONSTANTS = {'phi': PHI, 'phi_squared': PHI ** 2, 'inverse_phi': 1 / PHI,
             'sqrt2': math.sqrt(2), 'sqrt3': math.sqrt(3), 'pi': math.pi,
             'e': math.e, 'zeta2': math.pi ** 2 / 6, 'plastic': plastic,
             'rational_13_8': 13 / 8, 'fixed_random': 3.5337708052156676}


def orientations(filename, level):
    data = json.loads((SOURCE / filename).read_text())
    patch = next(p for p in data['patches'] if p['iteration'] == level)
    return [math.radians(t['orientationDegrees']) for t in patch['tiles']]


ANGLES = {'hat': orientations('hat_patches.json', 2),
          'pinwheel': orientations('pinwheel_patches.json', 4),
          'square': [j * math.pi / 2 for j in range(8)]}


def weights(config, m):
    if config['name'] == 'uniform':
        return [1.0] * m
    return [1 + config['amplitude'] * math.cos(ANGLES[config['tiling']][j]
            + 2 * math.pi * (j + 1) / config['constant'] + .37) for j in range(m)]


def descend(energy, start):
    state, path = start, [start]
    while energy[state] > 0:
        value, bit = min((energy[state ^ (1 << bit)], bit) for bit in range(3))
        if value >= energy[state] - 1e-12:
            break
        state ^= 1 << bit
        path.append(state)
    assert len(set(path)) == len(path)
    return path


def clause_forbidding(state):
    return [-(j + 1) if (state >> j) & 1 else j + 1 for j in range(3)]


def verify_clauses(clauses, state):
    return all(any(bool((state >> (abs(lit) - 1)) & 1) == (lit > 0)
                   for lit in clause) for clause in clauses)


def adversarial(config):
    w = weights(config, 7)
    lightest = min(range(7), key=lambda j: w[j])
    order = [None] * 7
    order[lightest] = 0
    other = iter(range(1, 7))
    for j in range(7):
        if order[j] is None:
            order[j] = next(other)
    clauses = [clause_forbidding(s) for s in order]
    energy = [0.0] * 8
    for j, forbidden in enumerate(order):
        energy[forbidden] = w[j]
    assert [s for s in range(8) if verify_clauses(clauses, s)] == [7]
    assert descend(energy, 0) == [0]
    assert all(energy[0] <= energy[1 << bit] + 1e-12 for bit in range(3))
    return {'config': config['name'], 'weights': w, 'clause_order_forbidden_states': order,
            'clauses': clauses, 'energy': energy, 'start': 0, 'unique_solution': 7}


def breakout(example):
    order, w = example['clause_order_forbidden_states'], example['weights'].copy()
    state, boosts, path = 0, 0, [0]
    while state != 7 and boosts < 100:
        energy = [0.0] * 8
        for j, forbidden in enumerate(order):
            energy[forbidden] = w[j]
        value, bit = min((energy[state ^ (1 << bit)], bit) for bit in range(3))
        if value < energy[state] - 1e-12:
            state ^= 1 << bit
            path.append(state)
        else:
            w[order.index(state)] += 1
            boosts += 1
    return {'success': state == 7, 'path': path, 'boosts': boosts, 'final_weights': w}


def main():
    configs = [{'name': 'uniform'}]
    for name, value in CONSTANTS.items():
        for tiling in ANGLES:
            for amplitude in (.25, .75):
                configs.append({'name': f'{name}/{tiling}/{amplitude}', 'constant': value,
                                'tiling': tiling, 'amplitude': amplitude})
    stats, traps = [], []
    for config in configs:
        solved = moves = stuck = 0
        for mask in range(255):  # 255 SAT subsets, including the empty formula.
            forbidden = [s for s in range(8) if mask & (1 << s)]
            w = weights(config, len(forbidden))
            energy = [0.0] * 8
            for j, state in enumerate(forbidden):
                energy[state] = w[j]
            for start in range(8):
                path = descend(energy, start)
                success = energy[path[-1]] == 0
                solved += success
                stuck += not success
                moves += len(path) - 1
        stats.append({'config': config['name'], 'runs': 2040, 'solved': solved,
                      'stuck': stuck, 'success_rate': solved / 2040, 'flips': moves})
        traps.append(adversarial(config))
    example = next(t for t in traps if t['config'] == 'phi/hat/0.75')
    output = {'scope': 'All 256 subsets of eight distinct clauses using all three variables; 255 SAT, one UNSAT.',
              'configurations': configs, 'summary': stats, 'adversarial_traps': traps,
              'runs': len(configs) * 2040, 'formulas': 256,
              'example': example, 'example_breakout': breakout(example),
              'proof_scope': 'Any seven fixed positive clause-position weights, independent of clause content; arbitrary clause permutations; strict one-bit descent. Does not cover adaptive weights, sideways moves, or nonlocal moves.',
              'tests': 'Each positive weight template has an independently checked unique-SAT clause formula and a strict-greedy stall at 000. All recorded paths strictly descend with no repeated state.'}
    destination = ROOT / 'exhaustive_small.json'
    destination.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'configs': len(configs), 'runs': output['runs'], 'all_templates_have_trap': True,
                      'example_breakout': output['example_breakout']}, indent=2))


if __name__ == '__main__':
    main()
