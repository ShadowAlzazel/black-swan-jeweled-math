#!/usr/bin/env python3
"""Two selected SAT demonstrations of adaptive Hat-orientation clause weights.

Standard library only.  Clauses and the source orientation snapshot are embedded;
when the original Hat patch JSON is present its orientation data are verified.
"""
import json
import math
from pathlib import Path

N, MAX_UPDATES = 4, 100
PHI = (1 + math.sqrt(5)) / 2
ROTATIONS = [240, 240, 120, 120, 240, 240, 120, 120, 120, 120, 0, 0, 180,
             240, 180, 240, 180, 0, 300, 0, 300, 300, 240, 120, 60]
ADDRESSES = [[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3],
             [2, 0], [2, 1], [2, 2], [2, 3], [3, 0], [4, 0], [4, 1], [5, 0],
             [5, 1], [6, 0], [6, 1], [7, 0], [7, 1], [8, 0], [8, 1], [9, 0], [9, 1]]
CASES = {
    'single_plateau14': {'start': 14, 'clauses': [
        [-1, -3, 4], [-1, 3, -4], [-1, -3, -4], [1, -2, 4], [-1, 2, -3],
        [-1, -2, 3], [1, 2, -4], [1, -2, -4], [-1, -2, -4], [1, 2, 3],
        [-2, 3, 4], [1, -2, 3], [-1, -2, 4], [-1, -2, -3]]},
    'rich17': {'start': 0, 'clauses': [
        [-1, 2, -3], [-1, -2, -3], [1, -2, 4], [-1, 3, 4], [2, -3, -4],
        [2, 3, -4], [-1, -3, -4], [1, -2, 3], [1, 3, -4], [-1, 2, -4],
        [-1, 2, 3], [1, 2, -4], [-2, 3, -4], [1, -3, 4], [-1, 3, -4],
        [-1, -3, 4], [1, 2, 3]]}}


def violated(state, clauses):
    return [j for j, clause in enumerate(clauses)
            if not any(bool((state >> (abs(lit) - 1)) & 1) == (lit > 0) for lit in clause)]


def energy(state, clauses, weights):
    return sum(weights[j] for j in violated(state, clauses))


def run(case, method):
    clauses, state = case['clauses'], case['start']
    weights = [1.0] * len(clauses)
    updates = flips = 0
    events, path = [], [state]

    def record(kind, **extra):
        events.append(dict(state=state, stateLabel=''.join(str((state >> j) & 1) for j in range(N)),
            weights=weights.copy(), violatedClauses=violated(state, clauses), kind=kind,
            updateCount=updates, flipCount=flips, energy=energy(state, clauses, weights), **extra))

    record('start')
    while True:
        bad = violated(state, clauses)
        if not bad:
            reason = 'sat'
            break
        best, variable = min((energy(state ^ (1 << j), clauses, weights), j) for j in range(N))
        if best < energy(state, clauses, weights) - 1e-12:
            state ^= 1 << variable
            flips += 1
            path.append(state)
            record('flip', flippedVariable=variable)
        elif method == 'uniform' or updates >= MAX_UPDATES:
            reason = 'trap' if method == 'uniform' else 'update_limit'
            break
        else:
            phase = 2 * math.pi * updates / PHI ** 2
            delta = [0.0] * len(clauses)
            for j in bad:
                a = math.radians(ROTATIONS[j]) - phase
                boost = (1.0 if method == 'breakout' else 1 + 0.5 * math.cos(a)
                         if method == 'orientation_cos' else 0.25 + 2 * (1 + math.cos(a)) ** 2)
                delta[j] = boost
                weights[j] += boost
            updates += 1
            record('boost', updateIndex=updates - 1, updatePhase=phase, deltaWeights=delta)
    return dict(success=reason == 'sat', reason=reason, updates=updates, flips=flips,
                path=path, events=events, finalState=state)


def validate(case, result):
    clauses, events = case['clauses'], result['events']
    assert all((a ^ b).bit_count() == 1 for a, b in zip(result['path'], result['path'][1:]))
    assert result['success'] == (not violated(result['finalState'], clauses))
    for event in events:
        assert all(w > 0 for w in event['weights'])
        for state in range(1 << N):
            assert (energy(state, clauses, event['weights']) == 0) == (not violated(state, clauses))
    for previous, current in zip(events, events[1:]):
        if current['kind'] == 'flip':
            assert previous['weights'] == current['weights']
            assert (previous['state'] ^ current['state']).bit_count() == 1
            assert current['energy'] < previous['energy'] - 1e-12
        else:
            assert current['kind'] == 'boost' and previous['state'] == current['state']
            for j, (before, after) in enumerate(zip(previous['weights'], current['weights'])):
                assert (after > before) == (j in previous['violatedClauses'])


def main():
    source = Path(__file__).parent / 'tiling_research' / 'hat_patches.json'
    source_verified = False
    if source.exists():
        patch = next(p for p in json.loads(source.read_text())['patches'] if p['iteration'] == 1)
        assert [t['orientationDegrees'] for t in patch['tiles']] == ROTATIONS
        assert [t['address'] for t in patch['tiles']] == ADDRESSES
        assert [i for i, t in enumerate(patch['tiles']) if t['reflected']] == [3, 7, 11]
        source_verified = True
    output = {'configuration': {'n': N, 'maxUpdates': MAX_UPDATES, 'phi': PHI,
        'indexing': 'Clause/tile/variable indices are zero-based; signed literals are +/-1..4; integer bit j is variable x_(j+1).',
        'ties': 'Choose the smallest variable index among minimum-energy neighbors; identical across methods.',
        'phase': 'First trap uses updateIndex r=0; phase=2*pi*r/phi^2.'},
        'tiling': {'source': 'https://github.com/isohedral/hatviz', 'license': 'BSD-3-Clause',
            'patchIteration': 1, 'metatile': 'H', 'sourceFileVerified': source_verified,
            'role': 'Clause j uses rotation of tile j. Reflections are recorded but do not change the update rule. Geometry neither adds clauses nor certifies SAT.',
            'tiles': [dict(tileIndex=j, orientationDegrees=angle, reflected=j in (3, 7, 11), address=ADDRESSES[j])
                      for j, angle in enumerate(ROTATIONS)]},
        'methods': {'uniform': 'Strict descent; stop on a non-solution trap.',
            'breakout': 'At a trap add 1 to each currently violated clause weight.',
            'orientation_cos': 'At a trap add 1+0.5*cos(theta_j-2*pi*r/phi^2) to violated clause j only.',
            'orientation_341': 'At a trap add .25+g(theta_j-2*pi*r/phi^2) to violated clause j only, with g(a)=3+4*cos(a)+cos(2*a)=2*(1+cos(a))^2; compute the squared form.'},
        'interpretation': 'Two deliberately selected tiny examples, not a performance benchmark. Each method uses different boost magnitudes. Energy decreases only between weight updates; total energy need not decrease across updates. All positive weights preserve exactly the SAT zero set.',
        'cases': {}, 'tests': ['Exact SAT zero sets checked for all16 assignments at every recorded event.',
            'All flips change exactly one bit and strictly decrease the current fixed weighted energy.',
            'Every boost changes only violated clause weights, with strictly positive increments.',
            'All reported final success states satisfy every original clause.']}
    for name, case in CASES.items():
        case = dict(case, n=N, m=len(case['clauses']), orientationTileIndices=list(range(len(case['clauses']))))
        case['runs'] = {method: run(case, method) for method in output['methods']}
        for result in case['runs'].values():
            validate(case, result)
        output['cases'][name] = case
    for j in range(100):
        a = j / 7
        assert math.isclose(3 + 4 * math.cos(a) + math.cos(2 * a), 2 * (1 + math.cos(a)) ** 2, abs_tol=1e-12)
    target = Path(__file__).with_name('tiling_adaptive_results.json')
    target.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({name: {method: {k: run[k] for k in ('success', 'updates', 'flips', 'path')}
          for method, run in case['runs'].items()} for name, case in output['cases'].items()}, indent=2))


if __name__ == '__main__':
    main()
