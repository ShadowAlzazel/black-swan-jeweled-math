#!/usr/bin/env python3
"""Complete capped DPLL, creative search guidance, and exact XOR recognition.

Search guidance never deletes a clause or identifies distinct Boolean states.
An interrupted exact search reports UNKNOWN, never UNSAT.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import random
import statistics
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SEED = 20261127
NODE_CAP = 20000
FRONTEND_CAP = 256
CONFIGS = [
    dict(name='occurrence_dpll', schedule=None, frontend=None, xor=False),
    dict(name='random_order_dpll', schedule=None, frontend=None, xor=False),
    dict(name='walksat_then_dpll', schedule=None, frontend='walksat', xor=False),
    dict(name='xor_then_dpll', schedule=None, frontend=None, xor=True),
] + [dict(name=name, schedule=name, frontend='creative', xor=name in ('all_ideas', 'shuffled_all_ideas'))
     for name in ('golden_hat', 'silver_pinwheel', 'zeta_quasicrystal',
                  'feigenbaum_mobius', 'collatz_fourier', 'all_ideas', 'shuffled_all_ideas')]


def seed_for(*parts):
    raw = ':'.join(map(str, (SEED,) + parts)).encode()
    return int.from_bytes(hashlib.blake2b(raw, digest_size=8).digest(), 'big')


def satisfied(state, clauses):
    return all(any(bool((state >> (abs(lit) - 1)) & 1) == (lit > 0)
                   for lit in clause) for clause in clauses)


def to_three_cnf(n, clauses):
    """Existentially equivalent chain reduction; retain every auxiliary mapping."""
    output, reductions = [], []
    for source_index, raw in enumerate(clauses):
        clause = list(raw)
        if len(clause) <= 3:
            output.append(clause)
            continue
        auxiliary = list(range(n + 1, n + len(clause) - 2))
        n += len(auxiliary)
        chain = [[clause[0], clause[1], auxiliary[0]]]
        for j in range(len(auxiliary) - 1):
            chain.append([-auxiliary[j], clause[j + 2], auxiliary[j + 1]])
        chain.append([-auxiliary[-1], clause[-2], clause[-1]])
        output.extend(chain)
        reductions.append(dict(source_clause_index=source_index, raw_clause=clause,
                               aux_variables=auxiliary, output_clauses=chain))
    return n, output, reductions


def encode(clauses):
    encoded = []
    for j, clause in enumerate(clauses):
        positive = {lit for lit in clause if lit > 0}
        negative = {-lit for lit in clause if lit < 0}
        if positive & negative:  # A tautology imposes no assignment restriction.
            continue
        encoded.append((j, sum(1 << (v-1) for v in positive),
                        sum(1 << (v-1) for v in negative)))
    return encoded


class NodeLimit(Exception):
    pass


def dpll(n, clauses, weights=None, priors=None, preferred_bits=None,
         node_cap=NODE_CAP, random_order=False, seed=0):
    """Exact binary partition with unit propagation; the sole pruning is conflict."""
    weights = weights or [1.] * len(clauses)
    rng = random.Random(seed)
    if priors is not None:
        low, high = min(priors, default=0), max(priors, default=0)
        priors = [(p-low)/(high-low) if high > low else .5 for p in priors]
    if random_order:
        order = list(range(n)); rng.shuffle(order)
        priors = [0.] * n
        for rank, variable in enumerate(order):
            priors[variable] = n-rank
        preferred_bits = [rng.randrange(2) for _ in range(n)]
    stats = dict(dpll_nodes=0, propagations=0, clause_visits=0, conflicts=0)

    def search(active, ones, assigned):
        if stats['dpll_nodes'] >= node_cap:
            raise NodeLimit
        stats['dpll_nodes'] += 1
        while True:
            remaining, unit_ones, unit_zeros = [], 0, 0
            for j, pos, neg in active:
                stats['clause_visits'] += 1
                if (pos & ones) or (neg & assigned & ~ones):
                    continue
                free = (pos | neg) & ~assigned
                if not free:
                    stats['conflicts'] += 1
                    return None
                remaining.append((j, pos, neg))
                if free & (free - 1) == 0:
                    if pos & free:
                        unit_ones |= free
                    else:
                        unit_zeros |= free
            if unit_ones & unit_zeros:
                stats['conflicts'] += 1
                return None
            units = unit_ones | unit_zeros
            active = remaining
            if not active:
                return ones
            if not units:
                break
            stats['propagations'] += units.bit_count()
            assigned |= units
            ones = (ones | unit_ones) & ~unit_zeros
        positive, negative = [0.] * n, [0.] * n
        appearing = 0
        for j, pos, neg in active:
            free = (pos | neg) & ~assigned
            appearing |= free
            while free:
                bit = free & -free; variable = bit.bit_length()-1
                (positive if pos & bit else negative)[variable] += weights[j]
                free ^= bit
        variables = [j for j in range(n) if appearing >> j & 1]
        if random_order:
            variable = max(variables, key=lambda j: (priors[j], -j))
        else:
            variable = max(variables, key=lambda j: ((positive[j]+negative[j]) *
                (1 if priors is None else .75+.5*priors[j]), -j))
        preferred = (preferred_bits[variable] if preferred_bits is not None
                     else int(positive[variable] >= negative[variable]))
        bit = 1 << variable
        for value in (preferred, 1-preferred):
            answer = search(active, (ones | bit) if value else (ones & ~bit), assigned | bit)
            if answer is not None:
                return answer
        return None

    try:
        witness = search(encode(clauses), 0, 0)
        status = 'SAT' if witness is not None else 'UNSAT'
    except NodeLimit:
        witness, status = None, 'UNKNOWN'
    if status == 'SAT':
        assert satisfied(witness, clauses)
    return dict(status=status, witness=witness, **stats)


def violations(state, clauses):
    return [not any(bool((state >> (abs(lit)-1)) & 1) == (lit > 0)
                    for lit in clause) for clause in clauses]


def frontend(n, clauses, schedule, seed, mode, diffuse=False):
    rng = random.Random(seed)
    state = rng.randrange(1 << n)
    weights = list(schedule.get('clause_weights', [1.] * len(clauses)))
    angles = schedule.get('boost_angles', [0.] * len(clauses))
    constant = schedule.get('phase_constant', (1+math.sqrt(5))/2)
    neighbors = schedule.get('clause_neighbors')
    if diffuse and neighbors is None:
        sets = [set(map(abs, c)) for c in clauses]
        neighbors = [[k for k in range(len(clauses)) if k != j and sets[j] & sets[k]]
                     for j in range(len(clauses))]
    evaluations = flips = updates = 0

    def evaluate(s):
        nonlocal evaluations
        evaluations += 1
        return violations(s, clauses)

    bad = evaluate(state)
    while any(bad):
        indices = [j for j, violated in enumerate(bad) if violated]
        if any(not clauses[j] for j in indices):
            break  # An empty clause is left for the exact solver to certify.
        clause = clauses[rng.choice(indices)]
        variables = sorted(set(abs(lit)-1 for lit in clause))
        if evaluations + len(variables) > FRONTEND_CAP:
            break
        candidates = [state ^ (1 << j) for j in variables]
        candidate_bad = [evaluate(s) for s in candidates]
        if mode == 'walksat':
            breaks = [sum(a and not b for a,b in zip(c, bad)) for c in candidate_bad]
            minimum = min(breaks)
            chosen = (rng.choice([j for j,b in enumerate(breaks) if b == minimum])
                      if minimum == 0 or rng.random() >= .5 else rng.randrange(len(variables)))
        else:
            current_energy = sum(w for w,b in zip(weights,bad) if b)
            scores = [sum(w for w,b in zip(weights,c) if b) for c in candidate_bad]
            minimum = min(scores)
            if rng.random() < .2:
                chosen = rng.randrange(len(variables))
            elif minimum < current_energy - 1e-12:
                chosen = rng.choice([j for j,s in enumerate(scores) if abs(s-minimum) <= 1e-12])
            else:
                phase = 2*math.pi*updates/constant
                raw = {j:.25 + 2*(1+math.cos(angles[j]-phase))**2 for j in indices}
                average = sum(raw.values())/len(raw)
                if diffuse:
                    weights = [.85*w + .15*sum(weights[k] for k in [j]+neighbors[j]) /
                               (len(neighbors[j])+1) for j,w in enumerate(weights)]
                for j in indices:
                    weights[j] += raw[j]/average
                updates += 1
                assert all(w > 0 and math.isfinite(w) for w in weights)
                scores = [sum(w for w,b in zip(weights,c) if b) for c in candidate_bad]
                minimum = min(scores)
                chosen = rng.choice([j for j,s in enumerate(scores) if abs(s-minimum) <= 1e-12])
        state, bad = candidates[chosen], candidate_bad[chosen]
        flips += 1
    return dict(state=state, success=not any(bad), weights=weights,
                frontend_evaluations=evaluations, frontend_flips=flips,
                weight_updates=updates, frontend_clause_checks=evaluations*len(clauses))


def random_formula(n, m, rng, planted):
    witness = rng.randrange(1 << n) if planted else None
    clauses, seen = [], set()
    while len(clauses) < m:
        clause = tuple((v+1)*rng.choice((-1,1)) for v in sorted(rng.sample(range(n),3)))
        if clause in seen or (planted and not satisfied(witness,[clause])):
            continue
        seen.add(clause);clauses.append(list(clause))
    rng.shuffle(clauses)
    return clauses,witness


def connected(vertices, edges):
    adjacency=[[] for _ in range(vertices)]
    for a,b in edges:adjacency[a].append(b);adjacency[b].append(a)
    seen={0};stack=[0]
    while stack:
        for b in adjacency[stack.pop()]:
            if b not in seen:seen.add(b);stack.append(b)
    return len(seen)==vertices


def cubic_graph(k, rng, mix=False):
    edges={tuple(sorted((layer*k+j,layer*k+(j+1)%k))) for layer in range(2) for j in range(k)}
    edges|={(j,j+k) for j in range(k)}
    if mix:
        for _ in range(80*k):
            first,second=rng.sample(sorted(edges),2);a,b=first;c,d=second
            if len({a,b,c,d})<4:continue
            new={tuple(sorted((a,c))),tuple(sorted((b,d)))}
            if len(new)<2 or any(e in edges for e in new):continue
            candidate=(edges-{first,second})|new
            if connected(2*k,candidate):edges=candidate
    assert connected(2*k,edges)
    assert all(sum(v in edge for edge in edges)==3 for v in range(2*k))
    return sorted(edges)


def tseitin_case(k, rng, unsat, mix=False):
    edges=cubic_graph(k,rng,mix);n=len(edges);witness=rng.randrange(1<<n)
    equations=[];clauses=[]
    for vertex in range(2*k):
        variables=[j for j,e in enumerate(edges) if vertex in e]
        charge=sum(witness>>j&1 for j in variables)%2
        if unsat and vertex==0:charge^=1
        equations.append(dict(variables=variables,rhs=charge))
        for bits in itertools.product((0,1),repeat=3):
            if sum(bits)%2!=charge:
                clauses.append([-(j+1) if bit else j+1 for j,bit in zip(variables,bits)])
    rng.shuffle(clauses)
    return n,clauses,None if unsat else witness,dict(vertices=2*k,edges=edges,equations=equations,
            charge_sum_mod2=sum(e['rhs'] for e in equations)%2,connected=True)


def pigeonhole(pigeons,holes):
    var=lambda p,h:p*holes+h+1
    raw=[[var(p,h) for h in range(holes)] for p in range(pigeons)]
    raw += [[-var(p,a),-var(p,b)] for p in range(pigeons) for a,b in itertools.combinations(range(holes),2)]
    raw += [[-var(a,h),-var(b,h)] for h in range(holes) for a,b in itertools.combinations(range(pigeons),2)]
    n,clauses,reductions=to_three_cnf(pigeons*holes,raw)
    return n,clauses,dict(pigeons=pigeons,holes=holes,raw_n=pigeons*holes,raw_clauses=raw,reductions=reductions)


def generate_dataset():
    cases=[]
    def add(name,n,clauses,family,split,expected=None,witness=None,construction=None):
        identifier=len(cases)
        assert all(0<len(c)<=3 and len(set(map(abs,c)))==len(c) for c in clauses)
        assert all(1<=abs(lit)<=n for c in clauses for lit in c)
        if witness is not None:assert satisfied(witness,clauses)
        cases.append(dict(id=identifier,name=name,n=n,m=len(clauses),clauses=clauses,family=family,
            split=split,expected_status=expected,expected_witness=witness,construction=construction or {}))
    for n in (12,18,24,32):
        for planted in (True,False):
            for replicate in range(6):
                rng=random.Random(seed_for('formula',n,planted,replicate))
                clauses,witness=random_formula(n,round(4.2*n),rng,planted)
                status='SAT' if planted else None
                if not planted:
                    label=dpll(n,clauses,node_cap=200000)
                    status,witness=label['status'],label['witness']
                add(f'{"planted" if planted else "random"}_{n}_{replicate}',n,clauses,
                    'planted' if planted else 'random','train' if replicate%2==0 else 'test',status,witness,
                    {'replicate':replicate,'labels_from':'plant' if planted else 'complete occurrence DPLL'})
    for mix,ks in ((False,(3,4,5,6)),(True,(8,10,12))):
        for k in ks:
            for unsat in (False,True):
                n,c,w,metadata=tseitin_case(k,random.Random(seed_for('tseitin',mix,k)),unsat,mix)
                add(f'xor_{"cubic" if mix else "prism"}_{2*k}_{"UNSAT" if unsat else "SAT"}',n,c,
                    'xor_cubic' if mix else 'xor_prism','train' if k in (3,5,8) else 'test',
                    'UNSAT' if unsat else 'SAT',w,metadata)
    for holes in (3,4,5):
        for pigeons in (holes,holes+1):
            n,c,metadata=pigeonhole(pigeons,holes)
            add(f'php_{pigeons}_{holes}',n,c,'pigeonhole','train' if holes==3 else 'test',
                'SAT' if pigeons<=holes else 'UNSAT',None,metadata)
    for n in (12,18,24,32):
        half=n//2;rng=random.Random(seed_for('components',n))
        left,a=random_formula(half,round(4.2*half),rng,True)
        right,b=random_formula(n-half,round(4.2*(n-half)),rng,True)
        shifted=[[(abs(lit)+half)*(1 if lit>0 else -1) for lit in c] for c in right]
        add(f'components_{n}',n,left+shifted,'components','train' if n in (12,24) else 'test',
            'SAT',a|(b<<half),{'partition':[half,n-half]})
    assert len(cases)==72
    return dict(seed=SEED,node_cap=NODE_CAP,frontend_eval_cap=FRONTEND_CAP,trials=2,
                configs=CONFIGS,instances=cases)


def solve_instance(instance, config, trial):
    n,clauses=instance['n'],instance['clauses'];rng_seed=seed_for('search',instance['id'],trial)
    schedule=dict(variable_scores=[.5]*n,preferred_bits=None,clause_weights=[1.]*len(clauses),metadata={})
    if config['schedule']:
        from creative_schedules import make_schedule
        schedule=make_schedule(n,clauses,config['schedule'],rng_seed)
        assert len(schedule['variable_scores'])==n and len(schedule['preferred_bits'])==n
        assert len(schedule['clause_weights'])==len(clauses)
        assert all(math.isfinite(w) and w>0 for w in schedule['clause_weights'])
    result=dict(instance=instance['id'],config=config['name'],trial=trial,status='UNKNOWN',witness=None,
                stage=None,frontend_evaluations=0,frontend_flips=0,weight_updates=0,frontend_clause_checks=0,
                dpll_nodes=0,propagations=0,clause_visits=0,conflicts=0,kernel=None)
    if config['xor']:
        from xor_kernel import analyze_cnf
        kernel=analyze_cnf(n,clauses)
        result['kernel']={key:kernel.get(key) for key in ('status','pure','rank','nullity','contradiction_combination')}
        result['kernel']['equation_count']=len(kernel.get('equations',[]))
        if kernel['status'] in ('SAT','UNSAT'):
            result.update(status=kernel['status'],witness=kernel.get('assignment'),stage='xor_kernel')
            if result['status']=='SAT':assert satisfied(result['witness'],clauses)
            return result
    weights=schedule['clause_weights'];preferred=schedule.get('preferred_bits')
    if config['frontend']:
        front=frontend(n,clauses,schedule,rng_seed,config['frontend'],
                       diffuse=config['name'] in ('all_ideas','shuffled_all_ideas'))
        for key in ('frontend_evaluations','frontend_flips','weight_updates','frontend_clause_checks'):
            result[key]=front[key]
        if front['success']:
            result.update(status='SAT',witness=front['state'],stage='frontend')
            assert satisfied(result['witness'],clauses)
            return result
        weights=front['weights'];preferred=[front['state']>>j&1 for j in range(n)]
    exact=dpll(n,clauses,weights=weights,priors=None if config['schedule'] is None else schedule['variable_scores'],
               preferred_bits=preferred,node_cap=NODE_CAP,random_order=config['name']=='random_order_dpll',seed=rng_seed)
    result.update(exact,stage='dpll')
    assert result['frontend_evaluations']<=FRONTEND_CAP and result['dpll_nodes']<=NODE_CAP
    if instance['expected_status'] in ('SAT','UNSAT') and result['status']!='UNKNOWN':
        assert result['status']==instance['expected_status'],(instance['name'],result)
    return result


def summarize(dataset,runs):
    instances={i['id']:i for i in dataset['instances']};answer={}
    for config in CONFIGS:
        selected=[r for r in runs if r['config']==config['name']];groups={}
        for name,records in [('all',selected)]+[(split,[r for r in selected if instances[r['instance']]['split']==split]) for split in ('train','test')]+[(f,[r for r in selected if instances[r['instance']]['family']==f]) for f in sorted({i['family'] for i in dataset['instances']})]:
            groups[name]=dict(runs=len(records),statuses=dict(Counter(r['status'] for r in records)),
                stages=dict(Counter(r['stage'] for r in records)),nodes=sum(r['dpll_nodes'] for r in records),
                median_nodes=statistics.median(r['dpll_nodes'] for r in records) if records else None,
                frontend_evaluations=sum(r['frontend_evaluations'] for r in records),
                propagations=sum(r['propagations'] for r in records),
                clause_visits=sum(r['clause_visits'] for r in records))
        answer[config['name']]=groups
    return answer


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate-only',action='store_true');parser.add_argument('--pilot',action='store_true')
    args=parser.parse_args();start=time.perf_counter();dataset=generate_dataset()
    (ROOT/'dataset.json').write_text(json.dumps(dataset,indent=2)+'\n')
    print('Dataset ready:',len(dataset['instances']),'cases',dict(Counter(i['expected_status'] for i in dataset['instances'])),flush=True)
    if args.generate_only:return
    runs=[]
    for config in CONFIGS:
        subset=dataset['instances'] if not args.pilot else dataset['instances'][:2]+dataset['instances'][54:56]+dataset['instances'][66:68]
        for instance in subset:
            for trial in range(2):runs.append(solve_instance(instance,config,trial))
        print('Completed',config['name'],len(runs),'runs',round(time.perf_counter()-start,2),'seconds',flush=True)
    output=dict(protocol='Complete unit-propagation DPLL with bounded guidance. Caps yield UNKNOWN. Geometric maps and weights guide branching, never merge states or delete constraints.',
        prespecified_candidate='all_ideas',configuration_count=len(CONFIGS),instance_count=len(dataset['instances']),
        run_count=len(runs),trials=2,node_cap=NODE_CAP,frontend_eval_cap=FRONTEND_CAP,
        dataset_counts=dict(Counter(i['expected_status'] for i in dataset['instances'])),
        summary=summarize(dataset,runs),runs=runs,elapsed_seconds=time.perf_counter()-start,
        interpretation=['No inference of a polynomial worst-case bound from finite tests.',
            'Kernel fast paths exploit recognized parity structure; general formulas retain exact DPLL fallback.',
            'Frontend candidate evaluations, exact nodes, clause visits, and kernel ranks are separate costs.',
            'Randomness changes bounded guidance and branch order, not the logical correctness of complete search.'])
    destination=ROOT/('pilot_results.json' if args.pilot else 'results.json')
    destination.write_text(json.dumps(output,indent=2)+'\n')
    failures=[r for r in runs if r['status']=='UNKNOWN']
    (ROOT/'unresolved_cases.json').write_text(json.dumps(failures,indent=2)+'\n')
    print('Wrote',destination,'unresolved runs',len(failures),flush=True)


if __name__=='__main__':main()
