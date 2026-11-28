"""Seeded geometric/spectral guidance. Guidance never identifies Boolean states.

All positions are finite features, not a replacement for the SAT constraints.
Only positive weights are supplied; branch scores are suggestions for exact DPLL.
"""
from __future__ import annotations

import cmath
import hashlib
import json
import math
import random
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PHI = (1 + math.sqrt(5)) / 2
SILVER = 1 + math.sqrt(2)
DELTA = 4.669201609102991
ALPHA = 2.502907875095893
C512 = -1.4011532908499239
GAMMA = (14.134725141734695,21.022039638771556,25.01085758014569,
         30.424876125859512,32.93506158773919,37.586178158825675,
         40.9187190121475,43.327073280915,48.00515088116716,49.7738324776723)
NAMES = ('golden_hat','silver_pinwheel','zeta_quasicrystal',
         'feigenbaum_mobius','collatz_fourier','all_ideas','shuffled_all_ideas')


def keyed_seed(seed, name):
    return int.from_bytes(hashlib.blake2b(f'{seed}:{name}'.encode(),digest_size=8).digest(),'big')


@lru_cache(maxsize=1)
def templates():
    h = json.loads((ROOT/'experiments/02_tiling_jewels/tiling_research/hat_patches.json').read_text())['patches'][2]
    p = json.loads((ROOT/'experiments/02_tiling_jewels/tiling_research/pinwheel_patches.json').read_text())['patches'][4]
    q = json.loads((ROOT/'experiments/01_initial_jewel/jewel_geometry.json').read_text())['patch']
    result = {}
    cached = json.loads((ROOT/'experiments/03_golden_sweep/tiling_templates.json').read_text())
    for name, patch in [('hat',h),('pinwheel',p)]:
        coords = [sum((complex(*v) for v in t['polygon']),0j)/len(t['polygon']) for t in patch['tiles']]
        edges = patch.get('adjacency',[])
        if not edges:
            edges = [(j,k) for j,row in enumerate(cached[name]['neighbors']) for k in row if j<k]
        result[name] = dict(coords=coords,angles=[math.radians(t['orientationDegrees']) for t in patch['tiles']],edges=edges)
    result['quasicrystal'] = dict(coords=[complex(*t['physical']) for t in q['points']],angles=[cmath.phase(complex(*t['physical'])) for t in q['points']],edges=q['edges'])
    return result


def normalize(values):
    if not values:
        return []
    low, high = min(values), max(values)
    return [(v-low)/(high-low) if high>low+1e-14 else .5 for v in values]


def unit_coords(values):
    center = sum(values,0j)/len(values)
    scale = max(abs(z-center) for z in values) or 1
    return [.8*(z-center)/scale for z in values]


def graph_modes(n, clauses, seed):
    """Two approximate low Laplacian modes by deflated symmetric diffusion."""
    adjacency = [dict() for _ in range(n)]
    for clause in clauses:
        variables = sorted({abs(lit)-1 for lit in clause})
        for i in variables:
            for j in variables:
                if i!=j:
                    adjacency[i][j] = adjacency[i].get(j,0)+1
    degrees = [sum(row.values()) for row in adjacency]
    divisor = max(degrees,default=0)+1
    rng = random.Random(keyed_seed(seed,'laplacian'))
    modes = []
    for _ in range(2):
        vector = [rng.uniform(-1,1) for _ in range(n)]
        for iteration in range(72):
            mean = sum(vector)/n
            vector = [v-mean for v in vector]
            for previous in modes:
                dot = sum(a*b for a,b in zip(vector,previous))
                vector = [a-dot*b for a,b in zip(vector,previous)]
            norm = math.sqrt(sum(v*v for v in vector))
            if norm<1e-14:
                vector = [rng.uniform(-1,1) for _ in range(n)]
                continue
            vector = [v/norm for v in vector]
            if iteration<71:
                vector = [vector[i]-(degrees[i]*vector[i]-sum(weight*vector[j] for j,weight in adjacency[i].items()))/divisor for i in range(n)]
        modes.append(vector)
    return modes


def primes(count):
    answer=[];candidate=2
    while len(answer)<count:
        if all(candidate%p for p in answer if p*p<=candidate):
            answer.append(candidate)
        candidate+=1
    return answer


def geometry(n, family, mapping):
    template=templates()[family]
    coords=unit_coords([template['coords'][j%len(template['coords'])] for j in range(n)])
    if mapping=='exp':
        return [cmath.exp(z) for z in coords]
    if mapping=='square':
        return [z*z for z in coords]
    a=.3+.1j
    return [(z/ALPHA-a)/(1-a.conjugate()*z/ALPHA) for z in coords]


def clause_graph(clauses, family, union_logical=False):
    m=len(clauses);template=templates()[family];size=len(template['coords'])
    graph=[set() for _ in clauses]
    for block in range(0,m,size):
        for u,v in template['edges']:
            if block+u<m and block+v<m:
                graph[block+u].add(block+v);graph[block+v].add(block+u)
    if union_logical:
        occurrences={}
        for j,clause in enumerate(clauses):
            for v in {abs(lit) for lit in clause}:
                occurrences.setdefault(v,[]).append(j)
        for group in occurrences.values():
            for j in group:
                graph[j].update(k for k in group if k!=j)
    return [sorted(row) for row in graph]


def make_schedule(n, clauses, name, seed):
    if name not in NAMES:
        raise ValueError(name)
    if n == 0 or any(len(clause) == 0 for clause in clauses):
        # Integration fix (2026-11-28): the geometric features divide by n and by clause
        # widths. Zero-variable inputs and empty clauses get a neutral schedule; the exact
        # DPLL component already decides them (an empty clause is UNSAT).
        return dict(variable_scores=[.5]*n, preferred_bits=[0]*n, clause_weights=[1.0]*len(clauses),
                    boost_angles=[0.0]*len(clauses), phase_constant=1.0,
                    clause_neighbors=[[] for _ in clauses], spectral_scores=[.5]*n,
                    metadata=dict(name=name, seed=seed, degenerate_input=True))
    if name=='shuffled_all_ideas':
        base=make_schedule(n,clauses,'all_ideas',seed)
        rng=random.Random(keyed_seed(seed,'matched-shuffle'))
        vp=list(range(n));cp=list(range(len(clauses)));rng.shuffle(vp);rng.shuffle(cp)
        inverse={old:new for new,old in enumerate(cp)}
        return dict(variable_scores=[base['variable_scores'][j] for j in vp],preferred_bits=[base['preferred_bits'][j] for j in vp],clause_weights=[base['clause_weights'][j] for j in cp],boost_angles=[base['boost_angles'][j] for j in cp],clause_neighbors=[[inverse[k] for k in base['clause_neighbors'][j]] for j in cp],phase_constant=base['phase_constant'],metadata=dict(base['metadata'],name=name,control='Same feature/weight multiset and isomorphic diffusion graph; attachment independently shuffled.'))
    rng=random.Random(keyed_seed(seed,'shared-phases'))
    phase=rng.random()*2*math.pi
    modes=graph_modes(n,clauses,seed)
    spectral=normalize(modes[0])
    p=primes(max(n,len(clauses),1))
    zeta=[sum(math.cos(g*math.log(p[j])+phase) for g in GAMMA)/len(GAMMA) for j in range(n)]
    orbit=[];value=0.0
    for _ in range(max(n,64)):
        value=value*value+C512;orbit.append(value)
    integer=1+keyed_seed(seed,'collatz')%100000
    parity=[]
    for _ in range(64):
        parity.append(integer&1)
        integer=(3*integer+1)//2 if integer&1 else integer//2
    fourier=[sum((2*bit-1)*cmath.exp(-2j*math.pi*k*j/len(parity)) for j,bit in enumerate(parity))/len(parity) for k in range(1,9)]
    hat=geometry(n,'hat','exp');pin=geometry(n,'pinwheel','square');quasi=geometry(n,'quasicrystal','exp')
    a=.3+.1j
    spectral_z=[complex(modes[0][j],modes[1][j])*math.sqrt(n)/ALPHA for j in range(n)]
    # Restrict the Möbius input to the unit disk before evaluating the denominator.
    spectral_z=unit_coords(spectral_z)
    mobius=[(z-a)/(1-a.conjugate()*z) for z in spectral_z]
    feature={
        'golden_hat':[math.cos(cmath.phase(z)+2*math.pi*j/PHI+phase) for j,z in enumerate(hat)],
        'silver_pinwheel':[math.cos(cmath.phase(z)+2*math.pi*j/SILVER+phase) for j,z in enumerate(pin)],
        'zeta_quasicrystal':[.6*zeta[j]+.4*math.cos(cmath.phase(quasi[j])+j/(math.pi**2/6)) for j in range(n)],
        'feigenbaum_mobius':[math.cos(cmath.phase(mobius[j])+2*math.pi*orbit[j]/ALPHA+j/DELTA) for j in range(n)],
        'collatz_fourier':[math.cos(cmath.phase(fourier[j%8])+2*math.pi*j/n)+(2*parity[j%64]-1)*.25 for j in range(n)]}
    if name=='all_ideas':
        score=[sum(values[j] for values in feature.values())/len(feature)+.35*(2*spectral[j]-1)+.025*rng.uniform(-1,1) for j in range(n)]
        positions=[(hat[j]+pin[j]+quasi[j]+mobius[j])/4 for j in range(n)]
        family='hat';constant=DELTA
    else:
        score=feature[name]
        family={'golden_hat':'hat','silver_pinwheel':'pinwheel','zeta_quasicrystal':'quasicrystal','feigenbaum_mobius':'hat','collatz_fourier':'pinwheel'}[name]
        constant={'golden_hat':PHI,'silver_pinwheel':SILVER,'zeta_quasicrystal':math.pi**2/6,'feigenbaum_mobius':DELTA,'collatz_fourier':math.e}[name]
        positions={'golden_hat':hat,'silver_pinwheel':pin,'zeta_quasicrystal':quasi,'feigenbaum_mobius':mobius,'collatz_fourier':pin}[name]
    normalized=normalize(score)
    weights=[];angles=[]
    for j,clause in enumerate(clauses):
        variables=[abs(lit)-1 for lit in clause]
        center=sum((positions[v] for v in variables),0j)/len(variables)
        angle=cmath.phase(center)+2*math.pi*(j+1)/constant+phase
        angles.append(angle)
        weights.append(1+.4*math.tanh(sum(score[v] for v in variables)/len(variables)))
    assert all(math.isfinite(w) and .59<w<1.41 for w in weights)
    return dict(variable_scores=normalized,preferred_bits=[int((score[j]>=0)^bool(parity[j%64])) for j in range(n)],clause_weights=weights,boost_angles=angles,phase_constant=constant,clause_neighbors=clause_graph(clauses,family,name=='all_ideas'),spectral_scores=spectral,metadata=dict(name=name,seed=seed,geometry_family=family,geometry_used_as='Features/ordering only; variable identities are never merged.',known_zero_count=len(GAMMA),collatz_steps=64,spectral_diffusion_iterations=72,constants=dict(phi=PHI,silver=SILVER,delta=DELTA,alpha=ALPHA,zeta2=math.pi**2/6),rh_role='Global interpretation only; the finite known ordinates do not require RH.'))
