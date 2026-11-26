#!/usr/bin/env python3
"""Reproducible broad SAT landscape sweep. No network or external packages."""
import collections,csv,gzip,hashlib,json,math,random,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SEED=20261125
BUDGET=2049
CONSTANTS={'phi':(1+math.sqrt(5))/2,'phi_squared':((1+math.sqrt(5))/2)**2,'phi_reciprocal':2/(1+math.sqrt(5)),
 'sqrt2':math.sqrt(2),'sqrt3':math.sqrt(3),'pi':math.pi,'e':math.e,'zeta2':math.pi**2/6,
 'plastic':1.324717957244746,'rational_13_8':13/8,'random_frequency':random.Random(SEED).uniform(1,4)}
METHODS={'greedy':0,'restarts':1,'breakout':2,'walksat':3,'static':4,'cos':5,'341':6,'diffusion':7}
def seed(*parts):return int.from_bytes(hashlib.blake2b((':'.join(map(str,(SEED,)+parts))).encode(),digest_size=8).digest(),'big')
def satisfied(state,clauses):return all(any(bool(state>>(abs(l)-1)&1)==(l>0) for l in c) for c in clauses)
def exact_solve(n,clauses,cap=2000000):
 masks=[(sum(1<<(l-1) for l in c if l>0),sum(1<<(-l-1) for l in c if l<0)) for c in clauses]
 nodes=0
 def visit(known,value):
  nonlocal nodes
  nodes+=1
  if nodes>cap:raise TimeoutError
  while True:
   unresolved=[];unit=0;sign=0
   for pos,neg in masks:
    if value&pos or (known&~value)&neg:continue
    left=(pos|neg)&~known
    if not left:return False,None
    if left&(left-1)==0:unit=left;sign=bool(left&pos);break
    unresolved.append((pos,neg,left))
   if not unit:break
   known|=unit;value=(value|unit) if sign else (value&~unit)
  if not unresolved:return True,value
  counts=[0.0]*n;polar=[0.0]*n
  for pos,neg,left in unresolved:
   importance=2.0**(-left.bit_count())
   while left:
    bit=left&-left;j=bit.bit_length()-1;left-=bit;counts[j]+=importance;polar[j]+=importance if bit&pos else -importance
  j=max(range(n),key=lambda j:counts[j]);bit=1<<j
  first=(value|bit) if polar[j]>=0 else (value&~bit)
  result=visit(known|bit,first)
  return result if result[0] else visit(known|bit,first^bit)
 try:
  sat,w=visit(0,0)
  if sat:assert satisfied(w,clauses)
  return ('SAT' if sat else 'UNSAT'),w,nodes
 except TimeoutError:return 'UNKNOWN',None,nodes
def make_formula(n,m,family,rng):
 planted=rng.randrange(1<<n) if family!='random' else None
 side=math.ceil(math.sqrt(n));coords=[(j%side,j//side) for j in range(n)]
 nearest=[sorted((k for k in range(n) if k!=j),key=lambda k:((coords[k][0]-coords[j][0])**2+(coords[k][1]-coords[j][1])**2,k))[:min(6,n-1)] for j in range(n)]
 clauses=[];seen=set()
 while len(clauses)<m:
  if family=='geometric':
   j=rng.randrange(n);variables=sorted([j]+rng.sample(nearest[j],2))
  else:variables=sorted(rng.sample(range(n),3))
  c=tuple((j+1)*rng.choice((-1,1)) for j in variables)
  if c in seen or (planted is not None and not any(bool(planted>>(abs(l)-1)&1)==(l>0) for l in c)):continue
  clauses.append(c);seen.add(c)
 rng.shuffle(clauses)
 return clauses,planted
def shared_edge(a,b,tol=1e-9):
 for u,v in zip(a,a[1:]+a[:1]):
  dx,dy=v[0]-u[0],v[1]-u[1];length=math.hypot(dx,dy)
  for x,y in zip(b,b[1:]+b[:1]):
   if abs(dx*(x[1]-u[1])-dy*(x[0]-u[0]))>tol*length or abs(dx*(y[1]-u[1])-dy*(y[0]-u[0]))>tol*length:continue
   t0=((x[0]-u[0])*dx+(x[1]-u[1])*dy)/length
   t1=((y[0]-u[0])*dx+(y[1]-u[1])*dy)/length
   if min(length,max(t0,t1))-max(0,min(t0,t1))>tol:return True
 return False
def tilings():
 cached=ROOT/'tiling_templates.json'
 if cached.exists():return json.loads(cached.read_text())
 source=ROOT/'tiling_research'
 if not source.exists():source=ROOT.parent/'tiling_research'
 h=json.loads((source/'hat_patches.json').read_text())['patches'][2]
 p=json.loads((source/'pinwheel_patches.json').read_text())['patches'][4]
 result={}
 for name,patch in [('hat',h),('pinwheel',p)]:
  count=len(patch['tiles']);neighbors=[[] for _ in range(count)]
  if name=='hat':edges=patch['adjacency']
  else:
   polygons=[t['polygon'] for t in patch['tiles']]
   boxes=[(min(q[0] for q in poly),max(q[0] for q in poly),min(q[1] for q in poly),max(q[1] for q in poly)) for poly in polygons]
   edges=[]
   for i in range(count):
    a=boxes[i]
    for j in range(i):
     b=boxes[j]
     if a[0]>b[1]+1e-9 or b[0]>a[1]+1e-9 or a[2]>b[3]+1e-9 or b[2]>a[3]+1e-9:continue
     if shared_edge(polygons[i],polygons[j]):edges.append([i,j])
  for i,j in edges:neighbors[i].append(j);neighbors[j].append(i)
  result[name]={'angles':[math.radians(t['orientationDegrees']) for t in patch['tiles']],
    'neighbors':neighbors,'centroids':[[sum(q[k] for q in t['polygon'])/len(t['polygon']) for k in (0,1)] for t in patch['tiles']],
    'polygons':[t['polygon'] for t in patch['tiles']],
    'source':'actual substitution patch; adjacency by positive-length shared polygon boundary'}
 side=13;count=side*side
 result['square']={'angles':[(j%4)*math.pi/2 for j in range(count)],
   'neighbors':[[k for k in (j-1,j+1,j-side,j+side) if 0<=k<count and (abs(k-j)==side or k//side==j//side)] for j in range(count)],
   'centroids':[[j%side,j//side] for j in range(count)],'source':'periodic square-grid control'}
 cached.write_text(json.dumps(result,indent=2)+'\n')
 return result
def generate():
 ts=tilings();configs=[]
 def add(name,method,tiling='square',constant=1,amplitude=0,cname='none'):
  configs.append(dict(id=len(configs),name=name,method=method,tiling=tiling,constant=constant,amplitude=amplitude,cname=cname))
 for m in ('greedy','restarts','breakout','walksat'):add(m,m)
 for cname,c in CONSTANTS.items():
  for tiling in ts:
   for method in ('static','cos'):
    for amplitude in (.25,.75):add(f'{method}:{cname}:{tiling}:a{amplitude}',method,tiling,c,amplitude,cname)
   for method in ('341','diffusion'):add(f'{method}:{cname}:{tiling}',method,tiling,c,0,cname)
 instances=[]
 for n in (12,18,24,32):
  for family in ('planted','random','geometric'):
   for rep in range(8):
    iid=len(instances);clauses,plant=make_formula(n,round(4.2*n),family,random.Random(seed('formula',n,family,rep)))
    label,witness,nodes=exact_solve(n,clauses)
    assert plant is None or label=='SAT'
    instances.append(dict(id=iid,n=n,m=len(clauses),clauses=clauses,family=family,replicate=rep,
      split='train' if rep%2==0 else 'test',sat_label=label,witness=witness,exact_nodes=nodes,
      planted_generation_only=plant,starts=[random.Random(seed('start',iid,j)).randrange(1<<n) for j in range(3)]))
    print('generated',iid,n,family,label,nodes,flush=True)
 dataset=dict(seed=SEED,budget_assignment_evaluations=BUDGET,max_flips=1000,max_weight_updates=100,
   constants=CONSTANTS,configurations=configs,instances=instances,tilings=ts,
   protocol='Replicate parity fixes training/test split before any heuristic is run. Every full assignment evaluation scans all clauses; adaptive updates recompute scores from cached exact current violations. All positive adaptive boosts normalize to mean1 over violated clauses. Diffusion before a boost uses .8I+.2P where P averages self and retained edge-sharing neighbors. Failed searches do not certify UNSAT. At the same cap methods consume different candidate counts per decision; no wall-clock speed claims.')
 (ROOT/'dataset.json').write_text(json.dumps(dataset,indent=2)+'\n')
 return dataset
def serialize(dataset,path):
 ts=dataset['tilings'];tnames=list(ts)
 with path.open('w') as f:
  f.write(f"{len(dataset['instances'])} {len(dataset['configurations'])} {len(ts)} {BUDGET}\n")
  for t in ts.values():
   f.write(str(len(t['angles']))+'\n'+' '.join(map(repr,t['angles']))+'\n')
   for row in t['neighbors']:f.write(str(len(row))+' '+' '.join(map(str,row))+'\n')
  for c in dataset['configurations']:f.write(f"{c['id']} {METHODS[c['method']]} {tnames.index(c['tiling'])} {c['constant']!r} {c['amplitude']!r}\n")
  for i in dataset['instances']:
   f.write(f"{i['id']} {i['n']} {i['m']}\n")
   for clause in i['clauses']:
    pos=sum(1<<(l-1) for l in clause if l>0);neg=sum(1<<(-l-1) for l in clause if l<0)
    f.write(f"{pos} {neg} "+' '.join(str(abs(l)-1) for l in clause)+'\n')
   f.write(' '.join(map(str,i['starts']))+'\n')
def summarize(dataset,runs):
 instances={i['id']:i for i in dataset['instances']};rows=list(csv.DictReader(runs.open()))
 for r in rows:
  for key in r:r[key]=int(r[key]) if key!='reason' else r[key]
  i=instances[r['instance']]
  assert r['success']==satisfied(r['final_state'],i['clauses'])
  bad=sum(not any(bool(r['final_state']>>(abs(l)-1)&1)==(l>0) for l in c) for c in i['clauses'])
  assert r['final_violations']==bad and r['evaluations']<=BUDGET
  assert not (i['sat_label']=='UNSAT' and r['success'])
 def stats(records):
  known=[r for r in records if instances[r['instance']]['sat_label']=='SAT']
  solved=[r for r in known if r['success']]
  return dict(runs=len(records),known_sat_runs=len(known),sat_successes=len(solved),sat_success_rate=len(solved)/len(known) if known else None,
   all_successes=sum(r['success'] for r in records),mean_assignment_evaluations=sum(r['evaluations'] for r in records)/len(records) if records else 0,
   mean_flips=sum(r['flips'] for r in records)/len(records) if records else 0,
   mean_updates=sum(r['updates'] for r in records)/len(records) if records else 0,
   reasons=dict(collections.Counter(r['reason'] for r in records)))
 summaries=[]
 for c in dataset['configurations']:
  cr=[r for r in rows if r['config']==c['id']]
  entry=dict(c)
  for split in ('train','test'):entry[split]=stats([r for r in cr if instances[r['instance']]['split']==split])
  entry['strata']={f'{split}:{family}:n{n}':stats([r for r in cr if instances[r['instance']]['split']==split and instances[r['instance']]['family']==family and instances[r['instance']]['n']==n]) for split in ('train','test') for family in ('planted','random','geometric') for n in (12,18,24,32)}
  summaries.append(entry)
 # Freeze selection by train SAT success, then training consumed evaluations, then fixed config id.
 noncontrols=[s for s in summaries if s['id']>=4]
 winner=min(noncontrols,key=lambda s:(-s['train']['sat_success_rate'],s['train']['mean_assignment_evaluations'],s['id']))
 summary=dict(protocol=dataset['protocol'],configuration_count=len(summaries),run_count=len(rows),instance_count=len(instances),
  exact_label_counts=dict(collections.Counter(i['sat_label'] for i in instances.values())),
  split_label_counts={split:dict(collections.Counter(i['sat_label'] for i in instances.values() if i['split']==split)) for split in ('train','test')},
  selection_rule='Highest training SAT success fraction; tie by lower training mean assignment evaluations, then smaller config id. Test outcomes are not used in selection.',
  frozen_training_winner=winner,configurations=summaries,
  validation={'independent_boolean_checks':len(rows),'all_eval_caps_pass':True,'no_unsat_successes':True})
 (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 with gzip.open(ROOT/'raw_runs.csv.gz','wt') as f:f.write(runs.read_text())
 print(json.dumps({'runs':len(rows),'labels':summary['exact_label_counts'],'winner':winner['name'],'winner_train':winner['train'],'winner_test':winner['test'],'controls':[{k:s[k] for k in ('name','train','test')} for s in summaries[:4]]},indent=2))
 return summary
def main():
 ROOT.mkdir(exist_ok=True)
 dataset=generate();serialize(dataset,ROOT/'engine_input.txt')
 print('Running',len(dataset['configurations'])*len(dataset['instances'])*3,'searches',flush=True)
 subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'sweep_engine.cpp'),'-o',str(ROOT/'sweep_engine')],check=True)
 subprocess.run([str(ROOT/'sweep_engine'),str(ROOT/'engine_input.txt'),str(ROOT/'raw_runs.csv')],check=True)
 summarize(dataset,ROOT/'raw_runs.csv')
if __name__=='__main__':main()
