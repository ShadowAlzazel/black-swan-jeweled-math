#!/usr/bin/env python3
"""Locked training winner under clause order, variable relabeling, and gauge flips."""
import copy,csv,gzip,json,random,subprocess
from pathlib import Path
from run_sweep import ROOT,seed,serialize,satisfied,BUDGET

def map_state(state,permutation,gauge):
 return sum((((state>>j)&1)^gauge[j])<<permutation[j] for j in range(len(permutation)))

def main():
 data=json.loads((ROOT/'dataset.json').read_text());main_summary=json.loads((ROOT/'summary.json').read_text())
 winner=main_summary['frozen_training_winner']['id']
 config_ids=[winner,2,3]
 original=[i for i in data['instances'] if i['split']=='test' and i['sat_label']=='SAT']
 configs=[c for c in data['configurations'] if c['id'] in config_ids]
 transformations=[];all_rows=[];records=[]
 if not (ROOT/'sweep_engine').exists():subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'sweep_engine.cpp'),'-o',str(ROOT/'sweep_engine')],check=True)
 for variant in range(4):
  transformed=[]
  for ins in original:
   item=copy.deepcopy(ins);rng=random.Random(seed('permutation_stress',variant,ins['id']))
   permutation=list(range(ins['n']));gauge=[0]*ins['n'];order=list(range(ins['m']))
   if variant:
    rng.shuffle(permutation);gauge=[rng.randrange(2) for _ in range(ins['n'])];rng.shuffle(order)
   clauses=[]
   for old_index in order:
    new=[]
    for literal in ins['clauses'][old_index]:
     j=abs(literal)-1;sign=1 if literal>0 else -1
     new.append((permutation[j]+1)*sign*(-1 if gauge[j] else 1))
    clauses.append(new)
   item['clauses']=clauses
   item['witness']=map_state(ins['witness'],permutation,gauge)
   item['starts']=[map_state(s,permutation,gauge) for s in ins['starts']]
   assert satisfied(item['witness'],clauses)
   for j in range(16):
    state=random.Random(seed('stress_equivalence',ins['id'],j)).randrange(1<<ins['n'])
    assert satisfied(state,ins['clauses'])==satisfied(map_state(state,permutation,gauge),clauses)
   transformed.append(item)
   transformations.append(dict(variant=variant,instance=ins['id'],variable_permutation=permutation,gauge_flips=gauge,clause_order=order))
  task=dict(data,instances=transformed,configurations=configs)
  input_path=ROOT/f'stress_input_{variant}.txt';run_path=ROOT/f'stress_runs_{variant}.csv'
  serialize(task,input_path)
  subprocess.run([str(ROOT/'sweep_engine'),str(input_path),str(run_path)],check=True)
  by_id={i['id']:i for i in transformed}
  rows=list(csv.DictReader(run_path.open()))
  for r in rows:
   for key in r:r[key]=int(r[key]) if key!='reason' else r[key]
   assert r['success']==satisfied(r['final_state'],by_id[r['instance']]['clauses'])
   assert r['evaluations']<=BUDGET
   r['variant']=variant
  all_rows.extend(rows)
  for c in configs:
   cr=[r for r in rows if r['config']==c['id']]
   records.append(dict(variant=variant,config=c['id'],name=c['name'],runs=len(cr),successes=sum(r['success'] for r in cr),
    success_rate=sum(r['success'] for r in cr)/len(cr),mean_evaluations=sum(r['evaluations'] for r in cr)/len(cr)))
  if variant==0:
   # Identity must reproduce the original engine outcomes exactly, including RNG.
   original_runs={(int(r['instance']),int(r['config']),int(r['start_index'])):r for r in csv.DictReader((ROOT/'raw_runs.csv').open())}
   for r in rows:
    old=original_runs[(r['instance'],r['config'],r['start_index'])]
    for key in old:assert str(r[key])==old[key],(key,r,old)
  input_path.unlink();run_path.unlink()
 summary={'frozen_training_winner':winner,'instance_count':len(original),'variants':4,'permuted_variants':3,'starts_per_instance':3,
  'run_count':len(all_rows),'selection':'Config chosen by original training only; no retuning after stress.',
  'protocol':'Identity plus three fresh seeded clause permutations, variable bijections, and Boolean gauge flips. Three original starts and a SAT witness transform exactly. Original instance IDs preserve RNG and phase seeds across variants. Identity runs reproduce raw main outcomes exactly.',
  'records':records,
  'aggregate_permuted':[{**{'config':c['id'],'name':c['name']},**{'runs':sum(r['config']==c['id'] and r['variant']>0 for r in all_rows),
   'successes':sum(r['success'] for r in all_rows if r['config']==c['id'] and r['variant']>0)}} for c in configs],
  'transformations':transformations}
 (ROOT/'permutation_stress_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 with gzip.open(ROOT/'permutation_stress_runs.csv.gz','wt',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(all_rows[0]));writer.writeheader();writer.writerows(all_rows)
 print(json.dumps({k:v for k,v in summary.items() if k!='transformations'},indent=2))

if __name__=='__main__':main()
