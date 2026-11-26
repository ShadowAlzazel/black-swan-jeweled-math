"""Exact affine pinwheel-style five-triangle dissection, independently implemented.

Conway--Radin basis: https://web.ma.utexas.edu/users/radin/papers/pinwheel.pdf
The explicit coordinates below were derived by partitioning (2+i)T, not copied
from the authors' software. Each child's canonical order is right angle,
long-leg endpoint, short-leg endpoint. Reflected children are retained.
"""
import json, math
from pathlib import Path

ROOT=Path(__file__).parent
BASE=[(0.,0.),(2.,0.),(0.,1.)]
INFLATED=[(0.,0.),(4.,2.),(-1.,2.)]
CHILDREN=[
 [(-1.,2.),(0.,0.),(0.,2.)],
 [(2.,1.),(0.,1.),(0.,0.)],
 [(2.,1.),(0.,2.),(0.,1.)],
 [(0.,2.),(2.,1.),(2.,2.)],
 [(4.,2.),(2.,2.),(2.,1.)],
]

def transform(t,p):
 a,b,c,d,e,f=t;x,y=p
 return [a*x+b*y+c,d*x+e*y+f]
def compose(a,b):
 return [a[0]*b[0]+a[1]*b[3],a[0]*b[1]+a[1]*b[4],a[0]*b[2]+a[1]*b[5]+a[2],
         a[3]*b[0]+a[4]*b[3],a[3]*b[1]+a[4]*b[4],a[3]*b[2]+a[4]*b[5]+a[5]]
def canonical(triangle):
 for right in range(3):
  o=triangle[right];others=[triangle[i] for i in range(3) if i!=right]
  ds=[math.dist(o,p)**2 for p in others]
  if abs(min(ds)-1)<1e-12 and abs(max(ds)-4)<1e-12:
   b=others[ds.index(max(ds))];c=others[ds.index(min(ds))]
   if abs((b[0]-o[0])*(c[0]-o[0])+(b[1]-o[1])*(c[1]-o[1]))<1e-12:
    return [o,b,c]
 raise AssertionError('Not a canonical 1:2:sqrt5 triangle')
def affine_from_canonical(points):
 o,b,c=points
 return [(b[0]-o[0])/2,c[0]-o[0],o[0],(b[1]-o[1])/2,c[1]-o[1],o[1]]
def area(triangle):
 a,b,c=triangle
 return abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def overlap(t,u):
 p=list(t)
 if cross(*u)<0:u=list(reversed(u))
 for a,b in zip(u,u[1:]+u[:1]):
  if not p:return 0
  result=[];previous=p[-1];cp=cross(a,b,previous)
  for current in p:
   cc=cross(a,b,current)
   if (cc>=-1e-10)!=(cp>=-1e-10):
    ratio=cp/(cp-cc)
    result.append([previous[0]+ratio*(current[0]-previous[0]),previous[1]+ratio*(current[1]-previous[1])])
   if cc>=-1e-10:result.append(current)
   previous,cp=current,cc
  p=result
 return abs(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(p,p[1:]+p[:1])))/2 if p else 0
def bbox(p):return [min(v[0] for v in p),max(v[0] for v in p),min(v[1] for v in p),max(v[1] for v in p)]

# A=[[2,-1],[1,2]], A^-1=[[2,1],[-1,2]]/5.
inverse_inflation=[.4,.2,0,-.2,.4,0]
rules=[compose(inverse_inflation,affine_from_canonical(canonical(t))) for t in CHILDREN]
assert abs(sum(area(t) for t in CHILDREN)-area(INFLATED))<1e-12
assert all(overlap(CHILDREN[i],CHILDREN[j])<1e-12 for i in range(5) for j in range(i))
patches=[];current=[([1,0,0,0,1,0],[])]
for level in range(5):
 tiles=[]
 for matrix,address in current:
  angle=round(math.degrees(math.atan2(matrix[3],matrix[0]))%360,9)%360
  tiles.append({'address':address,'transform':matrix,'polygon':[transform(matrix,p) for p in BASE],
    'orientationDegrees':angle,'reflected':matrix[0]*matrix[4]-matrix[1]*matrix[3]<0})
 assert len(tiles)==5**level
 assert abs(sum(area(t['polygon']) for t in tiles)-1)<1e-10
 assert all(abs(area(t['polygon'])-5**(-level))<1e-10 for t in tiles)
 maximum_overlap=0
 if level<=3:
  bounds=[bbox(t['polygon']) for t in tiles]
  for i,t in enumerate(tiles):
   a,b,c,d=bounds[i]
   for j in range(i):
    e,f,g,h=bounds[j]
    if b<e or f<a or d<g or h<c:continue
    overlap_area=overlap(t['polygon'],tiles[j]['polygon'])
    maximum_overlap=max(maximum_overlap,overlap_area)
    assert overlap_area<1e-9,(level,i,j,overlap_area)
 orientations=len({round(t['orientationDegrees'],6) for t in tiles})
 validation={'totalArea':1.,'triangleArea':5**(-level),'orientationClasses':orientations,
  'pairwiseOverlapChecked':level<=3,'maxPairOverlap':maximum_overlap if level<=3 else None,
  'constructionNote':'Every refinement is an affine copy of a validated, area-preserving, nonoverlapping five-triangle partition.'}
 patches.append({'iteration':level,'tiles':tiles,'validation':validation})
 print('level',level,'triangles',len(tiles),'orientation classes',orientations,'reflected',sum(t['reflected'] for t in tiles))
 current=[(compose(t,r),address+[j]) for t,address in current for j,r in enumerate(rules)]

output={'construction':'Independently implemented Conway--Radin pinwheel-style affine substitution',
 'primarySource':'https://web.ma.utexas.edu/users/radin/papers/pinwheel.pdf',
 'baseTriangle':BASE,'inflatedTriangle':INFLATED,'inflationMatrix':[[2,-1],[1,2]],
 'baseInflatedDissection':CHILDREN,'canonicalOrder':['right angle','long-leg endpoint','short-leg endpoint'],
 'orientationDefinition':'atan2 of the first affine column (direction of the long leg), modulo 360 degrees',
 'patches':patches}
(ROOT/'pinwheel_patches.json').write_text(json.dumps(output,separators=(',',':'))+'\n')
