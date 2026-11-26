"""Validate polygon disjointness and extract edge adjacency without dependencies."""
import json, math
from pathlib import Path

ROOT=Path(__file__).parent
EPS=1e-9
def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def area(p):return sum(p[i][0]*p[(i+1)%len(p)][1]-p[i][1]*p[(i+1)%len(p)][0] for i in range(len(p)))/2
def inside_triangle(p,a,b,c):return cross(a,b,p)>=-EPS and cross(b,c,p)>=-EPS and cross(c,a,p)>=-EPS
def triangulate(poly):
 p=list(poly)
 if area(p)<0:p.reverse()
 result=[]
 while len(p)>3:
  for i in range(len(p)):
   a,b,c=p[i-1],p[i],p[(i+1)%len(p)]
   if cross(a,b,c)<=EPS:continue
   if any(inside_triangle(v,a,b,c) for j,v in enumerate(p) if j not in ((i-1)%len(p),i,(i+1)%len(p))):continue
   result.append([a,b,c]);p.pop(i);break
  else:raise AssertionError('Ear clipping failed')
 result.append(p)
 return result
def clip(poly,a,b):
 if not poly:return []
 result=[]
 previous=poly[-1];cp=cross(a,b,previous)
 for current in poly:
  cc=cross(a,b,current)
  if (cc>=-EPS)!=(cp>=-EPS):
   t=cp/(cp-cc)
   result.append([previous[0]+t*(current[0]-previous[0]),previous[1]+t*(current[1]-previous[1])])
  if cc>=-EPS:result.append(current)
  previous,cp=current,cc
 return result
def intersection_area(t,u):
 p=t
 for a,b in zip(u,u[1:]+u[:1]):p=clip(p,a,b)
 return abs(area(p)) if len(p)>2 else 0
def shared_length(p,q):
 total=0
 for a,b in zip(p,p[1:]+p[:1]):
  dx,dy=b[0]-a[0],b[1]-a[1];length=math.hypot(dx,dy)
  for c,d in zip(q,q[1:]+q[:1]):
   if abs(cross(a,b,c))>EPS or abs(cross(a,b,d))>EPS:continue
   tc=((c[0]-a[0])*dx+(c[1]-a[1])*dy)/(length*length)
   td=((d[0]-a[0])*dx+(d[1]-a[1])*dy)/(length*length)
   total+=max(0,min(1,max(tc,td))-max(0,min(tc,td)))*length
 return total
data=json.loads((ROOT/'hat_patches.json').read_text())
for patch in data['patches']:
 polygons=[t['polygon'] for t in patch['tiles']]
 triangles=[triangulate(p) for p in polygons]
 bounds=[(min(v[0] for v in p),max(v[0] for v in p),min(v[1] for v in p),max(v[1] for v in p)) for p in polygons]
 expected=abs(area(polygons[0]));adjacency=[];maximum_overlap=0
 assert all(abs(abs(area(p))-expected)<1e-8 for p in polygons)
 for i,p in enumerate(polygons):
  a,b,c,d=bounds[i]
  for j in range(i):
   e,f,g,h=bounds[j]
   if b<e-EPS or f<a-EPS or d<g-EPS or h<c-EPS:continue
   overlap=sum(intersection_area(t,u) for t in triangles[i] for u in triangles[j])
   maximum_overlap=max(maximum_overlap,overlap)
   assert overlap<1e-8,(i,j,overlap)
   length=shared_length(p,polygons[j])
   if length>EPS:adjacency.append([j,i])
 connected={0}
 while True:
  nxt=connected|{v for a,b in adjacency if a in connected or b in connected for v in (a,b)}
  if nxt==connected:break
  connected=nxt
 assert len(connected)==len(polygons)
 patch['adjacency']=adjacency
 patch['validation']={'tileArea':expected,'totalArea':expected*len(polygons),'maxPairOverlap':maximum_overlap,'edgeConnected':True,'edgeCount':len(adjacency)}
 print('Validated',len(polygons),'hats,',len(adjacency),'adjacencies; max overlap',maximum_overlap)
(ROOT/'hat_patches.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
