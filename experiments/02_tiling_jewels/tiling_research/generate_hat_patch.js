// Derive finite valid hat patches from the authors' BSD-3-Clause hatviz.
// Source: https://github.com/isohedral/hatviz (see LICENSE alongside).
const fs = require('fs');
const vm = require('vm');
const path = require('path');
const source = ['geometry.js','hat.js'].map(name => fs.readFileSync(path.join(__dirname,name),'utf8')).join('\n');
const program = `
function gatherHatPolygons(node, T=ident, address=[], result=[]) {
  if (node instanceof HatTile) {
    const angle=((Math.atan2(T[3],T[0])*180/Math.PI)%360+360)%360;
    result.push({label:node.label, address, transform:T, reflected:T[0]*T[4]-T[1]*T[3]<0,
      orientationDegrees:(Math.round(angle*1e8)/1e8)%360,
      polygon:hat_outline.map(p=>{const q=transPt(T,p);return [q.x,q.y];})});
  } else {
    node.children.forEach((child,index)=>gatherHatPolygons(child.geom,mul(T,child.T),address.concat(index),result));
  }
  return result;
}
let family=[H_init,T_init,P_init,F_init];
const patches=[];
for(let iteration=0;iteration<=2;iteration++) {
 patches.push({iteration, metatile:'H', tiles:gatherHatPolygons(family[0])});
 if(iteration<2) family=constructMetatiles(constructPatch(...family));
}
JSON.stringify({source:'https://github.com/isohedral/hatviz', license:'BSD-3-Clause',
 geometrySha:'46e19f05545ed0152776c5fbeec74547932273f1',
 hatSha:'d8c1839bd4861c6f871eba83d6a42d010f1d9551', patches});
`;
const output=vm.runInNewContext(source+'\n'+program,{cos:Math.cos,sin:Math.sin,PI:Math.PI});
const destination=path.join(__dirname,'hat_patches.json');
fs.writeFileSync(destination,output+'\n');
for(const patch of JSON.parse(output).patches) console.log('iteration',patch.iteration,'hats',patch.tiles.length,'reflected',patch.tiles.filter(t=>t.reflected).length);
