#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
struct RNG { uint64_t x; uint64_t next(){uint64_t z=(x+=0x9e3779b97f4a7c15ULL);z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);} int pick(int n){return int(next()%uint64_t(n));} double unit(){return (next()>>11)*0x1.0p-53;} };
struct Clause {uint64_t pos,neg;array<int,3> var;};
struct Instance {int id,n,m;vector<Clause> clauses;array<uint64_t,3> starts;};
struct Tiling {vector<double> angle;vector<vector<int>> neighbors;};
struct Config {int id,method,tiling;double constant,amplitude;};
struct Eval {double score;int bad;vector<unsigned char> violation;};
int main(int argc,char**argv){
 if(argc!=3){cerr<<"usage: engine input.txt runs.csv\n";return 2;}
 ifstream in(argv[1]);ofstream out(argv[2]);int ni,nc,nt,budget;in>>ni>>nc>>nt>>budget;
 vector<Tiling> tilings(nt);
 for(auto &t:tilings){int l;in>>l;t.angle.resize(l);t.neighbors.resize(l);for(double&a:t.angle)in>>a;for(auto&row:t.neighbors){int k;in>>k;row.resize(k);for(int&v:row)in>>v;}}
 vector<Config> configs(nc);for(auto&c:configs)in>>c.id>>c.method>>c.tiling>>c.constant>>c.amplitude;
 vector<Instance> instances(ni);for(auto&i:instances){in>>i.id>>i.n>>i.m;i.clauses.resize(i.m);for(auto&c:i.clauses)in>>c.pos>>c.neg>>c.var[0]>>c.var[1]>>c.var[2];for(auto&s:i.starts)in>>s;}
 if(!in){cerr<<"invalid serialized dataset\n";return 3;}
 out<<"instance,config,start_index,success,reason,evaluations,flips,updates,restarts,final_state,final_violations\n";
 constexpr double PI=3.1415926535897932384626433832795;
 long long total=0;
 for(const auto &config:configs){
  for(const auto &ins:instances)for(int startidx=0;startidx<3;startidx++){
   RNG rng{0xfeed20261125ULL ^ (uint64_t(ins.id+1)*0x123456789abcdefULL) ^ (uint64_t(config.id+1)*0x9e3779b97f4a7c15ULL) ^ uint64_t(startidx+1)*0xa0761d6478bd642fULL};
   // Phase is shared across methods/constants for a given instance/start.
   RNG phaseRng{0x7adf0312ULL ^ uint64_t(ins.id+1)*0x9e3779b97f4a7c15ULL ^ uint64_t(startidx+1)*0x8ebc6af09c88c6e3ULL};
   const double phase0=2*PI*phaseRng.unit();const auto&t=tilings[config.tiling];
   vector<double> weights(ins.m,1.0);vector<double> angle(ins.m);
   for(int j=0;j<ins.m;j++)angle[j]=t.angle[j]+2*PI*(j+1)/config.constant+phase0;
   if(config.method==4)for(int j=0;j<ins.m;j++)weights[j]=1+config.amplitude*cos(angle[j]);
   int evaluations=0,flips=0,updates=0,restarts=0;uint64_t state=ins.starts[startidx];string reason="budget";
   auto evaluate=[&](uint64_t s){Eval e{0,0,vector<unsigned char>(ins.m)};evaluations++;for(int j=0;j<ins.m;j++){const auto&c=ins.clauses[j];bool bad=((s&c.pos)==0&&(s&c.neg)==c.neg);e.violation[j]=bad;if(bad){e.score+=weights[j];e.bad++;}}return e;};
   Eval current=evaluate(state);
   while(true){
    if(current.bad==0){reason="sat";break;}
    if(flips>=1000){reason="flip_limit";break;}
    vector<int> variables;
    if(config.method==3){vector<int> bad;for(int j=0;j<ins.m;j++)if(current.violation[j])bad.push_back(j);auto&c=ins.clauses[bad[rng.pick(bad.size())]];variables.assign(c.var.begin(),c.var.end());}
    else{variables.resize(ins.n);iota(variables.begin(),variables.end(),0);}
    if(evaluations+int(variables.size())>budget)break;
    vector<Eval> candidates;candidates.reserve(variables.size());for(int v:variables)candidates.push_back(evaluate(state^(1ULL<<v)));
    int chosen=-1;
    if(config.method==3){
     vector<int> breaks(variables.size());for(int k=0;k<int(variables.size());k++)for(int j=0;j<ins.m;j++)breaks[k]+=candidates[k].violation[j]&&!current.violation[j];
     int best=*min_element(breaks.begin(),breaks.end());vector<int> ties;for(int k=0;k<int(breaks.size());k++)if(breaks[k]==best)ties.push_back(k);
     if(best==0||rng.unit()>=0.5)chosen=ties[rng.pick(ties.size())];else chosen=rng.pick(variables.size());
    }else{
     double best=1e300;vector<int> ties;
     for(int k=0;k<int(candidates.size());k++){double s=candidates[k].score;if(s<best-1e-10){best=s;ties={k};}else if(abs(s-best)<=1e-10)ties.push_back(k);}
     if(best<current.score-1e-10)chosen=ties[rng.pick(ties.size())];
    }
    if(chosen>=0){state^=1ULL<<variables[chosen];current=move(candidates[chosen]);flips++;continue;}
    if(config.method==0||config.method==4){reason="trap";break;}
    if(config.method==1){if(evaluations>=budget)break;state=rng.next()&((1ULL<<ins.n)-1);current=evaluate(state);restarts++;continue;}
    if(updates>=100){reason="update_limit";break;}
    vector<double> boost(ins.m,0);double mean=0;
    for(int j=0;j<ins.m;j++)if(current.violation[j]){
     const double a=angle[j]-2*PI*updates/config.constant;
     double value=1;
     if(config.method==5)value=1+config.amplitude*cos(a);
     if(config.method==6||config.method==7)value=0.25+2*pow(1+cos(a),2);
     boost[j]=value;mean+=value;
    }
    mean/=current.bad;
    if(config.method==7){
     vector<double> next=weights;
     for(int j=0;j<ins.m;j++){double s=weights[j];int count=1;for(int neighbor:t.neighbors[j])if(neighbor<ins.m&&neighbor!=j){s+=weights[neighbor];count++;}next[j]=0.8*weights[j]+0.2*s/count;}
     weights.swap(next);
    }
    for(int j=0;j<ins.m;j++)weights[j]+=boost[j]/mean;
    current.score=0;for(int j=0;j<ins.m;j++)if(current.violation[j])current.score+=weights[j];
    updates++;
   }
   out<<ins.id<<','<<config.id<<','<<startidx<<','<<(current.bad==0)<<','<<reason<<','<<evaluations<<','<<flips<<','<<updates<<','<<restarts<<','<<state<<','<<current.bad<<'\n';
   total++;
  }
  out.flush();if((config.id+1)%20==0)cerr<<"completed "<<config.id+1<<" / "<<configs.size()<<" configurations; "<<total<<" runs\n";
 }
 cerr<<"complete: "<<total<<" runs\n";
}
