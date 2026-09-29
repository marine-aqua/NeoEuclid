#include "neo/problems/parabola.hpp"

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <complex>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>

namespace neo::parabola {
namespace {
constexpr double eps=1e-8, pi=3.14159265358979323846;
struct V{double x{},y{};}; struct T{double a{},b{},c{};};
enum class CK{Line,Circle}; enum class PK{Given,Pair,Parabola};
struct Signature {std::uint64_t first{},second{};bool operator==(const Signature&o)const{return first==o.first&&second==o.second;}};
struct P{std::vector<V> v;std::string r;PK kind{PK::Given};int ca{-1},cb{-1},branch{};Signature hash{};};
struct C{CK kind;std::vector<T> v;std::string r;char op{'G'};int pa{-1},pb{-1},pc{-1};std::uint64_t hash{};int branch{};};
struct SignatureHash {std::size_t operator()(const Signature&s)const{return static_cast<std::size_t>(s.first^(s.second+0x9e3779b97f4a7c15ULL+(s.first<<6)+(s.first>>2)));}};
struct StepNode{std::shared_ptr<const StepNode>parent;std::string recipe;std::size_t depth{};};
struct S{std::vector<P> p;std::vector<C> c;std::shared_ptr<const StepNode>history;int cost{};std::string mask;std::uint64_t hash_xor{},hash_sum{};bool used_circle_parabola_point{},used_bisector{};double beam_score{};engine::CostPolicy cost_policy{engine::CostPolicy::Atomic};};
struct Task{S state;std::string mask;};
struct SetShard { std::mutex mutex; std::unordered_set<Signature,SignatureHash> values; };
using Roots=std::vector<std::vector<V>>;
using SharedRoots=std::shared_ptr<const Roots>;
struct GeometryShard {
    std::mutex mutex;
    std::unordered_map<std::uint64_t,SharedRoots> values;
};
struct StateCache {
    std::vector<std::unique_ptr<SetShard>> shards;
    std::size_t per_shard{};
    std::atomic_size_t hits{0}, entries{0};
    explicit StateCache(std::size_t capacity) {
        constexpr std::size_t count=256;
        per_shard=capacity/count;
        for(std::size_t i=0;i<count;++i)shards.push_back(std::make_unique<SetShard>());
    }
    bool insert(Signature value) {
        if(per_shard==0)return true;
        auto& shard=*shards[SignatureHash{}(value)%shards.size()];
        std::lock_guard<std::mutex> lock(shard.mutex);
        if(shard.values.count(value)){++hits;return false;}
        if(shard.values.size()<per_shard){shard.values.insert(value);++entries;}
        return true;
    }
};
struct GeometryCache {
    std::vector<std::unique_ptr<GeometryShard>> parabola_shards,pair_shards;
    std::size_t per_shard{};
    std::atomic_size_t parabola_hits{0},pair_hits{0},entries{0};
    explicit GeometryCache(std::size_t capacity) {
        constexpr std::size_t count=128;
        per_shard=capacity/(2*count);
        for(std::size_t i=0;i<count;++i){parabola_shards.push_back(std::make_unique<GeometryShard>());pair_shards.push_back(std::make_unique<GeometryShard>());}
    }
};
struct LocalGeometryCache {
    struct Generations {
        std::unordered_map<std::uint64_t,SharedRoots> young,old;
        std::size_t capacity{32768};
        SharedRoots find(std::uint64_t key) const {auto it=young.find(key);if(it!=young.end())return it->second;it=old.find(key);return it==old.end()?SharedRoots{}:it->second;}
        void insert(std::uint64_t key,SharedRoots value){if(young.size()>=capacity){old=std::move(young);young.clear();young.reserve(capacity);}young.insert_or_assign(key,std::move(value));}
    };
    Generations parabola,pair;
    std::size_t parabola_hits{},pair_hits{};
};
struct Shared {
    const SearchConfig& cfg; bool tangent{}; std::chrono::steady_clock::time_point start;
    std::atomic<bool> stop{false}; std::atomic_size_t expanded{0},generated{0},duplicates{0},dense_rejections{0},local_parabola_hits{0},local_pair_hits{0},mitm_records{0},mitm_matches{0};
    StateCache states; GeometryCache geometry;
    std::mutex found_mutex; SearchReport report;
    Shared(const SearchConfig& c,bool t):cfg(c),tangent(t),start(std::chrono::steady_clock::now()),states(c.state_cache_entries),geometry(c.geometry_cache_entries){}
};
double dist(V a,V b){return std::hypot(a.x-b.x,a.y-b.y);} 
std::optional<T> line(V p,V q){double a=p.y-q.y,b=q.x-p.x,n=std::hypot(a,b);if(n<eps)return{};a/=n;b/=n;double c=-(a*p.x+b*p.y);if(a<0||(std::abs(a)<eps&&b<0)){a=-a;b=-b;c=-c;}return T{a,b,c};}
std::string q(double x){return std::to_string(std::llround(x*1e6));}
std::string key(const P&p){std::string s;for(auto x:p.v)s+=q(x.x)+","+q(x.y)+";";return s;}
std::string key(const C&c){std::string s=c.kind==CK::Line?"L":"C";for(auto x:c.v)s+=q(x.a)+","+q(x.b)+","+q(x.c)+";";return s;}
std::uint64_t fingerprint(const std::string& value){return static_cast<std::uint64_t>(std::hash<std::string>{}(value));}
std::uint64_t mix64(std::uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
std::int64_t quantized(double value){return std::llround(value*1e6);}
void append_hash(Signature& hash,std::int64_t value){const auto v=static_cast<std::uint64_t>(value);hash.first=mix64(hash.first^(v+0x9e3779b97f4a7c15ULL));hash.second+=mix64(v+hash.second+0xd6e8feb86659fd93ULL);}
Signature numeric_key(const P&point){Signature hash{0x243f6a8885a308d3ULL,0x13198a2e03707344ULL};for(const auto value:point.v){append_hash(hash,quantized(value.x));append_hash(hash,quantized(value.y));}return hash;}
Signature numeric_key(const C&curve){Signature hash{curve.kind==CK::Line?0xa4093822299f31d0ULL:0x082efa98ec4e6c89ULL,0x452821e638d01377ULL};for(const auto value:curve.v){append_hash(hash,quantized(value.a));append_hash(hash,quantized(value.b));append_hash(hash,quantized(value.c));}return hash;}
bool same_numeric_curve(const C&left,const C&right){if(left.kind!=right.kind||left.v.size()!=right.v.size())return false;for(std::size_t i=0;i<left.v.size();++i)if(quantized(left.v[i].a)!=quantized(right.v[i].a)||quantized(left.v[i].b)!=quantized(right.v[i].b)||quantized(left.v[i].c)!=quantized(right.v[i].c))return false;return true;}
Signature state_signature(const S&s,const std::string&suffix={}){auto h=fingerprint(suffix);if(s.used_circle_parabola_point)h^=0xd6e8feb86659fd93ULL;if(s.used_bisector)h^=0xa0761d6478bd642fULL;return{mix64(s.hash_xor^h),mix64(s.hash_sum+h*0x9e3779b97f4a7c15ULL)};}
std::size_t step_count(const S&s){return s.history?s.history->depth:0;}
std::vector<std::string> materialize_steps(const S&s){std::vector<std::string>result(step_count(s));for(auto node=s.history;node;node=node->parent)result[node->depth-1]=node->recipe;return result;}
std::vector<V> meet(T x,T y,CK a,CK b){
 if(a==CK::Line&&b==CK::Line){double d=x.a*y.b-y.a*x.b;if(std::abs(d)<eps)return{};return{{(x.b*y.c-y.b*x.c)/d,(x.c*y.a-y.c*x.a)/d}};}
 if(a==CK::Circle&&b==CK::Line){std::swap(x,y);std::swap(a,b);} 
 if(a==CK::Line){double sd=x.a*y.a+x.b*y.b+x.c;if(std::abs(sd)>y.c+eps)return{};double fx=y.a-sd*x.a,fy=y.b-sd*x.b,h2=y.c*y.c-sd*sd;if(h2<-eps)return{};double h=std::sqrt(std::max(0.,h2));if(h<eps)return{{fx,fy}};return{{fx-x.b*h,fy+x.a*h},{fx+x.b*h,fy-x.a*h}};}
 double dx=y.a-x.a,dy=y.b-x.b,d=std::hypot(dx,dy);if(d<eps||d>x.c+y.c+eps||d<std::abs(x.c-y.c)-eps)return{};double u=(x.c*x.c-y.c*y.c+d*d)/(2*d),h2=x.c*x.c-u*u;if(h2<-eps)return{};double h=std::sqrt(std::max(0.,h2)),mx=x.a+u*dx/d,my=x.b+u*dy/d;if(h<eps)return{{mx,my}};return{{mx-dy*h/d,my+dx*h/d},{mx+dy*h/d,my-dx*h/d}};
}
std::vector<V> parabola(const C&c,std::size_t i){
 std::vector<V>o;auto z=c.v[i];
 if(c.kind==CK::Line){ // a*y^2/4+b*y+c=0
  double A=z.a/4,B=z.b,D=z.c;if(std::abs(A)<eps){if(std::abs(B)>eps){double y=-D/B;o.push_back({y*y/4,y});}return o;}double disc=B*B-4*A*D;if(disc<-eps)return o;double h=std::sqrt(std::max(0.,disc));for(double y:{(-B-h)/(2*A),(-B+h)/(2*A)})if(o.empty()||std::abs(y-o[0].y)>1e-7)o.push_back({y*y/4,y});
 }else{ // y^4 + A*y^2 + B*y + C = 0 (scaled by 16).
  std::vector<double>co={1.,0.,16.-8.*z.a,-32.*z.b,16.*(z.a*z.a+z.b*z.b-z.c*z.c)};
  auto poly=[&](std::complex<double>x){std::complex<double>v=0.;for(double a:co)v=v*x+a;return v;};
  double radius=1.;for(std::size_t j=1;j<co.size();++j)radius=std::max(radius,1+std::abs(co[j]));
  std::vector<std::complex<double>>rr;for(int j=0;j<4;++j)rr.push_back(std::polar(radius,2*pi*j/4+.173));
  for(int it=0;it<160;++it){double movement=0;for(int j=0;j<4;++j){std::complex<double>den=1.;for(int k=0;k<4;++k)if(k!=j)den*=rr[j]-rr[k];if(std::abs(den)<1e-20)den={1e-20,1e-20};auto delta=poly(rr[j])/den;rr[j]-=delta;movement=std::max(movement,std::abs(delta));}if(movement<1e-13)break;}
  for(auto root:rr)if(std::abs(root.imag())<2e-6){double y=root.real();for(int it=0;it<10;++it){double f=((((co[0]*y+co[1])*y+co[2])*y+co[3])*y+co[4]);double df=((4*co[0]*y+3*co[1])*y+2*co[2])*y+co[3];if(std::abs(df)<1e-13)break;y-=f/df;}V p{y*y/4,y};if(std::none_of(o.begin(),o.end(),[&](V q){return dist(p,q)<1e-5;}))o.push_back(p);}
 }
 std::sort(o.begin(),o.end(),[](V a,V b){return a.y<b.y;});return o;
}
std::vector<std::vector<V>> stable(const std::vector<std::vector<V>>& roots){
 if(roots.empty()||roots[0].empty())return{};std::vector<std::vector<V>>tr;for(auto p:roots[0])tr.push_back({p});
 for(std::size_t s=1;s<roots.size();++s){std::set<int>used;std::vector<std::vector<V>>next;for(auto t:tr){int bi=-1;double bd=1e9;for(int j=0;j<(int)roots[s].size();++j)if(!used.count(j)){auto a=t.back(),b=roots[s][j];double d=dist({a.x/(1+std::abs(a.x)),a.y/(1+std::abs(a.y))},{b.x/(1+std::abs(b.x)),b.y/(1+std::abs(b.y))});if(d<bd){bd=d;bi=j;}}if(bi>=0&&bd<.9){used.insert(bi);t.push_back(roots[s][bi]);next.push_back(std::move(t));}}tr=std::move(next);if(tr.empty())break;}return tr;
}
SharedRoots cached_parabola(const C& curve,GeometryCache* cache,LocalGeometryCache* local){
 const auto hash=curve.hash?curve.hash:fingerprint(key(curve));
 if(local)if(auto result=local->parabola.find(hash)){++local->parabola_hits;return result;}
 if(cache&&cache->per_shard){auto& shard=*cache->parabola_shards[hash%cache->parabola_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);auto it=shard.values.find(hash);if(it!=shard.values.end()){++cache->parabola_hits;auto result=it->second;if(local)local->parabola.insert(hash,result);return result;}}
 std::vector<std::vector<V>> roots;for(std::size_t i=0;i<curve.v.size();++i)roots.push_back(parabola(curve,i));SharedRoots result=std::make_shared<const Roots>(stable(roots));
 if(cache&&cache->per_shard){auto& shard=*cache->parabola_shards[hash%cache->parabola_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);if(shard.values.size()<cache->per_shard){auto [it,inserted]=shard.values.emplace(hash,result);if(inserted)++cache->entries;else result=it->second;}}
 if(local)local->parabola.insert(hash,result);
 return result;
}
SharedRoots cached_pair(const C& first,const C& second,GeometryCache* cache,LocalGeometryCache* local){
 auto a=first.hash?first.hash:fingerprint(key(first)),b=second.hash?second.hash:fingerprint(key(second));if(b<a)std::swap(a,b);const auto hash=mix64(a^(mix64(b)+0x9e3779b97f4a7c15ULL));
 if(local)if(auto result=local->pair.find(hash)){++local->pair_hits;return result;}
 if(cache&&cache->per_shard){auto& shard=*cache->pair_shards[hash%cache->pair_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);auto it=shard.values.find(hash);if(it!=shard.values.end()){++cache->pair_hits;auto result=it->second;if(local)local->pair.insert(hash,result);return result;}}
 std::vector<std::vector<V>>roots;for(std::size_t i=0;i<first.v.size();++i)roots.push_back(meet(first.v[i],second.v[i],first.kind,second.kind));SharedRoots result=std::make_shared<const Roots>(stable(roots));
 if(cache&&cache->per_shard){auto& shard=*cache->pair_shards[hash%cache->pair_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);if(shard.values.size()<cache->per_shard){auto [it,inserted]=shard.values.emplace(hash,result);if(inserted)++cache->entries;else result=it->second;}}
 if(local)local->pair.insert(hash,result);
 return result;
}
S initial(const std::vector<double>&deg,bool vertex_mode=false){S s;auto vals=[&](V x){return std::vector<V>(deg.size(),x);};s.p.push_back({vals({0,0}),"O"});s.p.push_back({{},"P"});s.p.push_back({vals({1,0}),"F"});s.p.push_back({vals({-1,0}),"D[directrix-foot]"});for(double d:deg){double u=vertex_mode?1/std::tan(d*pi/180):1/std::tan(d*pi/360);s.p[1].v.push_back(vertex_mode?V{4*u*u,4*u}:V{u*u,2*u});}for(auto&point:s.p)point.hash=numeric_key(point);
 auto given=[&](T t,std::string r){s.c.push_back({CK::Line,std::vector<T>(deg.size(),t),std::move(r)});};given({0,1,0},"x-axis");given({1,0,0},"y-axis");given({1,0,1},"directrix");C side{CK::Line,{},"angle-side"};V vertex=vertex_mode?V{0,0}:V{1,0};for(auto p:s.p[1].v)side.v.push_back(*line(vertex,p));s.c.push_back(std::move(side));for(auto&c:s.c){c.hash=numeric_key(c).first;s.hash_xor^=mix64(c.hash);s.hash_sum+=mix64(c.hash+0x9e3779b97f4a7c15ULL);}return s;}
engine::OperationKind operation_kind(char op){if(op=='L')return engine::OperationKind::LineThrough;if(op=='C')return engine::OperationKind::CircleCenterThrough;if(op=='N')return engine::OperationKind::PerpendicularThrough;if(op=='B')return engine::OperationKind::PerpendicularBisector;if(op=='P')return engine::OperationKind::ParallelThrough;if(op=='A')return engine::OperationKind::AngleBisector;return engine::OperationKind::TransferCircle;}
int operation_cost(char op,engine::CostPolicy policy=engine::CostPolicy::Atomic){return engine::operation_cost(operation_kind(op),policy);}
S configured_initial(const std::vector<double>&degrees,bool vertex_mode,engine::CostPolicy policy){auto s=initial(degrees,vertex_mode);s.cost_policy=policy;return s;}
bool is_circle_parabola_point(const S&s,int index){return index>=0&&index<(int)s.p.size()&&s.p[index].kind==PK::Parabola&&s.p[index].ca>=0&&s.p[index].ca<(int)s.c.size()&&s.c[s.p[index].ca].kind==CK::Circle;}
bool operation_uses_p(const S&s,char op,int ia,int ib,int ic=-1){if(op=='A')return false;if(op=='N'||op=='P')return ia>=0&&ia<(int)s.p.size()&&s.p[ia].r=="P";return (ia>=0&&ia<(int)s.p.size()&&s.p[ia].r=="P")||(ib>=0&&ib<(int)s.p.size()&&s.p[ib].r=="P")||(ic>=0&&ic<(int)s.p.size()&&s.p[ic].r=="P");}
std::optional<S> finish_add(const S&s,C n,int cost,std::size_t maxp,GeometryCache* cache,LocalGeometryCache* local){
 n.hash=numeric_key(n).first;for(auto&c:s.c)if(c.hash==n.hash&&same_numeric_curve(c,n))return{};const int ni=(int)s.c.size();std::unordered_set<Signature,SignatureHash>known;known.reserve(s.p.size()*2+8);for(auto&p:s.p)known.insert(p.hash);std::vector<P>added;
 auto absorb=[&](const std::vector<std::vector<V>>& branches,PK pk,int other){int branch=0;for(auto&v:branches){P p{v,"",pk,ni,other,branch};p.r=(pk==PK::Parabola?"intersect("+n.r+", parabola)":"intersect("+n.r+", "+s.c[other].r+")")+"#"+std::to_string(branch);p.hash=numeric_key(p);if(known.insert(p.hash).second)added.push_back(std::move(p));++branch;}};
 absorb(*cached_parabola(n,cache,local),PK::Parabola,-1);
 for(int j=0;j<(int)s.c.size();++j){absorb(*cached_pair(n,s.c[j],cache,local),PK::Pair,j);if(s.p.size()+added.size()>maxp)return{};}
 S z=s;z.p.reserve(s.p.size()+added.size());for(auto&point:added)z.p.push_back(std::move(point));
 bool used=s.used_circle_parabola_point;if(n.op=='L'||n.op=='C'||n.op=='B')used=used||is_circle_parabola_point(s,n.pa)||is_circle_parabola_point(s,n.pb);else if(n.op=='N'||n.op=='P')used=used||is_circle_parabola_point(s,n.pa);z.used_circle_parabola_point=used;z.used_bisector=s.used_bisector||n.op=='B';z.hash_xor^=mix64(n.hash);z.hash_sum+=mix64(n.hash+0x9e3779b97f4a7c15ULL);char op=n.op;z.c.push_back(std::move(n));z.history=std::make_shared<StepNode>(StepNode{s.history,z.c.back().r,step_count(s)+1});z.cost+=cost;z.mask.push_back(op);return z;
}
std::optional<S> add(const S&s,char op,int ia,int ib,std::size_t maxp,GeometryCache* cache=nullptr,LocalGeometryCache* local=nullptr,int macro_branch=0,int ic=-1){
 C n;n.op=op;n.pa=ia;n.pb=ib;n.pc=ic;n.kind=(op=='C'||op=='T')?CK::Circle:CK::Line;n.branch=macro_branch;
 if(op=='L'||op=='C'){n.r=(op=='L'?"line(":"circle(")+s.p[ia].r+(op=='L'?", ":"; ")+s.p[ib].r+")";for(std::size_t k=0;k<s.p[ia].v.size();++k){auto a=s.p[ia].v[k],b=s.p[ib].v[k];if(op=='L'){auto l=line(a,b);if(!l)return{};n.v.push_back(*l);}else{double r=dist(a,b);if(r<eps||!std::isfinite(r))return{};n.v.push_back({a.x,a.y,r});}}}
 else if(op=='T'){if(ic<0||ic>=(int)s.p.size())return{};n.r="transfer-circle("+s.p[ia].r+", "+s.p[ib].r+"; "+s.p[ic].r+")[cost=5]";for(std::size_t k=0;k<s.p[ia].v.size();++k){double radius=dist(s.p[ia].v[k],s.p[ib].v[k]);auto center=s.p[ic].v[k];if(radius<eps||!std::isfinite(radius))return{};n.v.push_back({center.x,center.y,radius});}}
 else if(op=='B'){n.r="perpendicular-bisector("+s.p[ia].r+", "+s.p[ib].r+")";for(std::size_t k=0;k<s.p[ia].v.size();++k){auto a=s.p[ia].v[k],b=s.p[ib].v[k];double x=b.x-a.x,y=b.y-a.y,norm=std::hypot(x,y);if(norm<eps)return{};x/=norm;y/=norm;double c=-(x*(a.x+b.x)/2+y*(a.y+b.y)/2);if(x<0||(std::abs(x)<eps&&y<0)){x=-x;y=-y;c=-c;}n.v.push_back({x,y,c});}}
 else if(op=='N'||op=='P'){if(ib<0||ib>=(int)s.c.size()||s.c[ib].kind!=CK::Line)return{};n.r=std::string(op=='N'?"perpendicular-through(":"parallel-through(")+s.p[ia].r+", "+s.c[ib].r+")";for(std::size_t k=0;k<s.p[ia].v.size();++k){V p=s.p[ia].v[k];T l=s.c[ib].v[k];double a=op=='N'?-l.b:l.a,b=op=='N'?l.a:l.b,c=-(a*p.x+b*p.y);if(a<0||(std::abs(a)<eps&&b<0)){a=-a;b=-b;c=-c;}n.v.push_back({a,b,c});}}
 else if(op=='A'){if(ia<0||ib<0||ia>=(int)s.c.size()||ib>=(int)s.c.size()||s.c[ia].kind!=CK::Line||s.c[ib].kind!=CK::Line)return{};n.r="angle-bisector("+s.c[ia].r+", "+s.c[ib].r+")#"+std::to_string(macro_branch);double sign=macro_branch==0?1.:-1.;for(std::size_t k=0;k<s.c[ia].v.size();++k){auto x=s.c[ia].v[k],y=s.c[ib].v[k];double a=x.a+sign*y.a,b=x.b+sign*y.b,c=x.c+sign*y.c,norm=std::hypot(a,b);if(norm<eps)return{};a/=norm;b/=norm;c/=norm;if(a<0||(std::abs(a)<eps&&b<0)){a=-a;b=-b;c=-c;}n.v.push_back({a,b,c});}}
 else return{};
 return finish_add(s,std::move(n),operation_cost(op,s.cost_policy),maxp,cache,local);
}

std::optional<S> five_step_e_prefix(const SearchConfig&cfg,GeometryCache*cache,LocalGeometryCache*local){
 S state=configured_initial(cfg.search_degrees,true,engine::CostPolicy::Euclidea);
 auto point_matching=[&](auto expected){
  for(int index=0;index<(int)state.p.size();++index){bool match=true;for(std::size_t k=0;k<state.p[index].v.size();++k)if(dist(state.p[index].v[k],expected(k))>2e-6){match=false;break;}if(match)return index;}return -1;
 };
 auto append=[&](char op,int a,int b){auto next=add(state,op,a,b,cfg.max_points,cache,local);if(!next)return false;state=std::move(*next);return true;};
 // Gamma = circle(F;O), Delta = circle(O;(2,0)).
 if(!append('C',2,0))return{};
 int x2=point_matching([](std::size_t){return V{2,0};});if(x2<0||!append('C',0,x2))return{};
 // Q is the opposite-ray intersection of Delta with the given angle side.
 int qpoint=point_matching([&](std::size_t k){double a=cfg.search_degrees[k]*pi/180;return V{-2*std::cos(a),-2*std::sin(a)};});
 if(qpoint<0||!append('C',qpoint,0))return{};
 const int gamma=4,omega=6;
 int x=point_matching([&](std::size_t k){double a=cfg.search_degrees[k]*pi/180;return V{-4*std::cos(a),0};});
 int y=point_matching([&](std::size_t k){double a=cfg.search_degrees[k]*pi/180;return V{0,-4*std::sin(a)};});
 int r=point_matching([&](std::size_t k){double a=cfg.search_degrees[k]*pi/180;return V{-4*std::cos(a),-4*std::sin(a)};});
 int j=-1;for(int i=0;i<(int)state.p.size();++i)if(state.p[i].kind==PK::Pair&&((state.p[i].ca==gamma&&state.p[i].cb==omega)||(state.p[i].ca==omega&&state.p[i].cb==gamma))){bool origin=true;for(auto v:state.p[i].v)origin=origin&&dist(v,{0,0})<2e-6;if(!origin){j=i;break;}}
 if(x<0||y<0||r<0||j<0||!append('L',x,j)||!append('L',y,r))return{};
 // The prefix remains in history and geometry, while max_cost measures only
 // the extra suffix requested by the caller.
 state.cost=0;state.mask.clear();return state;
}
std::optional<S> k2_initial(const std::vector<double>&degrees,std::size_t maxp,GeometryCache*cache,LocalGeometryCache*local){S s=initial(degrees,true);auto first=add(s,'C',3,0,maxp,cache,local);if(!first)return{};s=std::move(*first);int a2=-1;for(int i=0;i<(int)s.p.size();++i){bool ok=true;for(auto p:s.p[i].v)if(std::abs(p.x+2.)>2e-6||std::abs(p.y)>2e-6){ok=false;break;}if(ok){a2=i;break;}}if(a2<0)return{};auto second=add(s,'L',1,a2,maxp,cache,local);if(!second)return{};return second;}
std::vector<S> expand_operation(const S&s,char op,std::size_t maxp,GeometryCache*cache,LocalGeometryCache*local,bool require_first_p=false){std::vector<S>out;auto push=[&](int a,int b,int branch=0,int c=-1){if(require_first_p&&!s.history&&!operation_uses_p(s,op,a,b,c))return;auto child=add(s,op,a,b,maxp,cache,local,branch,c);if(child)out.push_back(std::move(*child));};if(op=='L'||op=='B'){for(int i=0;i<(int)s.p.size();++i)for(int j=i+1;j<(int)s.p.size();++j)push(i,j);}else if(op=='C'){for(int i=0;i<(int)s.p.size();++i)for(int j=0;j<(int)s.p.size();++j)if(i!=j)push(i,j);}else if(op=='T'){for(int i=0;i<(int)s.p.size();++i)for(int j=i+1;j<(int)s.p.size();++j)for(int k=0;k<(int)s.p.size();++k)if(k!=i&&k!=j)push(i,j,0,k);}else if(op=='N'||op=='P'){for(int i=0;i<(int)s.p.size();++i)for(int j=0;j<(int)s.c.size();++j)if(s.c[j].kind==CK::Line)push(i,j);}else if(op=='A'){for(int i=0;i<(int)s.c.size();++i)if(s.c[i].kind==CK::Line)for(int j=i+1;j<(int)s.c.size();++j)if(s.c[j].kind==CK::Line){push(i,j,0);push(i,j,1);}}return out;}
double target_angle(double d,bool tangent){return (tangent?d/2:d/3)*pi/180;}
V goal_vertex(const S&s,std::size_t sample,bool tangent,bool vertex_mode){return tangent?s.p[1].v[sample]:(vertex_mode?V{0,0}:V{1,0});}
bool target_line(const S&s,const C&c,bool tangent,bool vertex_mode){if(c.kind!=CK::Line)return false;for(std::size_t i=0;i<c.v.size();++i){V v=goal_vertex(s,i,tangent,vertex_mode);double t=target_angle((i==0?0:0),tangent);(void)t;auto l=c.v[i]; // direction test is supplied separately below
 if(std::abs(l.a*v.x+l.b*v.y+l.c)>2e-6)return false;}return true;}
bool line_goal(const S&s,const C&c,const std::vector<double>&deg,bool tangent,bool vertex_mode){if(!target_line(s,c,tangent,vertex_mode))return false;for(std::size_t i=0;i<deg.size();++i){double t=target_angle(deg[i],tangent);auto l=c.v[i];if(std::abs(l.a*std::cos(t)+l.b*std::sin(t))>2e-6)return false;}return true;}
bool point_goal(const S&s,const P&p,const std::vector<double>&deg,bool tangent,bool vertex_mode){for(std::size_t i=0;i<deg.size();++i){V v=goal_vertex(s,i,tangent,vertex_mode);double t=target_angle(deg[i],tangent),dx=p.v[i].x-v.x,dy=p.v[i].y-v.y,radius=std::hypot(dx,dy),forward=dx*std::cos(t)+dy*std::sin(t),cross=-std::sin(t)*dx+std::cos(t)*dy;if(radius<=1e-5||std::abs(forward)<=1e-5||std::abs(cross)>2e-7*std::max(1.,radius))return false;}return true;}
bool center_goal(const S&s,const P&p){for(std::size_t i=0;i<p.v.size();++i){const auto input=s.p[1].v[i];V wanted{8.+4.5*input.x,4.*input.y};if(dist(p.v[i],wanted)>3e-6*std::max(1.,std::hypot(wanted.x,wanted.y)))return false;}return true;}
V e_target(const S&s,std::size_t i){const auto input=s.p[1].v[i];return {8.,-4.*input.y/std::hypot(input.x,input.y)};}
bool e_goal(const S&s,const P&p){for(std::size_t i=0;i<p.v.size();++i){const auto wanted=e_target(s,i);if(dist(p.v[i],wanted)>3e-6*std::max(1.,std::hypot(wanted.x,wanted.y)))return false;}return true;}
bool k2_circle_goal(const S&s,const C&curve){if(curve.kind!=CK::Circle)return false;for(std::size_t i=0;i<curve.v.size();++i){double t=s.p[1].v[i].y/2.;if(std::abs(t)<eps)return false;double qv=2./t,h=3.5+4.5*qv*qv,k=2.*qv,r=std::sqrt(h*h+qv*qv);auto actual=curve.v[i];double scale=std::max({1.,std::abs(h),std::abs(k),r});if(std::abs(actual.a-h)>3e-6*scale||std::abs(actual.b-k)>3e-6*scale||std::abs(actual.c-r)>3e-6*scale)return false;}return true;}
bool input_dependent(const std::string&recipe){return recipe.find('P')!=std::string::npos||recipe.find("angle-side")!=std::string::npos;}
std::optional<std::pair<std::string,int>> terminal(const S&s,const std::vector<double>&d,bool tangent,int maxcost,bool vertex_mode=false,bool require_circle_parabola=false,bool target_center=false,bool require_bisector=false,bool target_k2=false,bool target_e=false){if(require_bisector&&!s.used_bisector)return{};if(target_e){for(const auto&p:s.p)if(input_dependent(p.r)&&e_goal(s,p))return{{"constructed-E-point",s.cost}};return{};}if(target_k2){for(const auto&c:s.c)if(input_dependent(c.r)&&k2_circle_goal(s,c))return{{"k2-triple-angle-circle",s.cost}};return{};}if(target_center){for(const auto&p:s.p)if(input_dependent(p.r)&&center_goal(s,p))return{{"constructed-circle-center",s.cost}};return{};}for(auto&c:s.c)if((!require_circle_parabola||s.used_circle_parabola_point)&&input_dependent(c.r)&&line_goal(s,c,d,tangent,vertex_mode))return{{"direct-target-line",s.cost}};if(s.cost+1<=maxcost)for(int i=0;i<(int)s.p.size();++i){const auto&p=s.p[i];if((!require_circle_parabola||s.used_circle_parabola_point||is_circle_parabola_point(s,i))&&input_dependent(p.r)&&point_goal(s,p,d,tangent,vertex_mode))return{{p.kind==PK::Parabola?"curve-parabola-intersection+join":(p.kind==PK::Pair?"two-curves-intersection+join":"existing-point+join"),s.cost+1}};}return{};}
std::vector<std::string> masks(int n,const std::vector<std::string>&requested){if(!requested.empty())return requested;std::vector<std::string>o;for(int bits=0;bits<(1<<n);++bits){std::string m;for(int i=0;i<n;++i)m.push_back(bits&(1<<i)?'C':'L');o.push_back(m);}return o;}
bool expired(Shared&x){return x.stop||std::chrono::duration<double>(std::chrono::steady_clock::now()-x.start).count()>=x.cfg.time_limit_seconds||x.generated>=x.cfg.max_states;}
std::vector<std::vector<V>> target_ray_branches(const C&curve,const std::vector<double>&degrees,bool vertex_mode){std::vector<std::vector<V>>roots;V v=vertex_mode?V{0,0}:V{1,0};for(std::size_t i=0;i<degrees.size();++i){double t=target_angle(degrees[i],false);T ray=*line(v,{v.x+std::cos(t),v.y+std::sin(t)});std::vector<V>forward;for(auto p:meet(curve.v[i],ray,curve.kind,CK::Line))if((p.x-v.x)*std::cos(t)+(p.y-v.y)*std::sin(t)>1e-5)forward.push_back(p);roots.push_back(std::move(forward));}return stable(roots);}
std::vector<std::vector<V>> target_e_branches(const S&s,const C&curve){std::vector<V>branch;for(std::size_t i=0;i<curve.v.size();++i){const V p=e_target(s,i);const T c=curve.v[i];const double error=curve.kind==CK::Line?std::abs(c.a*p.x+c.b*p.y+c.c):std::abs(std::hypot(p.x-c.a,p.y-c.b)-c.c);if(error>3e-6*std::max(1.,std::hypot(p.x,p.y)))return{};branch.push_back(p);}return {std::move(branch)};}
Signature ray_signature(const std::vector<V>&branch,bool vertex_mode){Signature result{0x6a09e667f3bcc909ULL,0xbb67ae8584caa73bULL};const double vx=vertex_mode?0.:1.;for(auto point:branch){const double d=std::hypot(point.x-vx,point.y);append_hash(result,std::llround((d/(1+d))*1e6));}return result;}
bool same_branch(const std::vector<V>&a,const std::vector<V>&b){if(a.size()!=b.size())return false;for(std::size_t i=0;i<a.size();++i)if(dist(a[i],b[i])>2e-5*std::max({1.,dist(a[i],{1,0}),dist(b[i],{1,0})}))return false;return true;}
bool transverse(const C&a,const C&b,const std::vector<V>&points){if(a.kind==b.kind&&a.hash==b.hash)return false;for(std::size_t i=0;i<points.size();++i){V p=points[i];T x=a.v[i],y=b.v[i];double measure=0;if(a.kind==CK::Line&&b.kind==CK::Line)measure=std::abs(x.a*y.b-x.b*y.a);else if(a.kind==CK::Line)measure=std::abs(x.a*(p.y-y.b)-x.b*(p.x-y.a))/std::max(1.,y.c);else if(b.kind==CK::Line)measure=std::abs(y.a*(p.y-x.b)-y.b*(p.x-x.a))/std::max(1.,x.c);else measure=std::abs((p.x-x.a)*(p.y-y.b)-(p.y-x.b)*(p.x-y.a))/std::max(1.,x.c*y.c);if(measure<2e-6)return false;}return true;}
struct ArmRecord{std::size_t state_index{};std::vector<V>branch;Signature ray;};
struct ArmSet{std::vector<S>states;std::vector<ArmRecord>records;std::unordered_set<Signature,SignatureHash>record_seen,state_seen;};
Signature arm_identity(const S&state,const C&curve,Signature ray){append_hash(ray,static_cast<std::int64_t>(curve.hash));append_hash(ray,curve.op);for(int point_index:{curve.pa,curve.pb})if(point_index>=0&&point_index<(int)state.p.size()){append_hash(ray,static_cast<std::int64_t>(state.p[point_index].hash.first));append_hash(ray,static_cast<std::int64_t>(state.p[point_index].hash.second));}return ray;}
void enumerate_arm(const S&base,const std::string&arm,int pos,ArmSet&out,Shared&sh,LocalGeometryCache&geometry){if(expired(sh))return;if(pos==(int)arm.size()){const auto&curve=base.c.back();std::optional<std::size_t>state_index;for(auto branch:target_ray_branches(curve,sh.cfg.search_degrees,sh.cfg.angle_at_parabola_vertex)){auto ray=ray_signature(branch,sh.cfg.angle_at_parabola_vertex);auto identity=arm_identity(base,curve,ray);if(!out.record_seen.insert(identity).second)continue;if(!state_index){state_index=out.states.size();out.states.push_back(base);}out.records.push_back({*state_index,std::move(branch),ray});}return;}char op=arm[pos];for(int i=0;i<(int)base.p.size();++i)for(int j=0;j<(int)base.p.size();++j){if(i==j||(op=='L'&&j<i))continue;auto child=add(base,op,i,j,sh.cfg.max_points,&sh.geometry,&geometry);++sh.generated;if(!child)continue;auto sig=state_signature(*child,arm.substr(pos+1));if(!out.state_seen.insert(sig).second){++sh.duplicates;continue;}enumerate_arm(*child,arm,pos+1,out,sh,geometry);if(expired(sh))return;}}
const std::vector<char>& euclidea_mitm_operations(){static const std::vector<char>ops{'L','C','B','N','P','A'};return ops;}
void enumerate_dynamic_arm(const S&base,int remaining,ArmSet&out,Shared&sh,LocalGeometryCache&geometry){if(expired(sh))return;if(remaining==0){const auto&curve=base.c.back();std::optional<std::size_t>state_index;for(auto branch:target_e_branches(base,curve)){auto signature=ray_signature(branch,true);auto identity=arm_identity(base,curve,signature);if(!out.record_seen.insert(identity).second)continue;if(!state_index){state_index=out.states.size();out.states.push_back(base);}out.records.push_back({*state_index,std::move(branch),signature});}return;}for(char op:euclidea_mitm_operations())for(auto&child:expand_operation(base,op,sh.cfg.max_points,&sh.geometry,&geometry)){++sh.generated;if(expired(sh))return;auto signature=state_signature(child,std::to_string(remaining-1));if(!out.state_seen.insert(signature).second){++sh.duplicates;continue;}enumerate_dynamic_arm(child,remaining-1,out,sh,geometry);}}
void enumerate_schema_arm(const S&base,const std::string&schema,std::size_t pos,ArmSet&out,Shared&sh,LocalGeometryCache&geometry){if(expired(sh))return;if(pos==schema.size()){const auto&curve=base.c.back();std::optional<std::size_t>state_index;for(auto branch:target_e_branches(base,curve)){auto signature=ray_signature(branch,true);auto identity=arm_identity(base,curve,signature);if(!out.record_seen.insert(identity).second)continue;if(!state_index){state_index=out.states.size();out.states.push_back(base);}out.records.push_back({*state_index,std::move(branch),signature});}return;}for(auto&child:expand_operation(base,schema[pos],sh.cfg.max_points,&sh.geometry,&geometry)){++sh.generated;if(expired(sh))return;auto signature=state_signature(child,schema.substr(pos+1));if(!out.state_seen.insert(signature).second){++sh.duplicates;continue;}enumerate_schema_arm(child,schema,pos+1,out,sh,geometry);}}
std::optional<S> merge_arms(const S&base,const S&left,const S&right,Shared&sh,LocalGeometryCache&geometry){S merged=base;auto append=[&](const S&source){auto point_index=[&](int source_index){if(source_index<0||source_index>=(int)source.p.size())return -1;for(int i=0;i<(int)merged.p.size();++i)if(merged.p[i].r==source.p[source_index].r)return i;return -1;};auto curve_index=[&](int source_index){if(source_index<0||source_index>=(int)source.c.size())return -1;for(int i=0;i<(int)merged.c.size();++i)if(merged.c[i].r==source.c[source_index].r)return i;return -1;};for(std::size_t ci=base.c.size();ci<source.c.size();++ci){const auto&curve=source.c[ci];bool exists=false;for(const auto&known:merged.c)if(known.hash==curve.hash&&key(known)==key(curve)){exists=true;break;}if(exists)continue;int pa=-1,pb=-1,pc=-1;if(curve.op=='N'||curve.op=='P'){pa=point_index(curve.pa);pb=curve_index(curve.pb);}else if(curve.op=='A'){pa=curve_index(curve.pa);pb=curve_index(curve.pb);}else{pa=point_index(curve.pa);pb=point_index(curve.pb);if(curve.op=='T')pc=point_index(curve.pc);}if(pa<0||pb<0||(curve.op=='T'&&pc<0))return false;auto next=add(merged,curve.op,pa,pb,sh.cfg.max_points,&sh.geometry,&geometry,curve.branch,pc);if(!next)return false;merged=std::move(*next);}return true;};if(!append(left)||!append(right))return{};return merged;}
// Replay on a dense grid by re-running the exact operation choices. Branches are
// reconstructed geometrically; a candidate must independently satisfy every dense sample.
bool dense_verify(const S&hit,bool tangent,int final_cost,const SearchConfig&cfg){
 std::vector<double>d;double high=cfg.validation_max_degrees>0?cfg.validation_max_degrees:(cfg.angle_at_parabola_vertex?85.:175.);for(int i=0;i<cfg.validation_samples;++i)d.push_back(5.+(high-5.)*i/(cfg.validation_samples-1));S s=configured_initial(d,cfg.angle_at_parabola_vertex,cfg.cost_policy);
 for(std::size_t step=0;step<step_count(hit);++step){const C&hc=hit.c[4+step];int ia=-1,ib=-1,ic=-1;if(hc.op=='L'||hc.op=='C'||hc.op=='B'||hc.op=='T'){for(int i=0;i<(int)s.p.size();++i){if(s.p[i].r==hit.p[hc.pa].r)ia=i;if(s.p[i].r==hit.p[hc.pb].r)ib=i;if(hc.op=='T'&&s.p[i].r==hit.p[hc.pc].r)ic=i;}}else if(hc.op=='N'||hc.op=='P'){for(int i=0;i<(int)s.p.size();++i)if(s.p[i].r==hit.p[hc.pa].r)ia=i;for(int i=0;i<(int)s.c.size();++i)if(s.c[i].r==hit.c[hc.pb].r)ib=i;}else if(hc.op=='A'){for(int i=0;i<(int)s.c.size();++i){if(s.c[i].r==hit.c[hc.pa].r)ia=i;if(s.c[i].r==hit.c[hc.pb].r)ib=i;}}
  if(ia<0||ib<0||(hc.op=='T'&&ic<0))return false;auto z=add(s,hc.op,ia,ib,cfg.max_points*4,nullptr,nullptr,hc.branch,ic);if(!z)return false;s=std::move(*z);
  if(cfg.retain_five_step_e_prefix&&step==4){s.cost=0;s.mask.clear();}
 }
 return terminal(s,d,tangent,final_cost,cfg.angle_at_parabola_vertex,cfg.require_circle_parabola_use,cfg.target_circle_center,cfg.require_bisector_use,cfg.target_k2_circle,cfg.target_e_point).has_value();
}
void dfs(S s,const std::string&mask,int pos,Shared&sh,LocalGeometryCache&geometry){if(expired(sh))return;sh.expanded++;if(auto t=terminal(s,sh.cfg.search_degrees,sh.tangent,sh.cfg.max_cost,sh.cfg.angle_at_parabola_vertex,sh.cfg.require_circle_parabola_use,sh.cfg.target_circle_center,sh.cfg.require_bisector_use,sh.cfg.target_k2_circle)){if(dense_verify(s,sh.tangent,t->second,sh.cfg)){std::lock_guard<std::mutex>g(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=t->second;sh.report.steps=materialize_steps(s);if(t->second>s.cost)sh.report.steps.push_back(std::string("line(")+(sh.tangent?"P":(sh.cfg.angle_at_parabola_vertex?"O":"F")) +", target-point)");sh.report.terminal_schema=t->first;sh.report.mask=s.mask;sh.stop=true;}return;}sh.dense_rejections++;}
 if(pos>=(int)mask.size()||s.cost>=sh.cfg.max_cost)return;char op=mask[pos];if(s.cost+operation_cost(op)>sh.cfg.max_cost)return;std::unordered_set<Signature,SignatureHash>local;for(auto&z:expand_operation(s,op,sh.cfg.max_points,&sh.geometry,&geometry,sh.cfg.require_first_step_uses_p)){sh.generated++;auto signature=state_signature(z,mask.substr(pos+1));if(!local.insert(signature).second){sh.duplicates++;continue;}if(!sh.states.insert(signature)){sh.duplicates++;continue;}dfs(std::move(z),mask,pos+1,sh,geometry);if(expired(sh))return;}}

double beam_score(const S& state, const SearchConfig& cfg) {
 double best=1e30;
 if(cfg.target_e_point){for(const auto&point:state.p){if(!input_dependent(point.r))continue;double worst=0.;for(std::size_t i=0;i<point.v.size();++i){const auto wanted=e_target(state,i);worst=std::max(worst,dist(point.v[i],wanted)/std::max(1.,std::hypot(wanted.x,wanted.y)));}best=std::min(best,worst);}return best+1e-6*state.p.size();}
 for(const auto& point:state.p){if(!input_dependent(point.r))continue;double worst=0.;bool valid=true;for(std::size_t i=0;i<cfg.search_degrees.size();++i){const V vertex=cfg.angle_at_parabola_vertex?V{0,0}:V{1,0};const double angle=target_angle(cfg.search_degrees[i],false),dx=point.v[i].x-vertex.x,dy=point.v[i].y-vertex.y,r=std::hypot(dx,dy),along=dx*std::cos(angle)+dy*std::sin(angle),cross=std::abs(-std::sin(angle)*dx+std::cos(angle)*dy)/std::max(1.,r);if(!std::isfinite(cross)||r<1e-8){valid=false;break;}worst=std::max(worst,cross+(along<=0?1.:0.));}if(valid)best=std::min(best,worst);}
 for(const auto& curve:state.c){if(curve.kind!=CK::Line||!input_dependent(curve.r))continue;double worst=0.;for(std::size_t i=0;i<cfg.search_degrees.size();++i){const V vertex=cfg.angle_at_parabola_vertex?V{0,0}:V{1,0};const double angle=target_angle(cfg.search_degrees[i],false);const auto line=curve.v[i];worst=std::max(worst,std::abs(line.a*vertex.x+line.b*vertex.y+line.c)+std::abs(line.a*std::cos(angle)+line.b*std::sin(angle)));}best=std::min(best,worst);}
 return best+1e-6*state.p.size();
}

struct BoundedBeam {
 struct Entry{S state;std::uint64_t serial{};};
 struct Ranked{double score{};Signature signature{};std::uint64_t serial{};};
 struct WorseFirst{bool operator()(const Ranked&a,const Ranked&b)const{return a.score<b.score;}};
 std::size_t capacity{};std::uint64_t next_serial{};
 std::unordered_map<Signature,Entry,SignatureHash> states;
 std::priority_queue<Ranked,std::vector<Ranked>,WorseFirst> worst;
 explicit BoundedBeam(std::size_t cap=0):capacity(cap){states.reserve(cap*2+1);}
 void clean(){while(!worst.empty()){auto it=states.find(worst.top().signature);if(it!=states.end()&&it->second.serial==worst.top().serial)break;worst.pop();}}
 bool insert(Signature signature,S state,const SearchConfig&cfg){if(states.find(signature)!=states.end())return false;state.beam_score=beam_score(state,cfg);clean();if(capacity==0)return true;if(states.size()>=capacity){if(!worst.empty()&&state.beam_score>=worst.top().score)return true;auto doomed=worst.top();worst.pop();states.erase(doomed.signature);}const auto serial=++next_serial;worst.push({state.beam_score,signature,serial});states.emplace(signature,Entry{std::move(state),serial});return true;}
 std::vector<S> take(){std::vector<S>result;result.reserve(states.size());for(auto&entry:states)result.push_back(std::move(entry.second.state));states.clear();while(!worst.empty())worst.pop();return result;}
};

std::vector<char> beam_operations(const S& state,const SearchConfig& cfg){if(cfg.masks.empty()){if(cfg.cost_policy!=engine::CostPolicy::Euclidea)return{'L','C'};return cfg.allow_transfer_circle?std::vector<char>{'L','C','B','N','P','A','T'}:std::vector<char>{'L','C','B','N','P','A'};}std::vector<char> result;for(const auto& mask:cfg.masks)if(mask.size()>state.mask.size()&&mask.compare(0,state.mask.size(),state.mask)==0&&std::find(result.begin(),result.end(),mask[state.mask.size()])==result.end())result.push_back(mask[state.mask.size()]);return result;}

SearchReport run_beam(const SearchConfig& cfg){Shared sh(cfg,false);LocalGeometryCache geometry;std::vector<std::vector<S>>buckets(static_cast<std::size_t>(cfg.max_cost+1));if(cfg.retain_five_step_e_prefix){auto prefix=five_step_e_prefix(cfg,&sh.geometry,&geometry);if(!prefix)throw std::runtime_error("failed to reconstruct five-step E prefix");buckets[0].push_back(std::move(*prefix));}else buckets[0].push_back(configured_initial(cfg.search_degrees,cfg.angle_at_parabola_vertex,cfg.cost_policy));
 for(int cost=0;cost<=cfg.max_cost&&!expired(sh);++cost){auto&states=buckets[static_cast<std::size_t>(cost)];for(auto&state:states)state.beam_score=beam_score(state,cfg);auto by_score=[](const S&a,const S&b){return a.beam_score<b.beam_score;};if(states.size()>cfg.beam_width){std::nth_element(states.begin(),states.begin()+cfg.beam_width,states.end(),by_score);states.resize(cfg.beam_width);}std::sort(states.begin(),states.end(),by_score);std::vector<BoundedBeam>generated;generated.reserve(static_cast<std::size_t>(cfg.max_cost+1));for(int i=0;i<=cfg.max_cost;++i)generated.emplace_back(cfg.beam_width);
  for(const auto&state:states){if(expired(sh))break;++sh.expanded;if(auto hit=terminal(state,cfg.search_degrees,false,cfg.max_cost,cfg.angle_at_parabola_vertex,cfg.require_circle_parabola_use,cfg.target_circle_center,cfg.require_bisector_use,cfg.target_k2_circle,cfg.target_e_point)){if(dense_verify(state,false,hit->second,cfg)){sh.report.found=sh.report.densely_verified=true;sh.report.cost=hit->second;sh.report.steps=materialize_steps(state);if(hit->second>state.cost)sh.report.steps.push_back(std::string("line(")+(cfg.angle_at_parabola_vertex?"O":"F")+", target-point)");sh.report.terminal_schema=hit->first;sh.report.mask=state.mask;sh.stop=true;break;}++sh.dense_rejections;}
   for(char op:beam_operations(state,cfg)){if(state.cost+operation_cost(op,cfg.cost_policy)>cfg.max_cost)continue;for(auto&child:expand_operation(state,op,cfg.max_points,&sh.geometry,&geometry,cfg.require_first_step_uses_p)){++sh.generated;if(expired(sh))break;auto signature=state_signature(child);if(generated[static_cast<std::size_t>(child.cost)].insert(signature,std::move(child),cfg))continue;++sh.duplicates;}}}
  for(int next=cost+1;next<=cfg.max_cost;++next){auto selected=generated[static_cast<std::size_t>(next)].take();auto&bucket=buckets[static_cast<std::size_t>(next)];bucket.insert(bucket.end(),std::make_move_iterator(selected.begin()),std::make_move_iterator(selected.end()));}}
 sh.local_parabola_hits=geometry.parabola_hits;sh.local_pair_hits=geometry.pair_hits;sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_entries=sh.states.entries;sh.report.state_cache_hits=sh.states.hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;}

SearchReport run_exhaustive_three_after_e(const SearchConfig&input){
 SearchConfig cfg=input;cfg.max_cost=3;cfg.max_states=std::numeric_limits<std::size_t>::max();cfg.time_limit_seconds=1e12;cfg.retain_five_step_e_prefix=true;cfg.allow_transfer_circle=false;cfg.angle_at_parabola_vertex=true;cfg.cost_policy=engine::CostPolicy::Euclidea;
 Shared sh(cfg,false);LocalGeometryCache root_geometry;auto root=five_step_e_prefix(cfg,&sh.geometry,&root_geometry);if(!root)throw std::runtime_error("failed to reconstruct five-step E prefix");
 std::deque<S>queue;std::unordered_set<Signature,SignatureHash>first_seen;
 for(char op:euclidea_mitm_operations())for(auto&child:expand_operation(*root,op,cfg.max_points,&sh.geometry,&root_geometry)){++sh.generated;auto signature=state_signature(child,"2");if(first_seen.insert(signature).second&&sh.states.insert(signature))queue.push_back(std::move(child));else ++sh.duplicates;}
 sh.report.prefixes=queue.size();sh.local_parabola_hits+=root_geometry.parabola_hits;sh.local_pair_hits+=root_geometry.pair_hits;std::mutex queue_mutex;
 const int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>workers;
 for(int worker=0;worker<nt;++worker)workers.emplace_back([&]{LocalGeometryCache geometry;std::function<void(S&&,int)>visit=[&](S&&state,int depth){if(sh.stop)return;++sh.expanded;if(auto goal=terminal(state,cfg.search_degrees,false,3,true)){if(dense_verify(state,false,goal->second,cfg)){std::lock_guard<std::mutex>lock(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=goal->second;sh.report.steps=materialize_steps(state);if(goal->second>state.cost)sh.report.steps.push_back("line(O, target-point)");sh.report.terminal_schema=goal->first;sh.report.mask=state.mask;sh.stop=true;}return;}++sh.dense_rejections;}if(depth==3)return;for(char op:euclidea_mitm_operations())for(auto&child:expand_operation(state,op,cfg.max_points,&sh.geometry,&geometry)){++sh.generated;if(sh.stop)return;auto signature=state_signature(child,std::to_string(2-depth));if(!sh.states.insert(signature)){++sh.duplicates;continue;}visit(std::move(child),depth+1);if(sh.stop)return;}};for(;;){S state;{std::lock_guard<std::mutex>lock(queue_mutex);if(queue.empty()||sh.stop)break;state=std::move(queue.front());queue.pop_front();}visit(std::move(state),1);}sh.local_parabola_hits+=geometry.parabola_hits;sh.local_pair_hits+=geometry.pair_hits;});
 for(auto&worker:workers)worker.join();sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_hits=sh.states.hits;sh.report.state_cache_entries=sh.states.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;
}

SearchReport run_exhaustive_cost4_e_sequential(const SearchConfig&input){
 SearchConfig cfg=input;cfg.max_cost=4;cfg.max_states=std::numeric_limits<std::size_t>::max();cfg.target_e_point=true;cfg.angle_at_parabola_vertex=true;cfg.cost_policy=engine::CostPolicy::Euclidea;cfg.allow_transfer_circle=false;
 Shared sh(cfg,false);LocalGeometryCache root_geometry;S root=configured_initial(cfg.search_degrees,true,cfg.cost_policy);std::deque<S>queue;std::unordered_set<Signature,SignatureHash>first_seen;
 for(char op:euclidea_mitm_operations())for(auto&child:expand_operation(root,op,cfg.max_points,&sh.geometry,&root_geometry)){++sh.generated;auto signature=state_signature(child,"3");if(first_seen.insert(signature).second&&sh.states.insert(signature))queue.push_back(std::move(child));else ++sh.duplicates;}
 sh.report.prefixes=queue.size();sh.local_parabola_hits+=root_geometry.parabola_hits;sh.local_pair_hits+=root_geometry.pair_hits;std::mutex queue_mutex;
 const int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>workers;
 for(int worker=0;worker<nt;++worker)workers.emplace_back([&]{LocalGeometryCache geometry;std::function<void(S&&,int)>visit=[&](S&&state,int depth){if(expired(sh))return;++sh.expanded;if(auto goal=terminal(state,cfg.search_degrees,false,4,true,false,false,false,false,true)){if(dense_verify(state,false,goal->second,cfg)){std::lock_guard<std::mutex>lock(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=goal->second;sh.report.steps=materialize_steps(state);sh.report.terminal_schema=goal->first;sh.report.mask=state.mask;sh.stop=true;}return;}++sh.dense_rejections;}if(depth==4)return;for(char op:euclidea_mitm_operations())for(auto&child:expand_operation(state,op,cfg.max_points,&sh.geometry,&geometry)){++sh.generated;if(expired(sh))return;auto signature=state_signature(child,std::to_string(4-depth-1));if(!sh.states.insert(signature)){++sh.duplicates;continue;}visit(std::move(child),depth+1);if(expired(sh))return;}};for(;;){S state;{std::lock_guard<std::mutex>lock(queue_mutex);if(queue.empty()||expired(sh))break;state=std::move(queue.front());queue.pop_front();}visit(std::move(state),1);}sh.local_parabola_hits+=geometry.parabola_hits;sh.local_pair_hits+=geometry.pair_hits;});
 for(auto&worker:workers)worker.join();sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_hits=sh.states.hits;sh.report.state_cache_entries=sh.states.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;
}

SearchReport run(const SearchConfig&cfg,bool tangent){Shared sh(cfg,tangent);auto ms=masks(cfg.max_cost,cfg.masks);std::deque<Task>q;
 LocalGeometryCache prefix_geometry;std::unordered_set<Signature,SignatureHash>prefix_seen;auto configured_base=cfg.target_k2_circle?k2_initial(cfg.search_degrees,cfg.max_points,&sh.geometry,&prefix_geometry):std::optional<S>{initial(cfg.search_degrees,cfg.angle_at_parabola_vertex)};if(!configured_base)throw std::runtime_error("failed to build k=2 reciprocal prefix");for(auto&m:ms){std::vector<S>layer{*configured_base};for(int depth=0;depth<std::min(cfg.prefix_depth,(int)m.size());++depth){std::vector<S>next;for(auto&s:layer){if(s.cost+operation_cost(m[depth])>cfg.max_cost)continue;for(auto&z:expand_operation(s,m[depth],cfg.max_points,&sh.geometry,&prefix_geometry,cfg.require_first_step_uses_p)){sh.generated++;auto signature=state_signature(z,m.substr(depth+1));if(prefix_seen.insert(signature).second)next.push_back(std::move(z));else sh.duplicates++;}}layer=std::move(next);}for(auto&s:layer){sh.states.insert(state_signature(s,m.substr(cfg.prefix_depth)));q.push_back({std::move(s),m});}}
 sh.local_parabola_hits+=prefix_geometry.parabola_hits;sh.local_pair_hits+=prefix_geometry.pair_hits;sh.report.prefixes=q.size();std::mutex qm;int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>w;for(int t=0;t<nt;++t)w.emplace_back([&]{LocalGeometryCache local_geometry;while(!expired(sh)){Task x;{std::lock_guard<std::mutex>g(qm);if(q.empty())break;x=std::move(q.front());q.pop_front();}dfs(std::move(x.state),x.mask,cfg.prefix_depth,sh,local_geometry);}sh.local_parabola_hits+=local_geometry.parabola_hits;sh.local_pair_hits+=local_geometry.pair_hits;});for(auto&t:w)t.join();sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_hits=sh.states.hits;sh.report.state_cache_entries=sh.states.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;}
SearchReport run_mitm(const SearchConfig&cfg){Shared sh(cfg,false);std::deque<Task>queue;LocalGeometryCache prefix_geometry;std::unordered_set<Signature,SignatureHash>prefix_seen;
 for(const auto&mask:cfg.masks){if((mask.size()<5||mask.size()>8)||mask.back()!='L')throw std::runtime_error("MITM mask must have five to eight operations and end in L");const int shared_depth=static_cast<int>(mask.size())-5;std::vector<S>layer{initial(cfg.search_degrees,cfg.angle_at_parabola_vertex)};for(int depth=0;depth<shared_depth&&!expired(sh);++depth){std::vector<S>next;for(const auto&state:layer){if(expired(sh))break;for(int i=0;i<(int)state.p.size()&&!expired(sh);++i)for(int j=0;j<(int)state.p.size();++j){char op=mask[depth];if(i==j||(op=='L'&&j<i))continue;auto child=add(state,op,i,j,cfg.max_points,&sh.geometry,&prefix_geometry);++sh.generated;if(!child)continue;auto sig=state_signature(*child,mask.substr(depth+1));if(prefix_seen.insert(sig).second)next.push_back(std::move(*child));else ++sh.duplicates;}}if(next.size()>cfg.beam_width){for(auto&state:next)state.beam_score=beam_score(state,cfg);auto by_score=[](const S&a,const S&b){return a.beam_score<b.beam_score;};std::nth_element(next.begin(),next.begin()+cfg.beam_width,next.end(),by_score);next.resize(cfg.beam_width);}layer=std::move(next);}for(auto&state:layer)queue.push_back({std::move(state),mask});}
 sh.report.prefixes=queue.size();sh.local_parabola_hits+=prefix_geometry.parabola_hits;sh.local_pair_hits+=prefix_geometry.pair_hits;std::mutex queue_mutex;int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>workers;
 for(int worker=0;worker<nt;++worker)workers.emplace_back([&]{LocalGeometryCache geometry;while(!expired(sh)){Task task;{std::lock_guard<std::mutex>lock(queue_mutex);if(queue.empty())break;task=std::move(queue.front());queue.pop_front();}++sh.expanded;const std::size_t shared_depth=task.mask.size()-5;const auto left_mask=task.mask.substr(shared_depth,2),right_mask=task.mask.substr(shared_depth+2,2);ArmSet left,right;enumerate_arm(task.state,left_mask,0,left,sh,geometry);enumerate_arm(task.state,right_mask,0,right,sh,geometry);
   sh.mitm_records+=left.records.size()+right.records.size();std::unordered_map<Signature,std::vector<std::size_t>,SignatureHash>index;index.reserve(left.records.size()*2+1);for(std::size_t i=0;i<left.records.size();++i)index[left.records[i].ray].push_back(i);for(const auto&rrecord:right.records){auto found=index.find(rrecord.ray);if(found==index.end())continue;auto&rstate=right.states[rrecord.state_index];const auto&rcurve=rstate.c.back();for(auto record_index:found->second){auto&lrecord=left.records[record_index];auto&lstate=left.states[lrecord.state_index];const auto&lcurve=lstate.c.back();if(!same_branch(lrecord.branch,rrecord.branch)||!transverse(lcurve,rcurve,rrecord.branch))continue;++sh.mitm_matches;auto merged=merge_arms(task.state,lstate,rstate,sh,geometry);if(!merged)continue;auto goal=terminal(*merged,cfg.search_degrees,false,cfg.max_cost,cfg.angle_at_parabola_vertex,cfg.require_circle_parabola_use);if(!goal||goal->second>cfg.max_cost)continue;if(!dense_verify(*merged,false,goal->second,cfg)){++sh.dense_rejections;std::lock_guard<std::mutex>lock(sh.found_mutex);if(sh.report.sample_only_steps.empty()){sh.report.sample_only_steps=materialize_steps(*merged);if(goal->second>merged->cost)sh.report.sample_only_steps.push_back(std::string("line(")+(cfg.angle_at_parabola_vertex?"O":"F")+", common target-ray intersection)");sh.report.sample_only_schema=goal->first;sh.report.sample_only_mask=task.mask;}continue;}std::lock_guard<std::mutex>lock(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=goal->second;sh.report.steps=materialize_steps(*merged);if(goal->second>merged->cost)sh.report.steps.push_back(std::string("line(")+(cfg.angle_at_parabola_vertex?"O":"F")+", common target-ray intersection)");sh.report.terminal_schema="MITM two-curves-intersection+final-L";sh.report.mask=task.mask;sh.stop=true;}break;}if(sh.stop)break;}}
   sh.local_parabola_hits+=geometry.parabola_hits;sh.local_pair_hits+=geometry.pair_hits;});for(auto&worker:workers)worker.join();sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_hits=sh.states.hits;sh.report.state_cache_entries=sh.states.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.mitm_records=sh.mitm_records;sh.report.mitm_matches=sh.mitm_matches;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;}

SearchReport run_dynamic_e_mitm(const SearchConfig&input){
 SearchConfig exhaustive_cfg=input;
 if(input.exhaustive_cost4_e_mitm){exhaustive_cfg.max_cost=4;exhaustive_cfg.max_states=std::numeric_limits<std::size_t>::max();exhaustive_cfg.time_limit_seconds=1e12;}
 const SearchConfig&cfg=exhaustive_cfg;
 Shared sh(cfg,false);LocalGeometryCache prefix_geometry;
 const int left_depth=cfg.max_cost==4?1:2,right_depth=cfg.max_cost==4?1:(cfg.max_cost%2?1:2);
 const int shared_depth=std::max(0,cfg.max_cost-left_depth-right_depth);
 const bool fixed_schema=!cfg.masks.empty();
 const std::string schema=fixed_schema?cfg.masks.front():std::string{};
 if(fixed_schema&&(int)schema.size()!=cfg.max_cost)throw std::runtime_error("dynamic E MITM schema length must equal max cost");
 std::vector<S>layer{configured_initial(cfg.search_degrees,true,cfg.cost_policy)};
 for(int depth=0;depth<shared_depth&&!expired(sh);++depth){const std::vector<char>ops=fixed_schema?std::vector<char>{schema[depth]}:euclidea_mitm_operations();if(cfg.exhaustive_cost4_e_mitm){std::vector<S>next;std::unordered_set<Signature,SignatureHash>seen;for(const auto&state:layer)for(char op:ops)for(auto&child:expand_operation(state,op,cfg.max_points,&sh.geometry,&prefix_geometry)){++sh.generated;auto signature=state_signature(child,std::to_string(shared_depth-depth-1));if(seen.insert(signature).second)next.push_back(std::move(child));else ++sh.duplicates;}layer=std::move(next);}else{BoundedBeam next(cfg.beam_width);for(const auto&state:layer)for(char op:ops)for(auto&child:expand_operation(state,op,cfg.max_points,&sh.geometry,&prefix_geometry)){++sh.generated;if(expired(sh))break;auto signature=state_signature(child,std::to_string(shared_depth-depth-1));if(!next.insert(signature,std::move(child),cfg))++sh.duplicates;}layer=next.take();}}
 sh.report.prefixes=layer.size();sh.local_parabola_hits+=prefix_geometry.parabola_hits;sh.local_pair_hits+=prefix_geometry.pair_hits;std::deque<S>queue;for(auto&state:layer)queue.push_back(std::move(state));std::mutex queue_mutex;const int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>workers;
 for(int worker=0;worker<nt;++worker)workers.emplace_back([&]{LocalGeometryCache geometry;while(!expired(sh)){S base;{std::lock_guard<std::mutex>lock(queue_mutex);if(queue.empty())break;base=std::move(queue.front());queue.pop_front();}++sh.expanded;ArmSet left,right;if(fixed_schema){enumerate_schema_arm(base,schema.substr(shared_depth,left_depth),0,left,sh,geometry);enumerate_schema_arm(base,schema.substr(shared_depth+left_depth,right_depth),0,right,sh,geometry);}else{enumerate_dynamic_arm(base,left_depth,left,sh,geometry);enumerate_dynamic_arm(base,right_depth,right,sh,geometry);}sh.mitm_records+=left.records.size()+right.records.size();std::unordered_map<Signature,std::vector<std::size_t>,SignatureHash>index;for(std::size_t i=0;i<left.records.size();++i)index[left.records[i].ray].push_back(i);for(const auto&rrecord:right.records){auto found=index.find(rrecord.ray);if(found==index.end())continue;const auto&rstate=right.states[rrecord.state_index];const auto&rcurve=rstate.c.back();for(auto record_index:found->second){const auto&lrecord=left.records[record_index];const auto&lstate=left.states[lrecord.state_index];const auto&lcurve=lstate.c.back();if(!same_branch(lrecord.branch,rrecord.branch)||!transverse(lcurve,rcurve,rrecord.branch))continue;++sh.mitm_matches;auto merged=merge_arms(base,lstate,rstate,sh,geometry);if(!merged)continue;auto goal=terminal(*merged,cfg.search_degrees,false,cfg.max_cost,true,false,false,false,false,true);if(!goal)continue;if(!dense_verify(*merged,false,goal->second,cfg)){++sh.dense_rejections;continue;}std::lock_guard<std::mutex>lock(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=goal->second;sh.report.steps=materialize_steps(*merged);sh.report.terminal_schema="dynamic-MITM two-curves-intersection";sh.report.mask=fixed_schema?schema:"dynamic:L/C/B/N/P/A";sh.stop=true;}break;}if(sh.stop)break;}}sh.local_parabola_hits+=geometry.parabola_hits;sh.local_pair_hits+=geometry.pair_hits;});
 for(auto&worker:workers)worker.join();sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.mitm_records=sh.mitm_records;sh.report.mitm_matches=sh.mitm_matches;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;
}

SearchReport run_exhaustive_e_311(const SearchConfig&input){
 SearchConfig cfg=input;cfg.max_states=std::numeric_limits<std::size_t>::max();cfg.time_limit_seconds=1e12;cfg.max_cost=5;cfg.target_e_point=true;cfg.angle_at_parabola_vertex=true;cfg.cost_policy=engine::CostPolicy::Euclidea;
 Shared sh(cfg,false);std::atomic_int next_schema{0};constexpr int schema_count=216;const int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>workers;
 for(int worker=0;worker<nt;++worker)workers.emplace_back([&]{LocalGeometryCache geometry;for(;;){const int number=next_schema.fetch_add(1);if(number>=schema_count||sh.stop)break;int code=number;std::string schema(3,'L');for(int i=0;i<3;++i){schema[i]=euclidea_mitm_operations()[code%6];code/=6;}std::array<std::unordered_set<Signature,SignatureHash>,4>seen;S root=configured_initial(cfg.search_degrees,true,cfg.cost_policy);seen[0].insert(state_signature(root,schema));std::function<void(const S&,int)>visit=[&](const S&state,int depth){if(sh.stop)return;if(depth==3){++sh.expanded;ArmSet arms;enumerate_dynamic_arm(state,1,arms,sh,geometry);sh.mitm_records+=arms.records.size();for(std::size_t i=0;i<arms.records.size()&&!sh.stop;++i)for(std::size_t j=i+1;j<arms.records.size();++j){const auto&a=arms.records[i];const auto&b=arms.records[j];const auto&as=arms.states[a.state_index];const auto&bs=arms.states[b.state_index];const auto&ac=as.c.back();const auto&bc=bs.c.back();if(!same_branch(a.branch,b.branch)||!transverse(ac,bc,a.branch))continue;++sh.mitm_matches;auto merged=merge_arms(state,as,bs,sh,geometry);if(!merged)continue;auto goal=terminal(*merged,cfg.search_degrees,false,5,true,false,false,false,false,true);if(!goal)continue;if(!dense_verify(*merged,false,goal->second,cfg)){++sh.dense_rejections;continue;}std::lock_guard<std::mutex>lock(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=goal->second;sh.report.steps=materialize_steps(*merged);sh.report.terminal_schema="exhaustive-3+1+1 two-curves-intersection";sh.report.mask=schema;sh.stop=true;}break;}return;}for(auto&child:expand_operation(state,schema[depth],cfg.max_points,&sh.geometry,&geometry)){++sh.generated;auto signature=state_signature(child,schema.substr(depth+1));if(!seen[depth+1].insert(signature).second){++sh.duplicates;continue;}visit(child,depth+1);if(sh.stop)return;}};visit(root,0);}sh.local_parabola_hits+=geometry.parabola_hits;sh.local_pair_hits+=geometry.pair_hits;});
 for(auto&worker:workers)worker.join();sh.report.prefixes=sh.expanded;sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.mitm_records=sh.mitm_records;sh.report.mitm_matches=sh.mitm_matches;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;
}
}
SearchReport search(const SearchConfig&c){return run(c,false);} SearchReport search_beam(const SearchConfig&c){if(c.exhaustive_three_after_e)return run_exhaustive_three_after_e(c);if(c.exhaustive_cost4_e_sequential)return run_exhaustive_cost4_e_sequential(c);return run_beam(c);} SearchReport search_mitm(const SearchConfig&c){if(c.exhaustive_311)return run_exhaustive_e_311(c);if(c.target_e_point&&c.cost_policy==engine::CostPolicy::Euclidea)return run_dynamic_e_mitm(c);return run_mitm(c);} SearchReport search_tangent(const SearchConfig&c){return run(c,true);}
}
