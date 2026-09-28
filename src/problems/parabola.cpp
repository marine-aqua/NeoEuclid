#include "neo/problems/parabola.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <complex>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
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
struct P{std::vector<V> v;std::string r;PK kind{PK::Given};int ca{-1},cb{-1},branch{};};
struct C{CK kind;std::vector<T> v;std::string r;char op{'G'};int pa{-1},pb{-1},pc{-1};std::uint64_t hash{};int branch{};};
struct Signature {std::uint64_t first{},second{};bool operator==(const Signature&o)const{return first==o.first&&second==o.second;}};
struct SignatureHash {std::size_t operator()(const Signature&s)const{return static_cast<std::size_t>(s.first^(s.second+0x9e3779b97f4a7c15ULL+(s.first<<6)+(s.first>>2)));}};
struct S{std::vector<P> p;std::vector<C> c;std::vector<std::string> steps;int cost{};std::string mask;std::uint64_t hash_xor{},hash_sum{};bool used_circle_parabola_point{},used_bisector{};};
struct Task{S state;std::string mask;};
struct SetShard { std::mutex mutex; std::unordered_set<Signature,SignatureHash> values; };
struct GeometryShard {
    std::mutex mutex;
    std::unordered_map<std::uint64_t,std::vector<std::vector<V>>> values;
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
    std::unordered_map<std::uint64_t,std::vector<std::vector<V>>> parabola,pair;
    std::size_t capacity{4096};
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
Signature state_signature(const S&s,const std::string&suffix={}){auto h=fingerprint(suffix);if(s.used_circle_parabola_point)h^=0xd6e8feb86659fd93ULL;if(s.used_bisector)h^=0xa0761d6478bd642fULL;return{mix64(s.hash_xor^h),mix64(s.hash_sum+h*0x9e3779b97f4a7c15ULL)};}
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
std::vector<std::vector<V>> cached_parabola(const C& curve,GeometryCache* cache,LocalGeometryCache* local){
 const auto hash=curve.hash?curve.hash:fingerprint(key(curve));
 if(local){auto it=local->parabola.find(hash);if(it!=local->parabola.end()){++local->parabola_hits;return it->second;}}
 if(cache&&cache->per_shard){auto& shard=*cache->parabola_shards[hash%cache->parabola_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);auto it=shard.values.find(hash);if(it!=shard.values.end()){++cache->parabola_hits;auto result=it->second;if(local){if(local->parabola.size()>=local->capacity)local->parabola.clear();local->parabola.emplace(hash,result);}return result;}}
 std::vector<std::vector<V>> roots;for(std::size_t i=0;i<curve.v.size();++i)roots.push_back(parabola(curve,i));auto result=stable(roots);
 if(cache&&cache->per_shard){auto& shard=*cache->parabola_shards[hash%cache->parabola_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);if(shard.values.size()<cache->per_shard){auto [it,inserted]=shard.values.emplace(hash,result);if(inserted)++cache->entries;else result=it->second;}}
 if(local){if(local->parabola.size()>=local->capacity)local->parabola.clear();local->parabola.emplace(hash,result);}
 return result;
}
std::vector<std::vector<V>> cached_pair(const C& first,const C& second,GeometryCache* cache,LocalGeometryCache* local){
 auto a=first.hash?first.hash:fingerprint(key(first)),b=second.hash?second.hash:fingerprint(key(second));if(b<a)std::swap(a,b);const auto hash=mix64(a^(mix64(b)+0x9e3779b97f4a7c15ULL));
 if(local){auto it=local->pair.find(hash);if(it!=local->pair.end()){++local->pair_hits;return it->second;}}
 if(cache&&cache->per_shard){auto& shard=*cache->pair_shards[hash%cache->pair_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);auto it=shard.values.find(hash);if(it!=shard.values.end()){++cache->pair_hits;auto result=it->second;if(local){if(local->pair.size()>=local->capacity)local->pair.clear();local->pair.emplace(hash,result);}return result;}}
 std::vector<std::vector<V>>roots;for(std::size_t i=0;i<first.v.size();++i)roots.push_back(meet(first.v[i],second.v[i],first.kind,second.kind));auto result=stable(roots);
 if(cache&&cache->per_shard){auto& shard=*cache->pair_shards[hash%cache->pair_shards.size()];std::lock_guard<std::mutex>lock(shard.mutex);if(shard.values.size()<cache->per_shard){auto [it,inserted]=shard.values.emplace(hash,result);if(inserted)++cache->entries;else result=it->second;}}
 if(local){if(local->pair.size()>=local->capacity)local->pair.clear();local->pair.emplace(hash,result);}
 return result;
}
S initial(const std::vector<double>&deg,bool vertex_mode=false){S s;auto vals=[&](V x){return std::vector<V>(deg.size(),x);};s.p.push_back({vals({0,0}),"O"});s.p.push_back({{},"P"});s.p.push_back({vals({1,0}),"F"});s.p.push_back({vals({-1,0}),"D[directrix-foot]"});for(double d:deg){double u=vertex_mode?1/std::tan(d*pi/180):1/std::tan(d*pi/360);s.p[1].v.push_back(vertex_mode?V{4*u*u,4*u}:V{u*u,2*u});}
 auto given=[&](T t,std::string r){s.c.push_back({CK::Line,std::vector<T>(deg.size(),t),std::move(r)});};given({0,1,0},"x-axis");given({1,0,0},"y-axis");given({1,0,1},"directrix");C side{CK::Line,{},"angle-side"};V vertex=vertex_mode?V{0,0}:V{1,0};for(auto p:s.p[1].v)side.v.push_back(*line(vertex,p));s.c.push_back(std::move(side));for(auto&c:s.c){c.hash=fingerprint(key(c));s.hash_xor^=mix64(c.hash);s.hash_sum+=mix64(c.hash+0x9e3779b97f4a7c15ULL);}return s;}
int operation_cost(char op){return op=='L'||op=='C'?1:(op=='N'||op=='B'?3:(op=='T'?5:4));}
bool is_circle_parabola_point(const S&s,int index){return index>=0&&index<(int)s.p.size()&&s.p[index].kind==PK::Parabola&&s.p[index].ca>=0&&s.p[index].ca<(int)s.c.size()&&s.c[s.p[index].ca].kind==CK::Circle;}
bool operation_uses_p(const S&s,char op,int ia,int ib,int ic=-1){if(op=='A')return false;if(op=='N'||op=='P')return ia>=0&&ia<(int)s.p.size()&&s.p[ia].r=="P";return (ia>=0&&ia<(int)s.p.size()&&s.p[ia].r=="P")||(ib>=0&&ib<(int)s.p.size()&&s.p[ib].r=="P")||(ic>=0&&ic<(int)s.p.size()&&s.p[ic].r=="P");}
std::optional<S> finish_add(const S&s,C n,int cost,std::size_t maxp,GeometryCache* cache,LocalGeometryCache* local){
 auto nk=key(n);n.hash=fingerprint(nk);for(auto&c:s.c)if(c.hash==n.hash&&key(c)==nk)return{};S z=s;int ni=(int)z.c.size();std::unordered_set<std::string>known;for(auto&p:z.p)known.insert(key(p));
 auto absorb=[&](std::vector<std::vector<V>> branches,PK pk,int other){int branch=0;for(auto&v:branches){P p{v,"",pk,ni,other,branch};p.r=(pk==PK::Parabola?"intersect("+n.r+", parabola)":"intersect("+n.r+", "+z.c[other].r+")")+"#"+std::to_string(branch);auto k=key(p);if(known.insert(k).second)z.p.push_back(std::move(p));++branch;}};
 absorb(cached_parabola(n,cache,local),PK::Parabola,-1);
 for(int j=0;j<(int)s.c.size();++j){absorb(cached_pair(n,s.c[j],cache,local),PK::Pair,j);if(z.p.size()>maxp)return{};}
 bool used=s.used_circle_parabola_point;if(n.op=='L'||n.op=='C'||n.op=='B')used=used||is_circle_parabola_point(s,n.pa)||is_circle_parabola_point(s,n.pb);else if(n.op=='N'||n.op=='P')used=used||is_circle_parabola_point(s,n.pa);z.used_circle_parabola_point=used;z.used_bisector=s.used_bisector||n.op=='B';z.hash_xor^=mix64(n.hash);z.hash_sum+=mix64(n.hash+0x9e3779b97f4a7c15ULL);char op=n.op;z.c.push_back(std::move(n));z.steps.push_back(z.c.back().r);z.cost+=cost;z.mask.push_back(op);return z;
}
std::optional<S> add(const S&s,char op,int ia,int ib,std::size_t maxp,GeometryCache* cache=nullptr,LocalGeometryCache* local=nullptr,int macro_branch=0,int ic=-1){
 C n;n.op=op;n.pa=ia;n.pb=ib;n.pc=ic;n.kind=(op=='C'||op=='T')?CK::Circle:CK::Line;n.branch=macro_branch;
 if(op=='L'||op=='C'){n.r=(op=='L'?"line(":"circle(")+s.p[ia].r+(op=='L'?", ":"; ")+s.p[ib].r+")";for(std::size_t k=0;k<s.p[ia].v.size();++k){auto a=s.p[ia].v[k],b=s.p[ib].v[k];if(op=='L'){auto l=line(a,b);if(!l)return{};n.v.push_back(*l);}else{double r=dist(a,b);if(r<eps||!std::isfinite(r))return{};n.v.push_back({a.x,a.y,r});}}}
 else if(op=='T'){if(ic<0||ic>=(int)s.p.size())return{};n.r="transfer-circle("+s.p[ia].r+", "+s.p[ib].r+"; "+s.p[ic].r+")[cost=5]";for(std::size_t k=0;k<s.p[ia].v.size();++k){double radius=dist(s.p[ia].v[k],s.p[ib].v[k]);auto center=s.p[ic].v[k];if(radius<eps||!std::isfinite(radius))return{};n.v.push_back({center.x,center.y,radius});}}
 else if(op=='B'){n.r="perpendicular-bisector("+s.p[ia].r+", "+s.p[ib].r+")";for(std::size_t k=0;k<s.p[ia].v.size();++k){auto a=s.p[ia].v[k],b=s.p[ib].v[k];double x=b.x-a.x,y=b.y-a.y,norm=std::hypot(x,y);if(norm<eps)return{};x/=norm;y/=norm;double c=-(x*(a.x+b.x)/2+y*(a.y+b.y)/2);if(x<0||(std::abs(x)<eps&&y<0)){x=-x;y=-y;c=-c;}n.v.push_back({x,y,c});}}
 else if(op=='N'||op=='P'){if(ib<0||ib>=(int)s.c.size()||s.c[ib].kind!=CK::Line)return{};n.r=std::string(op=='N'?"perpendicular-through(":"parallel-through(")+s.p[ia].r+", "+s.c[ib].r+")";for(std::size_t k=0;k<s.p[ia].v.size();++k){V p=s.p[ia].v[k];T l=s.c[ib].v[k];double a=op=='N'?-l.b:l.a,b=op=='N'?l.a:l.b,c=-(a*p.x+b*p.y);if(a<0||(std::abs(a)<eps&&b<0)){a=-a;b=-b;c=-c;}n.v.push_back({a,b,c});}}
 else if(op=='A'){if(ia<0||ib<0||ia>=(int)s.c.size()||ib>=(int)s.c.size()||s.c[ia].kind!=CK::Line||s.c[ib].kind!=CK::Line)return{};n.r="angle-bisector("+s.c[ia].r+", "+s.c[ib].r+")#"+std::to_string(macro_branch);double sign=macro_branch==0?1.:-1.;for(std::size_t k=0;k<s.c[ia].v.size();++k){auto x=s.c[ia].v[k],y=s.c[ib].v[k];double a=x.a+sign*y.a,b=x.b+sign*y.b,c=x.c+sign*y.c,norm=std::hypot(a,b);if(norm<eps)return{};a/=norm;b/=norm;c/=norm;if(a<0||(std::abs(a)<eps&&b<0)){a=-a;b=-b;c=-c;}n.v.push_back({a,b,c});}}
 else return{};
 return finish_add(s,std::move(n),operation_cost(op),maxp,cache,local);
}
std::optional<S> k2_initial(const std::vector<double>&degrees,std::size_t maxp,GeometryCache*cache,LocalGeometryCache*local){S s=initial(degrees,true);auto first=add(s,'C',3,0,maxp,cache,local);if(!first)return{};s=std::move(*first);int a2=-1;for(int i=0;i<(int)s.p.size();++i){bool ok=true;for(auto p:s.p[i].v)if(std::abs(p.x+2.)>2e-6||std::abs(p.y)>2e-6){ok=false;break;}if(ok){a2=i;break;}}if(a2<0)return{};auto second=add(s,'L',1,a2,maxp,cache,local);if(!second)return{};return second;}
std::vector<S> expand_operation(const S&s,char op,std::size_t maxp,GeometryCache*cache,LocalGeometryCache*local,bool require_first_p=false){std::vector<S>out;auto push=[&](int a,int b,int branch=0,int c=-1){if(require_first_p&&s.steps.empty()&&!operation_uses_p(s,op,a,b,c))return;auto child=add(s,op,a,b,maxp,cache,local,branch,c);if(child)out.push_back(std::move(*child));};if(op=='L'||op=='B'){for(int i=0;i<(int)s.p.size();++i)for(int j=i+1;j<(int)s.p.size();++j)push(i,j);}else if(op=='C'){for(int i=0;i<(int)s.p.size();++i)for(int j=0;j<(int)s.p.size();++j)if(i!=j)push(i,j);}else if(op=='T'){for(int i=0;i<(int)s.p.size();++i)for(int j=i+1;j<(int)s.p.size();++j)for(int k=0;k<(int)s.p.size();++k)if(k!=i&&k!=j)push(i,j,0,k);}else if(op=='N'||op=='P'){for(int i=0;i<(int)s.p.size();++i)for(int j=0;j<(int)s.c.size();++j)if(s.c[j].kind==CK::Line)push(i,j);}else if(op=='A'){for(int i=0;i<(int)s.c.size();++i)if(s.c[i].kind==CK::Line)for(int j=i+1;j<(int)s.c.size();++j)if(s.c[j].kind==CK::Line){push(i,j,0);push(i,j,1);}}return out;}
double target_angle(double d,bool tangent){return (tangent?d/2:d/3)*pi/180;}
V goal_vertex(const S&s,std::size_t sample,bool tangent,bool vertex_mode){return tangent?s.p[1].v[sample]:(vertex_mode?V{0,0}:V{1,0});}
bool target_line(const S&s,const C&c,bool tangent,bool vertex_mode){if(c.kind!=CK::Line)return false;for(std::size_t i=0;i<c.v.size();++i){V v=goal_vertex(s,i,tangent,vertex_mode);double t=target_angle((i==0?0:0),tangent);(void)t;auto l=c.v[i]; // direction test is supplied separately below
 if(std::abs(l.a*v.x+l.b*v.y+l.c)>2e-6)return false;}return true;}
bool line_goal(const S&s,const C&c,const std::vector<double>&deg,bool tangent,bool vertex_mode){if(!target_line(s,c,tangent,vertex_mode))return false;for(std::size_t i=0;i<deg.size();++i){double t=target_angle(deg[i],tangent);auto l=c.v[i];if(std::abs(l.a*std::cos(t)+l.b*std::sin(t))>2e-6)return false;}return true;}
bool point_goal(const S&s,const P&p,const std::vector<double>&deg,bool tangent,bool vertex_mode){for(std::size_t i=0;i<deg.size();++i){V v=goal_vertex(s,i,tangent,vertex_mode);double t=target_angle(deg[i],tangent),dx=p.v[i].x-v.x,dy=p.v[i].y-v.y,radius=std::hypot(dx,dy),forward=dx*std::cos(t)+dy*std::sin(t),cross=-std::sin(t)*dx+std::cos(t)*dy;if(radius<=1e-5||std::abs(forward)<=1e-5||std::abs(cross)>2e-7*std::max(1.,radius))return false;}return true;}
bool center_goal(const S&s,const P&p){for(std::size_t i=0;i<p.v.size();++i){const auto input=s.p[1].v[i];V wanted{8.+4.5*input.x,4.*input.y};if(dist(p.v[i],wanted)>3e-6*std::max(1.,std::hypot(wanted.x,wanted.y)))return false;}return true;}
bool k2_circle_goal(const S&s,const C&curve){if(curve.kind!=CK::Circle)return false;for(std::size_t i=0;i<curve.v.size();++i){double t=s.p[1].v[i].y/2.;if(std::abs(t)<eps)return false;double qv=2./t,h=3.5+4.5*qv*qv,k=2.*qv,r=std::sqrt(h*h+qv*qv);auto actual=curve.v[i];double scale=std::max({1.,std::abs(h),std::abs(k),r});if(std::abs(actual.a-h)>3e-6*scale||std::abs(actual.b-k)>3e-6*scale||std::abs(actual.c-r)>3e-6*scale)return false;}return true;}
bool input_dependent(const std::string&recipe){return recipe.find('P')!=std::string::npos||recipe.find("angle-side")!=std::string::npos;}
std::optional<std::pair<std::string,int>> terminal(const S&s,const std::vector<double>&d,bool tangent,int maxcost,bool vertex_mode=false,bool require_circle_parabola=false,bool target_center=false,bool require_bisector=false,bool target_k2=false){if(require_bisector&&!s.used_bisector)return{};if(target_k2){for(const auto&c:s.c)if(input_dependent(c.r)&&k2_circle_goal(s,c))return{{"k2-triple-angle-circle",s.cost}};return{};}if(target_center){for(const auto&p:s.p)if(input_dependent(p.r)&&center_goal(s,p))return{{"constructed-circle-center",s.cost}};return{};}for(auto&c:s.c)if((!require_circle_parabola||s.used_circle_parabola_point)&&input_dependent(c.r)&&line_goal(s,c,d,tangent,vertex_mode))return{{"direct-target-line",s.cost}};if(s.cost+1<=maxcost)for(int i=0;i<(int)s.p.size();++i){const auto&p=s.p[i];if((!require_circle_parabola||s.used_circle_parabola_point||is_circle_parabola_point(s,i))&&input_dependent(p.r)&&point_goal(s,p,d,tangent,vertex_mode))return{{p.kind==PK::Parabola?"curve-parabola-intersection+join":(p.kind==PK::Pair?"two-curves-intersection+join":"existing-point+join"),s.cost+1}};}return{};}
std::vector<std::string> masks(int n,const std::vector<std::string>&requested){if(!requested.empty())return requested;std::vector<std::string>o;for(int bits=0;bits<(1<<n);++bits){std::string m;for(int i=0;i<n;++i)m.push_back(bits&(1<<i)?'C':'L');o.push_back(m);}return o;}
bool expired(Shared&x){return x.stop||std::chrono::duration<double>(std::chrono::steady_clock::now()-x.start).count()>=x.cfg.time_limit_seconds||x.generated>=x.cfg.max_states;}
std::vector<std::vector<V>> target_ray_branches(const C&curve,const std::vector<double>&degrees,bool vertex_mode){std::vector<std::vector<V>>roots;V v=vertex_mode?V{0,0}:V{1,0};for(std::size_t i=0;i<degrees.size();++i){double t=target_angle(degrees[i],false);T ray=*line(v,{v.x+std::cos(t),v.y+std::sin(t)});std::vector<V>forward;for(auto p:meet(curve.v[i],ray,curve.kind,CK::Line))if((p.x-v.x)*std::cos(t)+(p.y-v.y)*std::sin(t)>1e-5)forward.push_back(p);roots.push_back(std::move(forward));}return stable(roots);}
std::string ray_key(const std::vector<V>&branch,bool vertex_mode){std::string out;double vx=vertex_mode?0.:1.;for(auto p:branch){double d=std::hypot(p.x-vx,p.y);out+=std::to_string(std::llround((d/(1+d))*1e6))+',';}return out;}
bool same_branch(const std::vector<V>&a,const std::vector<V>&b){if(a.size()!=b.size())return false;for(std::size_t i=0;i<a.size();++i)if(dist(a[i],b[i])>2e-5*std::max({1.,dist(a[i],{1,0}),dist(b[i],{1,0})}))return false;return true;}
bool transverse(const C&a,const C&b,const std::vector<V>&points){if(a.kind==b.kind&&a.hash==b.hash)return false;for(std::size_t i=0;i<points.size();++i){V p=points[i];T x=a.v[i],y=b.v[i];double measure=0;if(a.kind==CK::Line&&b.kind==CK::Line)measure=std::abs(x.a*y.b-x.b*y.a);else if(a.kind==CK::Line)measure=std::abs(x.a*(p.y-y.b)-x.b*(p.x-y.a))/std::max(1.,y.c);else if(b.kind==CK::Line)measure=std::abs(y.a*(p.y-x.b)-y.b*(p.x-x.a))/std::max(1.,x.c);else measure=std::abs((p.x-x.a)*(p.y-y.b)-(p.y-x.b)*(p.x-y.a))/std::max(1.,x.c*y.c);if(measure<2e-6)return false;}return true;}
void enumerate_arm(const S&base,const std::string&arm,int pos,std::vector<S>&out,Shared&sh,LocalGeometryCache&geometry,std::unordered_set<Signature,SignatureHash>&seen){if(expired(sh))return;if(pos==(int)arm.size()){out.push_back(base);return;}char op=arm[pos];for(int i=0;i<(int)base.p.size();++i)for(int j=0;j<(int)base.p.size();++j){if(i==j||(op=='L'&&j<i))continue;auto child=add(base,op,i,j,sh.cfg.max_points,&sh.geometry,&geometry);++sh.generated;if(!child)continue;auto sig=state_signature(*child,arm.substr(pos+1));if(!seen.insert(sig).second){++sh.duplicates;continue;}enumerate_arm(*child,arm,pos+1,out,sh,geometry,seen);if(expired(sh))return;}}
std::optional<S> merge_arms(const S&base,const S&left,const S&right,Shared&sh,LocalGeometryCache&geometry){S merged=base;auto append=[&](const S&source){for(std::size_t ci=base.c.size();ci<source.c.size();++ci){const auto&curve=source.c[ci];int pa=-1,pb=-1;const auto ka=key(source.p[curve.pa]),kb=key(source.p[curve.pb]);for(int i=0;i<(int)merged.p.size();++i){auto candidate=key(merged.p[i]);if(candidate==ka)pa=i;if(candidate==kb)pb=i;}if(pa<0||pb<0)return false;bool exists=false;for(const auto&known:merged.c)if(known.hash==curve.hash&&key(known)==key(curve)){exists=true;break;}if(exists)continue;auto next=add(merged,curve.op,pa,pb,sh.cfg.max_points,&sh.geometry,&geometry);if(!next)return false;merged=std::move(*next);}return true;};if(!append(left)||!append(right))return{};return merged;}
// Replay on a dense grid by re-running the exact operation choices. Branches are
// reconstructed geometrically; a candidate must independently satisfy every dense sample.
bool dense_verify(const S&hit,bool tangent,int final_cost,const SearchConfig&cfg){
 std::vector<double>d;double high=cfg.validation_max_degrees>0?cfg.validation_max_degrees:(cfg.angle_at_parabola_vertex?85.:175.);for(int i=0;i<cfg.validation_samples;++i)d.push_back(5.+(high-5.)*i/(cfg.validation_samples-1));S s=initial(d,cfg.angle_at_parabola_vertex);
 for(std::size_t step=0;step<hit.steps.size();++step){const C&hc=hit.c[4+step];int ia=-1,ib=-1,ic=-1;if(hc.op=='L'||hc.op=='C'||hc.op=='B'||hc.op=='T'){for(int i=0;i<(int)s.p.size();++i){if(s.p[i].r==hit.p[hc.pa].r)ia=i;if(s.p[i].r==hit.p[hc.pb].r)ib=i;if(hc.op=='T'&&s.p[i].r==hit.p[hc.pc].r)ic=i;}}else if(hc.op=='N'||hc.op=='P'){for(int i=0;i<(int)s.p.size();++i)if(s.p[i].r==hit.p[hc.pa].r)ia=i;for(int i=0;i<(int)s.c.size();++i)if(s.c[i].r==hit.c[hc.pb].r)ib=i;}else if(hc.op=='A'){for(int i=0;i<(int)s.c.size();++i){if(s.c[i].r==hit.c[hc.pa].r)ia=i;if(s.c[i].r==hit.c[hc.pb].r)ib=i;}}
  if(ia<0||ib<0||(hc.op=='T'&&ic<0))return false;auto z=add(s,hc.op,ia,ib,cfg.max_points*4,nullptr,nullptr,hc.branch,ic);if(!z)return false;s=std::move(*z);
 }
 return terminal(s,d,tangent,final_cost,cfg.angle_at_parabola_vertex,cfg.require_circle_parabola_use,cfg.target_circle_center,cfg.require_bisector_use,cfg.target_k2_circle).has_value();
}
void dfs(S s,const std::string&mask,int pos,Shared&sh,LocalGeometryCache&geometry){if(expired(sh))return;sh.expanded++;if(auto t=terminal(s,sh.cfg.search_degrees,sh.tangent,sh.cfg.max_cost,sh.cfg.angle_at_parabola_vertex,sh.cfg.require_circle_parabola_use,sh.cfg.target_circle_center,sh.cfg.require_bisector_use,sh.cfg.target_k2_circle)){if(dense_verify(s,sh.tangent,t->second,sh.cfg)){std::lock_guard<std::mutex>g(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=t->second;sh.report.steps=s.steps;if(t->second>s.cost)sh.report.steps.push_back(std::string("line(")+(sh.tangent?"P":(sh.cfg.angle_at_parabola_vertex?"O":"F")) +", target-point)");sh.report.terminal_schema=t->first;sh.report.mask=s.mask;sh.stop=true;}return;}sh.dense_rejections++;}
 if(pos>=(int)mask.size()||s.cost>=sh.cfg.max_cost)return;char op=mask[pos];if(s.cost+operation_cost(op)>sh.cfg.max_cost)return;std::unordered_set<Signature,SignatureHash>local;for(auto&z:expand_operation(s,op,sh.cfg.max_points,&sh.geometry,&geometry,sh.cfg.require_first_step_uses_p)){sh.generated++;auto signature=state_signature(z,mask.substr(pos+1));if(!local.insert(signature).second){sh.duplicates++;continue;}if(!sh.states.insert(signature)){sh.duplicates++;continue;}dfs(std::move(z),mask,pos+1,sh,geometry);if(expired(sh))return;}}

double beam_score(const S& state, const SearchConfig& cfg) {
 double best=1e30;
 for(const auto& point:state.p){if(!input_dependent(point.r))continue;double worst=0.;bool valid=true;for(std::size_t i=0;i<cfg.search_degrees.size();++i){const V vertex=cfg.angle_at_parabola_vertex?V{0,0}:V{1,0};const double angle=target_angle(cfg.search_degrees[i],false),dx=point.v[i].x-vertex.x,dy=point.v[i].y-vertex.y,r=std::hypot(dx,dy),along=dx*std::cos(angle)+dy*std::sin(angle),cross=std::abs(-std::sin(angle)*dx+std::cos(angle)*dy)/std::max(1.,r);if(!std::isfinite(cross)||r<1e-8){valid=false;break;}worst=std::max(worst,cross+(along<=0?1.:0.));}if(valid)best=std::min(best,worst);}
 for(const auto& curve:state.c){if(curve.kind!=CK::Line||!input_dependent(curve.r))continue;double worst=0.;for(std::size_t i=0;i<cfg.search_degrees.size();++i){const V vertex=cfg.angle_at_parabola_vertex?V{0,0}:V{1,0};const double angle=target_angle(cfg.search_degrees[i],false);const auto line=curve.v[i];worst=std::max(worst,std::abs(line.a*vertex.x+line.b*vertex.y+line.c)+std::abs(line.a*std::cos(angle)+line.b*std::sin(angle)));}best=std::min(best,worst);}
 return best+1e-6*state.p.size();
}

std::vector<char> beam_operations(const S& state,const SearchConfig& cfg){if(cfg.masks.empty())return{'L','C'};std::vector<char> result;for(const auto& mask:cfg.masks)if(mask.size()>state.mask.size()&&mask.compare(0,state.mask.size(),state.mask)==0&&std::find(result.begin(),result.end(),mask[state.mask.size()])==result.end())result.push_back(mask[state.mask.size()]);return result;}

SearchReport run_beam(const SearchConfig& cfg){Shared sh(cfg,false);LocalGeometryCache geometry;std::vector<std::vector<S>>buckets(static_cast<std::size_t>(cfg.max_cost+1));buckets[0].push_back(initial(cfg.search_degrees,cfg.angle_at_parabola_vertex));
 for(int cost=0;cost<=cfg.max_cost&&!expired(sh);++cost){auto&states=buckets[static_cast<std::size_t>(cost)];std::sort(states.begin(),states.end(),[&](const S&a,const S&b){return beam_score(a,cfg)<beam_score(b,cfg);});if(states.size()>cfg.beam_width)states.resize(cfg.beam_width);std::vector<std::unordered_map<Signature,S,SignatureHash>>generated(static_cast<std::size_t>(cfg.max_cost+1));
  for(const auto&state:states){if(expired(sh))break;++sh.expanded;if(auto hit=terminal(state,cfg.search_degrees,false,cfg.max_cost,cfg.angle_at_parabola_vertex,cfg.require_circle_parabola_use,cfg.target_circle_center,cfg.require_bisector_use,cfg.target_k2_circle)){if(dense_verify(state,false,hit->second,cfg)){sh.report.found=sh.report.densely_verified=true;sh.report.cost=hit->second;sh.report.steps=state.steps;if(hit->second>state.cost)sh.report.steps.push_back(std::string("line(")+(cfg.angle_at_parabola_vertex?"O":"F")+", target-point)");sh.report.terminal_schema=hit->first;sh.report.mask=state.mask;sh.stop=true;break;}++sh.dense_rejections;}
   for(char op:beam_operations(state,cfg)){if(state.cost+operation_cost(op)>cfg.max_cost)continue;for(auto&child:expand_operation(state,op,cfg.max_points,&sh.geometry,&geometry,cfg.require_first_step_uses_p)){++sh.generated;if(expired(sh))break;auto signature=state_signature(child);auto&level=generated[static_cast<std::size_t>(child.cost)];if(level.emplace(signature,std::move(child)).second)continue;++sh.duplicates;}}}
  for(int next=cost+1;next<=cfg.max_cost;++next)for(auto&entry:generated[static_cast<std::size_t>(next)])buckets[static_cast<std::size_t>(next)].push_back(std::move(entry.second));}
 sh.local_parabola_hits=geometry.parabola_hits;sh.local_pair_hits=geometry.pair_hits;sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_entries=sh.states.entries;sh.report.state_cache_hits=sh.states.hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;}

SearchReport run(const SearchConfig&cfg,bool tangent){Shared sh(cfg,tangent);auto ms=masks(cfg.max_cost,cfg.masks);std::deque<Task>q;
 LocalGeometryCache prefix_geometry;std::unordered_set<Signature,SignatureHash>prefix_seen;auto configured_base=cfg.target_k2_circle?k2_initial(cfg.search_degrees,cfg.max_points,&sh.geometry,&prefix_geometry):std::optional<S>{initial(cfg.search_degrees,cfg.angle_at_parabola_vertex)};if(!configured_base)throw std::runtime_error("failed to build k=2 reciprocal prefix");for(auto&m:ms){std::vector<S>layer{*configured_base};for(int depth=0;depth<std::min(cfg.prefix_depth,(int)m.size());++depth){std::vector<S>next;for(auto&s:layer){if(s.cost+operation_cost(m[depth])>cfg.max_cost)continue;for(auto&z:expand_operation(s,m[depth],cfg.max_points,&sh.geometry,&prefix_geometry,cfg.require_first_step_uses_p)){sh.generated++;auto signature=state_signature(z,m.substr(depth+1));if(prefix_seen.insert(signature).second)next.push_back(std::move(z));else sh.duplicates++;}}layer=std::move(next);}for(auto&s:layer){sh.states.insert(state_signature(s,m.substr(cfg.prefix_depth)));q.push_back({std::move(s),m});}}
 sh.local_parabola_hits+=prefix_geometry.parabola_hits;sh.local_pair_hits+=prefix_geometry.pair_hits;sh.report.prefixes=q.size();std::mutex qm;int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>w;for(int t=0;t<nt;++t)w.emplace_back([&]{LocalGeometryCache local_geometry;while(!expired(sh)){Task x;{std::lock_guard<std::mutex>g(qm);if(q.empty())break;x=std::move(q.front());q.pop_front();}dfs(std::move(x.state),x.mask,cfg.prefix_depth,sh,local_geometry);}sh.local_parabola_hits+=local_geometry.parabola_hits;sh.local_pair_hits+=local_geometry.pair_hits;});for(auto&t:w)t.join();sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_hits=sh.states.hits;sh.report.state_cache_entries=sh.states.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;}
SearchReport run_mitm(const SearchConfig&cfg){Shared sh(cfg,false);std::deque<Task>queue;LocalGeometryCache prefix_geometry;std::unordered_set<Signature,SignatureHash>prefix_seen;
 for(const auto&mask:cfg.masks){if(mask.size()!=7||mask.back()!='L')throw std::runtime_error("MITM mask must have seven operations and end in L");std::vector<S>layer{initial(cfg.search_degrees,cfg.angle_at_parabola_vertex)};for(int depth=0;depth<2;++depth){std::vector<S>next;for(const auto&state:layer)for(int i=0;i<(int)state.p.size();++i)for(int j=0;j<(int)state.p.size();++j){char op=mask[depth];if(i==j||(op=='L'&&j<i))continue;auto child=add(state,op,i,j,cfg.max_points,&sh.geometry,&prefix_geometry);++sh.generated;if(!child)continue;auto sig=state_signature(*child,mask.substr(depth+1));if(prefix_seen.insert(sig).second)next.push_back(std::move(*child));else ++sh.duplicates;}layer=std::move(next);}for(auto&state:layer)queue.push_back({std::move(state),mask});}
 sh.report.prefixes=queue.size();sh.local_parabola_hits+=prefix_geometry.parabola_hits;sh.local_pair_hits+=prefix_geometry.pair_hits;std::mutex queue_mutex;int nt=cfg.threads>0?cfg.threads:(int)std::max(1u,std::thread::hardware_concurrency());std::vector<std::thread>workers;
 for(int worker=0;worker<nt;++worker)workers.emplace_back([&]{LocalGeometryCache geometry;while(!expired(sh)){Task task;{std::lock_guard<std::mutex>lock(queue_mutex);if(queue.empty())break;task=std::move(queue.front());queue.pop_front();}++sh.expanded;const auto left_mask=task.mask.substr(2,2),right_mask=task.mask.substr(4,2);std::vector<S>left,right;std::unordered_set<Signature,SignatureHash>left_seen,right_seen;enumerate_arm(task.state,left_mask,0,left,sh,geometry,left_seen);enumerate_arm(task.state,right_mask,0,right,sh,geometry,right_seen);sh.mitm_records+=left.size()+right.size();
   std::unordered_map<std::string,std::vector<std::pair<std::size_t,std::vector<V>>>>index;for(std::size_t i=0;i<left.size();++i){const auto&curve=left[i].c.back();for(auto branch:target_ray_branches(curve,cfg.search_degrees,cfg.angle_at_parabola_vertex))index[ray_key(branch,cfg.angle_at_parabola_vertex)].push_back({i,std::move(branch)});}for(auto&rstate:right){const auto&rcurve=rstate.c.back();for(auto rbranch:target_ray_branches(rcurve,cfg.search_degrees,cfg.angle_at_parabola_vertex)){auto found=index.find(ray_key(rbranch,cfg.angle_at_parabola_vertex));if(found==index.end())continue;for(auto&entry:found->second){auto&lstate=left[entry.first];const auto&lcurve=lstate.c.back();if(!same_branch(entry.second,rbranch)||!transverse(lcurve,rcurve,rbranch))continue;++sh.mitm_matches;auto merged=merge_arms(task.state,lstate,rstate,sh,geometry);if(!merged)continue;auto goal=terminal(*merged,cfg.search_degrees,false,7,cfg.angle_at_parabola_vertex,cfg.require_circle_parabola_use);if(!goal||goal->second>7)continue;if(!dense_verify(*merged,false,goal->second,cfg)){++sh.dense_rejections;continue;}std::lock_guard<std::mutex>lock(sh.found_mutex);if(!sh.stop){sh.report.found=sh.report.densely_verified=true;sh.report.cost=goal->second;sh.report.steps=merged->steps;if(goal->second>merged->cost)sh.report.steps.push_back(std::string("line(")+(cfg.angle_at_parabola_vertex?"O":"F")+", common target-ray intersection)");sh.report.terminal_schema="MITM two-curves-intersection+final-L";sh.report.mask=task.mask;sh.stop=true;}break;}if(sh.stop)break;}if(sh.stop)break;}}
   sh.local_parabola_hits+=geometry.parabola_hits;sh.local_pair_hits+=geometry.pair_hits;});for(auto&worker:workers)worker.join();sh.report.expanded=sh.expanded;sh.report.generated=sh.generated;sh.report.duplicates=sh.duplicates;sh.report.dense_rejections=sh.dense_rejections;sh.report.state_cache_hits=sh.states.hits;sh.report.state_cache_entries=sh.states.entries;sh.report.parabola_cache_hits=sh.geometry.parabola_hits;sh.report.pair_cache_hits=sh.geometry.pair_hits;sh.report.local_parabola_cache_hits=sh.local_parabola_hits;sh.report.local_pair_cache_hits=sh.local_pair_hits;sh.report.geometry_cache_entries=sh.geometry.entries;sh.report.mitm_records=sh.mitm_records;sh.report.mitm_matches=sh.mitm_matches;sh.report.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sh.start).count();return sh.report;}
}
SearchReport search(const SearchConfig&c){return run(c,false);} SearchReport search_beam(const SearchConfig&c){return run_beam(c);} SearchReport search_mitm(const SearchConfig&c){return run_mitm(c);} SearchReport search_tangent(const SearchConfig&c){return run(c,true);}
}
