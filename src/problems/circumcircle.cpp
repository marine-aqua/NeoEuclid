#include "neo/problems/circumcircle.hpp"

#include "neo/state.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace neo::circumcircle {
namespace {

using Clock = std::chrono::steady_clock;

Point circumcenter(Point a, Point b, Point c) {
    const double d=2*(a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y));
    const double aa=a.x*a.x+a.y*a.y,bb=b.x*b.x+b.y*b.y,cc=c.x*c.x+c.y*c.y;
    return {(aa*(b.y-c.y)+bb*(c.y-a.y)+cc*(a.y-b.y))/d,
            (aa*(c.x-b.x)+bb*(a.x-c.x)+cc*(b.x-a.x))/d,0};
}

bool contains(const Curve& curve, Point p) {
    if (const auto* l=std::get_if<Line>(&curve))
        return std::abs(l->a*p.x+l->b*p.y+l->c)<2e-7;
    const auto& c=std::get<Circle>(curve);
    return std::abs(std::hypot(p.x-c.x,p.y-c.y)-c.radius)<2e-7*std::max(1.,c.radius);
}

bool transverse(const Curve& a,const Curve& b,Point p) {
    if(const auto* x=std::get_if<Line>(&a)){
        if(const auto* y=std::get_if<Line>(&b))return std::abs(x->a*y->b-x->b*y->a)>1e-7;
        const auto& y=std::get<Circle>(b);return std::abs(x->a*(p.y-y.y)-x->b*(p.x-y.x))>1e-7;
    }
    if(const auto* y=std::get_if<Line>(&b)){const auto& x=std::get<Circle>(a);return std::abs(y->a*(p.y-x.y)-y->b*(p.x-x.x))>1e-7;}
    const auto& x=std::get<Circle>(a);const auto& y=std::get<Circle>(b);
    return std::abs((p.x-x.x)*(p.y-y.y)-(p.y-x.y)*(p.x-y.x))>1e-7;
}

bool target_circle(const Curve& curve,Point center,Point vertex){
    const auto* c=std::get_if<Circle>(&curve);if(!c)return false;
    return std::hypot(c->x-center.x,c->y-center.y)<2e-7&&std::abs(c->radius-std::hypot(vertex.x-center.x,vertex.y-center.y))<2e-7;
}

State initial(){
    State s;
    s.recipes={{OperationKind::Given,0,0,"A"},{OperationKind::Given,0,0,"B"},{OperationKind::Given,0,0,"C"}};
    s.points={{0.17,-0.31,0},{4.29,0.83,1},{1.11,3.76,2}};
    for(auto pair:std::vector<std::pair<int,int>>{{0,1},{1,2},{2,0}}){
        auto id=static_cast<std::uint32_t>(s.recipes.size());
        s.recipes.push_back({OperationKind::Given,0,0,pair.first==0&&pair.second==1?"AB-line":(pair.first==1?"BC-line":"CA-line")});
        auto l=*line_through(s.points[pair.first],s.points[pair.second],id);s.curves.push_back(l);
    }
    return s;
}

struct Record{Curve curve;int cost{};std::string recipe;};

} // namespace

Report search_mitm(double seconds,std::size_t max_states,std::size_t max_points){
    const auto start=Clock::now();Report report;State root=initial();Point goal=circumcenter(root.points[0],root.points[1],root.points[2]);
    ExpansionOptions options{max_points,true,8};std::unordered_map<std::string,int>seen;seen[state_key(root,7)]=0;
    std::vector<Record> records;std::unordered_map<std::string,int>best_curve;
    constexpr std::size_t seen_limit=1'000'000;
    constexpr std::size_t record_limit=250'000;
    auto expired=[&]{return std::chrono::duration<double>(Clock::now()-start).count()>=seconds||report.generated>=max_states;};
    std::function<void(const State&)> dfs;
    dfs=[&](const State& state){
        if(expired()||report.found||state.cost>=5)return;
        auto candidates=candidate_curves(state,true);
        std::stable_sort(candidates.begin(),candidates.end(),[](const Candidate&a,const Candidate&b){return a.cost>b.cost;});
        for(const auto& candidate:candidates){
            if(state.cost+candidate.cost>5)continue;++report.generated;
            auto child=add_curve(state,candidate,options);if(!child)continue;
            const auto& made=child->curves.back();
            if(target_circle(made,goal,root.points[0])){report.found=true;report.cost=child->cost;report.left=describe_recipe(*child,std::visit([](auto x){return x.recipe;},made));return;}
            if(contains(made,goal)){
                const auto recipe=describe_recipe(*child,std::visit([](auto x){return x.recipe;},made));
                const auto ck=std::to_string(curve_key(made,7).kind)+":"+std::to_string(curve_key(made,7).first)+","+std::to_string(curve_key(made,7).second)+","+std::to_string(curve_key(made,7).third);
                auto known=best_curve.find(ck);
                if(known==best_curve.end()||child->cost<records[known->second].cost){
                    Record current{made,child->cost,recipe};
                    for(const auto& other:records)if(current.cost+other.cost+1<=6&&transverse(current.curve,other.curve,goal)){
                        report.found=true;report.cost=current.cost+other.cost+1;report.left=other.recipe;report.right=current.recipe;break;
                    }
                    if(report.found)return;
                    if(records.size()<record_limit){best_curve[ck]=static_cast<int>(records.size());records.push_back(std::move(current));}
                }
            }
            auto sk=state_key(*child,7);auto it=seen.find(sk);
            if(it!=seen.end()&&it->second<=child->cost)continue;
            if(seen.size()<seen_limit)seen[sk]=child->cost;
            dfs(*child);
            if(expired()||report.found)return;
        }
    };
    dfs(root);
    report.states=seen.size();report.target_curves=records.size();report.elapsed_seconds=std::chrono::duration<double>(Clock::now()-start).count();return report;
}

} // namespace neo::circumcircle
