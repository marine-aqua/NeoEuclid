#include "neo/problems/obtuse.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace neo::obtuse {
namespace {

constexpr double kEps = 1e-8;
constexpr double kPi = 3.14159265358979323846;

struct Vec2 { double x{}, y{}; };
struct Triple { double a{}, b{}, c{}; };
enum class CurveKind { Line, Circle };
enum class OpKind { Line, Circle, PerpendicularBisector, Perpendicular, Parallel, AngleBisector };

struct Point {
    std::vector<Vec2> values;
    std::string recipe;
};

struct Curve {
    CurveKind kind{CurveKind::Line};
    std::vector<Triple> values;
    std::string recipe;
    std::vector<std::string> input_points;
};

struct Candidate {
    Curve curve;
    OpKind kind{OpKind::Line};
    int cost{1};
};

struct State {
    std::vector<Point> points;
    std::vector<Curve> curves;
    std::vector<std::string> steps;
    std::vector<std::string> point_definitions;
    int cost{3};
    bool used_hyperbola_point{false};
    std::unordered_set<std::string> hyperbola_points;
};

double distance(Vec2 p, Vec2 q) { return std::hypot(p.x - q.x, p.y - q.y); }

std::optional<Triple> canonical_line(Vec2 p, Vec2 q) {
    double a = p.y - q.y, b = q.x - p.x;
    const double norm = std::hypot(a, b);
    if (norm < kEps) return std::nullopt;
    a /= norm; b /= norm;
    double c = -(a * p.x + b * p.y);
    if (a < -kEps || (std::abs(a) <= kEps && b < 0)) { a = -a; b = -b; c = -c; }
    return Triple{a, b, c};
}

std::string point_key(const Point& point, int digits = 7) {
    const double scale = std::pow(10.0, digits);
    std::ostringstream out;
    for (const auto& p : point.values)
        out << std::llround(p.x * scale) << ',' << std::llround(p.y * scale) << ';';
    return out.str();
}

std::string curve_key(const Curve& curve, int digits = 7) {
    const double scale = std::pow(10.0, digits);
    std::ostringstream out;
    out << (curve.kind == CurveKind::Line ? 'L' : 'C');
    for (const auto& v : curve.values)
        out << std::llround(v.a * scale) << ',' << std::llround(v.b * scale)
            << ',' << std::llround(v.c * scale) << ';';
    return out.str();
}

std::string state_key(const State& state) {
    std::vector<std::string> keys;
    keys.reserve(state.curves.size());
    for (const auto& curve : state.curves) keys.push_back(curve_key(curve, 6));
    std::sort(keys.begin(), keys.end());
    return std::accumulate(keys.begin(), keys.end(), std::string{});
}

std::optional<Curve> make_line(const Point& p, const Point& q) {
    Curve curve{CurveKind::Line, {}, "line(" + p.recipe + ", " + q.recipe + ")",
                {point_key(p), point_key(q)}};
    for (std::size_t i = 0; i < p.values.size(); ++i) {
        auto line = canonical_line(p.values[i], q.values[i]);
        if (!line) return std::nullopt;
        curve.values.push_back(*line);
    }
    return curve;
}

std::optional<Curve> make_circle(const Point& center, const Point& through) {
    Curve curve{CurveKind::Circle, {}, "circle(" + center.recipe + "; " + through.recipe + ")",
                {point_key(center), point_key(through)}};
    for (std::size_t i = 0; i < center.values.size(); ++i) {
        const double radius = distance(center.values[i], through.values[i]);
        if (radius < kEps || !std::isfinite(radius)) return std::nullopt;
        curve.values.push_back({center.values[i].x, center.values[i].y, radius});
    }
    return curve;
}

std::optional<Curve> perpendicular_bisector(const Point& p, const Point& q) {
    Curve curve{CurveKind::Line, {}, "perp-bisector(" + p.recipe + ", " + q.recipe + ")",
                {point_key(p), point_key(q)}};
    for (std::size_t i = 0; i < p.values.size(); ++i) {
        const auto pv = p.values[i], qv = q.values[i];
        double a = qv.x - pv.x, b = qv.y - pv.y;
        const double norm = std::hypot(a, b);
        if (norm < kEps) return std::nullopt;
        a /= norm; b /= norm;
        double c = -(a * (pv.x + qv.x) / 2.0 + b * (pv.y + qv.y) / 2.0);
        if (a < -kEps || (std::abs(a) <= kEps && b < 0)) { a=-a; b=-b; c=-c; }
        curve.values.push_back({a,b,c});
    }
    return curve;
}

std::optional<Curve> perpendicular_through(const Point& p, const Curve& line) {
    if (line.kind != CurveKind::Line) return std::nullopt;
    Curve curve{CurveKind::Line, {}, "perpendicular-through(" + p.recipe + ", " + line.recipe + ")",
                {point_key(p)}};
    for (std::size_t i=0; i<p.values.size(); ++i) {
        double a=-line.values[i].b, b=line.values[i].a;
        double c=-(a*p.values[i].x+b*p.values[i].y);
        if (a < -kEps || (std::abs(a)<=kEps && b<0)) {a=-a;b=-b;c=-c;}
        curve.values.push_back({a,b,c});
    }
    return curve;
}

std::optional<Curve> parallel_through(const Point& p, const Curve& line) {
    if (line.kind != CurveKind::Line) return std::nullopt;
    Curve curve{CurveKind::Line, {}, "parallel-through(" + p.recipe + ", " + line.recipe + ")",
                {point_key(p)}};
    for (std::size_t i=0; i<p.values.size(); ++i) {
        const double a=line.values[i].a,b=line.values[i].b;
        curve.values.push_back({a,b,-(a*p.values[i].x+b*p.values[i].y)});
    }
    return curve;
}

std::vector<Curve> angle_bisectors(const Curve& first, const Curve& second) {
    if (first.kind != CurveKind::Line || second.kind != CurveKind::Line) return {};
    std::vector<Curve> result(2);
    for (int branch=0;branch<2;++branch) {
        const double sign=branch==0?1.0:-1.0;
        result[branch].kind=CurveKind::Line;
        result[branch].recipe="angle-bisector("+first.recipe+", "+second.recipe+")#"+std::to_string(branch);
        for (std::size_t i=0;i<first.values.size();++i) {
            double a=first.values[i].a+sign*second.values[i].a;
            double b=first.values[i].b+sign*second.values[i].b;
            double c=first.values[i].c+sign*second.values[i].c;
            const double norm=std::hypot(a,b);
            if(norm<kEps)return {};
            a/=norm;b/=norm;c/=norm;
            if(a < -kEps || (std::abs(a)<=kEps&&b<0)){a=-a;b=-b;c=-c;}
            result[branch].values.push_back({a,b,c});
        }
    }
    return result;
}

std::vector<Vec2> line_line(Triple x, Triple y) {
    const double d=x.a*y.b-y.a*x.b;
    if(std::abs(d)<kEps)return {};
    return {{(x.b*y.c-y.b*x.c)/d,(x.c*y.a-y.c*x.a)/d}};
}

std::vector<Vec2> line_circle(Triple l, Triple c) {
    const double signed_distance=l.a*c.a+l.b*c.b+l.c;
    if(std::abs(signed_distance)>c.c+kEps)return {};
    const double fx=c.a-signed_distance*l.a,fy=c.b-signed_distance*l.b;
    const double delta=std::sqrt(std::max(0.0,c.c*c.c-signed_distance*signed_distance));
    if(delta<kEps)return {{fx,fy}};
    return {{fx-l.b*delta,fy+l.a*delta},{fx+l.b*delta,fy-l.a*delta}};
}

std::vector<Vec2> circle_circle(Triple x, Triple y) {
    const double dx=y.a-x.a,dy=y.b-x.b,d=std::hypot(dx,dy);
    if(d<kEps||d>x.c+y.c+kEps||d<std::abs(x.c-y.c)-kEps)return {};
    const double along=(x.c*x.c-y.c*y.c+d*d)/(2*d);
    const double h2=x.c*x.c-along*along;
    if(h2 < -kEps)return {};
    const double h=std::sqrt(std::max(0.0,h2));
    const double mx=x.a+along*dx/d,my=x.b+along*dy/d;
    if(h<kEps)return {{mx,my}};
    return {{mx-dy*h/d,my+dx*h/d},{mx+dy*h/d,my-dx*h/d}};
}

std::vector<Vec2> pair_intersections(const Curve& x,const Curve& y,std::size_t sample) {
    if(x.kind==CurveKind::Line&&y.kind==CurveKind::Line)return line_line(x.values[sample],y.values[sample]);
    if(x.kind==CurveKind::Line)return line_circle(x.values[sample],y.values[sample]);
    if(y.kind==CurveKind::Line)return line_circle(y.values[sample],x.values[sample]);
    return circle_circle(x.values[sample],y.values[sample]);
}

std::vector<Vec2> line_hyperbola(Triple line) {
    std::vector<Vec2> out;
    if(std::abs(line.b)>kEps){
        // a*x^2+c*x+b=0
        if(std::abs(line.a)<kEps){
            if(std::abs(line.c)>kEps){const double x=-line.b/line.c;
                if(std::abs(x)>kEps)out.push_back({x,1/x});}
        }else{
            const double disc=line.c*line.c-4*line.a*line.b;
            if(disc>=-kEps){
                const double root=std::sqrt(std::max(0.0,disc));
                for(double x:{(-line.c-root)/(2*line.a),(-line.c+root)/(2*line.a)})
                    if(std::abs(x)>kEps)out.push_back({x,1/x});
                if(out.size()==2&&distance(out[0],out[1])<1e-7)out.resize(1);
            }
        }
    }else if(std::abs(line.a)>kEps){
        const double x=-line.c/line.a;if(std::abs(x)>kEps)out.push_back({x,1/x});
    }
    return out;
}

std::complex<double> polynomial(const std::vector<double>& c,std::complex<double> z){
    std::complex<double> value=0.0;for(double coefficient:c)value=value*z+coefficient;return value;
}

std::vector<Vec2> circle_hyperbola(Triple circle) {
    std::vector<double> c={1.0,-2*circle.a,circle.a*circle.a+circle.b*circle.b-circle.c*circle.c,-2*circle.b,1.0};
    const auto largest=*std::max_element(c.begin()+1,c.end(),[](double a,double b){return std::abs(a)<std::abs(b);});
    const double radius=1.0+std::abs(largest);
    std::vector<std::complex<double>> roots;
    for(int i=0;i<4;++i)roots.push_back(std::polar(std::abs(radius),2*kPi*i/4.0+0.173));
    for(int iteration=0;iteration<120;++iteration){
        double movement=0;
        for(int i=0;i<4;++i){
            std::complex<double> denom=1.0;
            for(int j=0;j<4;++j)if(i!=j)denom*=roots[i]-roots[j];
            if(std::abs(denom)<1e-20)denom={1e-20,1e-20};
            const auto delta=polynomial(c,roots[i])/denom;roots[i]-=delta;movement=std::max(movement,std::abs(delta));
        }
        if(movement<1e-13)break;
    }
    std::vector<Vec2> out;
    for(auto root:roots)if(std::abs(root.imag())<2e-7&&std::abs(root.real())>kEps){
        double x=root.real();
        for(int iteration=0;iteration<10;++iteration){
            const double value=((((c[0]*x+c[1])*x+c[2])*x+c[3])*x+c[4]);
            const double derivative=((4*c[0]*x+3*c[1])*x+2*c[2])*x+c[3];
            if(std::abs(derivative)<1e-14)break;
            const double next=x-value/derivative;if(!std::isfinite(next))break;x=next;
        }
        Vec2 point{x,1/x};
        if(std::none_of(out.begin(),out.end(),[&](Vec2 p){return distance(p,point)<1e-6;}))out.push_back(point);
    }
    return out;
}

std::vector<Point> stable_branches(std::vector<std::vector<Vec2>> roots,const std::string& recipe){
    if(roots.empty()||roots[0].empty())return {};
    for(auto& sample:roots)for(auto& point:sample){
        if(std::hypot(point.x,point.y)<1e-7){point={0,0};continue;}
        if(std::abs(point.x)<1e-12)point.x=0;
        if(std::abs(point.y)<1e-12)point.y=0;
    }
    for(auto& sample:roots){
        std::vector<Vec2> unique;
        for(const auto point:sample)if(std::none_of(unique.begin(),unique.end(),[&](Vec2 old){return distance(old,point)<1e-7;}))unique.push_back(point);
        sample=std::move(unique);
    }
    const std::size_t count=roots[0].size();
    if(std::any_of(roots.begin(),roots.end(),[&](const auto& v){return v.size()!=count;}))return {};
    for(auto& values:roots)std::sort(values.begin(),values.end(),[](Vec2 a,Vec2 b){
        return std::tuple{std::atan2(a.y,a.x),a.x,a.y}<std::tuple{std::atan2(b.y,b.x),b.x,b.y};});
    std::vector<Point> result(count);
    for(std::size_t branch=0;branch<count;++branch){
        result[branch].recipe=recipe+"#"+std::to_string(branch);
        for(const auto& sample:roots)result[branch].values.push_back(sample[branch]);
    }
    return result;
}

std::vector<Point> hyperbola_branches(const Curve& curve){
    static std::unordered_map<std::string,std::vector<Point>> cache;
    const auto key=curve_key(curve,10);
    if(const auto found=cache.find(key);found!=cache.end()){
        auto result=found->second;for(std::size_t i=0;i<result.size();++i)result[i].recipe="intersect("+curve.recipe+", H)#"+std::to_string(i);return result;
    }
    std::vector<std::vector<Vec2>> roots;
    for(const auto value:curve.values)roots.push_back(curve.kind==CurveKind::Line?line_hyperbola(value):circle_hyperbola(value));
    auto result=stable_branches(std::move(roots),"intersect("+curve.recipe+", H)");
    if(cache.size()>=20000) cache.clear();
    auto cached=result;for(auto& point:cached)point.recipe.clear();
    cache.emplace(key,std::move(cached));
    return result;
}

std::vector<Point> pair_branches(const Curve& first,const Curve& second){
    static std::unordered_map<std::string,std::vector<Point>> cache;
    auto first_key=curve_key(first,10),second_key=curve_key(second,10);if(second_key<first_key)std::swap(first_key,second_key);
    const auto key=first_key+'|'+second_key;
    if(const auto found=cache.find(key);found!=cache.end()){
        auto result=found->second;for(std::size_t i=0;i<result.size();++i)result[i].recipe="intersect("+first.recipe+", "+second.recipe+")#"+std::to_string(i);return result;
    }
    std::vector<std::vector<Vec2>> roots;
    for(std::size_t i=0;i<first.values.size();++i)roots.push_back(pair_intersections(first,second,i));
    auto result=stable_branches(std::move(roots),"intersect("+first.recipe+", "+second.recipe+")");
    if(cache.size()>=20000) cache.clear();
    auto cached=result;for(auto& point:cached)point.recipe.clear();
    cache.emplace(key,std::move(cached));
    return result;
}

std::optional<State> add_curve(const State& parent,const Candidate& candidate,std::size_t max_points){
    const auto key=curve_key(candidate.curve);
    if(std::any_of(parent.curves.begin(),parent.curves.end(),[&](const Curve& c){return curve_key(c)==key;}))return std::nullopt;
    State child=parent;
    Curve new_curve=candidate.curve;
    const std::string expression=new_curve.recipe;
    new_curve.recipe="c"+std::to_string(child.curves.size());
    std::unordered_set<std::string> known;
    for(const auto& point:child.points)known.insert(point_key(point));
    auto record=[&](Point point,bool on_hyperbola){
        const auto pk=point_key(point);
        if(on_hyperbola)child.hyperbola_points.insert(pk);
        if(!known.insert(pk).second)return true;
        if(child.points.size()>=max_points)return false;
        const std::string definition=point.recipe;
        point.recipe="p"+std::to_string(child.points.size());
        child.point_definitions.push_back(point.recipe+" = "+definition);
        child.points.push_back(std::move(point));return true;
    };
    for(const auto& point:hyperbola_branches(new_curve))if(!record(point,true))return std::nullopt;
    for(const auto& old:parent.curves)for(const auto& point:pair_branches(new_curve,old))if(!record(point,false))return std::nullopt;
    for(const auto& input:candidate.curve.input_points)if(child.hyperbola_points.count(input))child.used_hyperbola_point=true;
    child.curves.push_back(std::move(new_curve));child.steps.push_back(child.curves.back().recipe+" = "+expression);child.cost+=candidate.cost;
    return child;
}

std::vector<Candidate> candidates(const State& state,int remaining){
    std::vector<Candidate> result;
    for(std::size_t i=0;i<state.points.size();++i)for(std::size_t j=i+1;j<state.points.size();++j){
        if(auto c=make_line(state.points[i],state.points[j]))result.push_back({*c,OpKind::Line,1});
        if(auto c=make_circle(state.points[i],state.points[j]))result.push_back({*c,OpKind::Circle,1});
        if(auto c=make_circle(state.points[j],state.points[i]))result.push_back({*c,OpKind::Circle,1});
        if(remaining>=3)if(auto c=perpendicular_bisector(state.points[i],state.points[j]))result.push_back({*c,OpKind::PerpendicularBisector,3});
    }
    if(remaining>=3)for(const auto& point:state.points)for(const auto& curve:state.curves)if(curve.kind==CurveKind::Line){
        if(auto c=perpendicular_through(point,curve))result.push_back({*c,OpKind::Perpendicular,3});
        if(remaining>=4)if(auto c=parallel_through(point,curve))result.push_back({*c,OpKind::Parallel,4});
    }
    if(remaining>=4)for(std::size_t i=0;i<state.curves.size();++i)for(std::size_t j=i+1;j<state.curves.size();++j)
        for(auto& c:angle_bisectors(state.curves[i],state.curves[j]))result.push_back({std::move(c),OpKind::AngleBisector,4});
    std::unordered_set<std::string> seen;std::vector<Candidate> unique;unique.reserve(result.size());
    for(auto& candidate:result){const auto key=curve_key(candidate.curve,7)+'/'+std::to_string(candidate.cost);
        if(seen.insert(key).second)unique.push_back(std::move(candidate));}
    return unique;
}

std::optional<std::string> target_line(const State& state,const std::vector<double>& degrees){
    for(const auto& curve:state.curves)if(curve.kind==CurveKind::Line){
        bool ok=true;for(std::size_t i=0;i<degrees.size();++i){const double t=degrees[i]*kPi/540.0;const auto l=curve.values[i];
            ok&=std::abs(l.c)+std::abs(l.a*std::cos(t)+l.b*std::sin(t))<2e-7;}
        if(ok)return curve.recipe;
    }return std::nullopt;
}

std::optional<std::string> target_point(const State& state,const std::vector<double>& degrees){
    for(const auto& point:state.points){bool ok=true;for(std::size_t i=0;i<degrees.size();++i){
        const double t=degrees[i]*kPi/540.0,dx=std::cos(t),dy=std::sin(t);const auto p=point.values[i];
        ok&=p.x*dx+p.y*dy>kEps&&std::abs(-dy*p.x+dx*p.y)<2e-7;}if(ok)return point.recipe;
    }return std::nullopt;
}

double best_error(const State& state,const std::vector<double>& degrees){
    double best=std::numeric_limits<double>::infinity();
    for(const auto& point:state.points){
        double total=0;bool valid=true;
        for(std::size_t i=0;i<degrees.size();++i){const auto value=point.values[i];
            if(std::hypot(value.x,value.y)<=kEps){valid=false;break;}
            const double target=degrees[i]*kPi/540.0;
            total+=std::abs(std::remainder(std::atan2(value.y,value.x)-target,2*kPi));
        }
        if(valid)best=std::min(best,total);
    }
    return best;
}

State initial_state(const std::vector<double>& degrees){
    State state;Point o,p,a,b;o.recipe="O";p.recipe="P[given]";a.recipe="A[prefix]";b.recipe="B[prefix-reflection]";
    Curve x{CurveKind::Line,{},"x-axis",{}},y{CurveKind::Line,{},"y-axis",{}},side{CurveKind::Line,{},"angle-side",{}},reflected{CurveKind::Line,{},"l",{}};
    const std::vector<double> radii={1.0,std::sqrt(2.0),kPi/2.0,std::sqrt(3.0)/2.0,1.271,0.937,1.618,1.113};
    for(std::size_t i=0;i<degrees.size();++i){const double t=degrees[i]*kPi/180.0,r=radii[i%radii.size()];Vec2 pv{r*std::cos(t),r*std::sin(t)};
        o.values.push_back({0,0});p.values.push_back(pv);a.values.push_back({0,r});b.values.push_back({-pv.x,pv.y});
        x.values.push_back({0,1,0});y.values.push_back({1,0,0});side.values.push_back(*canonical_line({0,0},pv));reflected.values.push_back(*canonical_line({0,0},{-pv.x,pv.y}));}
    state.points={o,p,a,b};state.curves={x,y,side,reflected};state.steps={"l = hidden y-reflection prefix (cost=3)"};
    state.point_definitions={"O = given origin","P = given point on angle side","A = prefix axis point","B = prefix reflected point"};
    std::unordered_set<std::string> known;for(const auto& point:state.points)known.insert(point_key(point));
    for(auto point:hyperbola_branches(reflected)){state.hyperbola_points.insert(point_key(point));if(known.insert(point_key(point)).second){const auto definition=point.recipe;point.recipe="p"+std::to_string(state.points.size());state.point_definitions.push_back(point.recipe+" = "+definition);state.points.push_back(std::move(point));}}
    return state;
}

State raw_initial_state(const std::vector<double>& degrees){
    State state;state.cost=0;
    Point origin,point;origin.recipe="O";point.recipe="P[given]";
    Curve x_axis{CurveKind::Line,{},"x-axis",{}},y_axis{CurveKind::Line,{},"y-axis",{}},
          angle_side{CurveKind::Line,{},"angle-side",{}};
    const std::vector<double> radii={1.0,std::sqrt(2.0),kPi/2.0,std::sqrt(3.0)/2.0,1.271,0.937,1.618,1.113,1.414};
    for(std::size_t i=0;i<degrees.size();++i){const double theta=degrees[i]*kPi/180.0,radius=radii[i%radii.size()];
        const Vec2 p{radius*std::cos(theta),radius*std::sin(theta)};
        origin.values.push_back({0,0});point.values.push_back(p);
        x_axis.values.push_back({0,1,0});y_axis.values.push_back({1,0,0});
        angle_side.values.push_back(*canonical_line({0,0},p));
    }
    state.points={origin,point};state.curves={x_axis,y_axis,angle_side};
    state.point_definitions={"O = given origin","P = arbitrary given point on angle side"};
    return state;
}

std::vector<Candidate> weighted_order(std::vector<Candidate> input,std::mt19937_64& rng){
    std::vector<std::vector<Candidate>> groups(6);
    for(auto& candidate:input)groups[static_cast<std::size_t>(candidate.kind)].push_back(std::move(candidate));
    for(auto& group:groups)std::shuffle(group.begin(),group.end(),rng);
    const double weights[]={28,28,12,12,12,8};std::vector<Candidate> output;output.reserve(input.size());
    while(output.size()<input.size()){
        std::vector<double> active(6);for(std::size_t i=0;i<6;++i)if(!groups[i].empty())active[i]=weights[i];
        std::discrete_distribution<std::size_t> choose(active.begin(),active.end());const auto index=choose(rng);
        output.push_back(std::move(groups[index].back()));groups[index].pop_back();
    }return output;
}

struct BranchRecord {
    Curve curve;
    int extension_cost{};
    std::vector<std::string> steps;
    std::vector<std::string> point_definitions;
    std::optional<std::string> direct_final_step;
    int direct_total_cost{};
};

std::vector<std::vector<double>> target_ray_branches(const Curve& curve,
                                                     const std::vector<double>& degrees){
    std::vector<std::vector<double>> per_sample;
    for(std::size_t i=0;i<degrees.size();++i){
        const double target=degrees[i]*kPi/540.0,dx=std::cos(target),dy=std::sin(target);
        const auto value=curve.values[i];std::vector<double> roots;
        if(curve.kind==CurveKind::Line){
            const double denominator=value.a*dx+value.b*dy;
            if(std::abs(denominator)>kEps){const double distance=-value.c/denominator;
                if(distance>kEps&&std::isfinite(distance))roots.push_back(distance);}
        }else{
            const double projection=value.a*dx+value.b*dy;
            const double perpendicular2=value.a*value.a+value.b*value.b-projection*projection;
            const double delta2=value.c*value.c-perpendicular2;
            if(delta2>=-kEps){const double delta=std::sqrt(std::max(0.0,delta2));
                for(double distance:{projection-delta,projection+delta})
                    if(distance>kEps&&std::isfinite(distance))roots.push_back(distance);
                std::sort(roots.begin(),roots.end());
                if(roots.size()==2&&std::abs(roots[0]-roots[1])<1e-8)roots.resize(1);
            }
        }
        per_sample.push_back(std::move(roots));
    }
    if(per_sample.empty()||per_sample[0].empty())return {};
    const auto count=per_sample[0].size();
    if(std::any_of(per_sample.begin(),per_sample.end(),[&](const auto& roots){return roots.size()!=count;}))return {};
    std::vector<std::vector<double>> branches(count);
    for(std::size_t branch=0;branch<count;++branch)for(const auto& roots:per_sample)branches[branch].push_back(roots[branch]);
    return branches;
}

std::string ray_fingerprint(const std::vector<double>& distances){
    std::ostringstream out;for(double distance:distances)out<<std::llround(distance*1e4)<<',';return out.str();
}

bool same_ray_point(const std::vector<double>& first,const std::vector<double>& second){
    if(first.size()!=second.size())return false;
    for(std::size_t i=0;i<first.size();++i)if(std::abs(first[i]-second[i])>3e-6)return false;
    return true;
}

bool same_geometric_curve(const Curve& first,const Curve& second){
    if(first.kind!=second.kind||first.values.size()!=second.values.size())return false;
    for(std::size_t i=0;i<first.values.size();++i){const auto a=first.values[i],b=second.values[i];
        const double scale=1+std::max({std::abs(a.a),std::abs(a.b),std::abs(a.c),std::abs(b.a),std::abs(b.b),std::abs(b.c)});
        if(std::max({std::abs(a.a-b.a),std::abs(a.b-b.b),std::abs(a.c-b.c)})>1e-5*scale)return false;}
    return true;
}

bool transverse_target_meet(const Curve& first,const Curve& second,
                            const std::vector<double>& distances,const std::vector<double>& degrees){
    for(std::size_t i=0;i<degrees.size();++i){const double target=degrees[i]*kPi/540.0;
        const Vec2 point{distances[i]*std::cos(target),distances[i]*std::sin(target)};
        Vec2 first_normal,second_normal;const auto a=first.values[i],b=second.values[i];
        if(first.kind==CurveKind::Line)first_normal={a.a,a.b};
        else{const double norm=std::hypot(point.x-a.a,point.y-a.b);if(norm<kEps)return false;first_normal={(point.x-a.a)/norm,(point.y-a.b)/norm};}
        if(second.kind==CurveKind::Line)second_normal={b.a,b.b};
        else{const double norm=std::hypot(point.x-b.a,point.y-b.b);if(norm<kEps)return false;second_normal={(point.x-b.a)/norm,(point.y-b.b)/norm};}
        if(std::abs(first_normal.x*second_normal.y-first_normal.y*second_normal.x)<1e-5)return false;
    }
    return true;
}

std::optional<State> random_shared_state(const State& root,int desired_cost,std::size_t max_points,
                                         std::mt19937_64& rng,std::size_t& generated){
    State state=root;
    while(state.cost<desired_cost){
        auto options=weighted_order(candidates(state,desired_cost-state.cost),rng);bool advanced=false;
        for(auto& candidate:options){if(state.cost+candidate.cost>desired_cost)continue;
            ++generated;auto child=add_curve(state,candidate,max_points);if(!child)continue;
            state=std::move(*child);advanced=true;break;}
        if(!advanced)return std::nullopt;
    }
    return state;
}

std::vector<BranchRecord> enumerate_branch_records(const State& base,const MeetConfig& config,
                                                   std::mt19937_64& rng,std::size_t& expanded,
                                                   std::size_t& generated,std::size_t& duplicates,
                                                   const std::chrono::steady_clock::time_point& started){
    struct Frame{State state;std::vector<Candidate> candidates;std::size_t index{};bool prepared{};};
    std::vector<Frame> stack{{base,{},0,false}};std::unordered_map<std::string,int> seen;
    seen[state_key(base)+(base.used_hyperbola_point?"/1":"/0")]=base.cost;
    std::unordered_map<std::string,std::size_t> best_record;
    std::vector<BranchRecord> records;
    const std::size_t generated_at_start=generated;
    auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();};
    while(!stack.empty()&&records.size()<config.branch_records_per_base&&
          generated-generated_at_start<config.branch_attempts_per_base&&
          elapsed()<config.time_limit_seconds){
        auto& frame=stack.back();const int remaining=config.max_cost-frame.state.cost-2;
        if(remaining<1){stack.pop_back();continue;}
        if(!frame.prepared){frame.candidates=weighted_order(candidates(frame.state,remaining),rng);frame.prepared=true;++expanded;}
        if(frame.index>=frame.candidates.size()){stack.pop_back();continue;}
        Candidate candidate=std::move(frame.candidates[frame.index++]);if(candidate.cost>remaining)continue;
        ++generated;auto child=add_curve(frame.state,candidate,config.max_points);if(!child)continue;
        const auto skey=state_key(*child)+(child->used_hyperbola_point?"/1":"/0");
        if(auto old=seen.find(skey);old!=seen.end()&&old->second<=child->cost){++duplicates;continue;}seen[skey]=child->cost;
        BranchRecord record;record.curve=child->curves.back();record.extension_cost=child->cost-base.cost;
        record.steps.assign(child->steps.begin()+static_cast<std::ptrdiff_t>(base.steps.size()),child->steps.end());
        record.point_definitions.assign(
            child->point_definitions.begin()+static_cast<std::ptrdiff_t>(base.point_definitions.size()),
            child->point_definitions.end());
        if(target_line(*child,config.sample_degrees)){
            record.direct_final_step="target line drawn";record.direct_total_cost=child->cost;
        }else if(child->cost+1<=config.max_cost){if(auto hit=target_point(*child,config.sample_degrees)){
            record.direct_final_step="line(O, "+*hit+")";record.direct_total_cost=child->cost+1;}}
        const auto ckey=curve_key(record.curve,7);
        auto previous=best_record.find(ckey);
        if(previous==best_record.end()){best_record.emplace(ckey,records.size());records.push_back(std::move(record));}
        else if(record.extension_cost<records[previous->second].extension_cost)records[previous->second]=std::move(record);
        stack.push_back({std::move(*child),{},0,false});
    }
    return records;
}

} // namespace

ModelDiagnostics diagnose_model(const std::vector<double>& sample_degrees){
    const auto state=initial_state(sample_degrees);
    ModelDiagnostics diagnostics;
    diagnostics.prefix_cost=state.cost;
    diagnostics.visible_curves=state.curves.size();
    diagnostics.visible_points=state.points.size();
    diagnostics.stable_hyperbola_points=state.hyperbola_points.size();
    diagnostics.helper_circles_hidden=std::none_of(state.curves.begin(),state.curves.end(),[](const Curve& curve){return curve.kind==CurveKind::Circle;});

    Curve unstable{CurveKind::Circle,{},"unstable-test",{}};
    // The unit circle meets the x-axis in all samples, except the final circle,
    // which is moved away. Equal-count stability must reject the whole branch.
    Curve axis{CurveKind::Line,{},"axis",{}};
    for(std::size_t i=0;i<sample_degrees.size();++i){
        unstable.values.push_back(i+1==sample_degrees.size()?Triple{0,3,1}:Triple{0,0,1});
        axis.values.push_back({0,1,0});
    }
    diagnostics.unstable_branch_rejected=pair_branches(unstable,axis).empty();
    diagnostics.quartic_solver_verified=circle_hyperbola({0,0,2}).size()==4;
    const auto available=candidates(state,4);
    bool perpendicular=false,parallel=false,bisector=false,angle=false;
    for(const auto& candidate:available){
        perpendicular|=candidate.kind==OpKind::Perpendicular&&candidate.cost==3;
        parallel|=candidate.kind==OpKind::Parallel&&candidate.cost==4;
        bisector|=candidate.kind==OpKind::PerpendicularBisector&&candidate.cost==3;
        angle|=candidate.kind==OpKind::AngleBisector&&candidate.cost==4;
    }
    diagnostics.macro_costs_verified=perpendicular&&parallel&&bisector&&angle;
    for(const auto& first:available){
        if(first.cost!=1)continue;
        auto child=add_curve(state,first,40);if(!child||child->cost!=4)continue;
        for(const auto& second:candidates(*child,1)){
            if(second.cost!=1)continue;
            auto grandchild=add_curve(*child,second,40);
            if(grandchild&&grandchild->cost==5){diagnostics.existing_objects_reused_without_recharge=true;break;}
        }
        if(diagnostics.existing_objects_reused_without_recharge)break;
    }
    Curve join_first{CurveKind::Line,{},"join-a",{}},join_second{CurveKind::Line,{},"join-b",{}};
    for(double degree:sample_degrees){const double target=degree*kPi/540.0;const Vec2 hit{2*std::cos(target),2*std::sin(target)};
        join_first.values.push_back(*canonical_line(hit,{hit.x+1,hit.y}));
        join_second.values.push_back(*canonical_line(hit,{hit.x,hit.y+1}));}
    const auto first_branches=target_ray_branches(join_first,sample_degrees);
    const auto second_branches=target_ray_branches(join_second,sample_degrees);
    diagnostics.target_ray_join_verified=first_branches.size()==1&&second_branches.size()==1&&
        ray_fingerprint(first_branches[0])==ray_fingerprint(second_branches[0])&&
        same_ray_point(first_branches[0],second_branches[0]);
    Curve perturbed=join_first;for(auto& value:perturbed.values)value.c+=1e-8;
    diagnostics.near_duplicate_curve_rejected=same_geometric_curve(join_first,perturbed);
    diagnostics.transverse_join_verified=transverse_target_meet(join_first,join_second,first_branches[0],sample_degrees)&&
        !transverse_target_meet(join_first,perturbed,first_branches[0],sample_degrees);
    if(auto prefix_circle=make_circle(state.points[2],state.points[0])){
        const auto branches=pair_branches(*prefix_circle,state.curves[2]);
    diagnostics.origin_branch_verified=std::any_of(branches.begin(),branches.end(),[](const Point& point){
            return std::all_of(point.values.begin(),point.values.end(),[](Vec2 value){return value.x==0&&value.y==0;});});
    }
    const auto raw=raw_initial_state(sample_degrees);
    diagnostics.raw_initial_state_verified=raw.cost==0&&raw.points.size()==2&&raw.curves.size()==3&&
        std::none_of(raw.curves.begin(),raw.curves.end(),[](const Curve& curve){return curve.recipe=="l";});
    return diagnostics;
}

SearchReport search(const SearchConfig& input){
    SearchConfig config=input;if(config.sample_degrees.empty())config.sample_degrees={91,103,117,129,141,153,167,179};
    SearchReport report;const auto started=std::chrono::steady_clock::now();double best=std::numeric_limits<double>::infinity();
    auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();};
    while(elapsed()<config.time_limit_seconds){
        std::mt19937_64 rng(config.seed+report.restarts*1000003ULL);
        State root=config.reflected_prefix?initial_state(config.sample_degrees):raw_initial_state(config.sample_degrees);
        struct Frame{State state;std::vector<Candidate> candidates;std::size_t index{};bool prepared{};};
        std::vector<Frame> stack{{std::move(root),{},0,false}};std::unordered_map<std::string,int> seen;seen[state_key(stack[0].state)+"/0"]=3;
        std::size_t local_generated=0,local_terminals=0;
        while(!stack.empty()&&local_generated<config.expansions_per_restart&&local_terminals<config.terminal_states_per_restart&&elapsed()<config.time_limit_seconds){
            auto& frame=stack.back();
            if(frame.state.cost>=config.max_cost){++local_terminals;stack.pop_back();continue;}
            if(!frame.prepared){frame.candidates=weighted_order(candidates(frame.state,config.max_cost-frame.state.cost),rng);
                frame.prepared=true;++report.expanded;}
            if(frame.index>=frame.candidates.size()){stack.pop_back();continue;}
            Candidate candidate=std::move(frame.candidates[frame.index++]);if(frame.state.cost+candidate.cost>config.max_cost)continue;
            auto child=add_curve(frame.state,candidate,config.max_points);++local_generated;++report.generated;if(!child)continue;
            const auto key=state_key(*child)+(child->used_hyperbola_point?"/1":"/0");auto old=seen.find(key);if(old!=seen.end()&&old->second<=child->cost){++report.duplicates;continue;}seen[key]=child->cost;
            const double error=best_error(*child,config.sample_degrees);if(error<best){best=error;report.best_steps=child->steps;report.best_point_definitions=child->point_definitions;std::ostringstream text;text<<"best angular-error-sum="<<error*180/kPi<<" deg at cost="<<child->cost;report.best_description=text.str();}
            if(auto hit=target_line(*child,config.sample_degrees)){report.found=true;report.cost=child->cost;report.steps=child->steps;for(const auto& definition:child->point_definitions)report.steps.push_back("[free] "+definition);report.elapsed_seconds=elapsed();return report;}
            if(child->cost+1<=config.max_cost)if(auto hit=target_point(*child,config.sample_degrees)){report.found=true;report.cost=child->cost+1;report.steps=child->steps;report.steps.push_back("line(O, "+*hit+")");for(const auto& definition:child->point_definitions)report.steps.push_back("[free] "+definition);report.elapsed_seconds=elapsed();return report;}
            const int remaining=config.max_cost-child->cost;if(!child->used_hyperbola_point){if(child->hyperbola_points.empty()&&remaining<2)continue;if(!child->hyperbola_points.empty()&&remaining<1)continue;}
            stack.push_back({std::move(*child),{},0,false});
        }
        ++report.restarts;if(config.verbose)std::cout<<"restart="<<report.restarts<<" generated="<<report.generated<<" expanded="<<report.expanded<<" elapsed="<<elapsed()<<"s "<<report.best_description<<'\n';
    }
    report.elapsed_seconds=elapsed();return report;
}

SearchReport meet_in_the_middle(const MeetConfig& input){
    MeetConfig config=input;if(config.sample_degrees.empty())config.sample_degrees={91,103,117,129,141,153,167,179};
    const int root_cost=config.reflected_prefix?3:0;
    config.shared_cost_min=std::max(root_cost,config.shared_cost_min);
    config.shared_cost_max=std::min(config.shared_cost_max,config.max_cost-3);
    SearchReport report;const auto started=std::chrono::steady_clock::now();
    auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();};
    State root=config.reflected_prefix?initial_state(config.sample_degrees):raw_initial_state(config.sample_degrees);
    while(elapsed()<config.time_limit_seconds){
        std::mt19937_64 rng(config.seed+report.restarts*1000003ULL);
        std::uniform_int_distribution<int> choose_cost(config.shared_cost_min,config.shared_cost_max);
        auto base=random_shared_state(root,choose_cost(rng),config.max_points,rng,report.generated);
        ++report.restarts;if(!base)continue;
        if(auto line=target_line(*base,config.sample_degrees)){report.found=true;report.cost=base->cost;report.steps=base->steps;report.elapsed_seconds=elapsed();return report;}
        auto records=enumerate_branch_records(*base,config,rng,report.expanded,report.generated,report.duplicates,started);
        std::unordered_map<std::string,std::vector<std::pair<std::size_t,std::vector<double>>>> index;
        for(std::size_t current=0;current<records.size();++current){auto& record=records[current];
            if(record.direct_final_step){report.found=true;report.cost=record.direct_total_cost;report.steps=base->steps;
                for(const auto& step:record.steps)report.steps.push_back("[branch] "+step);
                report.steps.push_back(*record.direct_final_step);report.elapsed_seconds=elapsed();return report;}
            for(const auto& distances:target_ray_branches(record.curve,config.sample_degrees)){
                const auto fingerprint=ray_fingerprint(distances);auto& previous=index[fingerprint];
                for(const auto& [other_index,other_distances]:previous){const auto& other=records[other_index];
                    const int total=base->cost+record.extension_cost+other.extension_cost+1;
                    if(total>config.max_cost||same_geometric_curve(record.curve,other.curve)||
                       !same_ray_point(distances,other_distances)||
                       !transverse_target_meet(record.curve,other.curve,distances,config.sample_degrees))continue;
                    report.found=true;report.cost=total;report.steps=base->steps;
                    for(const auto& definition:base->point_definitions)report.steps.push_back("[shared free] "+definition);
                    for(const auto& step:other.steps)report.steps.push_back("[left] "+step);
                    for(const auto& definition:other.point_definitions)report.steps.push_back("[left free] "+definition);
                    for(const auto& step:record.steps)report.steps.push_back("[right] "+step);
                    for(const auto& definition:record.point_definitions)report.steps.push_back("[right free] "+definition);
                    report.steps.push_back("line(O, intersection(left-final-curve, right-final-curve))");
                    report.elapsed_seconds=elapsed();return report;
                }
                previous.push_back({current,distances});
            }
        }
        if(config.verbose)std::cout<<"mitm base="<<report.restarts<<" records="<<records.size()
            <<" generated="<<report.generated<<" elapsed="<<elapsed()<<"s\n";
    }
    report.elapsed_seconds=elapsed();return report;
}

} // namespace neo::obtuse
