#pragma once

// A search with a geometric goal.
//
// The walks of the engine look for a figure that satisfies a condition, and
// they measure the condition by how small it is.  That is the wrong measure
// when the question is about a configuration: a small residual says that a
// figure comes close to the target, not that it owns the target, and a search
// that chases the residual walks towards clever approximations instead of
// towards the configuration that was asked for.
//
// This search asks a different question.  The caller says which point, or
// which circle, or which line, the figure has to own, and the search looks for
// the cheapest way to get it.  The test is geometric: the point has to be the
// one.
//
// The search walks the same tree as the exhaustive walks, with one filter
// that keeps it affordable: a curve that cannot come near the target cannot
// put the target into the figure, so it is never built.

#include "neo/engine/fast_geometry.hpp"
#include "neo/engine/operation.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace neo::goal {

// The operations of the All model.  The classic model has only the two atomic
// ones, and the type of the step says which of them it is.
enum class Op : std::uint8_t {
    Legacy = 0,
    TransferCircle,
    PerpendicularBisector,
    PerpendicularThrough,
    ParallelThrough,
    AngleBisector,
    Intersection,
};

// One operation of a recipe: a line through two points of the figure, or a
// circle about one of them through the other, and in the All model the
// operations of that model as well.
//
// The places of a step are the places of the figure: the points and the curves
// have one list each, so `first` and `second` are points for the two atomic
// operations, a point and a curve for the perpendicular and the parallel, two
// curves for the bisector of an angle and for the crossing, and `third` is the
// third point of the rigid circle.
struct Step {
    std::uint8_t type{0};  // 0 = line, 1 = circle
    int first{0};
    int second{0};
    Op op{Op::Legacy};
    int third{0};
    int which{0};  // which of the crossings of the pair
};

// The figure of a search.
struct State {
    std::vector<fast::GPoint> points;
    std::vector<fast::GCurve> curves;
    std::vector<Step> steps;
    int cost{0};
};

// What the figure has to own.
struct Target {
    enum class Kind { Point, Circle, Line } kind{Kind::Point};
    fast::GPoint point{};   // for Kind::Point
    fast::GPoint centre{};  // for Kind::Circle
    double radius{0.0};
    // When this is set the figure has to own that point as well.  It matters
    // when the goal is the length of a segment that is then used as a radius:
    // the circle about P through X is only an operation of the model when P
    // is a point of the figure.
    bool also_own{false};
    fast::GPoint also{};
    // For Kind::Line: the line a x + b y + c = 0 with the unit normal, in the
    // same sign convention as the lines of the figures: the normal points to
    // the right, that is a > 0, or a = 0 and b > 0.
    double line_a{};
    double line_b{};
    double line_c{};
};

struct Config {
    int max_extra{3};            // the curves on top of the start
    double seconds{600.0};       // the time limit
    std::size_t max_points{64};  // the point budget of a state
    double tolerance{1e-9};      // how near a point has to be to the target
    // The filter of the search keeps the curves that can put the target on
    // the paper, and it is applied to every curve that is drawn.  That is
    // right when the target is a circle, and when the target is a point that
    // is one curve away from a known circle, because a curve that meets that
    // circle leaves behind a point from which the target is one curve away.
    // It is too strong for a point that is given outright: a curve that is
    // drawn before the last one need not pass through the target, it only
    // has to leave behind the two points that the last curve is built from.
    // With this flag the filter is applied to the last curve alone, and the
    // walk builds every curve before it.
    bool filter_only_last{false};
    // When this is true the search does not stop at the first figure that owns
    // the goal: it keeps walking the last level and writes down every figure
    // that owns it.  A question like "is there a road that does not go through
    // a tangency" needs all of the roads, not the first one the walk meets.
    bool collect_all{false};
    // The model of the construction.  In the All model every operation is one
    // step and a crossing is an operation like any other, so the filter of the
    // search is applied to the last operation alone in that model as well.
    engine::Mode mode{engine::Mode::Classic};
    bool verbose{true};
};

// The certificate of a construction: the places it makes, and the places
// among them that nothing uses.
//
// The recipe of a construction is a computation: every place of the figure is
// built from the places that came before it, and the goal is one place of the
// figure.  A place that is not on the way to the goal is a place that is
// drawn for nothing, for it can be left out of the recipe and the goal is
// still made.  So a shortest construction has none of them.
//
// The condition is a necessary one and not a sufficient one: a construction
// may have every one of its places on the way to the goal and still be
// shortened by replacing one of its steps by another one.
struct Certificate {
    // Whether every place of the construction is on the way to the goal.
    bool ok{false};
    // The places the construction makes.
    std::size_t places{0};
    // The places among them that nothing uses.
    std::size_t dead{0};
};

struct Report {
    bool found{false};
    int extra{0};                     // the curves on top of the start
    std::vector<Step> steps;          // the curves that were drawn
    std::vector<std::size_t> levels;  // the states of each level
    // the points of the figure that the hit leaves behind, so that the
    // construction can be written down and checked outside the engine
    std::vector<fast::GPoint> points;
    // the curves of the figure, for the same reason: the operations of the
    // All model that name a line need it.
    std::vector<fast::GCurve> curves;
    // the states that the point budget cut.  A state that is cut is a state
    // the search did not look at, so a search that found nothing is only a
    // proof when this number is zero.
    std::uint64_t states_cut{0};
    std::uint64_t states{0};
    std::uint64_t candidates{0};
    // Whether the clock stopped the walk.  A search that found nothing is a
    // proof only when nothing stopped it: not a state by the point budget, and
    // not the clock.
    bool out_of_time{false};
    double elapsed{0.0};
    // Every figure that owns the goal, when the search was asked to collect
    // them all.  The recipe of each one and the points it leaves behind.
    std::vector<std::vector<Step>> hits;
    std::vector<std::vector<fast::GPoint>> hit_points;
    std::vector<std::vector<fast::GCurve>> hit_curves;
    // The certificate of the construction the search found: every place of it
    // is on the way to the goal, so no step of it can simply be dropped.
    Certificate certificate;
    // The certificate of each of the figures that were collected.
    std::vector<Certificate> hit_certificates;
};

// The figure after one more curve: the step names two points of the figure
// that the new curve is drawn from, and the crossings of the new curve with
// the curves that are already there become points of the figure.
//
// In the All model the crossings are operations of their own and the figure
// does not get them for nothing, so the mode says whether they come along.
bool extend(State& state, const Step& step, std::size_t max_points,
            engine::Mode mode = engine::Mode::Classic);

// Builds the figure of a recipe: each step names two points of the figure
// that the step is drawn from, and the crossings of the new curve with the
// curves that are already there become points of the figure.
bool build(State& state, const std::vector<Step>& recipe, std::size_t max_points,
           engine::Mode mode = engine::Mode::Classic);

// The certificate of a recipe: the places it makes, and the places among them
// that nothing uses.  A place that is not on the way to the goal is a place
// that is drawn for nothing, so a shortest construction has none of them.
//
// The goal is looked for in the figure the recipe makes, so the state that is
// given is the state the walk starts from: in the classic model the crossings
// of the curves that are already there belong to it.
Certificate certify(const State& start, const std::vector<Step>& recipe,
                    const Target& target, double tolerance = 1e-9,
                    engine::Mode mode = engine::Mode::Classic);

// The search.  The start state is what the caller has already drawn.
Report search(const State& start, const Target& target, const Config& config);

// A readable account of one step of a figure, for the modules that print
// their results.  The caller names the points and the curves of the figure, so
// that every module can keep its own naming.
template <typename PointName, typename CurveName>
std::string describe_step(const Step& step, PointName point_name, CurveName curve_name) {
    switch (step.op) {
        case Op::Legacy:
            if (step.type == 0)
                return "line(" + point_name(step.first) + ", " + point_name(step.second) + ")";
            return "circle(" + point_name(step.first) + "; " + point_name(step.second) + ")";
        case Op::TransferCircle:
            return "transfer-circle(" + point_name(step.first) + "; " + point_name(step.second) +
                   ", " + point_name(step.third) + ")";
        case Op::PerpendicularBisector:
            return "perpendicular-bisector(" + point_name(step.first) + ", " +
                   point_name(step.second) + ")";
        case Op::PerpendicularThrough:
            return "perpendicular(" + point_name(step.first) + ", " + curve_name(step.second) + ")";
        case Op::ParallelThrough:
            return "parallel(" + point_name(step.first) + ", " + curve_name(step.second) + ")";
        case Op::AngleBisector:
            return "angle-bisector(" + curve_name(step.first) + ", " + curve_name(step.second) + ")";
        case Op::Intersection:
            return "intersect(" + curve_name(step.first) + ", " + curve_name(step.second) + ")";
    }
    return "operation";
}

}  // namespace neo::goal
