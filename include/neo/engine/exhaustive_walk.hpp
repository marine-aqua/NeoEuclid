#pragma once

// The exhaustive walk that the searches of the engine are built on.
//
// A walk starts from a root: some points that are given, some curves that are
// already drawn, and the points that those curves leave behind.  It then draws
// one more curve at a time, and every curve it may draw is a line through two
// known points, the circle about one known point through another, or the
// circle about that other point through the one.  Intersections are free, so a
// point of a figure is a crossing of two drawn curves, and the cost of a
// figure is its number of curves.
//
// The walk has two phases and they are apart on purpose.  The first phase
// materialises the layers of the figure, one cost at a time, and it is the
// phase that costs the memory.  The second phase materialises nothing: it
// takes the last layer and tries every curve that one more curve could be,
// and only the ones that pass near the target, which is not a cut of the
// search, because a figure that owns the target has to make it with the curve
// it draws last.
//
// The walk is complete for the budget it is given as long as no state was cut
// by the point budget, the state budget was not reached and the time was not
// up; the report says how many states were cut and whether the budget was
// reached, so that "not found" can be read for what it is.

#include "neo/engine/fast_geometry.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace neo::engine::exhaustive_walk {

using fast::GCurve;
using fast::GPoint;

struct Config {
    // The number of curves the figure may have, the curves of the root
    // included.  A root with four curves and a budget of eight may draw four
    // more.
    int max_cost{6};
    // The number of points a state may hold.  A state that would hold more is
    // cut, and the walk is not complete when that happens.
    std::size_t max_points{48};
    // The rule that the last two levels of the walk are filtered with.
    //
    // A figure owns the target only if two of its curves pass through it, so
    // the level before the last one has to hold a curve that passes through
    // the target, and the last curve has to pass through the target as well.
    // With this rule the walk materialises the layers only up to
    // `max_cost - 2` and filters the last two levels, and the level before the
    // last one is not materialised either: a state of it holds a curve
    // through the target only if the state below it held one, or if the curve
    // that was added passes through the target.  The layer of the last level
    // never has to be in memory, which is what lets a walk of seven curves run
    // on a machine that cannot hold the layer of six.
    //
    // The rule is a rule of the search and not a cut of it: a figure that owns
    // the target has two curves that pass through it, and the curve it draws
    // last is one of them.
    //
    // The curves that are given are part of the figure, so they count towards
    // the two curves through the target.  A given curve is not something the
    // construction draws, but the rule is about the figure and not about what
    // is drawn, and leaving the given curves out would cut away the figures
    // whose only curve through the target is a given one.
    bool target_curve_filter{false};
    // Whether the walk prints every figure of the last layer that holds a
    // curve through the target.  The figures are the whole of the work of the
    // last two levels and a reader may want to look at them, but a task whose
    // given curve passes through the target has such a figure in every state
    // of the layer, and the print is then only noise.
    bool print_target_figures{true};
    // The walk stops when it has visited this many states, so that a walk of a
    // task whose layers grow fast runs out of budget instead of out of memory.
    std::size_t max_states{2'000'000};
    // The walk stops when the layers it has built take more than this many
    // bytes, which is the budget that keeps the machine alive.
    std::size_t memory_bytes{12ull * 1024 * 1024 * 1024};
    double time_limit_seconds{900.0};
    // A curve that does not come this near the target cannot make it.
    double align_tolerance{1e-9};
    // The distance at which a point of a figure is the target.
    double confirm_tolerance{1e-11};
    // Whether the walk collects every figure that owns the target instead of
    // stopping at the first one.
    //
    // The walk stops at the first hit because that is all a search needs: the
    // question it answers is whether a figure of the budget owns the target.
    // The question "how many figures of the budget own it" needs the other
    // behaviour, and it is a different walk: the last two levels are then not
    // left half done, and the walk takes the time the whole of the last two
    // levels takes.
    bool collect_all_hits{false};
    bool verbose{true};
};

// A curve of a figure, with the two points that define it.  The coordinates
// are what a rebuild in more digits needs: a figure is described by the
// points that define its curves, and those points are the whole of it.
struct StepGeometry {
    std::string label;
    GPoint first{};
    GPoint second{};
};

struct Report {
    bool found{false};
    // The number of curves of the figure that owns the target, the curves of
    // the root included.
    int cost{0};
    // The description of every curve of that figure, in the order they were
    // drawn.
    std::vector<std::string> steps;
    // The same curves with the coordinates of the two points that define them.
    std::vector<StepGeometry> geometry;
    // The points of the figure that owns the target, in the order they were
    // found, with their coordinates.  The index of a point is its name.
    std::vector<GPoint> found_points;
    // Every figure that owns the target, when the walk was asked to collect
    // them all: the description of the curves of each of them, in the order
    // they were drawn.  The walk stops at the first hit otherwise, and then
    // this is empty.
    std::vector<std::vector<std::string>> all_hits;
    // The same figures with the coordinates of the two points that define each
    // of their curves, so that a reader can see what a figure is.
    std::vector<std::vector<StepGeometry>> all_hit_geometry;
    // The number of figures that own the target, whether they were collected
    // or not.
    std::size_t all_hit_count{0};
    // The number of curves of the cheapest of those figures.
    int all_hit_cost{0};
    // The number of states of every layer the walk built, the root first.
    std::vector<std::size_t> layer_sizes;
    std::size_t states{0};
    std::size_t candidates{0};
    // The states that the point budget cut, which would make the walk
    // incomplete.
    std::size_t states_cut{0};
    // The number of states of the last materialised layer that held a curve
    // through the target, and the number of children of those states that the
    // filter of the two curves looked at.  They are the whole of the work of
    // the last two levels.
    std::size_t target_curve_states{0};
    std::size_t target_curve_children{0};
    // The walk stopped at the state budget, at the memory budget or at the
    // time limit.
    bool budget_hit{false};
    bool memory_hit{false};
    bool time_limit_hit{false};
    // The peak of the estimated size of the layers that were held.
    std::size_t peak_bytes{0};
    double residual{};
    double elapsed_seconds{};
    std::string note;
};

// The root of a walk: the points that are given, the curves that are already
// drawn, and the names of the given points.  The points that the curves leave
// behind are not in here; `make_root` adds them.
struct Root {
    std::vector<GPoint> points;
    std::vector<GCurve> curves;
    std::vector<std::string> names;
    // The curves of the root with the two points that define them, so that a
    // figure that the walk finds can be rebuilt from the report alone.
    std::vector<StepGeometry> geometry;
};

// The root of a walk: the given points, the given curves, and the crossings of
// those curves with each other, which are points of the figure as well.
Root make_root(const std::vector<GPoint>& points, const std::vector<GCurve>& curves,
               const std::vector<std::string>& names);

// The walk.  The cost of the figure that the report describes is the number of
// curves in it, the curves of the root included.
Report run(const Root& root, GPoint target, const Config& config);

}  // namespace neo::engine::exhaustive_walk
