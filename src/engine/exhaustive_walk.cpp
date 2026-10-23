#include "neo/engine/exhaustive_walk.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace neo::engine::exhaustive_walk {
namespace {

using fast::fast_intersections;
using fast::make_circle;
using fast::make_line;
using fast::QuantPoint;
using fast::quantize_point;
using fast::same_curve;

// A step of a walk: the kind of curve and the two points that define it.
struct StepRecord {
    std::uint8_t type{0};  // 0 = line, 1 = circle
    int first{0};
    int second{0};
};

struct WalkState {
    std::vector<GPoint> points;
    std::vector<GCurve> curves;
    std::vector<StepRecord> steps;
    int cost{0};
};

// The key of a state is the set of its curves, not their order, so two ways of
// drawing the same curves are one state.
struct StateKey {
    std::uint64_t first{};
    std::uint64_t second{};
    bool operator==(const StateKey& other) const {
        return first == other.first && second == other.second;
    }
};

struct StateKeyHash {
    std::size_t operator()(const StateKey& key) const noexcept {
        return static_cast<std::size_t>(key.first ^ (key.second * 0x9e3779b97f4a7c15ULL));
    }
};

std::uint64_t quantize_scalar(double value) {
    return static_cast<std::uint64_t>(static_cast<std::int64_t>(std::nearbyint(value * 1e8)));
}

StateKey state_key(const WalkState& state) {
    std::uint64_t sum = 0xcbf29ce484222325ULL;
    std::uint64_t mixed = 0x9e3779b97f4a7c15ULL;
    for (const auto& curve : state.curves) {
        std::uint64_t hash = 1469598103934665603ULL;
        auto feed = [&hash](std::uint64_t value) {
            hash ^= value;
            hash *= 1099511628211ULL;
        };
        feed(static_cast<std::uint64_t>(curve.type) + 1);
        feed(quantize_scalar(curve.a));
        feed(quantize_scalar(curve.b));
        feed(quantize_scalar(curve.c));
        sum += hash;
        mixed ^= hash * 0x9e3779b97f4a7c15ULL;
    }
    return StateKey{sum, mixed};
}

// How near the curve comes to the point.  The normal of a line is a unit
// vector, so the value is the distance itself.
double distance_to_point(const GCurve& curve, GPoint p) {
    if (curve.type == 0) return std::abs(curve.a * p.x + curve.b * p.y + curve.c);
    return std::abs(std::hypot(p.x - curve.a, p.y - curve.b) - curve.c);
}

// Whether a point of a figure is the target.
bool owns_target(const WalkState& state, const GPoint& target, double tolerance) {
    for (const auto& point : state.points)
        if (std::hypot(point.x - target.x, point.y - target.y) <= tolerance) return true;
    return false;
}

// Whether the figure holds a curve that passes through the target.
//
// A figure owns the target only if two of its curves pass through it, so this
// is the property that the last two levels of the walk are filtered with: the
// level before the last one has to hold a curve through the target, and the
// last curve has to pass through the target as well.
//
// The curves of the root are part of the figure and are counted here.  A given
// curve is not something the construction draws, but the rule is about the
// figure and not about what is drawn: a figure whose only curve through the
// target is a given one owns the target as soon as one more curve is drawn
// through it.  Leaving the given curves out would cut such figures away, and
// the walk would miss them: the secant task is a case in point, for the
// circle that is given passes through the point the secant has to meet.
bool holds_target_curve(const WalkState& state, const GPoint& target, double tolerance) {
    for (const auto& curve : state.curves)
        if (distance_to_point(curve, target) <= tolerance) return true;
    return false;
}

// How near the target is to the points of a figure.
double residual_of(const WalkState& state, const GPoint& target) {
    double best = std::numeric_limits<double>::infinity();
    for (const auto& point : state.points)
        best = std::min(best, std::hypot(point.x - target.x, point.y - target.y));
    return best;
}

// The size of a state, as far as the memory budget can see it: two doubles for
// a point, a curve is a type and three doubles, a step is a type and two
// integers.
std::size_t bytes_of(const WalkState& state) {
    return state.points.size() * 16 + state.curves.size() * 32 + state.steps.size() * 12;
}

bool extend_state(const WalkState& parent, const GCurve& curve, const StepRecord& record,
                  WalkState& child, std::size_t max_points, bool& cut) {
    for (const auto& existing : parent.curves)
        if (same_curve(existing, curve)) return false;
    child.points = parent.points;
    child.curves = parent.curves;
    child.steps = parent.steps;
    child.cost = parent.cost + 1;

    std::vector<QuantPoint> known;
    known.reserve(parent.points.size() + 8);
    for (const auto& point : parent.points) known.push_back(quantize_point(point.x, point.y));
    for (const auto& old : parent.curves) {
        GPoint hits[2];
        const int count = fast_intersections(curve, old, hits);
        for (int i = 0; i < count; ++i) {
            const QuantPoint key = quantize_point(hits[i].x, hits[i].y);
            if (std::find(known.begin(), known.end(), key) != known.end()) continue;
            known.push_back(key);
            if (child.points.size() + 1 > max_points) {
                cut = true;
                return false;
            }
            child.points.push_back(hits[i]);
        }
    }
    child.curves.push_back(curve);
    child.steps.push_back(record);
    return true;
}

std::string name_of(int index, const std::vector<std::string>& names) {
    if (index < static_cast<int>(names.size())) return names[index];
    // The name of a point is its index, so that the names of the report and the
    // indices of the points of the figure are the same.
    return "P" + std::to_string(index);
}

std::string describe_step(const StepRecord& record, const std::vector<std::string>& names) {
    const std::string first = name_of(record.first, names);
    const std::string second = name_of(record.second, names);
    return record.type == 0 ? "line(" + first + ", " + second + ")"
                            : "circle(" + first + "; " + second + ")";
}

std::vector<std::string> describe_steps(const WalkState& state,
                                        const std::vector<std::string>& names) {
    std::vector<std::string> steps;
    for (const auto& record : state.steps) steps.push_back(describe_step(record, names));
    return steps;
}

// The same steps with the coordinates of the two points that define them.
std::vector<StepGeometry> geometry_of(const WalkState& state,
                                       const std::vector<std::string>& names) {
    std::vector<StepGeometry> geometry;
    geometry.reserve(state.steps.size());
    for (const auto& record : state.steps) {
        geometry.push_back({describe_step(record, names), state.points[record.first],
                            state.points[record.second]});
    }
    return geometry;
}

}  // namespace

Root make_root(const std::vector<GPoint>& points, const std::vector<GCurve>& curves,
               const std::vector<std::string>& names) {
    Root root;
    root.points = points;
    root.curves = curves;
    root.names = names;

    // The points that the curves of the root leave behind are points of the
    // figure, and the walk needs them.
    std::vector<QuantPoint> known;
    known.reserve(points.size() + 8);
    for (const auto& point : points) known.push_back(quantize_point(point.x, point.y));
    for (std::size_t i = 0; i < curves.size(); ++i) {
        for (std::size_t j = i + 1; j < curves.size(); ++j) {
            GPoint hits[2];
            const int count = fast_intersections(curves[i], curves[j], hits);
            for (int k = 0; k < count; ++k) {
                const QuantPoint key = quantize_point(hits[k].x, hits[k].y);
                if (std::find(known.begin(), known.end(), key) != known.end()) continue;
                known.push_back(key);
                root.points.push_back(hits[k]);
            }
        }
    }
    return root;
}

Report run(const Root& root, GPoint target, const Config& config) {
    const auto started = std::chrono::steady_clock::now();
    Report report;

    const int root_cost = static_cast<int>(root.curves.size());
    const int extra = config.max_cost - root_cost;
    if (extra <= 0) {
        report.note = "the root already has every curve the budget allows";
        return report;
    }

    // The root may own the target already.
    for (const auto& point : root.points) {
        if (std::hypot(point.x - target.x, point.y - target.y) <= config.confirm_tolerance) {
            report.found = true;
            report.cost = root_cost;
            report.residual = std::hypot(point.x - target.x, point.y - target.y);
            report.note = "the root owns the target";
            return report;
        }
    }

    WalkState initial;
    initial.points = root.points;
    initial.curves = root.curves;
    initial.cost = root_cost;
    std::vector<WalkState> layer{std::move(initial)};
    report.layer_sizes.push_back(1);

    std::unordered_set<StateKey, StateKeyHash> seen;
    seen.insert(state_key(layer.front()));

    const auto deadline = started + std::chrono::duration<double>(config.time_limit_seconds);
    bool out_of_time = false;
    bool out_of_budget = false;
    bool out_of_memory = false;

    // The number of bytes the layers that are held take, as far as the walk can
    // see it.
    std::size_t held = bytes_of(layer.front());
    report.peak_bytes = held;

    int layer_cost = root_cost;

    // The rule of the two curves needs a level below the last two, so it is
    // only used when there is one; with a budget of one curve over the root
    // the walk of one level is the plain one.
    const bool two_levels = config.target_curve_filter && extra >= 2;

    // The figures that own the target, by the set of their curves.
    //
    // A figure is the set of its curves and not the order they were drawn in,
    // the same way the key of a state is.  The last two levels of the walk can
    // reach the same figure in two ways -- the curve it draws last may be one
    // or the other of two curves -- and the count has to be of the figures and
    // not of the ways, or one figure is counted twice.
    std::unordered_set<StateKey, StateKeyHash> hit_keys;

    // The figure that is found, with the points that define its curves.
    const auto record_found = [&report, &root, &hit_keys, target, &config](const WalkState& state,
                                                                          int cost) {
        if (!hit_keys.insert(state_key(state)).second) return;
        ++report.all_hit_count;
        if (report.all_hit_count == 1 || cost < report.all_hit_cost)
            report.all_hit_cost = cost;
        if (config.collect_all_hits) {
            report.all_hits.push_back(describe_steps(state, root.names));
            report.all_hit_geometry.push_back(geometry_of(state, root.names));
            return;
        }
        report.found = true;
        report.cost = cost;
        report.steps = describe_steps(state, root.names);
        report.geometry = root.geometry;
        const auto more = geometry_of(state, root.names);
        report.geometry.insert(report.geometry.end(), more.begin(), more.end());
        report.found_points = state.points;
        report.residual = residual_of(state, target);
    };

    for (int depth = 1; depth < extra - (two_levels ? 1 : 0) && !report.found; ++depth) {
        std::vector<WalkState> next;
        for (std::size_t index = 0; index < layer.size(); ++index) {
            if (std::chrono::steady_clock::now() > deadline) {
                out_of_time = true;
                break;
            }
            if (report.states >= config.max_states) {
                out_of_budget = true;
                break;
            }
            // The parent is moved out of the layer and released when it is
            // done with, so that the two layers do not both sit in memory.
            WalkState state = std::move(layer[index]);
            const std::size_t points = state.points.size();
            for (std::size_t i = 0; i < points; ++i) {
                for (std::size_t j = i + 1; j < points; ++j) {
                    const GPoint first = state.points[i];
                    const GPoint second = state.points[j];
                    const GCurve curves[3] = {make_line(first, second),
                                              make_circle(first, second),
                                              make_circle(second, first)};
                    const StepRecord records[3] = {{0, static_cast<int>(i), static_cast<int>(j)},
                                                   {1, static_cast<int>(i), static_cast<int>(j)},
                                                   {1, static_cast<int>(j), static_cast<int>(i)}};
                    for (int c = 0; c < 3; ++c) {
                        WalkState child;
                        bool cut = false;
                        if (!extend_state(state, curves[c], records[c], child,
                                          config.max_points, cut)) {
                            if (cut) ++report.states_cut;
                            continue;
                        }
                        ++report.candidates;
                        const StateKey key = state_key(child);
                        if (!seen.insert(key).second) continue;
                        ++report.states;
                        for (const auto& point : child.points) {
                            if (std::hypot(point.x - target.x, point.y - target.y) <=
                                config.confirm_tolerance) {
                                record_found(child, root_cost + depth);
                                break;
                            }
                        }
                        if (report.found) break;
                        next.push_back(std::move(child));
                    }
                    if (report.found) break;
                }
                if (report.found) break;
            }
            if (report.found) break;
        }
        if (report.found) break;
        report.layer_sizes.push_back(next.size());
        layer = std::move(next);
        layer_cost = root_cost + depth;
        // The layer that is held is the one that is built; the one before it
        // has been released already.
        held = 0;
        for (const auto& state : layer) held += bytes_of(state);
        report.peak_bytes = std::max(report.peak_bytes, held);
        if (config.verbose) {
            std::cout << "cost " << layer_cost << ": " << layer.size() << " states, "
                      << report.states << " states in all, "
                      << (held / (1024 * 1024)) << " MB, "
                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
                             .count()
                      << "s" << std::endl;
        }
        if (out_of_time || out_of_budget) break;
        if (held > config.memory_bytes) {
            out_of_memory = true;
            break;
        }
    }

    // The second phase: one more curve on the last layer.  Nothing is
    // materialised here, so this costs no memory beyond the layer it walks
    // over, and it runs even when the state budget stopped the first phase.
    //
    // Only the curves that pass near the target are tried.  That is not a cut
    // of the search: a figure that owns the target point has to make it with
    // the curve it draws last, so that curve goes through the target.
    if (!report.found && !two_levels) {
        for (std::size_t index = 0; index < layer.size(); ++index) {
            if (std::chrono::steady_clock::now() > deadline) {
                out_of_time = true;
                break;
            }
            WalkState state = std::move(layer[index]);
            const std::size_t points = state.points.size();
            for (std::size_t i = 0; i < points; ++i) {
                for (std::size_t j = i + 1; j < points; ++j) {
                    const GPoint first = state.points[i];
                    const GPoint second = state.points[j];
                    const GCurve curves[3] = {make_line(first, second),
                                              make_circle(first, second),
                                              make_circle(second, first)};
                    const StepRecord records[3] = {{0, static_cast<int>(i), static_cast<int>(j)},
                                                   {1, static_cast<int>(i), static_cast<int>(j)},
                                                   {1, static_cast<int>(j), static_cast<int>(i)}};
                    for (int c = 0; c < 3; ++c) {
                        if (distance_to_point(curves[c], target) > config.align_tolerance) continue;
                        WalkState child;
                        bool cut = false;
                        if (!extend_state(state, curves[c], records[c], child,
                                          config.max_points, cut)) {
                            if (cut) ++report.states_cut;
                            continue;
                        }
                        for (const auto& point : child.points) {
                            if (std::hypot(point.x - target.x, point.y - target.y) >
                                config.confirm_tolerance)
                                continue;
                            if (config.collect_all_hits) {
                                record_found(child, layer_cost + 1);
                                break;
                            }
                            report.found = true;
                            report.cost = layer_cost + 1;
                            report.residual = std::hypot(point.x - target.x, point.y - target.y);
                            report.steps = describe_steps(child, root.names);
                            report.geometry = root.geometry;
                            const auto extra = geometry_of(child, root.names);
                            report.geometry.insert(report.geometry.end(), extra.begin(), extra.end());
                            report.found_points = child.points;
                            break;
                        }
                        if (report.found) break;
                    }
                    if (report.found) break;
                }
                if (report.found) break;
            }
            if (report.found) break;
        }
    }

    // The last two levels of the walk, filtered with the rule of the two
    // curves.
    //
    // A figure owns the target only if two of its curves pass through it, so
    // the level before the last one has to hold a curve that passes through the
    // target, and the last curve has to pass through the target as well.
    //
    // Neither of the two levels is materialised.  A state of the level before
    // the last one holds a curve through the target only if the state below it
    // held one, or if the curve that was added passes through the target, so
    // the states that have to be looked at are the children of the states that
    // hold such a curve, and the children of the curves through the target
    // themselves.  The layer of the last level is never in memory, which is
    // what lets a walk of seven curves run on a machine that cannot hold the
    // layer of six.
    if (!report.found && two_levels) {
        for (std::size_t index = 0; index < layer.size(); ++index) {
            if (std::chrono::steady_clock::now() > deadline) {
                out_of_time = true;
                break;
            }
            WalkState state = std::move(layer[index]);
            const bool state_holds = holds_target_curve(state, target, config.align_tolerance);
            if (state_holds) {
                ++report.target_curve_states;
                if (config.verbose && config.print_target_figures) {
                    // The figures that hold a curve through the target are the
                    // whole of the work of the last two levels, so they are
                    // printed: a reader can look at them.
                    std::cout << "a figure of " << layer_cost << " curves with a curve through the target:";
                    for (const auto& step : describe_steps(state, root.names)) std::cout << " " << step;
                    std::cout << "\n";
                }
            }
            const std::size_t points = state.points.size();
            for (std::size_t i = 0; i < points && !report.found; ++i) {
                for (std::size_t j = i + 1; j < points && !report.found; ++j) {
                    const GPoint first = state.points[i];
                    const GPoint second = state.points[j];
                    const GCurve curves[3] = {make_line(first, second),
                                              make_circle(first, second),
                                              make_circle(second, first)};
                    const StepRecord records[3] = {{0, static_cast<int>(i), static_cast<int>(j)},
                                                   {1, static_cast<int>(i), static_cast<int>(j)},
                                                   {1, static_cast<int>(j), static_cast<int>(i)}};
                    for (int c = 0; c < 3 && !report.found; ++c) {
                        const bool passes =
                            distance_to_point(curves[c], target) <= config.align_tolerance;
                        // The child holds a curve through the target when the
                        // state below it held one, or when the curve that is
                        // added passes through the target.  A child that holds
                        // none of them can never own the target, whatever is
                        // drawn on it, so it is not looked at.
                        if (!state_holds && !passes) continue;
                        WalkState child;
                        bool cut = false;
                        if (!extend_state(state, curves[c], records[c], child,
                                          config.max_points, cut)) {
                            if (cut) ++report.states_cut;
                            continue;
                        }
                        ++report.candidates;
                        ++report.target_curve_children;
                        // The child may own the target already: it holds a
                        // curve through the target, and the curve that was
                        // added is one too.
                        if (owns_target(child, target, config.confirm_tolerance)) {
                            record_found(child, layer_cost + 1);
                            if (!config.collect_all_hits) break;
                            continue;
                        }
                        // The last curve.  It has to pass through the target,
                        // and the child already holds a curve that passes
                        // through it, so the figure owns the target.
                        const std::size_t child_points = child.points.size();
                        for (std::size_t u = 0; u < child_points && !report.found; ++u) {
                            for (std::size_t v = u + 1; v < child_points && !report.found; ++v) {
                                const GPoint one = child.points[u];
                                const GPoint two = child.points[v];
                                const GCurve last[3] = {make_line(one, two),
                                                        make_circle(one, two),
                                                        make_circle(two, one)};
                                const StepRecord last_records[3] = {
                                    {0, static_cast<int>(u), static_cast<int>(v)},
                                    {1, static_cast<int>(u), static_cast<int>(v)},
                                    {1, static_cast<int>(v), static_cast<int>(u)}};
                                for (int k = 0; k < 3; ++k) {
                                    if (distance_to_point(last[k], target) >
                                        config.align_tolerance)
                                        continue;
                                    WalkState grandchild;
                                    bool cut2 = false;
                                    if (!extend_state(child, last[k], last_records[k], grandchild,
                                                      config.max_points, cut2)) {
                                        if (cut2) ++report.states_cut;
                                        continue;
                                    }
                                    ++report.candidates;
                                    if (!owns_target(grandchild, target, config.confirm_tolerance))
                                        continue;
                                    record_found(grandchild, layer_cost + 2);
                                    if (!config.collect_all_hits) break;
                                    continue;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    report.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    report.budget_hit = out_of_budget;
    report.memory_hit = out_of_memory;
    report.time_limit_hit = out_of_time;
    if (report.found) {
        std::ostringstream note;
        note << "the target is made by " << report.cost << " curves";
        report.note = note.str();
    } else if (config.collect_all_hits && report.all_hit_count) {
        std::ostringstream note;
        note << report.all_hit_count << " figures own the target, and they were all collected";
        report.note = note.str();
    } else if (out_of_budget) {
        report.note =
            "the walk stopped at the state budget, so the layers it did not reach are not covered";
    } else if (out_of_memory) {
        report.note =
            "the walk stopped at the memory budget, so the layers it did not reach are not covered";
    } else if (out_of_time) {
        report.note =
            "the walk ran out of time, so the layers it did not reach are not covered";
    } else if (report.states_cut) {
        report.note = "the point budget cut some states, so the walk is not complete";
    } else {
        report.note = "no construction of the target within the budget";
        if (two_levels)
            report.note +=
                "; the last two levels were filtered with the rule of the two curves, "
                "which is a rule of the search and not a cut of it";
    }
    return report;
}

}  // namespace neo::engine::exhaustive_walk
