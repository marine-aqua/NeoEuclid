#include "neo/engine/goal_search.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <array>
#include <functional>
#include <string>
#include <unordered_set>
#include <vector>

namespace neo::goal {
namespace {

using fast::fast_intersections;

// The counters of the walk.  They are written down when the environment asks
// for them, and the walk may be told to leave a family of operations out, so
// that the price of one family can be told from the price of the others.
struct Counters {
    // one for each family of operations, in the order of `Family`
    unsigned long long candidates[8]{};
    // the calls of the kernel of the crossings
    unsigned long long intersections{};
    // the figures that were built
    unsigned long long extensions{};
    // the places that were copied into them
    unsigned long long copies{};
    // the curves that the figure had already
    unsigned long long duplicates{};
    // the crossings that were taken
    unsigned long long crossings{};
    // the figures the walk has been in already
    unsigned long long duplicate_states{};
    // the candidates that were made at the last level, where the filter of the
    // search decides, and the ones that were made anywhere else
    unsigned long long candidates_last{};
    unsigned long long candidates_other{};
    unsigned long long extensions_last{};
    // the points and the curves of the widest figure
    unsigned long long widest_points{};
    unsigned long long widest_curves{};
    // the figures that own the goal and leave a place of the construction
    // unused, which a shortest construction never does
    unsigned long long unused_places{};
};

enum Family {
    FamilyLine = 0,
    FamilyCircle,
    FamilyBisector,
    FamilyTransfer,
    FamilyPerpendicular,
    FamilyParallel,
    FamilyAngle,
    FamilyCrossing,
    FamilyCount
};

Counters counters;

// The families the walk is told to leave out.  The name of a family is
// looked for as a whole word, for `perpendicular` is a part of
// `perpendicular-bisector`.
int skip_mask() {
    static const int mask = [] {
        const char* text = std::getenv("NEO_GOAL_SKIP");
        if (text == nullptr) return 0;
        const std::string names(text);
        static const char* const words[FamilyCount] = {
            "line", "circle", "bisector", "transfer",
            "perpendicular", "parallel", "angle", "crossing"};
        int mask = 0;
        for (int family = 0; family < FamilyCount; ++family) {
            const std::string word = words[family];
            std::size_t at = names.find(word);
            while (at != std::string::npos) {
                const bool left = at == 0 || names[at - 1] == ' ' || names[at - 1] == ',';
                const std::size_t end = at + word.size();
                const bool right = end == names.size() || names[end] == ' ' || names[end] == ',';
                if (left && right) {
                    mask |= 1 << family;
                    break;
                }
                at = names.find(word, at + 1);
            }
        }
        return mask;
    }();
    return mask;
}

bool skipped(Family family, int mask) {
    return (mask & (1 << family)) != 0;
}

// Whether the filter of the last step is put before the candidate is built.
// The filter is the same either way; the question is only whether the
// candidates that cannot pass it are built at all.  It is on by default, for
// the candidates that the filter throws away are three quarters of the
// candidates of a walk, and the environment variable turns it off to measure
// what it is worth.
bool filter_first_enabled() {
    static const bool on = [] {
        const char* text = std::getenv("NEO_GOAL_FILTER_FIRST");
        return text == nullptr || std::string(text) != "0";
    }();
    return on;
}

// Whether the walk remembers the figures it has already been in.  Two ways
// of building the same figure lead to the same continuations, so the second
// one of them can be left out and the walk becomes a walk of a graph and not
// of a tree.  It is on by default, and the environment variable turns it off
// to measure what it is worth.
bool memo_enabled() {
    static const bool on = [] {
        const char* text = std::getenv("NEO_GOAL_MEMO");
        return text == nullptr || std::string(text) != "0";
    }();
    return on;
}

// The hash of a run of words.  The run is the place of a figure, and the hash
// is the way the figures are told apart.  The words are folded together with
// a multiply and a xor, and the result is put through the finaliser of
// SplitMix, for the words of a place are the places of a figure and they are
// not spread out the way a hash wants them.
//
// A hash that tells two places apart by accident costs nothing but a
// comparison of the two places, for the places are compared and not only
// their hashes.
inline std::uint64_t hash_words(const std::uint64_t* words, std::size_t count) {
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    for (std::size_t i = 0; i < count; ++i)
        hash = (hash ^ words[i]) * 0x9e3779b97f4a7c15ULL;
    hash ^= hash >> 30;
    hash *= 0xbf58476d1ce4e5b9ULL;
    hash ^= hash >> 27;
    hash *= 0x94d049bb133111ebULL;
    hash ^= hash >> 31;
    return hash;
}

// Is the place of `Width` words at `left` before the one at `right`?  The
// places are compared as unsigned words.  The order they are put in does not
// matter, for the key only has to tell the same figures apart.
template <std::size_t Width>
inline bool place_before(const std::uint64_t* left, const std::uint64_t* right) {
    for (std::size_t k = 0; k < Width; ++k) {
        if (left[k] != right[k]) return left[k] < right[k];
    }
    return false;
}

// The place of a point on the grid of the engine: the two numbers of the
// grid.
inline void quantise_point(double x, double y, std::uint64_t* place) {
    const fast::QuantPoint quantised = fast::quantize_point(x, y);
    place[0] = static_cast<std::uint64_t>(quantised.x);
    place[1] = static_cast<std::uint64_t>(quantised.y);
}

// The place of a curve on the grid of the engine: the kind of it, the normal
// of it and the offset of it.
inline void quantise_curve(const fast::GCurve& curve, std::uint64_t* place) {
    const fast::QuantPoint normal = fast::quantize_point(curve.a, curve.b);
    const fast::QuantPoint offset = fast::quantize_point(curve.c, 0.0);
    place[0] = curve.type;
    place[1] = static_cast<std::uint64_t>(normal.x);
    place[2] = static_cast<std::uint64_t>(normal.y);
    place[3] = static_cast<std::uint64_t>(offset.x);
}

// The place is put among the places that are already there, in the order of
// their places.  The places that are after it are moved to make room for it,
// and the place it was put at is given back, for the walk takes the place
// away again when the branch is done.
template <std::size_t Width>
std::size_t insert_place(std::vector<std::uint64_t>* places, const std::uint64_t* place) {
    const std::size_t count = places->size() / Width;
    std::size_t at = count;
    while (at > 0 && place_before<Width>(place, places->data() + (at - 1) * Width)) --at;
    places->insert(places->begin() + at * Width, place, place + Width);
    return at;
}

// Is the place among the places that are already there?  The places are in
// the order of their places, so the place is looked for with a search.
template <std::size_t Width>
bool find_place(const std::vector<std::uint64_t>& places, const std::uint64_t* place) {
    std::size_t low = 0;
    std::size_t high = places.size() / Width;
    while (low < high) {
        const std::size_t mid = low + (high - low) / 2;
        const std::uint64_t* const other = places.data() + mid * Width;
        if (place_before<Width>(place, other)) high = mid;
        else if (place_before<Width>(other, place)) low = mid + 1;
        else return true;
    }
    return false;
}

// The places of a figure: the points and the curves it has, in the order of
// their places, and the hash of them.
//
// The walk of the All model adds one place at a time and takes it away again
// when the branch is done, and the places that are already there do not move.
// So the places are kept in the order of their places, a curve is looked for
// with a search and not with a walk of the curves of the figure, and the hash
// of the figure is the sum of the hashes of its places.
class Figure {
public:
    // The figure that the walk starts from.
    void load(const State& state) {
        points_.clear();
        curves_.clear();
        point_at_.clear();
        curve_at_.clear();
        point_hash_ = 0;
        curve_hash_ = 0;
        // The givens are free: they are not operations of the construction,
        // so they are not places that a step is spent on and they do not have
        // to be on the way to the goal.
        for (const auto& point : state.points) {
            std::uint64_t place[2];
            quantise_point(point.x, point.y, place);
            point_at_.push_back(static_cast<std::uint32_t>(insert_place<2>(&points_, place)));
            point_hash_ += hash_words(place, 2);
        }
        for (const auto& curve : state.curves) {
            std::uint64_t place[4];
            quantise_curve(curve, place);
            curve_at_.push_back(static_cast<std::uint32_t>(insert_place<4>(&curves_, place)));
            curve_hash_ += hash_words(place, 4);
        }
    }

    // Is the curve already in the figure?
    bool has_curve(const std::uint64_t* place) const {
        return find_place<4>(curves_, place);
    }

    // A curve is added to the figure, and taken away again when the branch of
    // the walk is done.
    void push_curve(const std::uint64_t* place) {
        curve_at_.push_back(static_cast<std::uint32_t>(insert_place<4>(&curves_, place)));
        curve_hash_ += hash_words(place, 4);
    }
    void pop_curve() {
        const std::size_t at = curve_at_.back();
        curve_at_.pop_back();
        curve_hash_ -= hash_words(curves_.data() + at * 4, 4);
        curves_.erase(curves_.begin() + at * 4, curves_.begin() + at * 4 + 4);
    }

    // A point is added to the figure, and taken away again when the branch of
    // the walk is done.
    void push_point(const std::uint64_t* place) {
        point_at_.push_back(static_cast<std::uint32_t>(insert_place<2>(&points_, place)));
        point_hash_ += hash_words(place, 2);
    }
    void pop_point() {
        const std::size_t at = point_at_.back();
        point_at_.pop_back();
        point_hash_ -= hash_words(points_.data() + at * 2, 2);
        points_.erase(points_.begin() + at * 2, points_.begin() + at * 2 + 2);
    }

    std::size_t point_count() const { return points_.size() / 2; }
    std::size_t curve_count() const { return curves_.size() / 4; }
    const std::uint64_t* point_words() const { return points_.data(); }
    const std::uint64_t* curve_words() const { return curves_.data(); }

    // The hash of the figure.  It is the sum of the hashes of the places, and
    // the order of the places does not matter, so the hash is the same for
    // every way of building the figure.
    std::uint64_t hash() const {
        std::uint64_t hash = 0xcbf29ce484222325ULL;
        hash = (hash ^ point_count()) * 0x9e3779b97f4a7c15ULL;
        hash = (hash ^ curve_count()) * 0x9e3779b97f4a7c15ULL;
        hash = (hash ^ point_hash_) * 0x9e3779b97f4a7c15ULL;
        hash = (hash ^ curve_hash_) * 0x9e3779b97f4a7c15ULL;
        hash ^= hash >> 30;
        hash *= 0xbf58476d1ce4e5b9ULL;
        hash ^= hash >> 27;
        hash *= 0x94d049bb133111ebULL;
        hash ^= hash >> 31;
        return hash;
    }

private:
    std::vector<std::uint64_t> points_;
    std::vector<std::uint64_t> curves_;
    std::vector<std::uint32_t> point_at_;
    std::vector<std::uint32_t> curve_at_;
    std::uint64_t point_hash_{0};
    std::uint64_t curve_hash_{0};
};

// The figures the walk has already been in.  A figure is told apart by its
// key: the points and the curves it has, in the order of their places.
//
// The key is a run of words in the arena of the keys, and a figure is kept as
// the place of its key in that arena.  So no key is ever copied and no string
// is ever made.
class FigureSet {
public:
    struct Key {
        std::uint32_t offset;
        std::uint32_t length;
        std::uint64_t hash;
    };

    // The figures are emptied at every level of the deepening, for the walk of
    // a level is a walk of its own.
    void clear() {
        set_.clear();
        arena_.clear();
    }

    // The key of a figure is written into the arena.  The places of the figure
    // are already in the order of their places, so the key is a copy of them
    // and the hash of it is the hash the figure carries.
    Key build(const Figure& figure) {
        const std::size_t point_count = figure.point_count();
        const std::size_t curve_count = figure.curve_count();
        const std::size_t start = arena_.size();
        arena_.resize(start + 2 + 2 * point_count + 4 * curve_count);
        std::uint64_t* const words = arena_.data() + start;

        words[0] = point_count;
        if (point_count != 0)
            std::memcpy(words + 1, figure.point_words(), 2 * point_count * sizeof(std::uint64_t));
        words[1 + 2 * point_count] = curve_count;
        if (curve_count != 0)
            std::memcpy(words + 2 + 2 * point_count, figure.curve_words(),
                        4 * curve_count * sizeof(std::uint64_t));

        const std::size_t length = 2 + 2 * point_count + 4 * curve_count;
        return Key{static_cast<std::uint32_t>(start), static_cast<std::uint32_t>(length),
                   figure.hash()};
    }

    // The figure of a key is added, and it is told whether it was already
    // there.  A key that was already there is thrown away with its words.
    bool insert(const Key& key) {
        if (set_.insert(key).second) return true;
        arena_.resize(key.offset);
        return false;
    }

private:
    struct Hasher {
        std::size_t operator()(const Key& key) const noexcept {
            return static_cast<std::size_t>(key.hash);
        }
    };

    // Two figures are the same when their keys are the same, and not when
    // their hashes are: a hash that tells two keys apart by accident is
    // caught here, for the words of the two keys are compared.
    struct Equal {
        const std::vector<std::uint64_t>* arena;
        bool operator()(const Key& left, const Key& right) const {
            if (left.hash != right.hash) return false;
            if (left.length != right.length) return false;
            const std::uint64_t* const words = arena->data();
            for (std::uint32_t i = 0; i < left.length; ++i)
                if (words[left.offset + i] != words[right.offset + i]) return false;
            return true;
        }
    };

    std::vector<std::uint64_t> arena_;
    std::unordered_set<Key, Hasher, Equal> set_;

public:
    // The set of the figures is told apart by the arena of the keys, and the
    // arena is a member of this one.
    FigureSet() : set_(0, Hasher{}, Equal{&arena_}) {}
};

// Can a family of operations finish a figure whose goal is the target?  A
// curve operation gives the figure no point, so it can only finish a figure
// whose goal is a line, and only the families that draw a line can draw it.
bool family_finishes(Family family, const Target& target) {
    if (target.kind != Target::Kind::Line) return false;
    return family != FamilyCircle && family != FamilyTransfer;
}
using fast::GCurve;
using fast::GPoint;
using fast::make_circle;
using fast::make_line;
using fast::QuantPoint;
using fast::quantize_point;
using fast::same_curve;

// Does the curve meet the circle about the centre of the given radius?  The
// distance of the two centres is compared without a square root, for the same
// reason the circle of the fast kernel does it: `d` is in `[r - c - t, r + c + t]`
// when `d^2` is in the squares of those bounds.
bool curve_meets_circle(const GCurve& curve, GPoint centre, double radius, double tolerance) {
    if (curve.type == 0) {
        const double distance = std::abs(curve.a * centre.x + curve.b * centre.y + curve.c);
        return distance - radius <= tolerance;
    }
    const double dx = curve.a - centre.x;
    const double dy = curve.b - centre.y;
    const double squared = dx * dx + dy * dy;
    const double outer = radius + curve.c + tolerance;
    if (squared > outer * outer) return false;
    const double inner = radius - curve.c - tolerance;
    if (inner <= 0.0) return true;
    return squared >= inner * inner;
}

// Is the curve the line of the target?  Two lines are the same when their
// unit normals and their offsets agree, and the sign convention of the
// figures fixes the sign of the normal.
bool curve_is_target(const GCurve& curve, const Target& target, double tolerance) {
    if (target.kind != Target::Kind::Line) return false;
    if (curve.type != 0) return false;
    return std::abs(curve.a - target.line_a) <= tolerance &&
           std::abs(curve.b - target.line_b) <= tolerance &&
           std::abs(curve.c - target.line_c) <= tolerance;
}

// Can the curve come near the target at all?  A curve that cannot cannot put
// the target into the figure, so the search never builds it.  This is what
// keeps the walk affordable: for a point target almost nothing passes.
//
// A curve that carries the second half of the goal, the point the figure has
// to own as well, is allowed through as well: it may be the curve that makes
// the goal reachable.
bool curve_reaches_target(const GCurve& curve, const Target& target, double tolerance) {
    if (target.kind == Target::Kind::Line) return curve_is_target(curve, target, tolerance);
    if (target.also_own && fast::gcurve_contains(curve, target.also, tolerance)) return true;
    if (target.kind == Target::Kind::Point) {
        if (fast::gcurve_contains(curve, target.point, tolerance)) return true;
        // The circle is a second way in: a curve that meets it may leave behind
        // a point from which the goal is one curve away.
        if (target.radius > 0.0)
            return curve_meets_circle(curve, target.centre, target.radius, tolerance);
        return false;
    }
    return curve_meets_circle(curve, target.centre, target.radius, tolerance);
}

// Does the crossing put the target into the figure?  In the All model a
// crossing is an operation of its own, and the only way a crossing can finish
// a figure is by being the target itself.
bool crossing_reaches_target(const GPoint& point, const Target& target, double tolerance) {
    if (target.kind == Target::Kind::Point) {
        const double dx = point.x - target.point.x;
        const double dy = point.y - target.point.y;
        return dx * dx + dy * dy <= tolerance * tolerance;
    }
    if (target.kind == Target::Kind::Circle) {
        const double dx = point.x - target.centre.x;
        const double dy = point.y - target.centre.y;
        const double away = std::sqrt(dx * dx + dy * dy) - target.radius;
        return std::abs(away) <= tolerance;
    }
    return false;
}

// Does the curve hold the point?  The kernel of the fast geometry has the
// test: the distance of a point from a circle is compared without a square
// root, for the distance `d` and the radius `r` are near each other when
// `|d^2 - r^2|` is small, and `|d^2 - r^2|` is `|d - r| (d + r)`.
bool curve_holds(const GCurve& curve, GPoint point, double tolerance) {
    return fast::gcurve_contains(curve, point, tolerance);
}

// Is the point the target?  The distance is compared without a square root:
// the test is the same one, for the tolerance is not negative.
bool point_is_target(const GPoint& point, const Target& target, double tolerance) {
    if (target.kind == Target::Kind::Point) {
        const double dx = point.x - target.point.x;
        const double dy = point.y - target.point.y;
        return dx * dx + dy * dy <= tolerance * tolerance;
    }
    const double dx = point.x - target.centre.x;
    const double dy = point.y - target.centre.y;
    const double away = std::sqrt(dx * dx + dy * dy) - target.radius;
    return std::abs(away) <= tolerance;
}

// Does the figure own the whole goal?  The second point of the goal has to be
// there as well.
bool owns_target(const State& state, const Target& target, double tolerance) {
    if (target.kind == Target::Kind::Line) {
        // The goal is a line, so the figure has to own that line among its
        // curves.
        for (const auto& curve : state.curves)
            if (curve_is_target(curve, target, tolerance)) return true;
        return false;
    }
    if (target.also_own) {
        bool found = false;
        for (const auto& point : state.points) {
            const double dx = point.x - target.also.x;
            const double dy = point.y - target.also.y;
            if (dx * dx + dy * dy <= tolerance * tolerance) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    for (const auto& point : state.points)
        if (point_is_target(point, target, tolerance)) return true;
    return false;
}

// The scratch buffer of the crossings.  The extension of a figure is not
// recursive, so one buffer is enough for the whole walk, and the walk asks
// for one at every figure it builds.
std::vector<QuantPoint>& scratch_known() {
    static thread_local std::vector<QuantPoint> buffer;
    return buffer;
}

// The figure after one more curve, in place.  The crossings of the new curve
// with the curves that are already there become points of the figure, and in
// the All model they do not: a crossing is an operation of its own there.
//
// The figure is changed in place and not copied, for the walk builds a figure
// for every candidate it tries and a copy of a figure is three allocations
// and three copies of a vector.  The caller remembers where the points and
// the curves were and puts them back when the branch is done.
bool add_curve(State& state, const GCurve& curve, std::size_t max_points, bool& cut,
               bool add_crossings) {
    for (const auto& existing : state.curves)
        if (same_curve(existing, curve)) {
            ++counters.duplicates;
            return false;
        }
    ++counters.extensions;
    counters.copies += state.points.size() + state.curves.size();
    counters.widest_points = std::max<unsigned long long>(counters.widest_points,
                                                          state.points.size());
    counters.widest_curves = std::max<unsigned long long>(counters.widest_curves,
                                                          state.curves.size());

    const std::size_t mark = state.points.size();
    if (add_crossings) {
        std::vector<QuantPoint>& known = scratch_known();
        known.clear();
        for (const auto& point : state.points) known.push_back(quantize_point(point.x, point.y));
        for (const auto& old : state.curves) {
            GPoint hits[2];
            ++counters.intersections;
            const int count = fast_intersections(curve, old, hits);
            for (int i = 0; i < count; ++i) {
                const fast::QuantPoint key = quantize_point(hits[i].x, hits[i].y);
                if (std::find(known.begin(), known.end(), key) != known.end()) continue;
                known.push_back(key);
                if (state.points.size() + 1 > max_points) {
                    // The figure is left as it was, for the branch is given up.
                    state.points.resize(mark);
                    cut = true;
                    return false;
                }
                state.points.push_back(hits[i]);
            }
        }
    }
    state.curves.push_back(curve);
    return true;
}

// The figure after one more curve, as a copy of another one.  This is the way
// a recipe is replayed and not the way the walk goes: the walk builds the
// figures in place.
bool extend_state(const State& parent, const GCurve& curve, const Step& step, State& child,
                  std::size_t max_points, bool& cut, bool add_crossings,
                  bool keep_steps = true) {
    child.points = parent.points;
    child.curves = parent.curves;
    if (keep_steps) child.steps = parent.steps;
    child.cost = parent.cost + 1;
    if (!add_curve(child, curve, max_points, cut, add_crossings)) return false;
    if (keep_steps) child.steps.push_back(step);
    return true;
}

// The curve of one step of a recipe.
//
// The places of a step are of two kinds, and every operation says which kind
// each of them is: the two atomic operations and the bisector of a segment
// name points, the perpendicular and the parallel name a point and a curve,
// the bisector of an angle names two curves.  The indices are therefore
// checked against the number of points or against the number of curves, as
// the operation says, and not all of them against the number of points.
bool curve_of_step(const State& state, const Step& step, GCurve& curve) {
    const int points = static_cast<int>(state.points.size());
    const int curves = static_cast<int>(state.curves.size());
    if (step.first < 0) return false;
    switch (step.op) {
        case Op::Legacy:
            if (step.first >= points || step.second < 0 || step.second >= points) return false;
            curve = step.type == 0 ? make_line(state.points[step.first], state.points[step.second])
                                   : make_circle(state.points[step.first], state.points[step.second]);
            return true;
        case Op::TransferCircle: {
            if (step.first >= points || step.second < 0 || step.second >= points) return false;
            if (step.third < 0 || step.third >= points) return false;
            curve = fast::make_transfer_circle(state.points[step.first], state.points[step.second],
                                               state.points[step.third]);
            return true;
        }
        case Op::PerpendicularBisector:
            if (step.first >= points || step.second < 0 || step.second >= points) return false;
            curve = fast::make_perpendicular_bisector(state.points[step.first],
                                                       state.points[step.second]);
            return true;
        case Op::PerpendicularThrough: {
            if (step.first >= points) return false;
            if (step.second < 0 || step.second >= curves) return false;
            curve = fast::make_perpendicular_through(state.points[step.first],
                                                      state.curves[step.second]);
            return true;
        }
        case Op::ParallelThrough: {
            if (step.first >= points) return false;
            if (step.second < 0 || step.second >= curves) return false;
            curve = fast::make_parallel_through(state.points[step.first], state.curves[step.second]);
            return true;
        }
        case Op::AngleBisector: {
            if (step.first >= curves || step.second < 0 || step.second >= curves) return false;
            GCurve bisectors[2];
            const int number = fast::make_angle_bisectors(state.curves[step.first],
                                                          state.curves[step.second], bisectors);
            if (step.which < 0 || step.which >= number) return false;
            curve = bisectors[step.which];
            return true;
        }
        case Op::Intersection:
            return false;
    }
    return false;
}

// The crossing of one step of a recipe: the point that the two curves of the
// figure have in common.
bool crossing_of_step(const State& state, const Step& step, GPoint& point) {
    if (step.op != Op::Intersection) return false;
    if (step.first < 0 || step.first >= static_cast<int>(state.curves.size())) return false;
    if (step.second < 0 || step.second >= static_cast<int>(state.curves.size())) return false;
    GPoint hits[2];
    const int count = fast_intersections(state.curves[step.first], state.curves[step.second], hits);
    if (step.which < 0 || step.which >= count) return false;
    point = hits[step.which];
    return true;
}

}  // namespace

// The crossings of the curves that the caller has already drawn.  The model
// says that a crossing costs nothing, so the crossings of the givens belong
// to the figure as well: the walk of the search starts from a figure whose
// curves meet, and it must be allowed to use the points of those crossings
// for free.  The search does this by itself when it starts, and it does not
// do it in the All model, where a crossing is an operation of its own.
void seed_crossings(State& state, std::size_t max_points) {
    std::vector<QuantPoint> known;
    known.reserve(state.points.size() + 4);
    for (const auto& point : state.points) known.push_back(quantize_point(point.x, point.y));
    for (std::size_t i = 0; i < state.curves.size(); ++i) {
        for (std::size_t j = i + 1; j < state.curves.size(); ++j) {
            GPoint hits[2];
            const int count = fast_intersections(state.curves[i], state.curves[j], hits);
            for (int k = 0; k < count; ++k) {
                const fast::QuantPoint key = quantize_point(hits[k].x, hits[k].y);
                if (std::find(known.begin(), known.end(), key) != known.end()) continue;
                known.push_back(key);
                if (state.points.size() + 1 > max_points) return;
                state.points.push_back(hits[k]);
            }
        }
    }
}

bool extend(State& state, const Step& step, std::size_t max_points,
            engine::Mode mode) {
    if (step.op == Op::Intersection) {
        GPoint point;
        if (!crossing_of_step(state, step, point)) return false;
        for (const auto& known : state.points) {
            const double dx = known.x - point.x;
            const double dy = known.y - point.y;
            if (dx * dx + dy * dy < 1e-18) return false;
        }
        if (state.points.size() + 1 > max_points) return false;
        state.points.push_back(point);
        state.steps.push_back(step);
        state.cost += 1;
        return true;
    }
    GCurve curve;
    if (!curve_of_step(state, step, curve)) return false;
    State child;
    bool cut = false;
    if (!extend_state(state, curve, step, child, max_points, cut,
                      mode == engine::Mode::Classic))
        return false;
    state = std::move(child);
    return true;
}

bool build(State& state, const std::vector<Step>& recipe, std::size_t max_points,
           engine::Mode mode) {
    for (const auto& step : recipe)
        if (!extend(state, step, max_points, mode)) return false;
    return true;
}

// The certificate of a construction: the places it makes, and the places
// among them that nothing uses.
//
// The recipe of a construction is a computation.  Every place of the figure
// is built from the places that came before it, and the goal is one place of
// the figure.  The certificate walks the recipe backwards from the goal and
// counts the places of the construction that are not on the way to it.
//
// A place that is not on the way to the goal is a place that is drawn for
// nothing: it can be left out of the recipe and the goal is still made, so a
// shortest construction has none of them.  The condition is a necessary one
// and not a sufficient one: a construction may have every one of its places
// on the way to the goal and still be shortened by replacing one of its
// steps by another one, so the certificate does not prove that the
// construction is the shortest one.  The search proves that by its deepening,
// and the certificate says that the construction it found has no step that
// can simply be dropped.
//
// In the All model every place of the figure is an operation of the
// construction and it costs a step, so every one of them has to be on the
// way to the goal.  In the classic model the crossings are free, so only the
// curves are counted: a crossing that nothing uses costs nothing, but the
// curve that was drawn to make it does.
Certificate certify(const State& start, const std::vector<Step>& recipe,
                    const Target& target, double tolerance, engine::Mode mode) {
    Certificate out;
    // The places of the figure, and the places each of them is built from.
    // The givens come first and they are free.
    std::vector<std::vector<int>> from;
    std::vector<int> point_id;
    std::vector<int> curve_id;
    std::vector<int> drawn_points;
    std::vector<int> drawn_curves;
    const auto new_place = [&from](std::vector<int> parents) {
        const int id = static_cast<int>(from.size());
        from.push_back(std::move(parents));
        return id;
    };
    for (std::size_t i = 0; i < start.points.size(); ++i) point_id.push_back(new_place({}));
    for (std::size_t i = 0; i < start.curves.size(); ++i) curve_id.push_back(new_place({}));

    // The figure of the replay, for the crossings of the classic model.
    State figure = start;
    std::vector<QuantPoint> known;
    known.reserve(start.points.size() + 8);
    for (const auto& point : figure.points) known.push_back(quantize_point(point.x, point.y));

    for (const auto& step : recipe) {
        if (step.op == Op::Intersection) {
            const int id = new_place({curve_id[step.first], curve_id[step.second]});
            point_id.push_back(id);
            drawn_points.push_back(id);
            GPoint point;
            if (!crossing_of_step(figure, step, point)) break;
            figure.points.push_back(point);
            known.push_back(quantize_point(point.x, point.y));
            continue;
        }
        GCurve curve;
        if (!curve_of_step(figure, step, curve)) break;
        std::vector<int> parents;
        switch (step.op) {
            case Op::Legacy:
            case Op::PerpendicularBisector:
                parents = {point_id[step.first], point_id[step.second]};
                break;
            case Op::TransferCircle:
                parents = {point_id[step.first], point_id[step.second], point_id[step.third]};
                break;
            case Op::PerpendicularThrough:
            case Op::ParallelThrough:
                parents = {point_id[step.first], curve_id[step.second]};
                break;
            case Op::AngleBisector:
                parents = {curve_id[step.first], curve_id[step.second]};
                break;
            case Op::Intersection:
                break;
        }
        const int id = new_place(std::move(parents));
        curve_id.push_back(id);
        drawn_curves.push_back(id);
        if (mode == engine::Mode::Classic) {
            // The crossings of the new curve with the curves that are already
            // there come along for nothing, and each of them is built from the
            // two curves that meet there.
            for (std::size_t c = 0; c < figure.curves.size(); ++c) {
                GPoint hits[2];
                const int count = fast_intersections(curve, figure.curves[c], hits);
                for (int k = 0; k < count; ++k) {
                    const QuantPoint key = quantize_point(hits[k].x, hits[k].y);
                    if (std::find(known.begin(), known.end(), key) != known.end()) continue;
                    known.push_back(key);
                    figure.points.push_back(hits[k]);
                    const int point = new_place({id, curve_id[c]});
                    point_id.push_back(point);
                    drawn_points.push_back(point);
                }
            }
        }
        figure.curves.push_back(curve);
    }

    // The place the goal is: the last place of the figure that owns it.
    int goal = -1;
    if (target.kind == Target::Kind::Line) {
        for (std::size_t c = 0; c < figure.curves.size(); ++c)
            if (curve_is_target(figure.curves[c], target, tolerance)) goal = curve_id[c];
    } else {
        for (std::size_t p = 0; p < figure.points.size(); ++p)
            if (point_is_target(figure.points[p], target, tolerance)) goal = point_id[p];
    }

    // The places on the way to the goal.
    std::vector<char> on_the_way(from.size(), 0);
    std::vector<int> stack;
    if (goal >= 0) {
        on_the_way[goal] = 1;
        stack.push_back(goal);
    }
    while (!stack.empty()) {
        const int id = stack.back();
        stack.pop_back();
        for (const int parent : from[id]) {
            if (on_the_way[parent]) continue;
            on_the_way[parent] = 1;
            stack.push_back(parent);
        }
    }

    // The places of the construction that nothing uses.  In the classic model
    // only the curves are counted, for the crossings are free.
    const bool classic = mode == engine::Mode::Classic;
    out.places = classic ? drawn_curves.size() : drawn_curves.size() + drawn_points.size();
    out.dead = 0;
    for (const int id : drawn_curves)
        if (!on_the_way[id]) ++out.dead;
    if (!classic) {
        for (const int id : drawn_points)
            if (!on_the_way[id]) ++out.dead;
    }
    out.ok = goal >= 0 && out.dead == 0;
    return out;
}

Report search(const State& start, const Target& target, const Config& config) {
    Report report;
    const auto started = std::chrono::steady_clock::now();
    const auto deadline = started + std::chrono::duration<double>(config.seconds);
    bool out_of_time = false;

    // The crossings of the curves that are already there are free points, so
    // the figure the walk starts from owns them.  In the All model they are
    // operations like any other, and the walk has to pay for them.
    State root = start;
    if (config.mode == engine::Mode::Classic) seed_crossings(root, config.max_points);
    // The givens of the walk, kept apart from the figure of the walk, for the
    // walk builds its figure in place.  The certificate of an answer is
    // written from the givens.
    const State givens = root;

    // The places of the figure the walk is in.  The walk of the All model
    // adds one place at a time and takes it away again when the branch is
    // done, and the places that are already there do not move.
    Figure figure;
    figure.load(root);

    // The figures the walk has already been in.  It is emptied at every level
    // of the deepening, for the walk of a level is a walk of its own.
    FigureSet figures;

    // The families the walk is told to leave out.  The mask is read once and
    // not at every candidate, for the walk asks about a family very often.
    const int skip = skip_mask();

    // Whether the walk remembers the figures it has already been in.  It is
    // read once and not at every node, for the walk asks about it very often.
    const bool memo = memo_enabled();

    std::vector<Step> path;

    // The walk of the classic model: a curve is drawn from two points of the
    // figure, and the crossings it makes with the curves that are already
    // there come along for nothing.
    //
    // The figure is built in place and not copied, and the points and the
    // curves that the branch adds are taken away again when the branch is
    // done.
    std::function<bool(State&, int)> walk_classic = [&](State& state, int depth) -> bool {
        ++report.states;
        // The clock is read here and not only between the levels: a level of
        // the walk can hold billions of nodes, and a search that only looks
        // at the clock between the levels would run for hours past its budget.
        if ((report.states & 0xFFFF) == 0 && std::chrono::steady_clock::now() > deadline) {
            out_of_time = true;
            return true;
        }
        if (owns_target(state, target, config.tolerance)) {
            if (config.collect_all && depth == 0) {
                // One more figure that owns the goal: the walk keeps going so
                // that every one of them is written down.
                report.hits.push_back(path);
                report.hit_points.push_back(state.points);
                report.hit_curves.push_back(state.curves);
                report.hit_certificates.push_back(
                    certify(givens, path, target, config.tolerance, config.mode));
                if (!report.hit_certificates.back().ok) ++counters.unused_places;
                return false;
            }
            report.found = true;
            report.steps = path;
            report.extra = static_cast<int>(path.size());
            report.points = state.points;
            report.curves = state.curves;
            // The certificate of the construction: every place of it is on the
            // way to the goal, so no step of it can simply be dropped.
            report.certificate = certify(givens, path, target, config.tolerance, config.mode);
            if (!report.certificate.ok) ++counters.unused_places;
            return true;
        }
        if (depth == 0) return false;
        const int count = static_cast<int>(state.points.size());
        for (int i = 0; i < count; ++i) {
            for (int j = i + 1; j < count; ++j) {
                const GPoint first = state.points[i];
                const GPoint second = state.points[j];
                GCurve curves[3] = {make_line(first, second), make_circle(first, second),
                                    make_circle(second, first)};
                const Step steps[3] = {{0, i, j}, {1, i, j}, {1, j, i}};
                for (int k = 0; k < 3; ++k) {
                    if (k == 1 && curves[1].c <= kEpsilon) continue;
                    if (k == 2 && curves[2].c <= kEpsilon) continue;
                    ++report.candidates;
                    // The filter is applied to the last curve of the walk
                    // alone when the target is a point that is given
                    // outright: the curve before it only has to leave behind
                    // the points that the last curve is built from.
                    const bool last = depth == 1;
                    if (!(config.filter_only_last && !last) &&
                        !curve_reaches_target(curves[k], target, config.tolerance))
                        continue;
                    // The figure of the branch is the figure of the walk with
                    // one more curve, and it is taken back when the branch is
                    // done.
                    const std::size_t point_mark = state.points.size();
                    const std::size_t curve_mark = state.curves.size();
                    bool cut = false;
                    if (!add_curve(state, curves[k], config.max_points, cut, true)) {
                        // A state that the point budget cut is a state the walk
                        // did not look at, so the count of them says whether the
                        // search was complete.
                        if (cut) ++report.states_cut;
                        continue;
                    }
                    path.push_back(steps[k]);
                    const bool hit = walk_classic(state, depth - 1);
                    if (hit) return true;
                    path.pop_back();
                    state.points.resize(point_mark);
                    state.curves.resize(curve_mark);
                    if (out_of_time) return false;
                }
            }
        }
        return false;
    };

    // The walk of the All model: every operation of the model is one step, and
    // a crossing is an operation like any other, so the walk has to choose the
    // pair of curves it crosses and the crossing it takes.
    //
    // The figure is built in place and not copied, and the point or the curve
    // that the branch adds is taken away again when the branch is done.
    std::function<bool(State&, int)> walk_all = [&](State& state, int depth) -> bool {
        ++report.states;
        if ((report.states & 0xFFFF) == 0 && std::chrono::steady_clock::now() > deadline) {
            out_of_time = true;
            return true;
        }
        if (owns_target(state, target, config.tolerance)) {
            if (config.collect_all && depth == 0) {
                report.hits.push_back(path);
                report.hit_points.push_back(state.points);
                report.hit_curves.push_back(state.curves);
                report.hit_certificates.push_back(
                    certify(givens, path, target, config.tolerance, config.mode));
                if (!report.hit_certificates.back().ok) ++counters.unused_places;
                return false;
            }
            report.found = true;
            report.steps = path;
            report.extra = static_cast<int>(path.size());
            report.points = state.points;
            report.curves = state.curves;
            // The certificate of the construction: every place of it is on the
            // way to the goal, so no step of it can simply be dropped.
            report.certificate = certify(givens, path, target, config.tolerance, config.mode);
            if (!report.certificate.ok) ++counters.unused_places;
            return true;
        }
        if (depth == 0) return false;
        const bool last = depth == 1;
        const int count = static_cast<int>(state.points.size());

        // A figure the walk has already been in has the same continuations as
        // the first time, for the goal and the steps that are left do not
        // depend on the way the figure was built.
        if (memo) {
            const FigureSet::Key key = figures.build(figure);
            if (!figures.insert(key)) {
                ++counters.duplicate_states;
                return false;
            }
        }

        // The filter of the last step, put before the candidate is built.
        //
        // A curve operation gives the figure no point, so it can only finish a
        // figure whose goal is a line, and only a family that draws a line can
        // draw that line.  A crossing operation gives the figure a point and no
        // curve, so it can only finish a figure whose goal is a point or a
        // circle.  Every candidate of the wrong kind is built, tested and
        // thrown away when the filter is put after the building, and the last
        // step is where most of the candidates of a walk are made.
        const bool filter_first = filter_first_enabled();

        // The slack of a figure whose goal is a point.  The point has to come
        // out of a crossing, and a crossing needs two curves that hold it, so
        // a figure that has `held` curves through the point and `depth` steps
        // left cannot reach the goal when `held + depth < 2`.  The test is
        // cheap -- it does not build anything -- and it cuts the branches that
        // draw a curve that does not hold the point one step too early.
        if (filter_first && target.kind == Target::Kind::Point && depth > 0) {
            int held = 0;
            for (const auto& curve : state.curves)
                if (curve_holds(curve, target.point, config.tolerance)) ++held;
            if (held + depth < 2) return false;
        }
        // The walk stops as soon as one of its branches has the goal, and the
        // branches below have to be told, for a lambda cannot end the walk it
        // is in.
        bool stop = false;
        bool hit_found = false;

        // A curve operation.  In this model a curve does not bring its
        // crossings along, so the only curve operation that can finish a
        // figure is the line of the target itself.
        const auto draw = [&](const GCurve& curve, const Step& step, Family family) {
            if (stop) return;
            if (skipped(family, skip)) return;
            if (last && filter_first && !family_finishes(family, target)) return;
            ++report.candidates;
            ++counters.candidates[family];
            if (last) {
                ++counters.candidates_last;
                if (target.kind != Target::Kind::Line) return;
                if (!curve_reaches_target(curve, target, config.tolerance)) return;
            } else {
                ++counters.candidates_other;
            }
            // The place of the curve: the figure tells whether it has it
            // already, and the places of the figure are the key that tells
            // the figures apart.
            std::uint64_t place[4];
            quantise_curve(curve, place);
            if (figure.has_curve(place)) {
                ++counters.duplicates;
                return;
            }
            ++counters.extensions;
            counters.copies += state.points.size() + state.curves.size();
            counters.widest_points = std::max<unsigned long long>(counters.widest_points,
                                                                  state.points.size());
            counters.widest_curves = std::max<unsigned long long>(counters.widest_curves,
                                                                  state.curves.size());
            // The figure of the branch is the figure of the walk with one more
            // curve, and it is taken back when the branch is done.
            figure.push_curve(place);
            state.curves.push_back(curve);
            path.push_back(step);
            const bool hit = walk_all(state, depth - 1);
            if (hit) {
                hit_found = true;
                stop = true;
                return;
            }
            path.pop_back();
            figure.pop_curve();
            state.curves.pop_back();
            if (out_of_time) {
                stop = true;
                return;
            }
        };

        // A crossing operation: it costs a step and it gives one point.
        const auto cross = [&](const GPoint& point, const Step& step) {
            if (stop) return;
            if (skipped(FamilyCrossing, skip)) return;
            ++report.candidates;
            ++counters.candidates[FamilyCrossing];
            ++counters.crossings;
            if (last && !crossing_reaches_target(point, target, config.tolerance)) return;
            // A crossing the figure already owns is not an operation.
            for (const auto& known : state.points) {
                const double dx = known.x - point.x;
                const double dy = known.y - point.y;
                if (dx * dx + dy * dy < 1e-18) return;
            }
            if (state.points.size() + 1 > config.max_points) return;
            // The place of the crossing, for the places of the figure are the
            // key that tells the figures apart.
            std::uint64_t place[2];
            quantise_point(point.x, point.y, place);
            figure.push_point(place);
            state.points.push_back(point);
            path.push_back(step);
            const bool hit = walk_all(state, depth - 1);
            if (hit) {
                hit_found = true;
                stop = true;
                return;
            }
            path.pop_back();
            figure.pop_point();
            state.points.pop_back();
            if (out_of_time) {
                stop = true;
                return;
            }
        };

        // Whether a curve operation may be made at this level at all.  At the
        // last level of a figure whose goal is a point or a circle it may not:
        // a curve operation gives the figure no point, so it cannot finish
        // such a figure.  The whole family is then left out, and not only the
        // candidates of it, for the loops that make them are the price of the
        // level.
        const bool curves_may_finish =
            !(last && filter_first && target.kind != Target::Kind::Line);

        // The two atomic operations, the perpendicular bisector of a segment
        // and the rigid circle.
        if (curves_may_finish)
        for (int i = 0; i < count && !stop; ++i) {
            for (int j = i + 1; j < count && !stop; ++j) {
                const GPoint first = state.points[i];
                const GPoint second = state.points[j];
                draw(make_line(first, second), Step{0, i, j}, FamilyLine);
                draw(make_circle(first, second), Step{1, i, j}, FamilyCircle);
                draw(make_circle(second, first), Step{1, j, i}, FamilyCircle);
                draw(fast::make_perpendicular_bisector(first, second),
                     Step{0, i, j, Op::PerpendicularBisector}, FamilyBisector);
                for (int k = 0; k < count && !stop; ++k) {
                    if (k == i || k == j) continue;
                    draw(fast::make_transfer_circle(first, second, state.points[k]),
                         Step{1, i, j, Op::TransferCircle, k}, FamilyTransfer);
                }
            }
        }
        // The perpendicular and the parallel to a line of the figure.
        if (curves_may_finish)
        for (int i = 0; i < count && !stop; ++i) {
            for (std::size_t c = 0; c < state.curves.size() && !stop; ++c) {
                if (state.curves[c].type != 0) continue;
                draw(fast::make_perpendicular_through(state.points[i], state.curves[c]),
                     Step{0, i, static_cast<int>(c), Op::PerpendicularThrough},
                     FamilyPerpendicular);
                draw(fast::make_parallel_through(state.points[i], state.curves[c]),
                     Step{0, i, static_cast<int>(c), Op::ParallelThrough}, FamilyParallel);
            }
        }
        // The bisectors of the angles of the figure.
        if (curves_may_finish)
        for (std::size_t c = 0; c < state.curves.size() && !stop; ++c) {
            if (state.curves[c].type != 0) continue;
            for (std::size_t d = c + 1; d < state.curves.size() && !stop; ++d) {
                if (state.curves[d].type != 0) continue;
                GCurve bisectors[2];
                const int number = fast::make_angle_bisectors(state.curves[c], state.curves[d],
                                                              bisectors);
                for (int k = 0; k < number && !stop; ++k)
                    draw(bisectors[k], Step{0, static_cast<int>(c), static_cast<int>(d),
                                            Op::AngleBisector, 0, k}, FamilyAngle);
            }
        }
        // The crossings of the figure.
        //
        // The last step of a figure whose goal is a line cannot be a crossing,
        // for a crossing gives the figure a point and no curve.  The last step
        // of a figure whose goal is a point can only be the crossing of two
        // curves that hold the point, so the curves that hold it are looked
        // for once and only the pairs of them are crossed: the pairs of the
        // other curves are not crossed at all.
        std::vector<std::size_t> holders;
        const bool only_holders = last && filter_first && target.kind == Target::Kind::Point;
        if (only_holders) {
            holders.reserve(state.curves.size());
            for (std::size_t c = 0; c < state.curves.size(); ++c)
                if (curve_holds(state.curves[c], target.point, config.tolerance))
                    holders.push_back(c);
        }
        const std::size_t pairs = only_holders ? holders.size() : state.curves.size();
        for (std::size_t c = 0; c < pairs && !stop; ++c) {
            // The last step of a figure whose goal is a line cannot be a
            // crossing, for a crossing gives the figure a point and no curve.
            if (last && filter_first && target.kind == Target::Kind::Line) break;
            for (std::size_t d = c + 1; d < pairs && !stop; ++d) {
                const std::size_t first = only_holders ? holders[c] : c;
                const std::size_t second = only_holders ? holders[d] : d;
                GPoint hits[2];
                ++counters.intersections;
                const int number = fast_intersections(state.curves[first], state.curves[second], hits);
                for (int k = 0; k < number && !stop; ++k)
                    cross(hits[k], Step{0, static_cast<int>(first), static_cast<int>(second),
                                        Op::Intersection, 0, k});
            }
        }
        return hit_found;
    };

    // Iterative deepening: the first limit that gives a hit is the cheapest
    // one, and the memory of the walk stays the size of one path.
    for (int limit = 0; limit <= config.max_extra; ++limit) {
        path.clear();
        report.states = 0;
        report.candidates = 0;
        report.states_cut = 0;
        out_of_time = false;
        figures.clear();
        const bool found = config.mode == engine::Mode::All ? walk_all(root, limit)
                                                            : walk_classic(root, limit);
        report.levels.push_back(static_cast<std::size_t>(report.states));
        if (config.verbose)
            std::fprintf(stderr, "limit %d: %llu states, %llu candidates%s\n", limit,
                         static_cast<unsigned long long>(report.states),
                         static_cast<unsigned long long>(report.candidates),
                         out_of_time ? ", out of time" : "");
        if (found) break;
        // The hits of the first level that has any are the cheapest ones, so
        // the walk stops there when it was asked to collect all of them.
        if (config.collect_all && !report.hits.empty()) {
            report.found = true;
            report.extra = limit;
            break;
        }
        if (out_of_time) break;
    }

    report.out_of_time = out_of_time;
    const std::chrono::duration<double> since = std::chrono::steady_clock::now() - started;
    report.elapsed = since.count();
    if (std::getenv("NEO_GOAL_PROFILE") != nullptr) {
        static const char* const names[FamilyCount] = {
            "line", "circle", "bisector", "transfer",
            "perpendicular", "parallel", "angle", "crossing"};
        std::fprintf(stderr,
                     "profile: candidates");
        for (int family = 0; family < FamilyCount; ++family)
            std::fprintf(stderr, " %s=%llu", names[family], counters.candidates[family]);
        std::fprintf(stderr,
                     " | last=%llu other=%llu intersections=%llu extensions=%llu"
                     " copies=%llu duplicates=%llu crossings=%llu repeated=%llu"
                     " widest_points=%llu widest_curves=%llu"
                     " | states=%llu candidates=%llu unused=%llu elapsed=%.6f\n",
                     counters.candidates_last, counters.candidates_other,
                     counters.intersections, counters.extensions, counters.copies,
                     counters.duplicates, counters.crossings, counters.duplicate_states,
                     counters.widest_points,
                     counters.widest_curves,
                     static_cast<unsigned long long>(report.states),
                     static_cast<unsigned long long>(report.candidates),
                     counters.unused_places, report.elapsed);
    }
    return report;
}

}  // namespace neo::goal
