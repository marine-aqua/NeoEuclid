#pragma once

// A small geometry kernel for the exhaustive walks.
//
// The engine's own State carries recipes and labels, which is what a report
// needs and what a search does not: a search that visits millions of states
// only wants the points and the curves, and it wants them without allocating.
// This header holds that kernel: curves in a flat form, intersections without
// any allocation, and the quantisation of a coordinate that decides whether
// two points are the same point.

#include "neo/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace neo::fast {

// The quantisation of a point on the 1e-8 grid.  It is the same grid that the
// state key of the engine uses, so a state of the walk and a state of the
// engine merge in the same way.
using QuantPoint = QuantizedPoint;

inline QuantPoint quantize_point(double x, double y) {
    // nearbyint rounds ties to even, which is what Python's round(x, 8) does;
    // the two implementations therefore agree on the discretisation.
    return {static_cast<std::int64_t>(std::nearbyint(x * 1e8)),
            static_cast<std::int64_t>(std::nearbyint(y * 1e8))};
}

struct GPoint {
    double x{};
    double y{};
};

struct GCurve {
    std::uint8_t type{0};  // 0 = line (a x + b y + c = 0, unit normal), 1 = circle (centre a,b radius c)
    double a{};
    double b{};
    double c{};
};

// The length of a segment.  It is `sqrt` and not `hypot`: the arguments are
// the coordinates of a figure, so they cannot overflow, and `hypot` is a call
// of the library that costs about a hundred instructions for a length that
// every candidate of a walk needs.
inline double length_of(double dx, double dy) {
    return std::sqrt(dx * dx + dy * dy);
}

inline GCurve make_line(GPoint p, GPoint q) {
    double a = p.y - q.y;
    double b = q.x - p.x;
    const double norm = length_of(a, b);
    a /= norm;
    b /= norm;
    double c = -(a * p.x + b * p.y);
    if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
        a = -a;
        b = -b;
        c = -c;
    }
    return GCurve{0, a, b, c};
}

inline GCurve make_circle(GPoint centre, GPoint through) {
    return GCurve{1, centre.x, centre.y,
                  length_of(centre.x - through.x, centre.y - through.y)};
}

// The operations of the All model that the kernel has to build: the
// perpendicular bisector of a segment, the perpendicular and the parallel to a
// line, the rigid circle, and the two bisectors of an angle.
inline GCurve make_perpendicular_bisector(GPoint p, GPoint q) {
    double a = q.x - p.x;
    double b = q.y - p.y;
    const double norm = length_of(a, b);
    if (norm < kEpsilon) return GCurve{0, 0.0, 0.0, 0.0};
    a /= norm; b /= norm;
    const double mx = (p.x + q.x) / 2.0;
    const double my = (p.y + q.y) / 2.0;
    double c = -(a * mx + b * my);
    if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
        a = -a; b = -b; c = -c;
    }
    return GCurve{0, a, b, c};
}

inline GCurve make_perpendicular_through(GPoint p, const GCurve& line) {
    double a = line.b;
    double b = -line.a;
    double c = -(a * p.x + b * p.y);
    if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
        a = -a; b = -b; c = -c;
    }
    return GCurve{0, a, b, c};
}

inline GCurve make_parallel_through(GPoint p, const GCurve& line) {
    double a = line.a;
    double b = line.b;
    double c = -(a * p.x + b * p.y);
    if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
        a = -a; b = -b; c = -c;
    }
    return GCurve{0, a, b, c};
}

// The rigid compass: the circle about a point whose radius is the length of a
// segment given elsewhere.
inline GCurve make_transfer_circle(GPoint centre, GPoint from, GPoint to) {
    return GCurve{1, centre.x, centre.y, length_of(from.x - to.x, from.y - to.y)};
}

// The two bisectors of the angle of two lines.  The number of them that the
// caller gets is written back, because two lines that are parallel have none.
inline int make_angle_bisectors(const GCurve& first, const GCurve& second, GCurve* out) {
    int count = 0;
    for (const double sign : {1.0, -1.0}) {
        double a = first.a + sign * second.a;
        double b = first.b + sign * second.b;
        double c = first.c + sign * second.c;
        const double norm = length_of(a, b);
        if (norm < kEpsilon) continue;
        a /= norm; b /= norm; c /= norm;
        if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
            a = -a; b = -b; c = -c;
        }
        out[count++] = GCurve{0, a, b, c};
    }
    return count;
}

// The crossings of two curves.  Returns how many there are, and writes them
// to out.  The order is the order of the curves: a line in front of a circle
// gives the foot of the perpendicular first and the other crossing second.
inline int fast_intersections(const GCurve& first, const GCurve& second, GPoint* out) {
    if (first.type == 0 && second.type == 0) {
        const double determinant = first.a * second.b - second.a * first.b;
        if (std::abs(determinant) < kEpsilon) return 0;
        out[0] = {(first.b * second.c - second.b * first.c) / determinant,
                  (first.c * second.a - second.c * first.a) / determinant};
        return 1;
    }
    if (first.type == 0) {
        const double distance = first.a * second.a + first.b * second.b + first.c;
        const double radius = second.c;
        if (std::abs(distance) > radius + kEpsilon) return 0;
        const double foot_x = second.a - distance * first.a;
        const double foot_y = second.b - distance * first.b;
        const double delta = std::sqrt(std::max(0.0, radius * radius - distance * distance));
        if (delta < kEpsilon) {
            out[0] = {foot_x, foot_y};
            return 1;
        }
        out[0] = {foot_x - first.b * delta, foot_y + first.a * delta};
        out[1] = {foot_x + first.b * delta, foot_y - first.a * delta};
        return 2;
    }
    if (second.type == 0) {
        const double distance = second.a * first.a + second.b * first.b + second.c;
        const double radius = first.c;
        if (std::abs(distance) > radius + kEpsilon) return 0;
        const double foot_x = first.a - distance * second.a;
        const double foot_y = first.b - distance * second.b;
        const double delta = std::sqrt(std::max(0.0, radius * radius - distance * distance));
        if (delta < kEpsilon) {
            out[0] = {foot_x, foot_y};
            return 1;
        }
        out[0] = {foot_x - second.b * delta, foot_y + second.a * delta};
        out[1] = {foot_x + second.b * delta, foot_y - second.a * delta};
        return 2;
    }
    const double dx = second.a - first.a;
    const double dy = second.b - first.b;
    const double distance = length_of(dx, dy);
    if (distance < kEpsilon || distance > first.c + second.c + kEpsilon ||
        distance < std::abs(first.c - second.c) - kEpsilon) return 0;
    const double along = (first.c * first.c - second.c * second.c + distance * distance) / (2.0 * distance);
    const double height_squared = first.c * first.c - along * along;
    if (height_squared < -kEpsilon) return 0;
    const double height = std::sqrt(std::max(0.0, height_squared));
    const double middle_x = first.a + along * dx / distance;
    const double middle_y = first.b + along * dy / distance;
    const double rx = -dy * height / distance;
    const double ry = dx * height / distance;
    out[0] = {middle_x + rx, middle_y + ry};
    out[1] = {middle_x - rx, middle_y - ry};
    return height < kEpsilon ? 1 : 2;
}

// Is a point on a circle, within the tolerance?  The distance `d` of the point
// from the centre and the radius `r` are near each other when `|d - r|` is
// small.  The test is made without a square root: `d` is in `[r - t, r + t]`
// when `d^2` is in `[(r - t)^2, (r + t)^2]`, and the lower bound is only there
// when the radius is bigger than the tolerance.  `squared` is `d^2`.
inline bool circle_holds(double squared, double radius, double tolerance) {
    const double outer = radius + tolerance;
    if (squared > outer * outer) return false;
    if (radius < tolerance) return true;
    const double inner = radius - tolerance;
    return squared >= inner * inner;
}

inline bool gcurve_contains(const GCurve& curve, GPoint target, double tolerance) {
    if (curve.type == 0) return std::abs(curve.a * target.x + curve.b * target.y + curve.c) <= tolerance;
    const double dx = target.x - curve.a;
    const double dy = target.y - curve.b;
    return circle_holds(dx * dx + dy * dy, curve.c, tolerance);
}

inline bool same_curve(const GCurve& first, const GCurve& second) {
    if (first.type != second.type) return false;
    if (first.type == 0)
        return quantize_point(first.a, first.b) == quantize_point(second.a, second.b) &&
               quantize_point(first.c, 0.0) == quantize_point(second.c, 0.0);
    return quantize_point(first.a, first.b) == quantize_point(second.a, second.b) &&
           quantize_point(first.c, 0.0) == quantize_point(second.c, 0.0);
}

}  // namespace neo::fast
