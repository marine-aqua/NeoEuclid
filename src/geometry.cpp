#include "neo/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <tuple>

namespace neo {
namespace {

std::int64_t quantize(double value, int digits) {
    return static_cast<std::int64_t>(std::llround(value * std::pow(10.0, digits)));
}

template <typename T>
void combine_hash(std::size_t& seed, const T& value) {
    seed ^= std::hash<T>{}(value) + 0x9e3779b9U + (seed << 6U) + (seed >> 2U);
}

std::vector<Point> line_line(const Line& first, const Line& second) {
    const double determinant = first.a * second.b - second.a * first.b;
    if (std::abs(determinant) < kEpsilon) return {};
    return {{(first.b * second.c - second.b * first.c) / determinant,
             (first.c * second.a - second.c * first.a) / determinant, 0}};
}

std::vector<Point> line_circle(const Line& line, const Circle& circle) {
    const double signed_distance = line.a * circle.x + line.b * circle.y + line.c;
    if (std::abs(signed_distance) > circle.radius + kEpsilon) return {};
    const double foot_x = circle.x - signed_distance * line.a;
    const double foot_y = circle.y - signed_distance * line.b;
    const double delta = std::sqrt(std::max(0.0, circle.radius * circle.radius -
                                                     signed_distance * signed_distance));
    if (delta < kEpsilon) return {{foot_x, foot_y, 0}};
    return {{foot_x - line.b * delta, foot_y + line.a * delta, 0},
            {foot_x + line.b * delta, foot_y - line.a * delta, 0}};
}

std::vector<Point> circle_circle(const Circle& first, const Circle& second) {
    const double dx = second.x - first.x;
    const double dy = second.y - first.y;
    const double distance = std::hypot(dx, dy);
    if (distance < kEpsilon || distance > first.radius + second.radius + kEpsilon ||
        distance < std::abs(first.radius - second.radius) - kEpsilon) return {};
    const double along = (first.radius * first.radius - second.radius * second.radius +
                          distance * distance) / (2.0 * distance);
    const double height_squared = first.radius * first.radius - along * along;
    if (height_squared < -kEpsilon) return {};
    const double height = std::sqrt(std::max(0.0, height_squared));
    const double middle_x = first.x + along * dx / distance;
    const double middle_y = first.y + along * dy / distance;
    if (height < kEpsilon) return {{middle_x, middle_y, 0}};
    const double rx = -dy * height / distance;
    const double ry = dx * height / distance;
    return {{middle_x + rx, middle_y + ry, 0},
            {middle_x - rx, middle_y - ry, 0}};
}

}  // namespace

bool QuantizedCurve::operator<(const QuantizedCurve& other) const {
    return std::tie(kind, first, second, third) <
           std::tie(other.kind, other.first, other.second, other.third);
}

std::size_t PointKeyHash::operator()(const QuantizedPoint& key) const noexcept {
    std::size_t seed = 0;
    combine_hash(seed, key.x); combine_hash(seed, key.y);
    return seed;
}

std::size_t CurveKeyHash::operator()(const QuantizedCurve& key) const noexcept {
    std::size_t seed = key.kind;
    combine_hash(seed, key.first); combine_hash(seed, key.second); combine_hash(seed, key.third);
    return seed;
}

std::optional<Line> line_through(const Point& first, const Point& second,
                                 std::uint32_t recipe) {
    double a = first.y - second.y;
    double b = second.x - first.x;
    const double norm = std::hypot(a, b);
    if (norm < kEpsilon) return std::nullopt;
    a /= norm; b /= norm;
    double c = -(a * first.x + b * first.y);
    if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
        a = -a; b = -b; c = -c;
    }
    return Line{a, b, c, recipe};
}

std::optional<Circle> circle_center_through(const Point& center,
                                            const Point& through,
                                            std::uint32_t recipe) {
    const double radius = std::hypot(center.x - through.x, center.y - through.y);
    if (radius < kEpsilon || !std::isfinite(radius)) return std::nullopt;
    return Circle{center.x, center.y, radius, recipe};
}

std::optional<Line> perpendicular_bisector(const Point& first,
                                           const Point& second,
                                           std::uint32_t recipe) {
    double a = second.x - first.x;
    double b = second.y - first.y;
    const double norm = std::hypot(a, b);
    if (norm < kEpsilon) return std::nullopt;
    a /= norm; b /= norm;
    const double mx = (first.x + second.x) / 2.0;
    const double my = (first.y + second.y) / 2.0;
    double c = -(a * mx + b * my);
    if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
        a = -a; b = -b; c = -c;
    }
    return Line{a, b, c, recipe};
}

std::vector<Line> angle_bisectors(const Line& first, const Line& second,
                                  std::uint32_t first_recipe,
                                  std::uint32_t second_recipe) {
    std::vector<Line> result;
    for (const auto [sign, recipe] : {std::pair{1.0, first_recipe},
                                      std::pair{-1.0, second_recipe}}) {
        double a = first.a + sign * second.a;
        double b = first.b + sign * second.b;
        double c = first.c + sign * second.c;
        const double norm = std::hypot(a, b);
        if (norm < kEpsilon) continue;
        a /= norm; b /= norm; c /= norm;
        if (a < -kEpsilon || (std::abs(a) <= kEpsilon && b < 0.0)) {
            a = -a; b = -b; c = -c;
        }
        result.push_back({a, b, c, recipe});
    }
    return result;
}

std::vector<Point> intersections(const Curve& first, const Curve& second) {
    if (const auto* line1 = std::get_if<Line>(&first)) {
        if (const auto* line2 = std::get_if<Line>(&second)) return line_line(*line1, *line2);
        return line_circle(*line1, std::get<Circle>(second));
    }
    if (const auto* line2 = std::get_if<Line>(&second)) {
        return line_circle(*line2, std::get<Circle>(first));
    }
    return circle_circle(std::get<Circle>(first), std::get<Circle>(second));
}

QuantizedPoint point_key(const Point& point, int digits) {
    return {quantize(point.x, digits), quantize(point.y, digits)};
}

QuantizedCurve curve_key(const Curve& curve, int digits) {
    if (const auto* line = std::get_if<Line>(&curve)) {
        return {0, quantize(line->a, digits), quantize(line->b, digits), quantize(line->c, digits)};
    }
    const auto& circle = std::get<Circle>(curve);
    return {1, quantize(circle.x, digits), quantize(circle.y, digits),
            quantize(circle.radius, digits)};
}

double line_direction_error(const Line& line, double target_radians) {
    return std::abs(line.a * std::cos(target_radians) +
                    line.b * std::sin(target_radians)) + std::abs(line.c);
}

}  // namespace neo
