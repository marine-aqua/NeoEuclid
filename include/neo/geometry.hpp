#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

namespace neo {

constexpr double kEpsilon = 1e-9;

struct Point {
    double x{};
    double y{};
    std::uint32_t recipe{};
};

struct Line {
    double a{};
    double b{};
    double c{};
    std::uint32_t recipe{};
};

struct Circle {
    double x{};
    double y{};
    double radius{};
    std::uint32_t recipe{};
};

using Curve = std::variant<Line, Circle>;

struct QuantizedPoint {
    std::int64_t x{};
    std::int64_t y{};
    bool operator==(const QuantizedPoint& other) const {
        return x == other.x && y == other.y;
    }
};

struct QuantizedCurve {
    std::uint8_t kind{};
    std::int64_t first{};
    std::int64_t second{};
    std::int64_t third{};
    bool operator==(const QuantizedCurve& other) const {
        return kind == other.kind && first == other.first &&
               second == other.second && third == other.third;
    }
    bool operator<(const QuantizedCurve& other) const;
};

struct PointKeyHash {
    std::size_t operator()(const QuantizedPoint& key) const noexcept;
};

struct CurveKeyHash {
    std::size_t operator()(const QuantizedCurve& key) const noexcept;
};

std::optional<Line> line_through(const Point& first, const Point& second,
                                 std::uint32_t recipe = 0);
std::optional<Circle> circle_center_through(const Point& center,
                                            const Point& through,
                                            std::uint32_t recipe = 0);
// The rigid compass: the circle about a point whose radius is the length of a
// segment given elsewhere.  The soft compass is `circle_center_through`.
std::optional<Circle> transfer_circle(const Point& center, const Point& from,
                                      const Point& to, std::uint32_t recipe = 0);
std::optional<Line> perpendicular_bisector(const Point& first,
                                           const Point& second,
                                           std::uint32_t recipe = 0);
// The perpendicular to a line through a point, and the parallel to a line
// through a point.
std::optional<Line> perpendicular_through(const Point& point, const Line& line,
                                          std::uint32_t recipe = 0);
std::optional<Line> parallel_through(const Point& point, const Line& line,
                                     std::uint32_t recipe = 0);
std::vector<Line> angle_bisectors(const Line& first, const Line& second,
                                  std::uint32_t first_recipe = 0,
                                  std::uint32_t second_recipe = 0);
std::vector<Point> intersections(const Curve& first, const Curve& second);
QuantizedPoint point_key(const Point& point, int digits = 8);
QuantizedCurve curve_key(const Curve& curve, int digits = 8);
double line_direction_error(const Line& line, double target_radians);

}  // namespace neo
