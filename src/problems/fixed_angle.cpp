#include "neo/problems/fixed_angle.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace neo::problems {
namespace {
constexpr double kPi = 3.14159265358979323846;
}

FixedAngle::FixedAngle(double target_degrees)
    : target_degrees_(target_degrees), target_radians_(target_degrees * kPi / 180.0) {}

State FixedAngle::initial_state() const {
    State state;
    state.recipes = {{OperationKind::Given, 0, 0, "O"},
                     {OperationKind::Given, 0, 0, "A"},
                     {OperationKind::Given, 0, 0, "baseline[given]"}};
    state.points = {{0.0, 0.0, 0}, {1.0, 0.0, 1}};
    state.curves = {Line{0.0, 1.0, 0.0, 2}};
    state.score = score(state);
    return state;
}

double FixedAngle::score(const State& state) const {
    double best = std::numeric_limits<double>::infinity();
    for (const auto& point : state.points) {
        if (std::hypot(point.x, point.y) <= kEpsilon || point.y < -kEpsilon) continue;
        best = std::min(best, std::abs(std::atan2(point.y, point.x) - target_radians_));
    }
    const double degrees = best * 180.0 / kPi;
    const double band = degrees < 1e-5 ? 0.0 : degrees < 0.1 ? 1.0 : degrees < 2.0 ? 2.0 : 3.0;
    return band * 100.0 + best;
}

bool FixedAngle::exact_target_check(const Point& point) const {
    const double radius = std::hypot(point.x, point.y);
    if (radius < kEpsilon || point.y <= 0.0) return false;
    const double cosine = point.x / radius;
    if (std::abs(target_degrees_ - 72.0) < 1e-12)
        return std::abs(4.0 * cosine * cosine + 2.0 * cosine - 1.0) < 2e-8 && cosine > 0.0;
    if (std::abs(target_degrees_ - 36.0) < 1e-12)
        return std::abs(4.0 * cosine * cosine - 2.0 * cosine - 1.0) < 2e-8 && cosine > 0.0;
    const double degrees = std::acos(std::clamp(cosine, -1.0, 1.0)) * 180.0 / kPi;
    return std::abs(degrees - target_degrees_) < 1e-7;
}

std::optional<Goal> FixedAngle::goal(const State& state, int max_cost) const {
    for (const auto& curve : state.curves) {
        if (const auto* line = std::get_if<Line>(&curve)) {
            if (line_direction_error(*line, target_radians_) < 2e-8)
                return Goal{"target line already drawn", state.cost};
        }
    }
    if (state.cost + 1 > max_cost) return std::nullopt;
    const double dx = std::cos(target_radians_);
    const double dy = std::sin(target_radians_);
    for (const auto& point : state.points) {
        const double along = point.x * dx + point.y * dy;
        const double across = std::abs(-point.x * dy + point.y * dx);
        if (along > kEpsilon && across < 2e-8 && exact_target_check(point)) {
            return Goal{"line(O, " + describe_recipe(state, point.recipe) + ")", state.cost + 1};
        }
    }
    return std::nullopt;
}

}  // namespace neo::problems

