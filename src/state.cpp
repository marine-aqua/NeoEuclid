#include "neo/state.hpp"

#include <algorithm>
#include <sstream>
#include <unordered_set>

namespace neo {
namespace {

std::uint32_t curve_recipe(const Curve& curve) {
    return std::visit([](const auto& value) { return value.recipe; }, curve);
}

Curve with_recipe(Curve curve, std::uint32_t recipe) {
    std::visit([recipe](auto& value) { value.recipe = recipe; }, curve);
    return curve;
}

std::string operation_name(OperationKind kind) {
    switch (kind) {
        case OperationKind::LineThrough: return "line";
        case OperationKind::CircleCenterThrough: return "circle";
        case OperationKind::PerpendicularBisector: return "perpendicular-bisector";
        case OperationKind::AngleBisector: return "angle-bisector";
        case OperationKind::Intersection: return "intersect";
        case OperationKind::Given: return "given";
    }
    return "operation";
}

}  // namespace

std::vector<Candidate> candidate_curves(const State& state, bool use_macros) {
    std::vector<Candidate> result;
    for (std::size_t i = 0; i < state.points.size(); ++i) {
        for (std::size_t j = i + 1; j < state.points.size(); ++j) {
            const auto& first = state.points[i];
            const auto& second = state.points[j];
            if (auto line = line_through(first, second))
                result.push_back({*line, 1, OperationKind::LineThrough, first.recipe, second.recipe});
            if (auto circle = circle_center_through(first, second))
                result.push_back({*circle, 1, OperationKind::CircleCenterThrough, first.recipe, second.recipe});
            if (auto circle = circle_center_through(second, first))
                result.push_back({*circle, 1, OperationKind::CircleCenterThrough, second.recipe, first.recipe});
            if (use_macros) {
                if (auto line = perpendicular_bisector(first, second))
                    result.push_back({*line, 3, OperationKind::PerpendicularBisector, first.recipe, second.recipe});
            }
        }
    }
    if (use_macros) {
        for (std::size_t i = 0; i < state.curves.size(); ++i) {
            const auto* first = std::get_if<Line>(&state.curves[i]);
            if (!first) continue;
            for (std::size_t j = i + 1; j < state.curves.size(); ++j) {
                const auto* second = std::get_if<Line>(&state.curves[j]);
                if (!second) continue;
                for (const auto& line : angle_bisectors(*first, *second))
                    result.push_back({line, 4, OperationKind::AngleBisector,
                                      first->recipe, second->recipe});
            }
        }
    }
    return result;
}

std::optional<State> add_curve(const State& parent, const Candidate& candidate,
                               const ExpansionOptions& options) {
    std::unordered_set<QuantizedCurve, CurveKeyHash> known_curves;
    for (const auto& curve : parent.curves) known_curves.insert(curve_key(curve, options.key_digits));
    if (known_curves.count(curve_key(candidate.curve, options.key_digits))) return std::nullopt;

    State child = parent;
    RecipeNode operation;
    operation.kind = candidate.operation;
    operation.first = candidate.first;
    operation.second = candidate.second;
    operation.label = operation_name(operation.kind);
    const std::uint32_t operation_id = static_cast<std::uint32_t>(child.recipes.size());
    child.recipes.push_back(std::move(operation));
    Curve new_curve = with_recipe(candidate.curve, operation_id);

    std::unordered_set<QuantizedPoint, PointKeyHash> known_points;
    for (const auto& point : child.points) known_points.insert(point_key(point, options.key_digits));
    for (const auto& old : parent.curves) {
        for (auto point : intersections(new_curve, old)) {
            if (!known_points.insert(point_key(point, options.key_digits)).second) continue;
            if (child.points.size() + 1 > options.max_points) return std::nullopt;
            const auto recipe_id = static_cast<std::uint32_t>(child.recipes.size());
            child.recipes.push_back({OperationKind::Intersection, operation_id,
                                     curve_recipe(old), ""});
            point.recipe = recipe_id;
            child.points.push_back(point);
        }
    }
    child.curves.push_back(new_curve);
    child.steps.push_back(operation_id);
    child.cost += candidate.cost;
    return child;
}

std::string state_key(const State& state, int digits) {
    std::vector<QuantizedCurve> keys;
    keys.reserve(state.curves.size());
    for (const auto& curve : state.curves) keys.push_back(curve_key(curve, digits));
    std::sort(keys.begin(), keys.end());
    std::ostringstream output;
    for (const auto& key : keys) {
        output << static_cast<int>(key.kind) << ':' << key.first << ',' << key.second
               << ',' << key.third << ';';
    }
    return output.str();
}

std::string describe_recipe(const State& state, std::uint32_t recipe_id) {
    if (recipe_id >= state.recipes.size()) return "?";
    const auto& node = state.recipes[recipe_id];
    if (node.kind == OperationKind::Given) return node.label;
    const auto first = describe_recipe(state, node.first);
    const auto second = describe_recipe(state, node.second);
    switch (node.kind) {
        case OperationKind::LineThrough: return "line(" + first + ", " + second + ")";
        case OperationKind::CircleCenterThrough: return "circle(" + first + "; " + second + ")";
        case OperationKind::PerpendicularBisector: return "perpendicular-bisector(" + first + ", " + second + ")";
        case OperationKind::AngleBisector: return "angle-bisector(" + first + ", " + second + ")";
        case OperationKind::Intersection: return "intersect(" + first + ", " + second + ")";
        case OperationKind::Given: break;
    }
    return "?";
}

}  // namespace neo
