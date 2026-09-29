#pragma once

#include "neo/engine/operation.hpp"
#include "neo/geometry.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace neo {

using OperationKind = engine::OperationKind;

struct RecipeNode {
    OperationKind kind{OperationKind::Given};
    std::uint32_t first{};
    std::uint32_t second{};
    std::string label;
};

struct Candidate {
    Curve curve;
    int cost{1};
    OperationKind operation{OperationKind::LineThrough};
    std::uint32_t first{};
    std::uint32_t second{};
};

struct State {
    std::vector<Point> points;
    std::vector<Curve> curves;
    std::vector<RecipeNode> recipes;
    std::vector<std::uint32_t> steps;
    int cost{};
    double score{};
};

struct ExpansionOptions {
    std::size_t max_points{28};
    bool use_macros{true};
    int key_digits{8};
    engine::CostPolicy cost_policy{engine::CostPolicy::Atomic};
};

std::vector<Candidate> candidate_curves(const State& state, bool use_macros,
                                        engine::CostPolicy cost_policy = engine::CostPolicy::Atomic);
std::optional<State> add_curve(const State& parent, const Candidate& candidate,
                               const ExpansionOptions& options);
std::string state_key(const State& state, int digits = 7);
std::string describe_recipe(const State& state, std::uint32_t recipe_id);

}  // namespace neo
