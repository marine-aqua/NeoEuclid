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
    // The rigid circle needs a third place: the two ends of the segment whose
    // length is the radius.  It comes after the label so that the places that
    // do not need it keep their shape.
    std::uint32_t third{};
};

struct Candidate {
    Curve curve;
    int cost{1};
    OperationKind operation{OperationKind::LineThrough};
    std::uint32_t first{};
    std::uint32_t second{};
    // The rigid circle needs three places: the centre, and the two ends of the
    // segment whose length is the radius.
    std::uint32_t third{};
};

// A crossing of two curves that the figure owns.  In the All model a crossing
// is an operation like any other: it costs a step, and the search has to say
// which pair of curves it crosses and which of the crossings of the pair it
// takes.
struct CrossingCandidate {
    Point point;
    int cost{1};
    std::uint32_t first{};      // the recipe of the first curve
    std::uint32_t second{};     // the recipe of the second curve
    std::size_t first_curve{};  // the place of the first curve in the figure
    std::size_t second_curve{};
    int which{};                // which of the crossings of the pair
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
    // The model of the construction.  In the All model the crossings are
    // operations that cost a step, and the figure gets them only when the
    // search asks for them.
    engine::Mode mode{engine::Mode::Classic};
};

std::vector<Candidate> candidate_curves(const State& state, bool use_macros,
                                        engine::CostPolicy cost_policy = engine::CostPolicy::Atomic);
// Every curve the figure can draw next, in the model of the given mode.  In the
// classic model this is the two atomic operations and the macros; in the All
// model it is every operation of the model, and each of them costs one step.
std::vector<Candidate> candidate_operations(const State& state, engine::Mode mode,
                                             bool use_macros = true);
// The crossings the figure can take, one operation each.  The classic model
// does not use them: it has them for nothing.
std::vector<CrossingCandidate> candidate_crossings(const State& state);
std::optional<State> add_curve(const State& parent, const Candidate& candidate,
                               const ExpansionOptions& options);
// The crossing of two curves, as an operation of its own.
std::optional<State> add_crossing(const State& parent, const CrossingCandidate& candidate,
                                   const ExpansionOptions& options);
std::string state_key(const State& state, int digits = 7);
std::string describe_recipe(const State& state, std::uint32_t recipe_id);

}  // namespace neo
