#pragma once

#include "neo/state.hpp"

#include <cstddef>
#include <optional>
#include <string>

namespace neo {

struct SearchConfig {
    int max_cost{6};
    std::size_t beam_width{500};
    std::size_t max_points{28};
    bool use_macros{true};
    bool verbose{true};
    engine::CostPolicy cost_policy{engine::CostPolicy::Atomic};
};

struct Goal {
    std::string final_step;
    int total_cost{};
};

struct SearchResult {
    State state;
    Goal goal;
    std::size_t expanded{};
    std::size_t generated{};
};

class SearchProblem {
public:
    virtual ~SearchProblem() = default;
    virtual State initial_state() const = 0;
    virtual double score(const State& state) const = 0;
    virtual std::optional<Goal> goal(const State& state,
                                     int max_cost) const = 0;
};

std::optional<SearchResult> beam_search(const SearchProblem& problem,
                                        const SearchConfig& config);

}  // namespace neo
