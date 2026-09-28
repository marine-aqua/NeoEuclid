#pragma once

#include "neo/search.hpp"

#include <optional>

namespace neo::problems {

struct FixedAngleMitmConfig {
    int max_cost{6};
    int shared_cost{3};
    int arm_cost{2};
    std::size_t prefix_beam{5000};
    std::size_t max_points{32};
    std::size_t max_states{2'000'000};
    double time_limit_seconds{60.0};
    bool use_macros{false};
};

std::optional<SearchResult> fixed_angle_mitm(double target_degrees,
                                             const FixedAngleMitmConfig& config);

}  // namespace neo::problems
