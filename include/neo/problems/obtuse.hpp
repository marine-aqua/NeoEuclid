#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace neo::obtuse {

struct SearchConfig {
    int max_cost{10};
    std::size_t max_points{40};
    std::size_t expansions_per_restart{20000};
    std::size_t terminal_states_per_restart{20000};
    std::uint64_t seed{20260926};
    double time_limit_seconds{60.0};
    bool verbose{true};
    bool reflected_prefix{true};
    std::vector<double> sample_degrees;
};

struct SearchReport {
    bool found{};
    int cost{};
    std::vector<std::string> steps;
    std::size_t restarts{};
    std::size_t expanded{};
    std::size_t generated{};
    std::size_t duplicates{};
    double elapsed_seconds{};
    std::string best_description;
    std::vector<std::string> best_steps;
    std::vector<std::string> best_point_definitions;
};

struct ModelDiagnostics {
    int prefix_cost{};
    std::size_t visible_curves{};
    std::size_t visible_points{};
    std::size_t stable_hyperbola_points{};
    bool helper_circles_hidden{};
    bool unstable_branch_rejected{};
    bool quartic_solver_verified{};
    bool macro_costs_verified{};
    bool existing_objects_reused_without_recharge{};
    bool target_ray_join_verified{};
    bool origin_branch_verified{};
    bool near_duplicate_curve_rejected{};
    bool transverse_join_verified{};
    bool raw_initial_state_verified{};
};

ModelDiagnostics diagnose_model(const std::vector<double>& sample_degrees);

// Randomised DFS for the hidden-circle, three-cost y-reflection prefix model.
// Only intersection branches present in every configured sample are retained.
SearchReport search(const SearchConfig& config);

struct MeetConfig : SearchConfig {
    int shared_cost_min{3};
    int shared_cost_max{7};
    std::size_t branch_records_per_base{1500};
    std::size_t branch_attempts_per_base{30000};
};

// Meet-in-the-middle search. Each trial constructs a shared prefix, then
// enumerates two independent descendants and joins curves that meet at the
// same positive point of the target ray in every sample.
SearchReport meet_in_the_middle(const MeetConfig& config);

}  // namespace neo::obtuse
