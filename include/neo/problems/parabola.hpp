#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace neo::parabola {

struct SearchConfig {
    int max_cost{6};
    int prefix_depth{2};
    int threads{0};
    std::size_t max_points{48};
    std::size_t max_states{2'000'000};
    std::size_t state_cache_entries{4'000'000};
    std::size_t geometry_cache_entries{200'000};
    double time_limit_seconds{60.0};
    std::uint64_t seed{20260927};
    // Search uses only these samples. A hit is replayed on validation_samples.
    std::vector<double> search_degrees{38.0, 103.0, 169.0};
    int validation_samples{81};
    bool verbose{true};
    bool angle_at_parabola_vertex{false};
    bool target_circle_center{false};
    bool require_bisector_use{false};
    bool target_k2_circle{false};
    double validation_max_degrees{0.0};
    bool require_first_step_uses_p{false};
    bool require_circle_parabola_use{false};
    std::vector<std::string> masks; // e.g. LLCCLC; empty means all L/C masks.
};

struct SearchReport {
    bool found{};
    bool densely_verified{};
    int cost{};
    std::vector<std::string> steps;
    std::string terminal_schema;
    std::string mask;
    std::size_t prefixes{};
    std::size_t expanded{};
    std::size_t generated{};
    std::size_t duplicates{};
    std::size_t dense_rejections{};
    std::size_t state_cache_hits{};
    std::size_t state_cache_entries{};
    std::size_t parabola_cache_hits{};
    std::size_t pair_cache_hits{};
    std::size_t local_parabola_cache_hits{};
    std::size_t local_pair_cache_hits{};
    std::size_t geometry_cache_entries{};
    std::size_t mitm_records{};
    std::size_t mitm_matches{};
    double elapsed_seconds{};
};

// y^2=4x, angle vertex F=(1,0), free axes, O, focus, directrix, parabola,
// and free P on the second side of the angle. The target is the internal
// alpha/3 ray. Intersections with every existing curve and with the parabola
// are free and only branches present in every search sample are retained.
SearchReport search(const SearchConfig& config);
SearchReport search_mitm(const SearchConfig& config);

// Small deterministic regression target: tangent at P. It exercises the same
// terminal schemas and should recover circle(F;P), then the line through P and
// the left x-axis intersection at total cost two.
SearchReport search_tangent(const SearchConfig& config);

} // namespace neo::parabola
