#include "neo/geometry.hpp"
#include "neo/problems/fixed_angle.hpp"
#include "neo/problems/obtuse.hpp"
#include "neo/problems/parabola.hpp"
#include "neo/search.hpp"
#include "neo/state.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    check(std::abs(actual - expected) <= tolerance, message);
}

void test_line_normalization() {
    const auto line = neo::line_through({0.0, 0.0, 0}, {2.0, 0.0, 1});
    check(line.has_value(), "line should exist");
    near(line->a, 0.0, 1e-12, "horizontal line a");
    near(line->b, 1.0, 1e-12, "horizontal line b");
    near(line->c, 0.0, 1e-12, "horizontal line c");
}

void test_line_circle_intersection() {
    const neo::Curve line = neo::Line{0.0, 1.0, 0.0, 0};
    const neo::Curve circle = neo::Circle{0.0, 0.0, 1.0, 0};
    const auto points = neo::intersections(line, circle);
    check(points.size() == 2, "line and circle should meet twice");
    near(std::abs(points[0].x), 1.0, 1e-12, "intersection x");
    near(points[0].y, 0.0, 1e-12, "intersection y");
}

void test_circle_circle_tangent() {
    const auto points = neo::intersections(neo::Circle{0.0, 0.0, 1.0, 0},
                                           neo::Circle{2.0, 0.0, 1.0, 0});
    check(points.size() == 1, "tangent circles should meet once");
    near(points[0].x, 1.0, 1e-12, "tangent x");
}

void test_quantized_keys() {
    check(neo::point_key({1.0, 2.0, 0}) == neo::point_key({1.0 + 1e-10, 2.0, 1}),
          "small numerical noise should merge");
}

void test_state_rejects_duplicate_curve() {
    neo::problems::FixedAngle problem(72.0);
    const auto state = problem.initial_state();
    const neo::Candidate duplicate{neo::Line{0.0, 1.0, 0.0, 0}, 1,
                                   neo::OperationKind::LineThrough, 0, 1};
    check(!neo::add_curve(state, duplicate, {}).has_value(), "duplicate curve must be rejected");
}

void test_small_search_misses_target() {
    neo::problems::FixedAngle problem(72.0);
    neo::SearchConfig config;
    config.max_cost = 1;
    config.beam_width = 20;
    config.verbose = false;
    check(!neo::beam_search(problem, config).has_value(), "72 degrees needs more than one step");
}

void test_known_72_degree_construction() {
    neo::problems::FixedAngle problem(72.0);
    neo::SearchConfig config;
    config.max_cost = 6;
    config.beam_width = 500;
    config.max_points = 28;
    config.use_macros = false;
    config.verbose = false;
    const auto result = neo::beam_search(problem, config);
    check(result.has_value(), "known 72-degree construction should be found");
    check(result->goal.total_cost == 6, "known construction should cost six");
    check(result->state.steps.size() == 5, "final line is reconstructed as the sixth step");
}

void test_obtuse_engine_smoke() {
    neo::obtuse::SearchConfig config;
    config.max_cost = 7;
    config.expansions_per_restart = 100;
    config.terminal_states_per_restart = 20;
    config.time_limit_seconds = 0.05;
    config.verbose = false;
    config.sample_degrees = {103.0, 120.0, 137.0, 161.0};
    const auto report = neo::obtuse::search(config);
    check(report.generated > 0, "obtuse engine should expand its prefixed state");
    check(report.elapsed_seconds >= 0.04, "obtuse engine should respect its time budget");
}

void test_obtuse_model_rules() {
    const auto diagnostics = neo::obtuse::diagnose_model({103.0, 120.0, 137.0, 161.0});
    check(diagnostics.prefix_cost == 3, "reflection prefix must cost three");
    check(diagnostics.visible_curves == 4, "only axes, angle side, and reflected line are visible");
    check(diagnostics.stable_hyperbola_points == 2, "reflected line should have two stable H intersections");
    check(diagnostics.helper_circles_hidden, "prefix helper circles must not be reusable");
    check(diagnostics.unstable_branch_rejected, "sample-specific intersections must be rejected");
    check(diagnostics.quartic_solver_verified, "circle-hyperbola quartic should return all real roots");
    check(diagnostics.macro_costs_verified, "macro costs must be 3, 4, 3, and 4");
    check(diagnostics.existing_objects_reused_without_recharge,
          "existing construction objects must be reusable without being charged again");
    check(diagnostics.target_ray_join_verified,
          "MITM curves must join at one coherent target-ray point across samples");
    check(diagnostics.origin_branch_verified,
          "the origin must remain one coherent branch across all angle samples");
    check(diagnostics.near_duplicate_curve_rejected,
          "MITM must not join numerically perturbed copies of the same curve");
    check(diagnostics.transverse_join_verified,
          "MITM must accept transverse target meets and reject coincident ones");
    check(diagnostics.raw_initial_state_verified,
          "raw mode must start at cost zero without a reflected line");
}

void test_obtuse_mitm_smoke() {
    neo::obtuse::MeetConfig config;
    config.max_cost = 9;
    config.shared_cost_min = 3;
    config.shared_cost_max = 5;
    config.branch_records_per_base = 50;
    config.time_limit_seconds = 0.08;
    config.verbose = false;
    config.sample_degrees = {103.0, 120.0, 137.0, 161.0};
    const auto report = neo::obtuse::meet_in_the_middle(config);
    check(report.restarts > 0, "MITM should construct at least one shared base");
    check(report.generated > 0, "MITM should enumerate branch candidates");
}

void test_parallel_parabola_terminal_schemas() {
    neo::parabola::SearchConfig config;
    config.max_cost = 2;
    config.prefix_depth = 1;
    config.threads = 2;
    config.time_limit_seconds = 2.0;
    config.max_states = 100000;
    config.verbose = false;
    const auto report = neo::parabola::search_tangent(config);
    check(report.found, "parameterized C++ engine should recover the cost-two tangent");
    check(report.densely_verified, "a sparse-sample hit must pass dense replay");
    check(report.cost == 2, "focus-circle tangent construction should cost two");
    check(report.terminal_schema == "two-curves-intersection+join",
          "circle/axis intersection followed by a join should use the pair terminal schema");
    check(report.prefixes > 0 && report.generated > 0,
          "layered prefixes and worker expansion should both run");
}

void test_parabola_mitm_smoke() {
    neo::parabola::SearchConfig config;
    config.max_cost = 7;
    config.threads = 2;
    config.time_limit_seconds = 1.0;
    config.max_states = 5000;
    config.state_cache_entries = 0;
    config.geometry_cache_entries = 2000;
    config.masks = {"CLCLCLL"};
    config.verbose = false;
    const auto report = neo::parabola::search_mitm(config);
    check(report.prefixes > 0, "MITM should generate shared prefixes");
    check(report.mitm_records > 0, "MITM should enumerate independent arm records");
}

void test_parabola_weighted_macros_smoke() {
    neo::parabola::SearchConfig config;
    config.max_cost = 4;
    config.prefix_depth = 1;
    config.threads = 2;
    config.time_limit_seconds = 1.0;
    config.max_states = 10000;
    config.masks = {"N", "B", "P", "A"};
    config.verbose = false;
    const auto report = neo::parabola::search(config);
    check(report.prefixes > 0 && report.generated > 0,
          "weighted perpendicular, bisector, parallel, and angle macros should generate candidates");
}

}  // namespace

int main() {
    try {
        test_line_normalization();
        test_line_circle_intersection();
        test_circle_circle_tangent();
        test_quantized_keys();
        test_state_rejects_duplicate_curve();
        test_small_search_misses_target();
        test_known_72_degree_construction();
        test_obtuse_engine_smoke();
        test_obtuse_model_rules();
        test_obtuse_mitm_smoke();
        test_parallel_parabola_terminal_schemas();
        test_parabola_mitm_smoke();
        test_parabola_weighted_macros_smoke();
        std::cout << "All tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failed: " << error.what() << '\n';
        return 1;
    }
}
