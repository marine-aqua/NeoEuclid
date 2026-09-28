#include "neo/problems/fixed_angle_mitm.hpp"

#include "neo/engine/strategies.hpp"
#include "neo/problems/fixed_angle.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace neo::problems {
namespace {

constexpr double kPi = 3.14159265358979323846;

struct CurveRecord {
    State state;
    std::size_t curve_index{};
    int extension_cost{};
};

class FixedAngleMitmPolicy {
public:
    using Key = std::int64_t;

    FixedAngleMitmPolicy(double target_degrees, const FixedAngleMitmConfig& config,
                         engine::StrategyStats& stats)
        : problem_(target_degrees), config_(config), stats_(stats),
          target_radians_(target_degrees * kPi / 180.0),
          target_line_{-std::sin(target_radians_), std::cos(target_radians_), 0.0, 0} {}

    std::vector<State> shared_prefixes() {
        std::vector<std::vector<State>> buckets(
            static_cast<std::size_t>(config_.shared_cost + 1));
        auto root = problem_.initial_state();
        buckets[0].push_back(std::move(root));

        for (int cost = 0; cost < config_.shared_cost; ++cost) {
            auto& states = buckets[static_cast<std::size_t>(cost)];
            std::sort(states.begin(), states.end(), [](const State& first, const State& second) {
                return first.score < second.score;
            });
            if (states.size() > config_.prefix_beam) states.resize(config_.prefix_beam);

            std::unordered_map<std::string, State> generated;
            for (const auto& state : states) {
                for (const auto& candidate : candidate_curves(state, config_.use_macros)) {
                    if (state.cost + candidate.cost > config_.shared_cost) continue;
                    ++stats_.generated;
                    auto child = add_curve(state, candidate,
                                           {config_.max_points, config_.use_macros, 8});
                    if (!child) continue;
                    child->score = problem_.score(*child);
                    const auto key = state_key(*child, 7);
                    auto known = generated.find(key);
                    if (known == generated.end() || child->score < known->second.score)
                        generated.insert_or_assign(key, std::move(*child));
                    else
                        ++stats_.duplicates;
                }
            }
            for (auto& item : generated) {
                auto& child = item.second;
                buckets[static_cast<std::size_t>(child.cost)].push_back(std::move(child));
            }
        }

        auto result = std::move(buckets[static_cast<std::size_t>(config_.shared_cost)]);
        std::sort(result.begin(), result.end(), [](const State& first, const State& second) {
            return first.score < second.score;
        });
        if (result.size() > config_.prefix_beam) result.resize(config_.prefix_beam);
        return result;
    }

    std::vector<CurveRecord> left_records(const State& base) {
        std::vector<CurveRecord> records;
        records.reserve(base.curves.size());
        for (std::size_t i = 0; i < base.curves.size(); ++i)
            records.push_back({base, i, 0});
        return records;
    }

    std::vector<CurveRecord> right_records(const State& base) {
        std::vector<CurveRecord> records;
        std::unordered_set<std::string> seen;

        auto enumerate = [&](auto&& self, const State& state, int spent) -> void {
            if (spent >= config_.arm_cost) return;
            for (const auto& candidate : candidate_curves(state, config_.use_macros)) {
                const int next_spent = spent + candidate.cost;
                if (next_spent > config_.arm_cost ||
                    base.cost + next_spent + 1 > config_.max_cost)
                    continue;
                ++stats_.generated;
                auto child = add_curve(state, candidate,
                                       {config_.max_points, config_.use_macros, 8});
                if (!child) continue;
                const auto key = state_key(*child, 7);
                if (!seen.insert(key).second) {
                    ++stats_.duplicates;
                    continue;
                }
                records.push_back({*child, child->curves.size() - 1, next_spent});
                self(self, *child, next_spent);
            }
        };
        enumerate(enumerate, base, 0);
        return records;
    }

    std::vector<Key> keys(const CurveRecord& record) const {
        std::vector<Key> result;
        const auto& curve = record.state.curves[record.curve_index];
        const double tx = std::cos(target_radians_);
        const double ty = std::sin(target_radians_);
        for (const auto& point : intersections(curve, target_line_)) {
            const double distance = point.x * tx + point.y * ty;
            const double residual = std::abs(point.x * ty - point.y * tx);
            if (distance <= kEpsilon || residual > 2e-8 * std::max(1.0, distance)) continue;
            result.push_back(static_cast<Key>(std::llround(distance * 1e7)));
        }
        return result;
    }

    std::optional<SearchResult> join(const State& base, const CurveRecord& left,
                                     const CurveRecord& right) const {
        if (left.extension_cost != 0 || right.extension_cost == 0) return std::nullopt;
        if (base.cost + right.extension_cost + 1 > config_.max_cost) return std::nullopt;

        const auto& left_curve = left.state.curves[left.curve_index];
        const auto& right_curve = right.state.curves[right.curve_index];
        const double tx = std::cos(target_radians_);
        const double ty = std::sin(target_radians_);
        for (const auto& point : intersections(left_curve, right_curve)) {
            const double distance = point.x * tx + point.y * ty;
            const double residual = std::abs(point.x * ty - point.y * tx);
            if (distance <= kEpsilon || residual > 2e-8 * std::max(1.0, distance)) continue;

            const auto wanted = point_key(point, 8);
            for (const auto& known : right.state.points) {
                if (!(point_key(known, 8) == wanted)) continue;
                const auto verified = problem_.goal(right.state, config_.max_cost);
                if (!verified || verified->total_cost != right.state.cost + 1)
                    continue;
                return SearchResult{right.state, *verified, stats_.expanded,
                                    stats_.generated};
            }
        }
        return std::nullopt;
    }

private:
    FixedAngle problem_;
    const FixedAngleMitmConfig& config_;
    engine::StrategyStats& stats_;
    double target_radians_;
    Line target_line_;
};

}  // namespace

std::optional<SearchResult> fixed_angle_mitm(double target_degrees,
                                             const FixedAngleMitmConfig& config) {
    engine::StrategyStats stats;
    FixedAngleMitmPolicy policy(target_degrees, config, stats);
    engine::StrategyLimits limits{config.max_cost, config.max_states,
                                  config.time_limit_seconds};
    auto result = engine::meet_in_the_middle<State, CurveRecord,
                                              FixedAngleMitmPolicy::Key,
                                              SearchResult>(limits, policy, stats);
    if (result) {
        result->expanded = stats.expanded;
        result->generated = stats.generated;
    }
    return result;
}

}  // namespace neo::problems
