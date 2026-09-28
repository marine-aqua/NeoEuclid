#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace neo::engine {

struct StrategyLimits {
    int max_cost{6};
    std::size_t max_states{2'000'000};
    double time_limit_seconds{60.0};
};

struct StrategyStats {
    std::size_t expanded{};
    std::size_t generated{};
    std::size_t duplicates{};
    double elapsed_seconds{};
};

class SearchBudget {
public:
    explicit SearchBudget(StrategyLimits limits)
        : limits_(limits), started_(std::chrono::steady_clock::now()) {}

    bool exhausted(std::size_t generated) const {
        return generated >= limits_.max_states || elapsed_seconds() >= limits_.time_limit_seconds;
    }

    double elapsed_seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - started_).count();
    }

    const StrategyLimits& limits() const noexcept { return limits_; }

private:
    StrategyLimits limits_;
    std::chrono::steady_clock::time_point started_;
};

// Policy contract:
//   optional<Result> goal(const State&)
//   vector<State> expand(const State&)
//   int cost(const State&)
//   string key(const State&)
template <class State, class Result, class Policy, class Random>
std::optional<Result> randomized_depth_first(State root,
                                             const StrategyLimits& limits,
                                             Policy& policy,
                                             Random& random,
                                             StrategyStats& stats) {
    SearchBudget budget(limits);
    std::unordered_set<std::string> seen;

    auto visit = [&](auto&& self, State state) -> std::optional<Result> {
        if (budget.exhausted(stats.generated)) return std::nullopt;
        if (auto result = policy.goal(state)) return result;
        if (policy.cost(state) >= limits.max_cost) return std::nullopt;

        ++stats.expanded;
        auto children = policy.expand(state);
        std::shuffle(children.begin(), children.end(), random);
        for (auto& child : children) {
            ++stats.generated;
            if (policy.cost(child) > limits.max_cost) continue;
            if (!seen.insert(policy.key(child)).second) {
                ++stats.duplicates;
                continue;
            }
            if (auto result = self(self, std::move(child))) return result;
            if (budget.exhausted(stats.generated)) break;
        }
        return std::nullopt;
    };

    seen.insert(policy.key(root));
    auto result = visit(visit, std::move(root));
    stats.elapsed_seconds = budget.elapsed_seconds();
    return result;
}

// Policy contract is the same except expand(state, operation) accepts one
// operation token from the mask. Weighted operation costs remain in State.
template <class State, class Result, class Policy>
std::optional<Result> masked_depth_first(State root,
                                         const std::string& mask,
                                         const StrategyLimits& limits,
                                         Policy& policy,
                                         StrategyStats& stats) {
    SearchBudget budget(limits);
    std::unordered_set<std::string> seen;

    auto visit = [&](auto&& self, State state, std::size_t position)
        -> std::optional<Result> {
        if (budget.exhausted(stats.generated)) return std::nullopt;
        if (auto result = policy.goal(state)) return result;
        if (position >= mask.size() || policy.cost(state) >= limits.max_cost)
            return std::nullopt;

        ++stats.expanded;
        for (auto& child : policy.expand(state, mask[position])) {
            ++stats.generated;
            if (policy.cost(child) > limits.max_cost) continue;
            const auto key = policy.key(child) + '|' + mask.substr(position + 1);
            if (!seen.insert(key).second) {
                ++stats.duplicates;
                continue;
            }
            if (auto result = self(self, std::move(child), position + 1)) return result;
            if (budget.exhausted(stats.generated)) break;
        }
        return std::nullopt;
    };

    seen.insert(policy.key(root) + '|' + mask);
    auto result = visit(visit, std::move(root), 0);
    stats.elapsed_seconds = budget.elapsed_seconds();
    return result;
}

// Policy contract:
//   vector<Base> shared_prefixes()
//   vector<Record> left_records(const Base&)
//   vector<Record> right_records(const Base&)
//   vector<Key> keys(const Record&)
//   optional<Result> join(const Base&, const Record&, const Record&)
// Key must be hashable. Numerical confirmation belongs in join(), never in the
// quantized index key alone.
template <class Base, class Record, class Key, class Result, class Policy>
std::optional<Result> meet_in_the_middle(const StrategyLimits& limits,
                                         Policy& policy,
                                         StrategyStats& stats) {
    SearchBudget budget(limits);
    for (auto& base : policy.shared_prefixes()) {
        if (budget.exhausted(stats.generated)) break;
        ++stats.expanded;
        auto left = policy.left_records(base);
        auto right = policy.right_records(base);
        stats.generated += left.size() + right.size();

        std::unordered_multimap<Key, std::size_t> index;
        for (std::size_t i = 0; i < left.size(); ++i)
            for (const auto& key : policy.keys(left[i])) index.emplace(key, i);

        for (const auto& right_record : right) {
            for (const auto& key : policy.keys(right_record)) {
                const auto range = index.equal_range(key);
                for (auto it = range.first; it != range.second; ++it)
                    if (auto result = policy.join(base, left[it->second], right_record)) {
                        stats.elapsed_seconds = budget.elapsed_seconds();
                        return result;
                    }
            }
        }
    }
    stats.elapsed_seconds = budget.elapsed_seconds();
    return std::nullopt;
}

}  // namespace neo::engine
