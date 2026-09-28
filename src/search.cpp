#include "neo/search.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace neo {

std::optional<SearchResult> beam_search(const SearchProblem& problem,
                                        const SearchConfig& config) {
    std::vector<std::vector<State>> buckets(static_cast<std::size_t>(config.max_cost + 1));
    auto initial = problem.initial_state();
    initial.score = problem.score(initial);
    buckets[0].push_back(std::move(initial));
    std::size_t total_expanded = 0;
    std::size_t total_generated = 0;
    const auto started = std::chrono::steady_clock::now();

    for (int cost = 0; cost <= config.max_cost; ++cost) {
        auto& states = buckets[static_cast<std::size_t>(cost)];
        std::sort(states.begin(), states.end(), [](const State& a, const State& b) {
            if (a.score != b.score) return a.score < b.score;
            return a.points.size() < b.points.size();
        });
        if (states.size() > config.beam_width) states.resize(config.beam_width);
        if (cost == config.max_cost) continue;
        std::unordered_map<std::string, State> generated;

        for (const auto& state : states) {
            ++total_expanded;
            std::unordered_set<QuantizedCurve, CurveKeyHash> local;
            for (const auto& candidate : candidate_curves(state, config.use_macros)) {
                if (state.cost + candidate.cost > config.max_cost) continue;
                if (!local.insert(curve_key(candidate.curve, 7)).second) continue;
                auto child = add_curve(state, candidate,
                                       {config.max_points, config.use_macros, 8});
                if (!child) continue;
                child->score = problem.score(*child);
                ++total_generated;
                if (const auto hit = problem.goal(*child, config.max_cost)) {
                    return SearchResult{std::move(*child), *hit, total_expanded, total_generated};
                }
                const auto key = state_key(*child, 7);
                const auto existing = generated.find(key);
                if (existing == generated.end() || child->score < existing->second.score) {
                    generated.insert_or_assign(key, std::move(*child));
                }
            }
        }
        for (auto& [key, child] : generated) {
            (void)key;
            buckets[static_cast<std::size_t>(child.cost)].push_back(std::move(child));
        }
        if (config.verbose) {
            const auto seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            std::cout << "cost=" << cost << " expanded=" << states.size()
                      << " generated=" << generated.size() << " elapsed="
                      << seconds << "s\n";
        }
    }
    return std::nullopt;
}

}  // namespace neo

