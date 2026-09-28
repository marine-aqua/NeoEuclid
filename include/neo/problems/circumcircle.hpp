#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace neo::circumcircle {

struct Report {
    bool found{};
    int cost{};
    std::string left;
    std::string right;
    std::size_t generated{};
    std::size_t states{};
    std::size_t target_curves{};
    double elapsed_seconds{};
};

Report search_mitm(double seconds, std::size_t max_states, std::size_t max_points);

}  // namespace neo::circumcircle
