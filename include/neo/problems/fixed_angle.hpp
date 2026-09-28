#pragma once

#include "neo/search.hpp"

namespace neo::problems {

class FixedAngle final : public SearchProblem {
public:
    explicit FixedAngle(double target_degrees);
    State initial_state() const override;
    double score(const State& state) const override;
    std::optional<Goal> goal(const State& state, int max_cost) const override;

private:
    bool exact_target_check(const Point& point) const;
    double target_degrees_;
    double target_radians_;
};

}  // namespace neo::problems

