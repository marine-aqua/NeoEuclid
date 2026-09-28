#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace neo::engine {

// Structure-of-arrays storage is intentional: each coordinate stream is
// contiguous so MSVC and clang-cl can vectorize the sample loops.
struct SampledPoints {
    std::vector<double> x;
    std::vector<double> y;

    std::size_t size() const noexcept { return x.size(); }
    bool consistent() const noexcept { return x.size() == y.size(); }
};

struct SampledDirections {
    std::vector<double> x;
    std::vector<double> y;

    std::size_t size() const noexcept { return x.size(); }
    bool consistent() const noexcept { return x.size() == y.size(); }
};

enum class HitClass {
    Miss,
    Candidate,
    CandidateWithDegeneracies,
    Invalid
};

struct HitSummary {
    HitClass classification{HitClass::Invalid};
    std::size_t samples{};
    std::size_t hits{};
    std::size_t degeneracies{};
    double maximum_normalized_residual{};
};

// Tests points against forward target rays. Degenerate samples are retained
// separately from genuine misses so an isolated special case does not erase a
// potentially useful invariant. All arrays use the same sample ordering.
HitSummary evaluate_target_rays(const SampledPoints& points,
                                const SampledPoints& vertices,
                                const SampledDirections& directions,
                                double relative_tolerance,
                                double minimum_radius,
                                const std::vector<std::uint8_t>& degenerate = {});

}  // namespace neo::engine
