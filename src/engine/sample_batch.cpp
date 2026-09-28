#include "neo/engine/sample_batch.hpp"

#include <limits>

namespace neo::engine {

HitSummary evaluate_target_rays(const SampledPoints& points,
                                const SampledPoints& vertices,
                                const SampledDirections& directions,
                                double relative_tolerance,
                                double minimum_radius,
                                const std::vector<std::uint8_t>& degenerate) {
    HitSummary result;
    if (!points.consistent() || !vertices.consistent() || !directions.consistent() ||
        points.size() != vertices.size() || points.size() != directions.size() ||
        (!degenerate.empty() && degenerate.size() != points.size()) ||
        relative_tolerance < 0.0 || minimum_radius < 0.0) {
        return result;
    }

    result.samples = points.size();
    if (result.samples == 0) return result;

    bool missed = false;
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (!degenerate.empty() && degenerate[i] != 0) {
            ++result.degeneracies;
            continue;
        }

        const double dx = points.x[i] - vertices.x[i];
        const double dy = points.y[i] - vertices.y[i];
        const double direction_norm = std::hypot(directions.x[i], directions.y[i]);
        const double radius = std::hypot(dx, dy);
        if (!std::isfinite(radius) || !std::isfinite(direction_norm) ||
            radius <= minimum_radius || direction_norm <= 0.0) {
            missed = true;
            continue;
        }

        const double along =
            (dx * directions.x[i] + dy * directions.y[i]) / direction_norm;
        const double cross =
            (dx * directions.y[i] - dy * directions.x[i]) / direction_norm;
        const double normalized = std::abs(cross) / std::max(1.0, radius);
        result.maximum_normalized_residual =
            std::max(result.maximum_normalized_residual, normalized);
        if (along > minimum_radius && normalized <= relative_tolerance) {
            ++result.hits;
        } else {
            missed = true;
        }
    }

    if (missed) {
        result.classification = HitClass::Miss;
    } else if (result.hits == 0) {
        result.classification = HitClass::Invalid;
    } else if (result.degeneracies != 0) {
        result.classification = HitClass::CandidateWithDegeneracies;
    } else {
        result.classification = HitClass::Candidate;
    }
    return result;
}

}  // namespace neo::engine
