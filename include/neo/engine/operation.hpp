#pragma once

namespace neo::engine {

enum class OperationKind {
    Given,
    LineThrough,
    CircleCenterThrough,
    PerpendicularBisector,
    PerpendicularThrough,
    ParallelThrough,
    AngleBisector,
    Intersection
};

enum class VisibilityPolicy {
    FinalOnly,
    SelectedHelpers,
    FullyExpanded
};

struct OperationSpec {
    int atomic_cost{};
    VisibilityPolicy visibility{VisibilityPolicy::FinalOnly};
};

constexpr OperationSpec operation_spec(OperationKind kind) noexcept {
    switch (kind) {
        case OperationKind::Given:
        case OperationKind::Intersection:
            return {0, VisibilityPolicy::FullyExpanded};
        case OperationKind::LineThrough:
        case OperationKind::CircleCenterThrough:
            return {1, VisibilityPolicy::FullyExpanded};
        case OperationKind::PerpendicularBisector:
            return {3, VisibilityPolicy::SelectedHelpers};
        case OperationKind::PerpendicularThrough:
            return {3, VisibilityPolicy::FinalOnly};
        case OperationKind::ParallelThrough:
        case OperationKind::AngleBisector:
            return {4, VisibilityPolicy::FinalOnly};
    }
    return {};
}

}  // namespace neo::engine
