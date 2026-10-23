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
    TransferCircle,
    Intersection
};

enum class CostPolicy { Atomic, Euclidea };

// The model of the construction.
//
//  * **Classic** is the model of the first day.  A curve is drawn with one of
//    the two atomic operations -- a line through two points, or a circle about
//    one point through another -- the macros cost what they cost in atomic
//    operations, and the crossing of two curves is free.
//
//  * **All** says that every operation of the model is one step: a line, a
//    circle about a point through another, the perpendicular bisector of a
//    segment, the bisector of an angle, the perpendicular to a line through a
//    point, the parallel to a line through a point, a circle whose radius is
//    carried over from a segment -- the rigid compass, where the old model has
//    only the soft one -- and the crossing of two curves, which the classic
//    model had for nothing.
//
// The two models differ in two places and in nothing else: the price of the
// macros, which the classic model pays in atomic operations, and the price of
// a crossing, which the classic model does not pay at all.  Everything else --
// the figure, the points, the curves, the recipe -- is the same.
enum class Mode { Classic, All };

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
        case OperationKind::TransferCircle:
            return {5, VisibilityPolicy::FinalOnly};
    }
    return {};
}

constexpr int operation_cost(OperationKind kind, CostPolicy policy) noexcept {
    if (kind == OperationKind::Given || kind == OperationKind::Intersection) return 0;
    return policy == CostPolicy::Euclidea ? 1 : operation_spec(kind).atomic_cost;
}

// The price of an operation in the two models.  In the All model every
// operation is one step, the crossing of two curves among them.
constexpr int operation_cost(OperationKind kind, Mode mode) noexcept {
    if (kind == OperationKind::Given) return 0;
    if (mode == Mode::All) return 1;
    if (kind == OperationKind::Intersection) return 0;
    return operation_spec(kind).atomic_cost;
}

constexpr const char* mode_name(Mode mode) noexcept {
    return mode == Mode::All ? "all" : "classic";
}

}  // namespace neo::engine
