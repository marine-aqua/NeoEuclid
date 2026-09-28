# NeoEuclidSearch

An independent C++17 refactor of the numerical straightedge-and-compass search
experiments. The existing Python implementation is retained as a reference and
is not modified by this project.

The code is separated into four layers:

- `geometry`: canonical lines/circles, intersections and numerical keys;
- `state`: construction state, candidate generation and incremental expansion;
- `search`: a problem-independent weighted beam-search engine;
- `problems`: initial states, ranking and goal checks for individual problems.

Recipes are represented by compact parent references. Human-readable text is
materialised only when a result is printed, avoiding recursively copied recipe
strings in the search frontier.

## Build and test

```powershell
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

## Run

```powershell
build\neo-euclid.exe fixed-angle --target 72 --max-cost 6 --beam 500
build\neo-euclid.exe fixed-angle-mitm --target 144 --max-cost 6 --shared-cost 3 --arm-cost 2
build\neo-euclid.exe obtuse --max-cost 10 --seconds 3600 --expansions 20000
build\neo-euclid.exe obtuse-mitm --max-cost 11 --seconds 300 --shared-min 3 --shared-max 7 --records 1500
build\neo-euclid.exe obtuse-mitm --raw --max-cost 7 --seconds 300 --shared-min 0 --shared-max 4
```

The obtuse search uses a charged three-step reflection prefix. Its two helper
circles are hidden, while their intended points and the reflected line remain
available. The coordinate axes and `xy=1` are free. Only intersection branches
that exist in every configured angle sample are retained.

`obtuse-mitm` constructs a reusable shared prefix, enumerates two independent
descendant branches, and indexes each final curve by where it crosses the
target trisector ray across all samples. Matching branches are charged as the
shared prefix once, both extensions, and one final line.

With `--raw`, no reflected line or reflection macro is supplied. The initial
state has cost zero and contains only `O`, arbitrary `P`, both axes, the obtuse
angle side, and implicit free access to `xy=1` intersections.

This is a numerical heuristic search, not a proof of constructibility or
minimality. A candidate should still be checked symbolically before it is used
as a mathematical result.

## Reusable search strategies

`neo/engine/strategies.hpp` provides shared execution budgets, counters, and
policy-based randomized DFS, mask-constrained DFS, and meet-in-the-middle
drivers. The existing `beam_search` is the corresponding problem-independent
beam strategy. Problem modules retain control of state layout, candidate
generation, pruning, matching keys, and numerical confirmation so their hot
paths can remain specialized.

MITM keys are only an approximate index. Every indexed match is passed back to
the problem policy for full numerical confirmation. `fixed-angle-mitm` uses an
asymmetric split: a shared construction prefix, a descendant arm, a reusable
curve from the prefix, and the charged final join.

## Parameterized parabola search

`parabola-search` implements the focal-vertex model `y^2=4x`, with free axes,
origin, focus `F=(1,0)`, directrix `x=-1`, and the supplied parabola point `P`.
It is a separate C++ search pipeline with:

- two or three sparse search samples (`--samples`), followed by an independent
  81-sample replay over 5–175 degrees;
- explicit line/circle operation masks (`--mask LCCL`); without masks all
  binary masks up to the cost bound are scheduled;
- breadth-wise shared-prefix generation and global canonical deduplication by
  state plus remaining mask;
- terminal schemas for an existing target line, an existing target-ray point,
  a paid-curve/free-parabola intersection, and a two-curve intersection,
  followed by the final join where needed;
- a dynamic multi-threaded prefix queue (`--threads`, `--prefix-depth`).

Explicit masks may also use weighted hidden-helper macros: `N` for a
perpendicular through a point (cost 3), `B` for a perpendicular bisector
(cost 3), `P` for a parallel through a point (cost 4), and `A` for either
angle bisector of two lines (cost 4). Only the resulting line is exposed to
later operations; the macro's internal construction objects are not reusable.

```powershell
build-integrated-ninja\neo-euclid.exe parabola-search --max-cost 6 --threads 8
build-integrated-ninja\neo-euclid.exe parabola-search --tangent --max-cost 2
build-integrated-ninja\neo-euclid.exe parabola-search --max-cost 6 --mask LCCLCC
build-integrated-ninja\neo-euclid.exe parabola-search --mitm --max-cost 7 --mask CLCLCLL
build-integrated-ninja\neo-euclid.exe parabola-search --vertex-angle --max-cost 7
build-integrated-ninja\neo-euclid.exe parabola-search --beam-search --beam 1000 --vertex-angle --max-cost 7
```

The `--tangent` target is a regression oracle: it must rediscover
`circle(F;P)` followed by joining `P` to the circle's left intersection with
the x-axis, at total cost two.

`--beam-search` selects a cost-layered beam over atomic line and circle
operations. It reuses the parameterized geometry caches, state signatures,
terminal schemas, and dense replay. The beam implementation is currently
single-threaded; `--threads` continues to apply to the prefix/DFS engine.

For a seven-operation MITM mask ending in `L`, the engine interprets the mask
as `shared[0:2] | left[2:4] | right[4:6] | final L`. Both arms are enumerated
independently from the same shared prefix and indexed by their stable
intersection with the target trisector ray. A match is accepted only for two
distinct transverse curves and is replayed on the dense validation grid.

`--vertex-angle` keeps the horizontal parabola `y^2=4x` but moves the angle
vertex from its focus to the parabola vertex `O`. It restricts the input to
acute angles, uses `P=(4 cot^2(alpha), 4 cot(alpha))` on the upper branch, and
validates hits over 5–85 degrees. The focus and directrix remain free givens.

`--structured-pruning` requires the first paid operation to use the supplied
point `P` directly and accepts a terminal construction only after a
circle–parabola intersection has been consumed by a later paid operation (or
is itself the terminal point). Merely generating such an intersection does
not satisfy the condition. This history bit is part of the global state key.
