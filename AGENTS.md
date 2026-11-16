# AGENTS.md

Guidance for coding agents (Claude Code, Codex and others) working in this repository. `CLAUDE.md` imports this file,
so it is the single source for both.

## What this is

Blast is a header-only C++17 library for fast trajectory optimization of robot manipulators: given a `Manipulator`, a
`Task` and constraints, it finds a time-optimal, collision-free B-spline trajectory. It accompanies the ICRA 2026 paper
*"Accelerating Trajectory Optimization by Exploiting B-Spline Gradient Structure"*. Python bindings (`blast-py`, import
`blast`) are built with nanobind.

Start with `readme.md` and `docs/guides/` (`concepts.md`, `constraints.md`, `objectives-and-guesses.md`) for the
mental model and the solver knobs.

## Building and testing

Developer build, as CI runs it (presets in `CMakePresets.json`):

```sh
cmake --preset ci-gcc            # also: ci-clang, ci (MSVC), ci-float, dev
cmake --build --preset ci-gcc
ctest --preset ci-gcc
```

One test (each test is registered as `<name>_test`):

```sh
ctest --preset ci-gcc -R test_decision_bounds --output-on-failure
```

Single precision with the in-house solver: `cmake --preset ci-float`. Benchmarks are opt-in:
`-DBLAST_BUILD_BENCHMARKS=ON`. Tests, examples and benchmarks are only built when Blast is the top-level project.

Python bindings and their tests (not run in CI, so run them yourself when you touch `python/`):

```sh
python -m venv .venv && source .venv/bin/activate
pip install . pytest
pytest tests/python
```

Documentation build: see `docs/about/contributing.md`.

## Layout

```
blast/blast              umbrella header (#include <blast>): integer aliases, MAX_JOINTS, MAX_CAPSULES, AVX2 guard
blast/blast_*.hpp        one aggregate header per area, included by the umbrella header
blast/math/              Array, Matrix, Vec3, Mat3 (AVX2 and FMA intrinsics)
blast/manipulator/       Manipulator, robot models (Kinova_Gen3, Kinova_Link6, UR5e)
blast/trajectory/        B-splines, 5th order polynomials
blast/optimization/      Optimization, constraints, objectives, initial guesses, in-house SQP
blast/world/             static obstacles, primitives, distances, GJK/EPA, broadphase
blast/gpu/               CUDA code; no CMake target builds it
blast/extern/            vendored: NLopt, VCL, Taskflow, nlohmann json
python/src/              nanobind bindings (bind_*.cpp); python/blast/ is the package
tests/                   Catch2 tests by area; tests/python/ for pytest; tests/benchmarks/
examples/                example_01 to example_04, all on the UR5e of examples/ur5e.hpp
```

## Conventions

- Format C++ with the root `.clang-format` (clang-format 17, as CI uses): 2-space indent, no column limit, aligned
  consecutive assignments and declarations, left pointer alignment. On a pull request that touches C++ files, the
  `clang-format` workflow reformats and pushes a `style: apply clang-format` commit to the branch: pull before pushing
  again.
- Write scalar literals in library code as `real`, for example `(real) 1`, never `1.0f` or a bare double literal where
  the value meets a `real`. The float build (`ci-float`) must compile and pass.
- Document public C++ API with Doxygen comments (`/** @brief ... */`). Give every Python binding a docstring in its
  `.def(...)` call. Both feed the documentation site.
- Mark functions that GPU code also calls with `blast_fn` (or `host_fn`, `device_fn`), as the existing headers do.
- Register a new test in the `TESTS` list of `tests/CMakeLists.txt` (target name, then source path); a test that is not
  listed is never built.
- Work on a branch named `<name>/<topic>` and open a pull request to `main`. Commit subjects are short, in sentence case,
  and say what changed.

## Things that are easy to get wrong

- `blast::real` is `double` by default and `float` when `BLAST_USE_DOUBLES=OFF`. Float requires
  `BLAST_USE_NATIVE_SQP=ON`, because the bundled NLopt is double only; the root `CMakeLists.txt` fails early otherwise.
  Changes to the optimizer must work with both solvers.
- The library requires AVX2 and FMA. `blast/blast` stops compilation with `#error` without them. Dev targets add
  `-march=native` or `/arch:AVX2` privately; the `install-test` CI job checks that an installed `Blast::blast` still
  compiles for a consumer that sets no flags. Do not move ISA flags into the public interface casually.
- `MAX_JOINTS` (8) and `MAX_CAPSULES` (7) are compile-time sizes of fixed arrays, overridable through `BLAST_MAX_JOINTS`
  and `BLAST_MAX_CAPSULES`. Loops must stop at `n_joints`, not at `MAX_JOINTS`.
- Do not modify `blast/extern/`. `NLOPT_LUKSAN` is forced off so that the compiled NLopt, and therefore Blast, stays
  outside the LGPL. Never re-enable it.
- The UR5e definition in `examples/ur5e.hpp` is duplicated, value for value, in `tests/python/conftest.py` and in the
  Python example of `readme.md`. Change all three together.
- Python trajectory arrays (`pos`, `vel`, `acc`) are Fortran-contiguous with shape `(n_joints, n_points)`. Inputs must
  be built with `blast.array(...)` or `dtype=blast.REAL_DTYPE` to avoid silent float32 and float64 mismatches.
- `world/test_gjk_epa.cpp` exists but is commented out in `tests/CMakeLists.txt`.

## Downstream users

`Dynamium-Lab/blast-playground` (the real-time Kinova drivers) embeds Blast as the submodule `extern/blast`. Its Gen3
driver calls `blast::dynamics`, `blast::forward_kinematics`, `blast::compute_5order_trajectory`, `blast::get_gen3_home`,
`blast::wrap2pi`, `blast::deg2rad` and `blast::rad2deg` in its 1 kHz control loop, and its MuJoCo digital twin is
validated against the Gen3 model of `blast/manipulator/Kinova_Gen3.hpp`. A change to those functions, to their
signatures, or to the Gen3 parameters changes the behaviour of a real arm: say so in the pull request.
