#pragma once

#include <blast>

namespace blast::test {

// Absolute tolerance for "exact" numeric checks in tests, scaled to the active
// precision of blast::real. Float carries ~7 significant digits (machine eps
// ~1.2e-7), so the 1e-9 absolute bound that is appropriate for double is simply
// unattainable in single precision. The quantities checked here are O(1)..O(1e2),
// so 1e-4 comfortably absorbs float rounding while still catching real errors.
#if BLAST_USE_DOUBLES
constexpr real abs_tol = 1e-9;
#else
constexpr real abs_tol = 1e-4;
#endif

// Looser tolerance for inherently approximate checks: finite-difference vs.
// analytical derivatives and accumulated basis-function sums. Float finite
// differencing with a tiny dt loses most of its significant digits to
// cancellation, so it needs a far looser bound than double. Double keeps the
// historical 1e-5 (== BLAST_EPSILON) so those checks stay as strict as before.
#if BLAST_USE_DOUBLES
constexpr real approx_tol = 1e-5;
#else
constexpr real approx_tol = 1e-1;
#endif

// A DETERMINISTIC initial guess: straight line from start to goal in control-point
// space, plus a fixed total time.
//
// Guess::random (the default) seeds a thread_local mt19937 from std::random_device,
// so any test asserting `result.success == true` is really asserting "this task
// solves from an arbitrary random guess", which is a statement about the solver,
// not about the behaviour under test. Measured on the UR5e fixtures, that fails
// ~9/100 -- and it was 0/200 before the accept gate became a single
// success_tolerance, so the flakiness is real and CI-visible, not hypothetical.
//
// Layout comes from Bspline::compute_control(): joint-major, each joint contributing
// its free interior control points (i = 3 .. n_ctrl-4) consecutively, then T last.
// Valid only for a fully specified task (no NaN boundary values) -- which is what
// Task::stop_to_stop produces.
inline Array straight_line_guess(const Optimization& opt, const Array& start,
                                 const Array& goal, real total_time = 2.0) {
  Array     x(opt.bspline.x_len(opt.task));
  const u32 n_free = opt.bspline.n_ctrl - 6;
  Assert(x.size == (u32) opt.manip.n_joints * n_free + 1);
  u32 k = 0;
  for (int j = 0; j < opt.manip.n_joints; j++)
    for (u32 i = 0; i < n_free; i++) {
      const real a = (real) (i + 1) / (real) (n_free + 1);
      x[k++]       = start[j] + a * (goal[j] - start[j]);
    }
  x.back() = total_time;
  return x;
}

} // namespace blast::test
