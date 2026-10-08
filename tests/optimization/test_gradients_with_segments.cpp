#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include <cmath>
#include <random>
#include "test_helper.hpp"

using namespace blast;

// Analytic constraint gradients against central finite differences, every constraint type on, for
// each segment-based method; the broadphase methods also against with_segments.
// Float differencing loses most significant digits to cancellation (as in test_helper.hpp), so
// it takes a larger step and a looser bound; the duration bug is still caught (~13% of entries).
#if BLAST_USE_DOUBLES
constexpr real finite_difference_step      = 1e-6;
constexpr real finite_difference_tolerance = 1e-3;
#else
constexpr real finite_difference_step      = 1e-3;
constexpr real finite_difference_tolerance = 1e-1;
#endif

using ConstraintsAndGradientsFunction = void (*)(const Array&, Optimization&, Array&, Matrix&);

static void check_gradients(ConstraintsAndGradientsFunction constraints_and_gradients, const char* method_name) {
  Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {-0.5825, -1.0700, 1.3892, -1.3196, 0.9927, 1.14};
  World world;
  world.add_sphere(Vec3{0.35, 0.25, 0.45}, 0.08);
  world.add_sphere(Vec3{-0.2, 0.45, 0.30}, 0.06);
  world.add_sphere(Vec3{0.1, -0.4, 0.6}, 0.07);
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, end));
  opt.world = world;
  create_static_bounding_volume_hierarchy(opt.world, opt.world.static_bounding_volume_hierarchy); // the broadphase methods read it; the caller builds it
  auto& selection    = opt.constraints;
  selection.position = selection.velocity = selection.acceleration = selection.torque = selection.tool_speed = true;
  selection.self_collisions = selection.external_collisions = true;
  initialize_optimization_with_segments(&opt);
  n_con_with_segments(&opt);

  const u32 n_variables = opt.bspline.x_len(opt.task), n_constraints = opt.constraints.n_constraints;
  std::mt19937                     random_engine(7);
  std::normal_distribution<double> noise(0, 0.3);
  int mismatched_control_points = 0, checked_control_points = 0;
  int mismatched_durations = 0, checked_durations = 0; // the duration column on its own
  for (int trial = 0; trial < 10; trial++) {
    Array x  = blast::guess_straight_line(&opt);
    x.back() = 1.0 + 0.2 * trial;
    for (u32 variable = 0; variable + 1 < n_variables; variable++)
      x[variable] += noise(random_engine);
    Array  constraints(n_constraints);
    Matrix gradient(n_variables, n_constraints);
    for (u32 entry = 0; entry < n_variables * n_constraints; entry++)
      gradient.data[entry] = 0;
    constraints_and_gradients(x, opt, constraints, gradient);
    if (constraints_and_gradients != constraints_and_gradients_with_segments) { // same values and gradient as with_segments
      Array  segments_constraints(n_constraints);
      Matrix segments_gradient(n_variables, n_constraints);
      for (u32 entry = 0; entry < n_variables * n_constraints; entry++)
        segments_gradient.data[entry] = 0;
      constraints_and_gradients_with_segments(x, opt, segments_constraints, segments_gradient);
      INFO(method_name);
      CHECK(is_close(constraints, segments_constraints));
      CHECK(is_close(gradient, segments_gradient));
    }
    for (u32 variable = 0; variable < n_variables; variable++) {
      const real step    = (variable + 1 == n_variables) ? finite_difference_step * x[variable] : finite_difference_step;
      Array      x_plus  = x;
      Array      x_minus = x;
      x_plus[variable] += step;
      x_minus[variable] -= step;
      Array  constraints_plus(n_constraints), constraints_minus(n_constraints);
      Matrix no_gradient;
      constraints_and_gradients(x_plus, opt, constraints_plus, no_gradient);
      constraints_and_gradients(x_minus, opt, constraints_minus, no_gradient);
      for (u32 row = 0; row < n_constraints; row++) {
        const real finite_difference = (constraints_plus[row] - constraints_minus[row]) / (2 * step);
        // Rows are a max over samples: where the argmax switches inside +/-step the finite
        // difference mixes two pieces, so allow a few, but a systematic error shows up in hundreds.
        const bool mismatch = std::abs(gradient(variable, row) - finite_difference) > finite_difference_tolerance * (std::abs(finite_difference) + finite_difference_tolerance);
        if (variable + 1 == n_variables)
          mismatched_durations += mismatch, checked_durations++;
        else
          mismatched_control_points += mismatch, checked_control_points++;
      }
    }
  }
  INFO(method_name << " mismatching entries: control points " << mismatched_control_points << "/" << checked_control_points
                   << ", duration " << mismatched_durations << "/" << checked_durations);
  CHECK(mismatched_control_points <= checked_control_points / 500);
  CHECK(mismatched_durations <= checked_durations / 100);
}

TEST_CASE("constraint gradients match finite differences", "[Optimization]") {
  check_gradients(constraints_and_gradients_with_segments, "with_segments");
  check_gradients(constraints_and_gradients_with_broadphase, "broadphase");
  check_gradients(constraints_and_gradients_with_double_broadphase, "double_broadphase");
}
