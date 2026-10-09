#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include <cmath>
#include <random>
#include "test_helper.hpp"

using namespace blast;

// External-collision gradients against central finite differences, world made of capsules only.
// Float differencing loses most significant digits to cancellation (as in test_helper.hpp), so
// it takes a larger step and a looser bound.
#if BLAST_USE_DOUBLES
constexpr real finite_difference_step      = 1e-6;
constexpr real finite_difference_tolerance = 1e-3;
#else
constexpr real finite_difference_step      = 1e-3;
constexpr real finite_difference_tolerance = 1e-1;
#endif

using ConstraintsAndGradientsFunction = void (*)(const Array&, Optimization&, Array&, Matrix&);

static World capsule_world() {
  World world;
  world.add_capsule({Vec3{0.35, 0.25, 0.35}, Vec3{0.35, 0.25, 0.55}, 0.06});
  world.add_capsule({Vec3{-0.3, 0.45, 0.30}, Vec3{-0.1, 0.45, 0.30}, 0.05});
  world.add_capsule({Vec3{0.1, -0.4, 0.5}, Vec3{0.1, -0.4, 0.7}, 0.07});
  return world;
}

static const Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
static const Array end   = {-0.5825, -1.0700, 1.3892, -1.3196, 0.9927, 1.14};

static void check_gradients(ConstraintsAndGradientsFunction constraints_and_gradients, const char* method_name) {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, end));
  opt.world = capsule_world();
  create_static_bounding_volume_hierarchy(opt.world, opt.world.static_bounding_volume_hierarchy); // the broadphase methods read it; the caller builds it
  auto& selection    = opt.constraints;
  selection.position = selection.velocity = selection.acceleration = selection.torque = selection.tool_speed = selection.self_collisions = false;
  selection.external_collisions                                                                                                          = true;
  initialize_optimization_with_segments(&opt);
  n_con_with_segments(&opt);

  const u32                        n_variables = opt.bspline.x_len(opt.task), n_constraints = opt.constraints.n_constraints;
  std::mt19937                     random_engine(7);
  std::normal_distribution<double> noise(0, 0.3);
  int                              mismatched = 0, checked = 0;
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
    for (u32 variable = 0; variable + 1 < n_variables; variable++) { // control points; the duration column has no collision term
      Array x_plus  = x;
      Array x_minus = x;
      x_plus[variable] += finite_difference_step;
      x_minus[variable] -= finite_difference_step;
      Array  constraints_plus(n_constraints), constraints_minus(n_constraints);
      Matrix no_gradient;
      constraints_and_gradients(x_plus, opt, constraints_plus, no_gradient);
      constraints_and_gradients(x_minus, opt, constraints_minus, no_gradient);
      for (u32 row = 0; row < n_constraints; row++) {
        const real finite_difference = (constraints_plus[row] - constraints_minus[row]) / (2 * finite_difference_step);
        // Rows are a max over samples: where the argmax switches inside the step the finite
        // difference mixes two pieces, so allow a few, but a systematic error shows up in hundreds.
        mismatched += std::abs(gradient(variable, row) - finite_difference) > finite_difference_tolerance * (std::abs(finite_difference) + finite_difference_tolerance);
        checked++;
      }
    }
  }
  INFO(method_name << " mismatching entries: " << mismatched << "/" << checked);
  CHECK(mismatched <= checked / 500);
}

TEST_CASE("collision gradients against world capsules match finite differences", "[Optimization]") {
  check_gradients(constraints_and_gradients_with_segments, "with_segments");
  check_gradients(constraints_and_gradients_with_broadphase, "broadphase");
  check_gradients(constraints_and_gradients_with_double_broadphase, "double_broadphase");
}

// The same for the point-based analytical-dynamics method, through its NLopt entry point.
TEST_CASE("analytical_dynamics collision gradients against world capsules match finite differences", "[Optimization]") {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, end));
  opt.world        = capsule_world();
  auto& selection  = opt.constraints;
  selection.torque = selection.tool_speed = selection.self_collisions = false;
  selection.position = selection.velocity = selection.acceleration = selection.external_collisions = true; // the gradient fill assumes PVA rows exist
  initialize_optimization(&opt);
  n_con(&opt);

  const u32                        n_variables = opt.bspline.x_len(opt.task), n_constraints = opt.constraints.n_constraints;
  std::mt19937                     random_engine(7);
  std::normal_distribution<double> noise(0, 0.3);
  int                              mismatched = 0, checked = 0;
  for (int trial = 0; trial < 10; trial++) {
    Array x  = blast::guess_straight_line(&opt);
    x.back() = 1.0 + 0.2 * trial;
    for (u32 variable = 0; variable + 1 < n_variables; variable++)
      x[variable] += noise(random_engine);
    Array constraints(n_constraints), gradient(n_constraints * n_variables); // NLopt layout: gradient[row * n_variables + variable]
    nlopt_constraints_with_analytical_dynamics(n_constraints, constraints.data, n_variables, x.data, gradient.data, &opt);
    for (u32 variable = 0; variable + 1 < n_variables; variable++) {
      Array x_plus  = x;
      Array x_minus = x;
      x_plus[variable] += finite_difference_step;
      x_minus[variable] -= finite_difference_step;
      Array constraints_plus(n_constraints), constraints_minus(n_constraints);
      nlopt_constraints_with_analytical_dynamics(n_constraints, constraints_plus.data, n_variables, x_plus.data, nullptr, &opt);
      nlopt_constraints_with_analytical_dynamics(n_constraints, constraints_minus.data, n_variables, x_minus.data, nullptr, &opt);
      for (u32 row = 0; row < n_constraints; row++) {
        const real finite_difference = (constraints_plus[row] - constraints_minus[row]) / (2 * finite_difference_step);
        mismatched += std::abs(gradient[row * n_variables + variable] - finite_difference) > finite_difference_tolerance * (std::abs(finite_difference) + finite_difference_tolerance);
        checked++;
      }
    }
  }
  INFO("analytical_dynamics mismatching entries: " << mismatched << "/" << checked);
  CHECK(mismatched <= checked / 500);
}
