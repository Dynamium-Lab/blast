#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>
#include "test_helper.hpp"

using namespace blast;

// Returns a representative stop-to-stop task for the UR5e.
inline blast::Task make_UR5e_task() {
  blast::Array start = {(real) 1.94822, (real) 0.473555, -(real) 0.0255247, -(real) 0.448375, (real) 0.370356, -(real) 3.12883};
  blast::Array end   = {(real) 2.5825, (real) 0.0700, -(real) 0.3892, (real) 0.3196, (real) 0.9927, -(real) 3.17328};
  return blast::Task::stop_to_stop(start, end);
}

Optimization make_base_optimization() {
  Manipulator manip = make_UR5e();

  int     n_pts  = 105;
  int     p      = 5;
  int     n_ctrl = 12;
  Bspline bspline(n_ctrl, n_pts, p, manip.n_joints);

  auto task = make_UR5e_task();

  return Optimization(manip, Task::stop_to_stop(Array(6, 0), Array(6, 1)), bspline);
}

TEST_CASE("time objective", "[Optimization][Objective]") {
  auto opt = make_base_optimization();

  opt.objective.time_weight = 1;
  Array x(opt.bspline.x_len(opt.task), 1);

  Array grad(x.size);

  SECTION("basic time objective") {
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 1);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == 0);
    }
    CHECK(grad.back() == Approx(1).epsilon(test::approx_tol));

    // Changing objective value
    x.back() = 2;

    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 2);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == 0);
    }
    CHECK(grad.back() == Approx(1).epsilon(test::approx_tol));

    // Changing time weight
    opt.objective.time_weight = 2;

    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 4);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == 0);
    }
    CHECK(grad.back() == Approx(2).epsilon(test::approx_tol));

    // Zero-ing time weight
    opt.objective.time_weight = 0;

    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 0);
    for (int i = 0; i < grad.size; i++) {
      CHECK(grad[i] == 0);
    }
  }
}

real custom_objective_velocity_sum(Optimization* opt, Array* x) {
  real result = 0;

  opt->bspline.compute_trajectory((*x), opt->task);
  for (int i = 0; i < opt->bspline.traj.vel.cols; i++) {
    for (int j = 0; j < opt->bspline.traj.vel.rows; j++) {
      result += std::abs(opt->bspline.traj.vel(j, i));
    }
  }

  return result;
}

real custom_objective_acceleration_sum(Optimization* opt, Array* x) {
  real result = 0;

  opt->bspline.compute_trajectory((*x), opt->task);
  for (int i = 0; i < opt->bspline.traj.acc.cols; i++) {
    for (int j = 0; j < opt->bspline.traj.acc.rows; j++) {
      result += std::abs(opt->bspline.traj.acc(j, i));
    }
  }

  return result;
}

TEST_CASE("Custom objective", "[Optimization][Objective]") {
  auto opt = make_base_optimization();

  opt.objective.time_weight = 0;
  opt.objective.add_custom_objective(custom_objective_velocity_sum, 1);

  Array x(opt.bspline.x_len(opt.task), 1);

  Array grad(x.size);

  SECTION("velocity objective") {
    real  original_objective_value = objective_function(x.size, x.data, grad.data, &opt);
    Array original_gradient_value  = grad;

    // Changing time value
    x.back() = 2;
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == (real) 0.5 * original_objective_value);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == (real) 0.5 * original_gradient_value[i]);
    }
    CHECK(grad.back() == Approx((real) 0.25 * original_gradient_value.back()).epsilon(test::approx_tol));

    x.back() = 1;

    // Changing time weight
    opt.objective.k_extra_objectives[0] = 2;
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 2 * original_objective_value);
    for (int i = 0; i < grad.size; i++) {
      CHECK(grad[i] == 2 * original_gradient_value[i]);
    }

    // Zero-ing custom weight
    opt.objective.k_extra_objectives[0] = 0;
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 0);
    for (int i = 0; i < grad.size; i++) {
      CHECK(grad[i] == 0);
    }
  }
}

Array custom_constraint_time(Optimization* opt, Array* x) {
  Array result(1);

  auto* max_time = (real*) opt->custom_data;

  return Array(1, x->back() - (*max_time));
}

TEST_CASE("Smoothness optimization using custom velocity objective", "[Optimization][Objective]") {
  auto opt = make_base_optimization();

  opt.objective.time_weight = 0;

  auto result1 = optimize(&opt);

  opt.guess.type      = opt.guess.custom;
  opt.guess.initial_x = result1.x;

  opt.objective.add_custom_objective(custom_objective_velocity_sum, 1);
  opt.objective.time_weight = 0;

  opt.constraints.add_global_constraint(custom_constraint_time, 1);

  real custom_value = result1.x.back();
  opt.custom_data   = &custom_value;

  opt.method = OptimizationMethod::baseline;

  auto result2 = optimize(&opt);

  real vel_sum_1 = custom_objective_velocity_sum(&opt, &result1.x);
  real vel_sum_2 = custom_objective_velocity_sum(&opt, &result2.x);

  real acc_sum_1 = custom_objective_acceleration_sum(&opt, &result1.x);
  real acc_sum_2 = custom_objective_acceleration_sum(&opt, &result2.x);

  // std::cout << "Result 1 velocity smoothness " << vel_sum_1 << std::endl;
  // std::cout << "Result 2 velocity smoothness " << vel_sum_2 << std::endl;

  // std::cout << "Result 1 acceleration smoothness " << acc_sum_1 << std::endl;
  // std::cout << "Result 2 acceleration smoothness " << acc_sum_2 << std::endl;

  // std::cout << "Result 1 time " << result1.x.back() << std::endl;
  // std::cout << "Result 2 time " << result2.x.back() << std::endl;

  CHECK(vel_sum_2 <= vel_sum_1);
}

// Sum of |vel| with its analytical gradient, for a stop-to-stop task (the pinned control points do not depend on T).
real custom_objective_velocity_sum_with_gradient(Optimization* opt, Array* x, Array* grad) {
  const real sum = custom_objective_velocity_sum(opt, x); // computes the trajectory
  if (grad) {
    const int  n_free = (int) opt->bspline.n_ctrl - 6;
    const real T      = x->back();
    for (int i = 0; i < opt->bspline.traj.vel.cols; i++)
      for (int j = 0; j < opt->bspline.traj.vel.rows; j++) {
        const real s = sign(opt->bspline.traj.vel(j, i));
        for (int k = 3; k < (int) opt->bspline.n_ctrl - 3; k++)
          (*grad)[j * n_free + k - 3] += s * opt->bspline.basis_v(k, i) / T;
      }
    grad->back() = -sum / T;
  }
  return sum;
}

TEST_CASE("Custom objective with analytical gradient", "[Optimization][Objective]") {
  auto analytical = make_base_optimization();
  analytical.objective.time_weight = (real) 0.5;
  analytical.objective.add_custom_objective_with_gradient(custom_objective_velocity_sum_with_gradient, 2);

  auto finite_differences = make_base_optimization();
  finite_differences.objective.time_weight = (real) 0.5;
  finite_differences.objective.add_custom_objective(custom_objective_velocity_sum, 2);

  Array x = blast::guess_straight_line(&analytical);
  for (u32 variable = 0; variable + 1 < x.size; variable++)
    x[variable] += (real) 0.1 * std::sin((real) variable);
  x.back() = (real) 1.5;

  Array      grad_analytical(x.size), grad_finite_differences(x.size);
  const real value_analytical         = objective_function(x.size, x.data, grad_analytical.data, &analytical);
  const real value_finite_differences = objective_function(x.size, x.data, grad_finite_differences.data, &finite_differences);
  CHECK(value_analytical == Approx(value_finite_differences).epsilon(test::approx_tol));
  CHECK(objective_function(x.size, x.data, nullptr, &analytical) == Approx(value_analytical).epsilon(test::approx_tol)); // value-only path
#if BLAST_USE_DOUBLES // in float, the forward differences of a sum over every sample are mostly rounding noise
  for (u32 variable = 0; variable < x.size; variable++) {
    INFO("variable " << variable);
    CHECK(grad_analytical[variable] == Approx(grad_finite_differences[variable]).epsilon(1e-3).margin(1e-3));
  }
#endif
}
