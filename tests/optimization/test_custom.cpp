#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;

// Returns a representative stop-to-stop task for the UR5e.
inline blast::Task make_UR5e_task() {
  blast::Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  blast::Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  return blast::Task::stop_to_stop(start, end);
}

Optimization make_base_optimization() {
  Manipulator manip = make_UR5e();

  int     n_pts  = 100;
  int     p      = 5;
  int     n_ctrl = 12;
  Bspline bspline(n_ctrl, n_pts, p, manip.n_joints);

  auto task = make_UR5e_task();

  return Optimization(manip, Task::stop_to_stop(Array(6, 0), Array(6, 1)), bspline);
}

TEST_CASE("time objective", "[Optimization][Objective]") {
  auto opt = make_base_optimization();

  opt.objective.time_weight = 1.0;
  Array x(opt.bspline.x_len(opt.task), 1);

  Array grad(x.size);

  SECTION("basic time objective") {
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 1.0);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == 0.0);
    }
    CHECK(grad.back() == Approx(1.0));

    // Changing objective value
    x.back() = 2;

    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 2.0);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == 0.0);
    }
    CHECK(grad.back() == Approx(1.0));

    // Changing time weight
    opt.objective.time_weight = 2;

    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 4.0);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == 0.0);
    }
    CHECK(grad.back() == Approx(2.0));

    // Zero-ing time weight
    opt.objective.time_weight = 0;

    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 0.0);
    for (int i = 0; i < grad.size; i++) {
      CHECK(grad[i] == 0.0);
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

  opt.objective.time_weight = 0.0;
  opt.objective.add_custom_objective(custom_objective_velocity_sum, 1.0);

  Array x(opt.bspline.x_len(opt.task), 1);

  Array grad(x.size);

  SECTION("velocity objective") {
    real  original_objective_value = objective_function(x.size, x.data, grad.data, &opt);
    Array original_gradient_value  = grad;

    // Changing time value
    x.back() = 2;
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 0.5 * original_objective_value);
    for (int i = 0; i < grad.size - 1; i++) {
      CHECK(grad[i] == 0.5 * original_gradient_value[i]);
    }
    CHECK(grad.back() == Approx(0.25 * original_gradient_value.back()));

    x.back() = 1;

    // Changing time weight
    opt.objective.k_extra_objectives[0] = 2;
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 2 * original_objective_value);
    for (int i = 0; i < grad.size; i++) {
      CHECK(grad[i] == 2 * original_gradient_value[i]);
    }

    // Zero-ing custom weight
    opt.objective.k_extra_objectives[0] = 0;
    CHECK(objective_function(x.size, x.data, grad.data, &opt) == 0.0);
    for (int i = 0; i < grad.size; i++) {
      CHECK(grad[i] == 0.0);
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

  opt.objective.time_weight = 0.0;

  auto result1 = optimize(&opt);

  opt.guess.type      = opt.guess.custom;
  opt.guess.initial_x = result1.x;

  opt.objective.add_custom_objective(custom_objective_velocity_sum, 1.0);
  opt.objective.time_weight = 0;

  opt.constraints.add_constraint(custom_constraint_time, 1);

  real custom_value = result1.x.back() + 0.0;
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
