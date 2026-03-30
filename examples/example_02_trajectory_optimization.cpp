// example_02_trajectory_optimization.cpp
//
// This example shows how to:
//   1. Define a point-to-point task (start and goal joint positions)
//   2. Configure a B-spline trajectory representation
//   3. Run trajectory optimization with kinematic constraints
//   4. Inspect the result
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_02_trajectory_optimization
// Run:
//   ./build/examples/example_02_trajectory_optimization

#define BLAST_USE_NATIVE_SQP  // use the built-in SQP solver (no extra install required)

#include <blast>
#include <iostream>
#include "ur5e.hpp"

int main() {
  using namespace blast;

  // -----------------------------------------------------------------------
  // Step 1 — Build the robot model and define the task.
  // A task is a (n_joints x 6) matrix encoding boundary conditions:
  //   column 0 → start joint positions  (rad)
  //   column 3 → goal  joint positions  (rad)
  //   columns 1,2,4,5 → boundary velocities/accelerations (zero = stop-to-stop)
  // -----------------------------------------------------------------------
  Manipulator ur5e = make_ur5e();
  Matrix      task = make_ur5e_task();

  // -----------------------------------------------------------------------
  // Step 2 — Create the Optimization problem.
  // Optimization bundles the robot, task, B-spline, constraints, and solver
  // settings into a single object.
  // -----------------------------------------------------------------------
  Optimization opt(ur5e, task);

  // -----------------------------------------------------------------------
  // Step 3 — Configure the B-spline trajectory representation.
  // Bspline(n_control_points, n_eval_points, degree, n_joints)
  //   n_control_points : resolution of the optimizer's decision variable
  //   n_eval_points    : number of samples used when evaluating constraints
  //   degree           : smoothness (5 → C4-continuous)
  // -----------------------------------------------------------------------
  opt.bspline = Bspline(16, 110, 5, ur5e.n_joints);

  // -----------------------------------------------------------------------
  // Step 4 — Select which constraints to enforce.
  // Start with the three kinematic constraints; collision is added in
  // example_03.
  // -----------------------------------------------------------------------
  opt.constraints.position     = true;  // enforce start and goal positions
  opt.constraints.velocity     = true;  // stay within joint velocity limits
  opt.constraints.acceleration = true;  // stay within joint acceleration limits

  // Allow up to a small relative violation before declaring failure
  // (0.01 = 1 % of the constraint bound).
  opt.success_tolerance = 0.01;

  // -----------------------------------------------------------------------
  // Step 5 — Configure the initial guess strategy.
  // Guess::random tries n_shot random starting points and picks the best
  // one as the initial solution for the optimizer.
  // -----------------------------------------------------------------------
  opt.guess.type   = Guess::random;
  opt.guess.n_random_shots = 50;

  // -----------------------------------------------------------------------
  // Step 6 — Run the optimizer.
  // -----------------------------------------------------------------------
  std::cout << "Running trajectory optimization...\n";
  Result result = optimize(&opt, OptimizationMethod::baseline);

  // -----------------------------------------------------------------------
  // Step 7 — Inspect the result.
  // -----------------------------------------------------------------------
  std::cout << "\n--- Result ---\n";
  std::cout << "Success:             " << (result.success ? "yes" : "no") << "\n";
  std::cout << "Compute time (ms):   " << result.compute_time << "\n";
  std::cout << "Function evaluations:" << result.num_eval << "\n";

  if (!result.x.is_empty()) {
    std::cout << "Trajectory duration (s): " << result.x.back() << "\n";
  }

  std::cout << "Max constraint violation: " << result.max_constraint_value << "\n";

  return 0;
}
