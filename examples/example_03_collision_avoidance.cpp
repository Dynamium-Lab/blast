// example_03_collision_avoidance.cpp
//
// This example shows how to:
//   1. Add obstacle geometry to a World
//   2. Enable self-collision and external-collision constraints
//   3. Use optimize() (default: with_segments method) which is better suited
//      for collision-constrained problems (it evaluates constraints at all
//      points within each B-spline segment and reduces them to a single
//      worst-case violation per segment, giving the optimizer fewer but
//      representative constraints over the full trajectory)
//   4. Handle the case where the optimizer does not find a solution and
//      retry with a fresh random guess
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_03_collision_avoidance
// Run:
//   ./build/examples/example_03_collision_avoidance

#define BLAST_USE_NATIVE_SQP  // use the built-in SQP solver

#include <blast>
#include <iostream>
#include "ur5e.hpp"

int main() {
  using namespace blast;

  // -----------------------------------------------------------------------
  // Step 1 — Robot and task (same as example_02).
  // -----------------------------------------------------------------------
  Manipulator ur5e = make_ur5e();
  Matrix      task = make_ur5e_task();

  // -----------------------------------------------------------------------
  // Step 2 — Build a collision World.
  // A World holds static obstacles: boxes, spheres, capsules, and their
  // dynamic (time-varying) equivalents. Here we add a single box obstacle
  // positioned in front of the robot.
  //
  // add_box(center, half-extents, rotation_matrix)
  //   center      : world-frame centre of the box (m)
  //   half-extents: half-width along each axis (m)
  //   rotation    : box orientation as a 3x3 rotation matrix
  // -----------------------------------------------------------------------
  World world;
  world.add_box(
    Vec3{0.4,  0.0, 0.6},           // centre: 40 cm in front, 60 cm high
    Vec3{0.05, 0.3, 0.3},           // half-extents: thin vertical slab
    Mat3{1, 0, 0,  0, 1, 0,  0, 0, 1} // upright, axis-aligned
  );

  // -----------------------------------------------------------------------
  // Step 3 — Create the Optimization problem and attach the world.
  // -----------------------------------------------------------------------
  Optimization opt(ur5e, task);
  opt.bspline = Bspline(16, 110, 5, ur5e.n_joints);
  opt.world   = world;

  // -----------------------------------------------------------------------
  // Step 4 — Enable all constraints including collisions.
  // -----------------------------------------------------------------------
  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.self_collisions     = true;  // avoid self-contact
  opt.constraints.external_collisions = true;  // avoid world obstacles

  opt.success_tolerance = 0.01;

  // -----------------------------------------------------------------------
  // Step 5 — Retry loop.
  // Collision-constrained problems are harder to solve; multiple random
  // starting points may be needed. We try up to max_attempts times.
  // -----------------------------------------------------------------------
  const int max_attempts = 5;
  Result    result(&opt); // initialise with a pointer to opt

  std::cout << "Running collision-aware trajectory optimization...\n";

  for (int attempt = 1; attempt <= max_attempts; ++attempt) {
    std::cout << "Attempt " << attempt << "/" << max_attempts << "... ";

    // Draw a fresh random initial guess for each attempt.
    opt.guess.type = Guess::random;
    opt.guess.n_random_shots = 30;

    // optimize() with with_segments (the default) reduces collision constraints
    // to one worst-case value per B-spline segment, keeping the problem tractable
    // while still covering the full trajectory.
    result = optimize(&opt);

    std::cout << (result.success ? "success" : "failed")
              << " (time: " << result.compute_time << " ms)\n";

    if (result.success)
      break;
  }

  // -----------------------------------------------------------------------
  // Step 6 — Report the final result.
  // -----------------------------------------------------------------------
  std::cout << "\n--- Final result ---\n";
  std::cout << "Success:                  " << (result.success ? "yes" : "no") << "\n";
  std::cout << "Compute time (ms):        " << result.compute_time << "\n";
  std::cout << "Function evaluations:     " << result.num_eval << "\n";
  std::cout << "Max constraint violation: " << result.max_constraint_value << "\n";

  if (!result.x.is_empty()) {
    std::cout << "Trajectory duration (s):  " << result.x.back() << "\n";
  }

  if (!result.success) {
    std::cout << "\nNote: no collision-free trajectory was found in "
              << max_attempts << " attempts. "
              << "Consider adjusting the obstacle, the task, or increasing max_attempts.\n";
  }

  return 0;
}
