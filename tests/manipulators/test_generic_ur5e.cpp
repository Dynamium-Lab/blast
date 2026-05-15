#define CATCH_CONFIG_MAIN

#include "blast_rush.h"
#include "catch2/catch.hpp"
#include "manipulator/UR5e.hpp"
#include "test_helper/test_functions.hpp"
#include "test_helper/test_helper.hpp"

using namespace blast;

TEST_CASE("UR5e forward_kinematics() test", "[Generic]") {

  GenericManipulator generic_manip = make_UR5e();

  //  @ home position
  Vec3  expected_tool_position = {0.0, -0.233, 1.08};
  Array test_position          = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};
  forward_kinematics(generic_manip, test_position);
  CHECK(is_close(generic_manip._p_j[generic_manip.joints], expected_tool_position, 1e-2));

  //  @ 2nd position
  Vec3  expected_tool_position_2 = {-0.1179, -0.436, 0.152};
  Array test_position_2          = {-1.54, -1.83, -2.28, -0.59, 1.60, 0.023};
  forward_kinematics(generic_manip, test_position_2);
  CHECK(is_close(generic_manip._p_j[generic_manip.joints], expected_tool_position_2, 1e-2));
}

TEST_CASE("UR5e optimize() trajectory with obstacles", "[Generic]") {
  using namespace blast;

  GenericManipulator generic_ur5e = make_UR5e();

  // --- Define Tasks ---
  std::vector<Array> positions(2);
  positions[0] = {2.0473, -0.3930, 1.2048, -2.3880, -1.6118, 1.8403}; // First touch point
  positions[1] = {0.7529, -0.1918, 0.8131, -2.2195, -1.6113, 1.8410}; // Second touch point

  Array top_screen_pos = {1.4610, -0.8622, 0.6568, -1.5820, -1.6106, 1.8282};

  std::vector<Matrix> tasks;
  for (u32 i = 0; i < positions.size() - 1; i++) {
    Matrix tmp(6, 6);

    for (u32 j = 0; j < generic_ur5e.joints; j++) {
      tmp(j, 0) = positions[i][j];
      tmp(j, 3) = positions[i + 1][j];
    }
    tasks.push_back(tmp);
  }

  // --- Define World ---
  World world;
  add_static_box({0.038, -0.676, 0.146}, {0.025, 0.35, 0.2}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // vertical plate
  add_static_box({0.038, -0.676, -0.404}, {0.6, 0.38, 0.35}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // table

  // --- Define Bspline ---
  Bspline bspline(16, 200, 5, generic_ur5e.joints);

  // --- Define Guess ---
  Guess guess; // using default guess
  guess.n_random_shots = 1000;

  // -- Define Constraints ---
  Constraints<GenericManipulator> constraints;
  constraints.position                = true;
  constraints.velocity                = true;
  constraints.acceleration            = false;
  constraints.torque                  = true;
  constraints.tool_speed              = true;
  constraints.self_collisions         = true;
  constraints.external_collisions     = true;
  constraints.n_collision_constraints = 100;
  constraints.n_collision_skip        = 2;

  // --- Define Objective ---
  Objective<GenericManipulator> objective;
  objective.time_weight = 1;

  real start_time = 0.0;

  for (u32 t = 0; t < tasks.size(); t++) {
    std::cout << "Task " << t + 1 << " of " << tasks.size() << std::endl;
    std::cout << "Trajectory starts at time t = " << start_time << " s" << std::endl;
    Optimization<GenericManipulator> opt(generic_ur5e, tasks[t]);
    opt.set_world(world);
    opt.set_bspline(bspline);
    opt.set_guess(guess);
    opt.set_constraints(constraints);
    opt.set_objective(objective);

    auto results = optimize(&opt, OptimizationMethod::baseline, 2);
    start_time += results.x[results.x.size - 1];
    std::cout << "Optimization done in " << results.compute_time / 1000.0 << std::endl;
    CHECK(results.success == true);
    CHECK(results.success_false == false);
  }
}
