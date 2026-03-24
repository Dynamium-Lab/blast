#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>
#include "test_helper/test_functions.hpp"
#include "test_helper/test_helper.hpp"

using namespace blast;

// GEN3 TESTS
TEST_CASE("Gen3 forward_kinematics() test", "[Generic]") {
  int  n_tests = 100;
  real epsilon = 1e-6;
  Gen3 expected_manip;
  auto generic_manip    = get_generic_gen3();
  expected_manip.base_position = generic_manip.base_position;
  expected_manip.base_rotation = generic_manip.base_rotation;

  // todo: remove?
  // setup_manip(&expected_manip);
  // setup_manip(&generic_manip);

  for (int i = 0; i < n_tests; i++) {
    Array test_position(7);
    fill_random(test_position, PI);
    forward_kinematics(generic_manip, test_position);
    forward_kinematics(expected_manip, test_position);

    CHECK(is_close(generic_manip._p_j, expected_manip._p_j, epsilon));
    CHECK(is_close(generic_manip._rotations, expected_manip._rotations, epsilon));
    CHECK(is_close(generic_manip._rotations_mult, expected_manip._rotations_mult, epsilon));
  }
}

TEST_CASE("Gen3 dynamics() test", "[Generic]") {
  int  n_tests = 100;
  real epsilon = 1e-6;
  Gen3 expected_manip;
  auto generic_manip    = get_generic_gen3();
  expected_manip.base_position = generic_manip.base_position;
  expected_manip.base_rotation = generic_manip.base_rotation;

  // todo: remove?
  // setup_manip(&expected_manip);
  // setup_manip(&generic_manip);

  for (int i = 0; i < n_tests; i++) {
    Matrix random_task(7, 6);
    fill_random(random_task, PI);
    Trajectory test_trajectory = compute_5order_trajectory(2.0, random_task);

    Matrix result_gen(7, test_trajectory.t.size);
    Matrix result_hc(7, test_trajectory.t.size);
    for (u32 p = 0; p < test_trajectory.t.size; p++) {
      auto tmp_pos = test_trajectory.pos.col(p);
      forward_kinematics(generic_manip, tmp_pos);
      forward_kinematics(expected_manip, tmp_pos);

      auto tmp_vel = test_trajectory.vel.col(p);
      auto tmp_acc = test_trajectory.acc.col(p);
      dynamics(generic_manip, tmp_vel, tmp_acc);
      dynamics(expected_manip, tmp_vel, tmp_acc);

      for (int j = 0; j < 7; j++) {
        result_gen(j, p) = generic_manip._efforts[j];
        result_hc(j, p)  = expected_manip._efforts[j];
      }
    }
    CHECK(is_close(result_gen, result_hc));
  }
}

// todo: collisions not working properly
TEST_CASE("Gen3 compute_capsules() test", "[Generic]") {
  int  n_tests = 100;
  Gen3 expected_manip;
  auto generic_manip    = get_generic_gen3();
  expected_manip.base_position = generic_manip.base_position;
  expected_manip.base_rotation = generic_manip.base_rotation;

  for (int i = 0; i < n_tests; i++) {
    Array test_position(7);
    fill_random(test_position, PI);

    forward_kinematics(generic_manip, test_position);
    forward_kinematics(expected_manip, test_position);

    generic_manip.compute_capsules();
    expected_manip.compute_capsules();
    auto generic_internal       = generic_manip.get_internal_collisions();
    auto expected_internal      = expected_manip.internal_collisions();
    auto generic_internal_dist  = min(generic_internal);
    auto expected_internal_dist = min(expected_internal);
    // CHECK(is_close(generic_internal_dist, expected_internal_dist, epsilon));

    // CHECK(is_close(generic_manip._capsule_list, expected_manip._capsule_list));

    // auto gen_self_collision_distances = generic_manip.internal_collisions();
    // auto expected_self_collision_distances = expected_manip.internal_collisions();

    // CHECK(is_close(gen_self_collision_distances, expected_self_collision_distances));
  }
}

TEST_CASE("Gen3 compute_constraints() test", "[Generic]") {
  int  n_tests = 100;
  Gen3 manip_hc;

  auto         opt_gen = get_generic_gen3_opt();
  Optimization opt_hc  = get_gen3_gen3_opt();
  opt_hc.manip.base_position  = opt_gen.manip.base_position;
  opt_hc.manip.base_rotation  = opt_gen.manip.base_rotation;

  opt_gen.constraints.self_collisions     = false;
  opt_gen.constraints.external_collisions = false;
  opt_hc.constraints.self_collisions      = false;
  opt_hc.constraints.external_collisions  = false;

  ncon(&opt_gen);
  ncon(&opt_hc);

  real x_len = opt_gen.bspline.x_len(opt_gen.task);

  Array x_test(x_len);
  Array result_gen(opt_gen.constraints.n_constraints);
  Array result_hc(opt_hc.constraints.n_constraints);
  for (int i = 0; i < n_tests; i++) {
    fill_random(x_test, PI);

    compute_constraints(result_gen.data, x_test, &opt_gen);
    compute_constraints(result_hc.data, x_test, &opt_hc);

    CHECK(is_close(result_gen, result_hc));
  }
}

TEST_CASE("Gen3 nlopt_constraints() test", "[Generic]") {
  int  n_tests = 100;
  Gen3 manip_hc;

  auto         opt_gen = get_generic_gen3_opt();
  Optimization opt_hc  = get_gen3_gen3_opt();
  opt_hc.manip.base_position  = opt_gen.manip.base_position;
  opt_hc.manip.base_rotation  = opt_gen.manip.base_rotation;

  opt_gen.constraints.self_collisions     = false;
  opt_gen.constraints.external_collisions = false;
  opt_hc.constraints.self_collisions      = false;
  opt_hc.constraints.external_collisions  = false;

  ncon(&opt_gen);
  ncon(&opt_hc);

  real x_len = opt_gen.bspline.x_len(opt_gen.task);

  Array x_test(x_len);
  Array constraints_value_gen(opt_gen.constraints.n_constraints);
  Array constraints_value_hc(opt_hc.constraints.n_constraints);
  Array grad_gen(opt_gen.constraints.n_constraints * x_len);
  Array grad_hc(opt_hc.constraints.n_constraints * x_len);
  for (int i = 0; i < n_tests; i++) {
    fill_random(x_test, PI);

    nlopt_constraints(opt_gen.constraints.n_constraints, constraints_value_gen.data, (unsigned int) x_len, x_test.data, grad_gen.data, (void*) (&opt_gen));
    nlopt_constraints(opt_hc.constraints.n_constraints, constraints_value_hc.data, x_len, x_test.data, grad_hc.data, &opt_hc);

    CHECK(is_close(constraints_value_gen, constraints_value_hc));
    CHECK(is_close(grad_gen, grad_hc));
  }
}

TEST_CASE("Gen3 compute_objective() test", "[objective]") {
  int  n_tests = 100;
  Gen3 manip_hc;

  auto         opt_gen = get_generic_gen3_opt();
  Optimization opt_hc  = get_gen3_gen3_opt();
  opt_hc.manip.base_position  = opt_gen.manip.base_position;
  opt_hc.manip.base_rotation  = opt_gen.manip.base_rotation;

  for (int i = 0; i < n_tests; i++) {
    Array x(37, 5 * (get_random() + 1));

    double obj_fun_gen = compute_objective(x, &opt_gen);
    double obj_fun_hc  = compute_objective(x, &opt_hc);

    CHECK(is_close(obj_fun_gen, obj_fun_hc));
  }
}

TEST_CASE("Gen3 objective_function() test", "[objective]") {
  int  n_tests = 100;
  Gen3 manip_hc;

  auto         opt_gen = get_generic_gen3_opt();
  Optimization opt_hc  = get_gen3_gen3_opt();
  opt_hc.manip.base_position  = opt_gen.manip.base_position;
  opt_hc.manip.base_rotation  = opt_gen.manip.base_rotation;

  real x_len = opt_gen.bspline.x_len(opt_gen.task);

  Array x_test(x_len);
  Array grad_gen(opt_gen.constraints.n_constraints * x_len);
  Array grad_hc(opt_hc.constraints.n_constraints * x_len);
  for (int i = 0; i < n_tests; i++) {
    fill_random(x_test, PI);

    double obj_fun_gen = objective_function(x_len, x_test.data, grad_gen.data, &opt_gen);
    double obj_fun_hc  = objective_function(x_len, x_test.data, grad_hc.data, &opt_hc);

    CHECK(is_close(obj_fun_gen, obj_fun_hc));
    CHECK(is_close(grad_gen, grad_hc));
  }
}

TEST_CASE("Gen3 optimize() test", "[Optimization][objective]") {
  int n_tests = 10;

  for (int i = 0; i < n_tests; i++) {
    auto         opt_gen = get_generic_gen3_opt();
    Optimization opt_hc  = get_gen3_gen3_opt();

    opt_gen.constraints.self_collisions     = false;
    opt_gen.constraints.external_collisions = false;
    opt_hc.constraints.self_collisions      = false;
    opt_hc.constraints.external_collisions  = false;

    World world;
    add_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // vertical plate (no coll)
    // add_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // table

    opt_gen.set_world(world);
    opt_hc.set_world(world);

    auto result_gen = optimize(&opt_gen);

    opt_hc.guess.type = Guess::custom;
    opt_hc.guess.initial_x   = (result_gen.x0);
    auto result_hc    = optimize(&opt_hc);

    CHECK(is_close(result_gen.x0, result_hc.x0));
    CHECK(is_close(result_gen.x, result_hc.x));
    CHECK(result_gen.success == result_hc.success);
    CHECK(result_gen.success_false == result_hc.success_false);
  }
}
