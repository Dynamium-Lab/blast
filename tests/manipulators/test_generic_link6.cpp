#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>
#include "test_helper/test_functions.hpp"
#include "test_helper/test_helper.hpp"


// From OG blast :
blast::Matrix jacobian(blast::Link6 manip, const blast::Array& pos) {
  blast::Array s(6);
  blast::Array c(6);
  blast::sincos(pos, s, c);
  blast::Mat3 Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1}; // base -> 0
  blast::Mat3 Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};  // 0    -> 1
  blast::Mat3 Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1}; // 1    -> 2
  blast::Mat3 Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};  // 2    -> 3
  blast::Mat3 Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};  // 3    -> 4
  blast::Mat3 Q6 = {0, s[5], c[5], 0, c[5], -s[5], -1, 0, 0};   // 4    -> 5

  // unit vectors in 1st reference
  blast::Vec3 e[6];
  auto        Q_tmp = Q1;
  e[0]              = Q_tmp * manip.joint_axes[0];
  e[1]              = (Q_tmp *= Q2) * manip.joint_axes[1];
  e[2]              = (Q_tmp *= Q3) * manip.joint_axes[2];
  e[3]              = (Q_tmp *= Q4) * manip.joint_axes[3];
  e[4]              = (Q_tmp *= Q5) * manip.joint_axes[4];
  e[5]              = (Q_tmp *= Q6) * manip.joint_axes[5];

  blast::Vec3 r[6];
  r[5] = manip.joint_offsets[5];
  r[4] = manip.joint_offsets[4] + Q6 * r[5];
  r[3] = manip.joint_offsets[3] + Q5 * r[4];
  r[2] = manip.joint_offsets[2] + Q4 * r[3];
  r[1] = manip.joint_offsets[1] + Q3 * r[2];
  r[0] = manip.joint_offsets[0] + Q2 * r[1];

  Q_tmp = Q1;
  r[0]  = (Q_tmp) *r[0];
  r[1]  = (Q_tmp *= Q2) * r[1];
  r[2]  = (Q_tmp *= Q3) * r[2];
  r[3]  = (Q_tmp *= Q4) * r[3];
  r[4]  = (Q_tmp *= Q5) * r[4];
  r[5]  = (Q_tmp *= Q6) * r[5];


  auto cr0 = cross(e[0], r[0]);
  auto cr1 = cross(e[1], r[1]);
  auto cr2 = cross(e[2], r[2]);
  auto cr3 = cross(e[3], r[3]);
  auto cr4 = cross(e[4], r[4]);
  auto cr5 = cross(e[5], r[5]);

  // jacobian matrix
  blast::Matrix J(6, 6);
  J(3, 0) = e[0].x;
  J(4, 0) = e[0].y;
  J(5, 0) = e[0].z;
  J(3, 1) = e[1].x;
  J(4, 1) = e[1].y;
  J(5, 1) = e[1].z;
  J(3, 2) = e[2].x;
  J(4, 2) = e[2].y;
  J(5, 2) = e[2].z;
  J(3, 3) = e[3].x;
  J(4, 3) = e[3].y;
  J(5, 3) = e[3].z;
  J(3, 4) = e[4].x;
  J(4, 4) = e[4].y;
  J(5, 4) = e[4].z;
  J(3, 5) = e[5].x;
  J(4, 5) = e[5].y;
  J(5, 5) = e[5].z;

  J(0, 0) = cr0.x;
  J(1, 0) = cr0.y;
  J(2, 0) = cr0.z;
  J(0, 1) = cr1.x;
  J(1, 1) = cr1.y;
  J(2, 1) = cr1.z;
  J(0, 2) = cr2.x;
  J(1, 2) = cr2.y;
  J(2, 2) = cr2.z;
  J(0, 3) = cr3.x;
  J(1, 3) = cr3.y;
  J(2, 3) = cr3.z;
  J(0, 4) = cr4.x;
  J(1, 4) = cr4.y;
  J(2, 4) = cr4.z;
  J(0, 5) = cr5.x;
  J(1, 5) = cr5.y;
  J(2, 5) = cr5.z;

  return J;
}

// LINK6 TESTS

TEST_CASE("Link6 jacobian() test", "[Generic]") {
  int  n_tests       = 100;
  auto generic_manip = blast::get_generic_Link6();

  for (int i = 0; i < n_tests; i++) {
    blast::Link6 expected_manip;
    blast::Array test_position(6);
    fill_random(test_position, blast::PI);
    forward_kinematics(generic_manip, test_position);
    blast::Matrix jac_gen = jacobian(generic_manip);
    blast::Matrix jac_exp = jacobian(expected_manip, test_position);

    CHECK(is_close(jac_gen, jac_exp));
  }
}

TEST_CASE("Link6 dynamics() test", "[Generic]") {
  int          n_tests = 100;
  blast::real  epsilon = 1e-6;
  blast::Link6 expected_manip;
  auto         generic_manip = blast::get_generic_Link6();

  for (int i = 0; i < n_tests; i++) {
    blast::Matrix random_task(6, 6);
    fill_random(random_task, blast::PI);
    auto test_trajectory = compute_5order_trajectory(2.0, random_task);

    blast::Matrix result_gen(6, test_trajectory.t.size);
    blast::Matrix result_hc(6, test_trajectory.t.size);
    for (blast::u32 p = 0; p < test_trajectory.t.size; p++) {
      auto tmp_pos = test_trajectory.pos.col(p);
      forward_kinematics(generic_manip, tmp_pos);
      forward_kinematics(expected_manip, tmp_pos);

      auto tmp_vel = test_trajectory.vel.col(p);
      auto tmp_acc = test_trajectory.acc.col(p);
      dynamics(generic_manip, tmp_vel, tmp_acc);
      dynamics(expected_manip, tmp_vel, tmp_acc);

      for (int j = 0; j < 6; j++) {
        result_gen(j, p) = generic_manip._efforts[j];
        result_hc(j, p)  = expected_manip._efforts[j];
      }
    }
    CHECK(is_close(result_gen, result_hc));
  }
}

TEST_CASE("Link6 forward_kinematics() test", "[Generic]") {
  int  n_tests       = 100;
  auto generic_manip = blast::get_generic_Link6();

  for (int i = 0; i < n_tests; i++) {
    blast::Link6 expected_manip;
    blast::Array test_position(6);
    // fill_random(test_position, PI);
    forward_kinematics(generic_manip, test_position);
    forward_kinematics(expected_manip, test_position);

    CHECK(is_close(generic_manip._p_j, expected_manip._p_j));
    CHECK(is_close(generic_manip._Q, expected_manip._Q));
    CHECK(is_close(generic_manip._Q_mult, expected_manip._Q_mult));
  }
}

TEST_CASE("Link6 compute_capsules() test", "[Generic]") {
  int                n_tests = 100;
  real               epsilon = 1e-6;
  blast::Link6       expected_manip;
  GenericManipulator generic_manip = blast::get_generic_Link6();

  for (int i = 0; i < n_tests; i++) {
    blast::Array test_position(6);
    fill_random(test_position, blast::PI);

    forward_kinematics(generic_manip, test_position);
    forward_kinematics(expected_manip, test_position);

    generic_manip.compute_capsules();
    expected_manip.compute_capsules();

    auto generic_internal_dist  = min(generic_manip.internal_collisions());
    auto expected_internal_dist = min(expected_manip.internal_collisions());
    CHECK(is_close(generic_internal_dist, expected_internal_dist, epsilon));
  }
}

TEST_CASE("Link6 compute_constraints() test", "[Generic]") {
  int          n_tests = 100;
  blast::Link6 manip_hc;

  auto                       opt_gen = get_generic_link6_opt();
  Optimization<blast::Link6> opt_hc  = get_hardcoded_link6_opt();

  ncon(&opt_gen);
  ncon(&opt_hc);

  real xlen = opt_gen.bspline.xlen(opt_gen.task);

  blast::Array x_test(xlen);
  blast::Array result_gen(opt_gen.constraints.n_constraints);
  blast::Array result_hc(opt_hc.constraints.n_constraints);
  for (int i = 0; i < n_tests; i++) {
    fill_random(x_test, blast::PI);

    compute_constraints(result_gen.data, x_test, &opt_gen);
    compute_constraints(result_hc.data, x_test, &opt_hc);

    CHECK(is_close(result_gen, result_hc));
  }
}

TEST_CASE("Link6 nlopt_constraints() test", "[Generic]") {
  int          n_tests = 100;
  blast::Link6 manip_hc;

  auto                       opt_gen = get_generic_link6_opt();
  Optimization<blast::Link6> opt_hc  = get_hardcoded_link6_opt();

  ncon(&opt_gen);
  ncon(&opt_hc);

  real xlen = opt_gen.bspline.xlen(opt_gen.task);

  blast::Array x_test(xlen);
  blast::Array constraints_value_gen(opt_gen.constraints.n_constraints);
  blast::Array constraints_value_hc(opt_hc.constraints.n_constraints);
  blast::Array grad_gen(opt_gen.constraints.n_constraints * xlen);
  blast::Array grad_hc(opt_hc.constraints.n_constraints * xlen);
  for (int i = 0; i < n_tests; i++) {
    fill_random(x_test, blast::PI);

    nlopt_constraints<GenericManipulator>((unsigned int) opt_gen.constraints.n_constraints, constraints_value_gen.data, (unsigned int) xlen, x_test.data, grad_gen.data, (void*) (&opt_gen));
    nlopt_constraints<blast::Link6>(opt_hc.constraints.n_constraints, constraints_value_hc.data, xlen, x_test.data, grad_hc.data, &opt_hc);

    CHECK(is_close(constraints_value_gen, constraints_value_hc));
    CHECK(is_close(grad_gen, grad_hc));
  }
}

TEST_CASE("Link6 compute_objective() test", "[objective]") {
  int          n_tests = 100;
  blast::Link6 manip_hc;

  auto                       opt_gen = get_generic_link6_opt();
  Optimization<blast::Link6> opt_hc  = get_hardcoded_link6_opt();

  for (int i = 0; i < n_tests; i++) {
    blast::Array x(37, 5 * (get_random() + 1));

    double obj_fun_gen = compute_objective(x, &opt_gen);
    double obj_fun_hc  = compute_objective(x, &opt_hc);

    CHECK(is_close(obj_fun_gen, obj_fun_hc));
  }
}

TEST_CASE("Link6 objective_function() test", "[objective]") {
  int          n_tests = 100;
  blast::Link6 manip_hc;

  auto                       opt_gen = get_generic_link6_opt();
  Optimization<blast::Link6> opt_hc  = get_hardcoded_link6_opt();

  real xlen = opt_gen.bspline.xlen(opt_gen.task);

  blast::Array x_test(xlen);
  blast::Array grad_gen(opt_gen.constraints.n_constraints * xlen);
  blast::Array grad_hc(opt_hc.constraints.n_constraints * xlen);
  for (int i = 0; i < n_tests; i++) {
    fill_random(x_test, blast::PI);

    double obj_fun_gen = objective_function<GenericManipulator>(xlen, x_test.data, grad_gen.data, &opt_gen);
    double obj_fun_hc  = objective_function<blast::Link6>(xlen, x_test.data, grad_hc.data, &opt_hc);

    CHECK(is_close(obj_fun_gen, obj_fun_hc));
    CHECK(is_close(grad_gen, grad_hc));
  }
}

TEST_CASE("Link6 optimize() test", "[Optimization]") {
  blast::Link6 manip_hc;
  int          n_tests = 1;

  auto                       opt_gen = get_generic_link6_opt();
  Optimization<blast::Link6> opt_hc  = get_hardcoded_link6_opt();

  opt_gen.constraints.self_collisions     = false;
  opt_gen.constraints.external_collisions = false;
  opt_hc.constraints.self_collisions      = false;
  opt_hc.constraints.external_collisions  = false;

  World world;
  add_static_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world);  // vertical plate (no coll)
  add_static_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // table

  opt_gen.set_world(world);
  opt_hc.set_world(world);

  std::vector<blast::Matrix> tasks = get_Link6_demo1_tasks();

  for (int j = 0; j < n_tests; j++) {
    for (int i = 0; i < tasks.size(); i++) {
      opt_gen.set_task(tasks[i]);
      opt_hc.set_task(tasks[i]);
      auto result_gen = optimize(&opt_gen, OptimizationMethod::baseline);

      opt_hc.guess.type = GuessType::custom;
      opt_hc.guess.initial_x   = result_gen.x0;
      auto result_hc    = optimize(&opt_hc, OptimizationMethod::baseline);

      CHECK(is_close(result_gen.x0, result_hc.x0));
      CHECK(is_close(result_gen.x, result_hc.x));
      if (!is_close(result_gen.x, result_hc.x)) {
        std::cout << "Task " << i << " was not the same x" << std::endl;
        std::cout << "Exit criteria = " << result_gen.nlopt_exit_criteria << std::endl;
        std::cout << "Exit criteria = " << result_hc.nlopt_exit_criteria << std::endl;
        std::cout << "num eval = " << result_gen.num_eval << std::endl;
        std::cout << "num eval = " << result_hc.num_eval << std::endl;
        print(result_gen.x - result_hc.x);
      }
      CHECK(result_gen.success == result_hc.success);
      CHECK(result_gen.success_false == result_hc.success_false);
    }
  }
}
