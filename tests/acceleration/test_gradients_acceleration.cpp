#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include "blast"
#include "catch2/catch.hpp"
#include "test_helper/test_functions.hpp"
#include "test_helper/test_helper.hpp"

using namespace blast;

std::vector<blast::Guess> read_array_from_csv(std::string file_path, blast::u32 xlen) {
  std::ifstream file(file_path);
  if (!file.is_open()) {
    std::cerr << "Error opening file" << std::endl;
    return {};
  }
  std::vector<blast::Guess> initial_guess_list;
  std::string               line;

  // Read each row from the CSV file
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string       token;
    blast::Array      arr(xlen);

    blast::u32 idx = 0;
    // Split the row using comma as the delimiter
    while (std::getline(ss, token, ',')) {
      try {
        // Convert token to double and add it to the Array's values
        arr[idx] = std::stod(token);
      } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid number found: " << token << std::endl;
      }
      idx++;
    }
    Guess guess_array(arr);
    initial_guess_list.push_back(guess_array);
  }
  file.close();
  return initial_guess_list;
}

TEST_CASE("test gradient upper and lower bounds nctrl = 16 npoints = 200", "[Acceleration]") {

  // find lb and ub for x, given p = 5
  u32 nctrl   = 16;
  u32 npoints = 200;
  u32 joints  = 6;

  Bspline bspline(nctrl, npoints, 5, joints);

  // print_to_csv(bspline.basis_p, "../../../tests/acceleration/basis_p.csv");
  // print_to_csv(bspline.basis_v, "../../../tests/acceleration/basis_v.csv");
  // print_to_csv(bspline.basis_a, "../../../tests/acceleration/basis_a.csv");

  for (u32 i = 0; i < nctrl; i++) {
    // we do not test the first and last points
    auto lb = bspline.lower_bounds[i];
    auto ub = bspline.upper_bounds[i];
    for (u32 k = lb; k <= ub; k++) {
      if (bspline.basis_p(i, k) == 0.0)
        std::cout << "i: " << i << " k: " << k << " basis_p: " << bspline.basis_p(i, k) << std::endl;
      CHECK(bspline.basis_p(i, k) != 0.0);
      CHECK(bspline.basis_v(i, k) != 0.0);
      CHECK(bspline.basis_a(i, k) != 0.0);
    }
    for (u32 k = 1; k < lb; k++) {
      if (bspline.basis_a(i, k) != 0.0)
        std::cout << "i: " << i << " k: " << k << " basis_a: " << bspline.basis_a(i, k) << std::endl;
      CHECK(is_close(bspline.basis_p(i, k), 0.0));
      CHECK(is_close(bspline.basis_v(i, k), 0.0));
      CHECK(is_close(bspline.basis_a(i, k), 0.0));
    }
    for (u32 k = ub + 1; k < npoints - 1; k++) {
      if (bspline.basis_p(i, k) != 0.0)
        std::cout << "i: " << i << " k: " << k << " basis_p: " << bspline.basis_p(i, k) << std::endl;
      CHECK(is_close(bspline.basis_p(i, k), 0.0));
      CHECK(is_close(bspline.basis_v(i, k), 0.0));
      CHECK(is_close(bspline.basis_a(i, k), 0.0));
    }
  }
}

TEST_CASE("test gradient upper and lower bounds random nctrl & npoints", "[Acceleration]") {
  u32 n_tests = 1000;

  for (u32 t = 0; t < n_tests; t++) {

    // find lb and ub for x, given p = 5
    u32 nctrl   = 10 + (int) std::abs(10 * get_random());
    u32 npoints = 200 + (int) std::abs(300 * get_random());
    u32 joints  = 6;

    Bspline bspline(nctrl, npoints, 5, joints);

    for (u32 i = 0; i < nctrl; i++) {
      auto lb = bspline.lower_bounds[i];
      auto ub = bspline.upper_bounds[i];
      for (u32 k = lb; k <= ub; k++) {
        // todo: make sure value is max/min
        bool b_p = bspline.basis_p(i, k) != 0.0 ? true : bspline.basis_p(i - 1, k) != 0.0 && bspline.basis_p(i + 1, k) != 0.0 ? true
                                                                                                                              : false; // remove outliers
        bool b_v = bspline.basis_v(i, k) != 0.0 ? true : bspline.basis_v(i - 1, k) != 0.0 && bspline.basis_v(i + 1, k) != 0.0 ? true
                                                                                                                              : false; // remove outliers
        bool b_a = bspline.basis_a(i, k) != 0.0 ? true : bspline.basis_a(i - 1, k) != 0.0 && bspline.basis_a(i + 1, k) != 0.0 ? true
                                                                                                                              : false; // remove outliers
        CHECK(b_p);
        CHECK(b_v);
        CHECK(b_a);
      }
      for (u32 k = 1; k < lb; k++) {
        if (bspline.basis_p(i, k) != 0.0 || bspline.basis_v(i, k) != 0.0 || bspline.basis_a(i, k) != 0.0)
          std::cout << "i: " << i << " k: " << k << " basis_p: " << bspline.basis_p(i, k) << " basis_v: " << bspline.basis_v(i, k) << " basis_a: " << bspline.basis_a(i, k) << std::endl;
        CHECK(is_close(bspline.basis_p(i, k), 0.0));
        CHECK(is_close(bspline.basis_v(i, k), 0.0));
        CHECK(is_close(bspline.basis_a(i, k), 0.0));
      }
      for (u32 k = ub + 1; k < npoints - 1; k++) {
        if (bspline.basis_p(i, k) != 0.0 || bspline.basis_v(i, k) != 0.0 || bspline.basis_a(i, k) != 0.0)
          std::cout << "i: " << i << " k: " << k << " basis_p: " << bspline.basis_p(i, k) << " basis_v: " << bspline.basis_v(i, k) << " basis_a: " << bspline.basis_a(i, k) << std::endl;
        CHECK(is_close(bspline.basis_p(i, k), 0.0));
        CHECK(is_close(bspline.basis_v(i, k), 0.0));
        CHECK(is_close(bspline.basis_a(i, k), 0.0));
      }
    }
  }
}

TEST_CASE("test gradient points with upper & lower bounds", "[Acceleration]") {
  using namespace blast;

  std::string input_file_initial_guess = "../../../examples/acceleration/files/link6_initial_guess.csv";

  // blast::begin_profile();
  // --- Define Manipulator ---
  auto link6 = get_generic_Link6();

  // --- Define Tasks ---
  std::vector<Array> positions(14);
  positions[0]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};                   // zero
  positions[1]  = deg2rad({-40.4, -26.9, 83.6, 1.5, 20.0, -42.2});  // w1 + 10
  positions[2]  = deg2rad({-40.2, -34.9, 96.4, 1.4, 41.7, -41.5});  // w1
  positions[3]  = deg2rad({-40.4, -26.9, 83.6, 1.5, 20.0, -42.2});  // w1 + 10
  positions[4]  = deg2rad({51.9, -13.6, 107.9, 3.6, 33.1, 51.2});   // wb1
  positions[5]  = deg2rad({-33.9, -35.6, 70.9, -0.7, 16.9, -31.8}); // w3 + 10
  positions[6]  = deg2rad({-34.8, -42.5, 81.6, 0.0, 33.3, -35.1});  // w3
  positions[7]  = deg2rad({-33.9, -35.6, 70.9, -0.7, 16.9, -31.8}); // w3 + 10
  positions[8]  = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});    // wb3
  positions[9]  = deg2rad({-20.0, -19.3, 95.9, -3.0, 23.5, -17.6}); // w4 + 10
  positions[10] = deg2rad({-20.5, -30.9, 105.2, 1.3, 43.3, -22.1}); // w4
  positions[11] = deg2rad({-20.0, -19.3, 95.9, -3.0, 23.5, -17.6}); // w4 + 10
  positions[12] = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});    // wb4
  positions[13] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};                   // zero

  std::vector<Matrix> tasks;
  for (u32 i = 0; i < positions.size() - 1; i++) {
    Matrix tmp(6, 6);

    for (u32 j = 0; j < link6.joints; j++) {
      tmp(j, 0) = positions[i][j];
      tmp(j, 3) = positions[i + 1][j];
    }
    tasks.push_back(tmp);
  }

  // --- Define Bspline ---
  Bspline bspline(16, 200, 5, 6);

  // --- Define Guess ---
  std::vector<Guess> guess_predefined;


  auto xlen        = bspline.xlen(tasks[0]);
  guess_predefined = read_array_from_csv(input_file_initial_guess, xlen);
  Assert(!guess_predefined.empty());

  // -- Define Constraints ---
  Constraints<GenericManipulator> constraints;
  constraints.show_info = true;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.torque       = true;

  constraints.tool_speed              = false;
  constraints.self_collisions         = false;
  constraints.external_collisions     = false;
  constraints.n_collision_constraints = 100;
  constraints.n_collision_skip        = 2;

  // --- Define Objective ---
  Objective<GenericManipulator> objective;
  objective.time_weight = 1;

  real start_time = 0.0;
  for (u32 t = 0; t < tasks.size(); t++) {
    Optimization<GenericManipulator> opt(link6, tasks[t]);

    opt.set_bspline(bspline);
    opt.set_guess(guess_predefined[t]);
    opt.set_constraints(constraints);
    opt.set_objective(objective);
    auto results = optimize(&opt, OptimizationMethod::baseline);

    Optimization<GenericManipulator> opt_acc(link6, tasks[t]);
    opt_acc.set_bspline(bspline);
    opt_acc.set_guess(guess_predefined[t]);
    opt_acc.set_constraints(constraints);
    opt_acc.set_objective(objective);

    auto result_acc = optimize_acc(&opt_acc);

    // CHECK(results.opt.constraints.grad_list.size() == result_acc.opt.constraints.grad_list.size());
    // CHECK(results.opt.constraints.constr_list.size() == result_acc.opt.constraints.constr_list.size());
    // CHECK(results.opt.constraints.x_list.size() == result_acc.opt.constraints.x_list.size());
    CHECK(is_close(results.x, result_acc.x));
    // for (u32 i = 0; i <results.opt.constraints.grad_list.size(); i++) {
    //     // u32 i = 10;
    //     for (u32 j = 0; j < results.opt.constraints.grad_list[i].rows; j++) {
    //         for (u32 k = 1; k < results.opt.constraints.grad_list[i].cols - 1; k++) { // skip first and last point
    //             CHECK(is_close(results.opt.constraints.grad_list[i](j, k), result_acc.opt.constraints.grad_list[i](j, k)));
    //             CHECK(is_close(results.opt.constraints.constr_list[i][k], result_acc.opt.constraints.constr_list[i][k]));
    //         }
    //     }
    // }
  }
}

TEST_CASE("test gradient acceleration accuracy", "[Acceleration]") {
  using namespace blast;

  std::string input_file_initial_guess = "../../../examples/acceleration/files/link6_initial_guess.csv";

  // blast::begin_profile();
  // --- Define Manipulator ---
  auto link6 = get_generic_Link6();

  // --- Define Tasks ---
  std::vector<Array> positions(14);
  positions[0]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};                   // zero
  positions[1]  = deg2rad({-40.4, -26.9, 83.6, 1.5, 20.0, -42.2});  // w1 + 10
  positions[2]  = deg2rad({-40.2, -34.9, 96.4, 1.4, 41.7, -41.5});  // w1
  positions[3]  = deg2rad({-40.4, -26.9, 83.6, 1.5, 20.0, -42.2});  // w1 + 10
  positions[4]  = deg2rad({51.9, -13.6, 107.9, 3.6, 33.1, 51.2});   // wb1
  positions[5]  = deg2rad({-33.9, -35.6, 70.9, -0.7, 16.9, -31.8}); // w3 + 10
  positions[6]  = deg2rad({-34.8, -42.5, 81.6, 0.0, 33.3, -35.1});  // w3
  positions[7]  = deg2rad({-33.9, -35.6, 70.9, -0.7, 16.9, -31.8}); // w3 + 10
  positions[8]  = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});    // wb3
  positions[9]  = deg2rad({-20.0, -19.3, 95.9, -3.0, 23.5, -17.6}); // w4 + 10
  positions[10] = deg2rad({-20.5, -30.9, 105.2, 1.3, 43.3, -22.1}); // w4
  positions[11] = deg2rad({-20.0, -19.3, 95.9, -3.0, 23.5, -17.6}); // w4 + 10
  positions[12] = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});    // wb4
  positions[13] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};                   // zero

  std::vector<Matrix> tasks;
  for (u32 i = 0; i < positions.size() - 1; i++) {
    Matrix tmp(6, 6);

    for (u32 j = 0; j < link6.joints; j++) {
      tmp(j, 0) = positions[i][j];
      tmp(j, 3) = positions[i + 1][j];
    }
    tasks.push_back(tmp);
  }

  // --- Define Bspline ---
  Bspline bspline(16, 200, 5, 6);

  // --- Define Guess ---
  std::vector<Guess> guess_predefined;


  auto xlen        = bspline.xlen(tasks[0]);
  guess_predefined = read_array_from_csv(input_file_initial_guess, xlen);
  Assert(!guess_predefined.empty());

  // -- Define Constraints ---
  Constraints<GenericManipulator> constraints;
  constraints.show_info = true;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.torque       = true;

  constraints.tool_speed              = false;
  constraints.self_collisions         = false;
  constraints.external_collisions     = false;
  constraints.n_collision_constraints = 100;
  constraints.n_collision_skip        = 2;

  // --- Define Objective ---
  Objective<GenericManipulator> objective;
  objective.time_weight = 1;

  real start_time = 0.0;
  for (u32 t = 0; t < tasks.size(); t++) {
    Optimization<GenericManipulator> opt(link6, tasks[t]);

    opt.set_bspline(bspline);
    opt.set_guess(guess_predefined[t]);
    opt.set_constraints(constraints);
    opt.set_objective(objective);
    auto results = optimize(&opt, OptimizationMethod::baseline);

    Optimization<GenericManipulator> opt_acc(link6, tasks[t]);
    opt_acc.set_bspline(bspline);
    opt_acc.set_guess(guess_predefined[t]);
    opt_acc.set_constraints(constraints);
    opt_acc.set_objective(objective);

    auto result_acc = optimize(&opt_acc, OptimizationMethod::with_analytical_pva);

    // CHECK(results.opt.constraints.grad_list.size() == result_acc.opt.constraints.grad_list.size());
    // CHECK(results.opt.constraints.constr_list.size() == result_acc.opt.constraints.constr_list.size());
    // CHECK(results.opt.constraints.x_list.size() == result_acc.opt.constraints.x_list.size());
    // CHECK(is_close(results.x, result_acc.x));
    // for (u32 i = 0; i <results.opt.constraints.grad_list.size(); i++) {
    u32 i = 0;
    for (u32 j = 0; j < results.opt.constraints.grad_list[i].rows; j++) {
      for (u32 k = 0; k < results.opt.constraints.grad_list[i].cols; k++) {
        // if (!is_close(results.opt.constraints.grad_list[i](j, k), result_acc.opt.constraints.grad_list[i](j, k), 0.01))
        //     std::cout << "i: " << i << " j: " << j << " k: " << k << " grad: " << results.opt.constraints.grad_list[i](j, k) << " grad_acc: " << result_acc.opt.constraints.grad_list[i](j, k) << std::endl;
        CHECK(is_close(results.opt.constraints.grad_list[i](j, k), result_acc.opt.constraints.grad_list[i](j, k), 0.01));
        CHECK(is_close(results.opt.constraints.constr_list[i][k], result_acc.opt.constraints.constr_list[i][k], 0.01));
      }
    }
    // }
  }
}
