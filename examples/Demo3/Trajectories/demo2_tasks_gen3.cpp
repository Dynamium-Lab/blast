
#define BLAST_TRACE_LEVEL 0
#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"

using namespace blast;

void fill_positions(Matrix& task, const Array& p1, const Array& p2) {
  Assert(task.rows == p1.size);
  Assert(task.rows == p2.size);
  for (int i = 0; i < task.rows; i++) {
    task(i, 0) = p1[i];
    task(i, 3) = p2[i];
  }
}

std::vector<Matrix> get_tasks() {
  std::vector<Matrix> result;

  // Define positions
  Array pos_0   = {0, 0, 0, 0, 0, 0, 0};
  Array pos_bin = wrap2pi(deg2rad({-44.7, 33.9, 162.4, 264.0, 84.6, 306.9, -26}));

  Array pos_5_10cm = wrap2pi(deg2rad({84.8, 75.3, 99.6, 301.9, 40.4, 312.1, 118.2}));
  Array pos_5      = wrap2pi(deg2rad({85.3, 86.7, 103.5, 298.6, 44.7, 322.0, 123.0}));

  Array pos_6_10cm = wrap2pi(deg2rad({66.0, 50.0, 126.2, 273.4, 41.3, 323.0, 92.7}));
  Array pos_6      = wrap2pi(deg2rad({62.7, 55.1, 134.0, 269.5, 47.0, 333.7, 88.5}));

  // Fill task and add to list
  Matrix task(7, 6);

  fill_positions(task, pos_0, pos_5_10cm);
  result.push_back(task);
  fill_positions(task, pos_5_10cm, pos_5);
  result.push_back(task);
  fill_positions(task, pos_5, pos_5_10cm);
  result.push_back(task);
  fill_positions(task, pos_5_10cm, pos_bin);
  result.push_back(task);

  fill_positions(task, pos_bin, pos_6_10cm);
  result.push_back(task);
  fill_positions(task, pos_6_10cm, pos_6);
  result.push_back(task);
  fill_positions(task, pos_6, pos_6_10cm);
  result.push_back(task);
  fill_positions(task, pos_6_10cm, pos_bin);
  result.push_back(task);

  fill_positions(task, pos_bin, pos_0);
  result.push_back(task);

  return result;
}

Matrix append_rows(const Matrix& m1, const Matrix& m2) {
  Assert(m1.cols == m2.cols);
  Matrix result(m1.rows + m2.rows, m1.cols);
  for (int i = 0; i < m1.rows; i++) {
    for (int j = 0; j < m1.cols; j++) {
      result(i, j) = m1(i, j);
    }
  }
  for (int i = 0; i < m2.rows; i++) {
    for (int j = 0; j < m2.cols; j++) {
      result(i + m1.rows, j) = m2(i, j);
    }
  }
  return result;
}

using nlohmann::json;
// helper: convert one Result to the JSON shape you’ve been using

void fill_array(const Result& res, const int task_id, json& array) {
  json row;
  row["task_idx"]                         = task_id;
  row["n_ctrl"]                           = res.opt->bspline.n_ctrl;
  row["n_points"]                         = res.opt->bspline.n_points;
  row["success"]                          = res.success;
  row["success_false"]                    = res.success_false;
  row["max_constraint_value"]             = res.max_constraint_value;
  row["max_constraint_idx"]               = res.max_constraint_idx;
  row["max_constraint_more_points_value"] = res.max_constraint_more_points_value;
  row["max_constraint_more_points_idx"]   = res.max_constraint_more_points_idx;
  row["nlopt_exit_criteria"]              = res.nlopt_exit_criteria;
  row["num_eval"]                         = res.num_eval;
  row["compute_time"]                     = res.compute_time;
  if (res.x.size)
    row["trajectory_time"] = res.x.back();

  // dump arrays x and x0
  row["x"] = json::array();
  for (u32 k = 0; k < res.x.size; ++k)
    row["x"].push_back(res.x[k]);

  row["x0"] = json::array();
  for (u32 k = 0; k < res.x0.size; ++k)
    row["x0"].push_back(res.x0[k]);

  array.push_back(std::move(row));
}

void print_to_json(const std::vector<Result>& res, const int n_tests, const std::string path) {
  json arr     = json::array();
  int  n_tasks = res.size() / n_tests;
  for (int i = 0; i < n_tasks; i++) {
    for (int j = 0; j < n_tests; j++) {
      fill_array(res[i], i, arr);
    }
  }

  std::ofstream f(path);
  f << "{\n\"results\": " << arr.dump(2) << "\n}\n";
}

int main() {
  // todo: confirm manip
  auto manip   = get_generic_gen3_fixed();
  manip.p_base = {1.28, 0.025, 0.73};
  manip.Q_base = {-1, 0, 0, 0, -1, 0, 0, 0, 1};
  auto world   = get_demo2_world();

  Bspline bspline(12, 70, 5, manip.n_joints);

  ConstraintSelection cons;
  cons.position                = true;
  cons.velocity                = true;
  cons.acceleration            = true;
  cons.torque                  = true;
  cons.tcp_speed               = true;
  cons.self_collisions         = true;
  cons.external_collisions     = true;
  cons.n_collision_constraints = 1;

  Guess guess;
  guess.type = Guess::custom;

  auto tasks = get_tasks();
  // tasks.resize(1);

  std::vector<Result> res;
  res.reserve(tasks.size());

  real total_time = 0.0;
  for (int t_id = 0; t_id < tasks.size(); t_id++) {
    cout << "--- Task " << t_id << " of " << tasks.size() << " ---" << std::endl;
    auto t = tasks[t_id];

    Optimization opt(manip, t);
    opt.world       = world;
    opt.bspline     = bspline;
    opt.constraints = cons;

    opt.guess = guess;

    opt.max_tries         = 1;
    opt.success_tolerance = 0.01;

    Result result(&opt);
    while (true) {
      opt.guess.x0        = random_array(opt.bspline.x_len(opt.task), 1);
      opt.guess.x0.back() = 0.5;
      result              = optimize_with_segments(&opt);
      tasks[t_id]         = result.opt->task;
      // result = optimize(&opt);
      if (result.success && !result.success_false) {
        res.push_back(result);
        break;
      }
    }
  }

  // print trajectory
  std::cout << "Printing trajectory..." << std::endl;
  std::vector<Trajectory> trajectories;
  trajectories.reserve(tasks.size());
  for (auto& results: res)
    trajectories.push_back(results.trajectory);

  print_to_csv(trajectories, "../../../examples/Demo3/Trajectories/trajectory_full_gen3.csv");

  std::cout << "Trajectory printed." << std::endl;

  std::cout << "Printing results..." << std::endl;
  print_to_json(res, 1, "../../../examples/Demo3/Trajectories/results_demo2_gen3.json");
  std::cout << std::endl;
  std::cout << "Results printed." << std::endl;

  return 0;
}
