
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
  Array pos_0 = {0, 0, 0, 0, 0, 0};

  Array pos_1      = deg2rad({-40.31, -50.85, 99.75, 1.88, 61.1, -43.1});
  Array pos_1_10cm = deg2rad({-40.445762634277344, -26.876392364501953, 83.60868835449219, 1.49383544921875, 19.951095581054688, -42.22943115234375});

  Array pos_2      = deg2rad({-23.56, -44.99, 125.77, -1.89, 78.84, -24.38});
  Array pos_2_10cm = deg2rad({-24.826316833496094, -9.07125473022461, 110.22100830078125, 1.2542877197265625, 30.40911865234375, -27.223480224609375});

  Array pos_3      = deg2rad({-35.7, -55.6, 86.55, 0.74, 53.04, -25.71});
  Array pos_3_10cm = deg2rad({-33.93860626220703, -35.56209945678711, 70.89104461669922, -0.6509857177734375, 16.881134033203125, -31.8890380859375});

  Array pos_4      = deg2rad({-20.64, -46.93, 112.33, 4.27, 70.27, -25.71});
  Array pos_4_10cm = deg2rad({-19.595504760742188, -19.332134246826172, 95.92303466796875, -3.042205810546875, 23.505020141601562, -17.61700439453125});

  Array pos_white_bin = deg2rad({62.90712356567383, -2.7591018676757812, 119.23771667480469, 8.18115234375, 28.954452514648438, 47.103912353515625});

  Array pos_black_bin = deg2rad({27.443552017211914, -54.93880081176758, 28.312477111816406, 1.807525634765625, -4.7853546142578125, 29.731155395507812});

  // Fill task and add to list
  // todo: change tasks
  Matrix task(6, 6);
  fill_positions(task, pos_0, pos_1_10cm);
  result.push_back(task);
  fill_positions(task, pos_1_10cm, pos_1);
  result.push_back(task);
  fill_positions(task, pos_1, pos_1_10cm);
  result.push_back(task);
  fill_positions(task, pos_1_10cm, pos_white_bin);
  result.push_back(task);

  fill_positions(task, pos_white_bin, pos_2_10cm);
  result.push_back(task);
  fill_positions(task, pos_2_10cm, pos_2);
  result.push_back(task);
  fill_positions(task, pos_2, pos_2_10cm);
  result.push_back(task);
  fill_positions(task, pos_2_10cm, pos_black_bin);
  result.push_back(task);

  fill_positions(task, pos_black_bin, pos_3_10cm);
  result.push_back(task);
  fill_positions(task, pos_3_10cm, pos_3);
  result.push_back(task);
  fill_positions(task, pos_3, pos_3_10cm);
  result.push_back(task);
  fill_positions(task, pos_3_10cm, pos_white_bin);
  result.push_back(task);

  fill_positions(task, pos_white_bin, pos_4_10cm);
  result.push_back(task);
  fill_positions(task, pos_4_10cm, pos_4);
  result.push_back(task);
  fill_positions(task, pos_4, pos_4_10cm);
  result.push_back(task);
  fill_positions(task, pos_4_10cm, pos_white_bin);
  result.push_back(task);

  fill_positions(task, pos_white_bin, pos_0);
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

template<typename T_manip>
inline void add_robot_obstacles(T_manip manip, const Matrix& trajectory, const real start_time, const real end_time, World* world) {
  std::vector<std::vector<Capsule>> capsules;

  ManipulatorTempData data;

  capsules.resize(manip._collision_model.size());
  for (u32 i = 0; i < manip._collision_model.size(); i++)
    capsules[i].reserve(trajectory.cols);

  Capsule temp_caps;
  for (u32 i = 0; i < trajectory.cols; i++) {
    forward_kinematics(manip, data, trajectory.col(0));
    compute_capsules(manip, data);
    for (u32 j = 0; j < manip.n_joints; j++) {
      capsules[j].push_back(data.capsule_list[j]);
    }
  }

  for (u32 i = 0; i < manip.n_joints; i++)
    world->add_dynamic_capsule(capsules[i], (u32) capsules[i].size(), start_time, end_time);
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
      fill_array(res[i + j], i, arr);
    }
  }

  std::ofstream f(path);
  f << "{\n\"results\": " << arr.dump(2) << "\n}\n";
}

int main() {
  auto manip   = get_generic_Link6_fixed();
  manip.p_base = {0, 0, 0.833}; // link6
  manip.Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  auto world   = get_demo2_world();

  auto inflated_gen3   = get_generic_gen3_fixed();
  inflated_gen3.p_base = {1.28, 0.025, 0.73};
  inflated_gen3.Q_base = {-1, 0, 0, 0, -1, 0, 0, 0, 1};
  for (auto& caps: inflated_gen3._collision_model) {
    caps.r += 0.05; // 5 cm buffer
  }

  auto trajectory = read_csv_trajectory_no_header("../../../examples/Demo3/Trajectories/trajectory_full_gen3.csv", ",");

  auto gen3_trajectory_pos = trajectory.pos;
  add_robot_obstacles(inflated_gen3, gen3_trajectory_pos, 0.0, gen3_trajectory_pos.cols / 1000, &world);

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

    opt.trajectory_start_time = total_time;

    Result result(&opt);
    while (true) {
      opt.guess.x0 = guess_random(opt.bspline, result.opt->task);
      // opt.guess.x0.back() = 0.5;
      result      = optimize_with_segments(&opt);
      tasks[t_id] = result.opt->task;
      // result = optimize(&opt);

      if (result.success && !result.success_false) {
        total_time += result.x.back();
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

  print_to_csv(trajectories, "../../../examples/Demo3/Trajectories/trajectory_full_link6.csv");
  std::cout << "Trajectory printed." << std::endl;

  std::cout << "Printing results..." << std::endl;
  print_to_json(res, 1, "../../../examples/Demo3/Trajectories/results_demo2_Link6.json");
  std::cout << std::endl;
  std::cout << "Results printed." << std::endl;

  return 0;
}
