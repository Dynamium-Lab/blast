
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
  Array pos_0 = {0, 0, 0, 0, 0, 0, 0};
  // Array pos_bin = wrap2pi(deg2rad(    {-44.7, 33.9,  162.4, 264.0, 84.6, 306.9, -26}));

  // Array pos_b2_10cm = wrap2pi(deg2rad({36.1,  79.9,  114.8, 339.7, 41.9, 321.6, 116.2}));
  // Array pos_b2 = wrap2pi(deg2rad(     {37.6,  89.1,  114.6, 335.7, 42.7, 332.0, 119.4}));

  // Array pos_b5_10cm = wrap2pi(deg2rad(     {73.3, 88.8,  98.4, 318.0, 43.4, 321.6, 142.4}));
  // Array pos_b5 = wrap2pi(deg2rad({72.4, 98.5,  98.5, 318.8, 42.0, 315.7, 135.7}));

  // Array pos_b6_10cm = wrap2pi(deg2rad({69.2,  59.7,  122.7, 274.6, 55.1, 352.1, 87.7}));
  // Array pos_b6 = wrap2pi(deg2rad(     {62.2,  64.9,  142.5, 270.5, 53.1, 347.5,  88.7}));

  Array pos_bin = wrap2pi(deg2rad(    {302.097, 34.659, 159.209, 264.607, 78.817, 294.022, 328.983}));

  // Array pos_b2_10cm = wrap2pi({ 1.1067,  1.4473,  1.5834, -0.7312,  0.6452, -0.7130,  2.3648 });
  // Array pos_b2 = wrap2pi({ 1.2636,  1.7191,  1.7191, -0.7191,  0.7330, -0.7732,  2.3684 });

  Array pos_b5_10cm = wrap2pi({ 0.6761,  0.9147,  1.2801, -0.6486,  0.4588, -0.5870,  2.3513 });
  Array pos_b5 = wrap2pi({ 0.8933,  1.1840,  1.4310, -0.6902,  0.5524, -0.6508,  2.3581 });

  Array pos_b6_10cm = wrap2pi({ 0.2450,  0.3815,  0.9775, -0.5638,  0.2701, -0.4574,  2.3381 });
  Array pos_b6 = wrap2pi({ 0.4605,  0.6483,  1.1282, -0.6064,  0.3646, -0.5223,  2.3447 });

  // Fill task and add to list
  Matrix task(7,6);
  // fill_positions(task, pos_0, pos_b2_10cm);
  // result.push_back(task);
  // fill_positions(task, pos_b2_10cm, pos_b2);
  // result.push_back(task);
  // fill_positions(task, pos_b2, pos_b2_10cm);
  // result.push_back(task);
  // fill_positions(task, pos_b2_10cm, pos_bin);
  // result.push_back(task);

  // fill_positions(task, pos_bin, pos_b5_10cm);
  // result.push_back(task);
  
  fill_positions(task, pos_0, pos_b5_10cm);
  result.push_back(task);
  fill_positions(task, pos_b5_10cm, pos_b5);
  result.push_back(task);
  fill_positions(task, pos_b5, pos_b5_10cm);
  result.push_back(task);
  fill_positions(task, pos_b5_10cm, pos_bin);
  result.push_back(task);
  
  fill_positions(task, pos_bin, pos_b6_10cm);
  result.push_back(task);
  fill_positions(task, pos_b6_10cm, pos_b6);
  result.push_back(task);
  fill_positions(task, pos_b6, pos_b6_10cm);
  result.push_back(task);
  fill_positions(task, pos_b6_10cm, pos_bin);
  result.push_back(task);

  fill_positions(task, pos_bin, pos_0);
  result.push_back(task);
  
  return result;
}

Matrix append_rows(const Matrix& m1, const Matrix& m2) {
  Assert(m1.cols == m2.cols);
  Matrix result(m1.rows+m2.rows, m1.cols);
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
  json arr = json::array();
  int n_tasks = res.size()/n_tests;
  for (int i = 0; i < n_tasks; i++) {
    for (int j = 0; j < n_tests; j++) {
      fill_array(res[i*n_tests + j], i, arr);
    }
  }
  
  std::ofstream f(path);
  f << "{\n\"results\": " << arr.dump(2) << "\n}\n";
}

int main() {
  int n_tests = 1;
  
  // todo: confirm manip
  auto manip = get_generic_gen3_fixed();
  manip.p_base = {1.27, 0.05, 0.7112};
  manip.Q_base = {-1, 0, 0, 0, -1, 0, 0, 0, 1};
  auto world = get_demo2_world();

  Bspline bspline(12, 70, 5, manip.n_joints);

  ConstraintSelection cons;
  cons.position            = true;
  cons.velocity            = true;
  cons.acceleration        = true;
  cons.torque              = true;
  cons.tcp_speed           = true;
  cons.self_collisions     = true;
  cons.external_collisions = true;
  cons.n_collision_constraints = 1;

  Guess guess;
  guess.type = Guess::custom;

  auto tasks = get_tasks();
  // tasks.resize(1);

  std::vector<Result> res;
  res.reserve(tasks.size()*n_tests);

  real total_time = 0.0;
  for (int t_id = 0; t_id < tasks.size(); t_id++) {
    cout << "--- Task " << t_id << " of " << tasks.size() << " ---" << std::endl;
    auto t = tasks[t_id];

    Optimization opt(manip, t);
    opt.world = world;
    opt.bspline = bspline;
    opt.constraints = cons;

    opt.guess = guess;

    opt.max_tries                           = 1;
    opt.success_tolerance                   = 0.01;

    Result result(&opt);
    for (int i = 0; i < n_tests; i++) {
      opt.guess.x0        = random_array(opt.bspline.x_len(opt.task), 1);
      opt.guess.x0.back() = 0.5;
      result = optimize_with_segments(&opt);
      tasks[t_id] = result.opt->task;
      // result = optimize(&opt);
      
      res.push_back(result);
    }
  }

  std::cout << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Results (compute and trajectory time is average)" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;

  real total_compute_time = 0.0;
  real total_trajectory_time = 0.0;
  real total_convergence = 0.0;
  for (int i = 0; i < tasks.size(); i++) {
    real compute_time = 0.0;
    real trajectory_time = 0.0;
    real convergence = 0.0;
    for (int j = 0; j < n_tests; j++) {
      if (res[i*n_tests + j].success && !res[i*n_tests + j].success_false) {
        compute_time += res[i*n_tests + j].compute_time;
        trajectory_time += res[i*n_tests + j].x.back();
        convergence ++;
      }
    }
    std::cout << "- Task " << i << " -" << std::endl;
    std::cout << " Convergence: " << convergence / n_tests << std::endl;
    std::cout << " Compute time: " << compute_time / n_tests << std::endl;
    std::cout << " Trajectory time: " << trajectory_time / n_tests << std::endl;

    total_compute_time += compute_time / n_tests;
    total_trajectory_time += trajectory_time / n_tests;
    total_convergence += convergence / n_tests;
  }

  std::cout << "- Total results -" << std::endl;
  std::cout << " Convergence: " << total_convergence / tasks.size() << std::endl;
  std::cout << " Compute time (avg): " << total_compute_time / tasks.size() << std::endl;
  std::cout << " Compute time (total): " << total_compute_time << std::endl;
  std::cout << " Trajectory time (avg): " << total_trajectory_time / tasks.size() << std::endl;
  std::cout << " Trajectory time (total): " << total_trajectory_time << std::endl;

  // print trajectory
  std::cout << std::endl << "Printing trajectory..." << std::endl;
  Matrix trajectory(0, manip.n_joints);
  for (int i = 0; i < tasks.size(); i++) {
    for (int j = 0; j < n_tests; j++) {
      if (res[i*n_tests + j].success && !res[i*n_tests + j].success_false) {
        int points_more = (int) std::ceil(res[i*n_tests + j].x.back() * 1000.0) + 1;

        Bspline bspline_val_more(res[i*n_tests + j].opt->bspline.n_ctrl, points_more, res[i*n_tests + j].opt->bspline.p, res[i*n_tests + j].opt->manip.n_joints); // todo: this is expensive
        bspline_val_more.compute_trajectory(res[i*n_tests + j].x, tasks[i]);

        trajectory = append_rows(trajectory, transpose(bspline_val_more.traj.pos));
        break;
      }
    }
  }
  print_to_csv(trajectory, "../../../examples/Demo3/Trajectories/trajectory_demo2_gen3.csv");
  std::cout << "Trajectory printed." << std::endl;

  std::cout << "Printing results..." << std::endl;
  print_to_json(res, n_tests, "../../../examples/Demo3/Trajectories/results_demo2_gen3.json");
  std::cout << std::endl;
  std::cout << "Results printed." << std::endl;

  return 0;
}
