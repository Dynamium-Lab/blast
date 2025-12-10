
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
  std::vector<Array> start_pos(10);
  start_pos[0] = {0.6962, -2.6857, 0.1043, -2.4137, 0.2672, 3.1071};
  start_pos[1] = {0.71916, -2.62426, 0.0011186, -2.3251, 0.372545, 2.89637};
  start_pos[2] = {-2.88535, 2.30619, -0.519746, 0.274299, 0.304614, 3.0678};
  start_pos[3] = {0.130399, -2.99611, -0.368365, -2.98505, 0.508639, 3.07793};
  start_pos[4] = {0.860489, -2.571, -0.0154261, -2.18895, 0.367129, 2.68166};
  start_pos[5] = {-2.3923, 2.1047, -0.0988, 1.0073, 0.6294, 2.3977};
  start_pos[6] = {-2.67318, 2.87552, 1.15243, 1.28614, 1.07142, 1.91098};
  start_pos[7] = {-2.52236, 1.6071, -1.40444, 0.630671, 0.105179, 3.04982};
  start_pos[8] = {0.5436, -1.4681, 2.3156, -2.4816, -0.5336, -2.7410};
  start_pos[9] = {-2.1372, 1.9866, -0.8635, -2.1164, 2.9868, -0.0970};

  std::vector<Array> end_pos(10);
  end_pos[0] = {0.3750, -0.4404, -0.0142, -0.4393, 0.3909, 0.1175};
  end_pos[1] = {-0.134401, -0.333873, 0.287012, 0.160065, 0.614795, -0.0852355};
  end_pos[2] = {0.6349, -0.300309, 0.896664, -1.13046, 0.825613, 1.05578};
  end_pos[3] = {0.151525, -0.861139, -1.03372, -0.15689, -0.170655, 0.0701045};
  end_pos[4] = {0.906918, -0.714292, 0.147075, -1.12048, 0.468276, 0.805861};
  end_pos[5] = {-2.84675, 1.59479, 1.26879, 2.8216, 0.310568, 0.129916};
  end_pos[6] = {0.776273, -1.65795, -1.38197, -0.809817, 0.192536, 0.227118};
  end_pos[7] = {0.02013, -0.0137358, 1.79172, 0.147243, 1.80794, -0.171793};
  end_pos[8] = {-0.0245301, 0.00536648, 1.5534, -2.45517, 1.60003, 2.46791};
  end_pos[9] = {-2.34874, 1.59395, 0.961949, 2.23711, 0.424871, 0.471586};

  // Fill task and add to list
  Matrix task(6, 6);
  for (auto& s: start_pos) {
    for (auto& e: end_pos) {
      fill_positions(task, s, e);
      result.push_back(task);
    }
  }

  return result;
}

inline World get_kitchen_open_doors() {
  World world;
  world.add_box(Vec3{0.427, -1.0847, 1.1974}, Vec3{0.36029, 0.06, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, 0.48791, 1.1974}, Vec3{0.36029, 0.06, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 1.5298}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 1.2074}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 0.86498}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.72729, -0.29838, 1.2074}, Vec3{0.06, 0.84629, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{-0.20911, 0.86968, 1.2009}, Vec3{0.032523, 0.42256, 0.36242}, Mat3{-0.70706, -0.70716, 0, -0.70716, -0.70706, 0, 0, 0, 1});
  world.add_box(Vec3{-0.20549, -1.4701, 1.2009}, Vec3{0.037648, 0.42257, 0.36244}, Mat3{-0.70706, 0.70716, 0, 0.70716, -0.70706, 0, 0, 0, 1});
  world.add_box(Vec3{-0.47232, -0.77127, -0.48346}, Vec3{0.023735, 0.41378, 0.41552}, Mat3{4.6327e-05, 0, 1, 0, 1, 0, 1, 0, 4.6327e-05});
  world.add_box(Vec3{0.35698, 0.425, -0.066059}, Vec3{0.41378, 0.8425, 0.43736}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -1.1546, -0.066059}, Vec3{0.41378, 0.03, 0.43736}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, 0.31136}, Vec3{0.41378, 0.38356, 0.06}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, -0.055}, Vec3{0.41378, 0.38356, 0.025}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, -0.47336}, Vec3{0.41378, 0.38356, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  return world;
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
      fill_array(res[i * n_tests + j], i, arr);
    }
  }

  std::ofstream f(path);
  f << "{\n\"results\": " << arr.dump(2) << "\n}\n";
}

int main() {
  int n_tests = 10;

  auto manip   = get_generic_Link6_fixed();
  manip.p_base = {-0.5, -0.3, 0.35}; // link6
  manip.Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  auto world   = get_kitchen_open_doors();

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
  res.reserve(tasks.size() * n_tests);

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
    for (int i = 0; i < n_tests; i++) {
      opt.guess.x0 = guess_random(opt.bspline, result.opt->task);
      // opt.guess.x0.back() = 0.5;
      result      = optimize_with_segments(&opt);
      tasks[t_id] = result.opt->task;
      // result = optimize(&opt);

      res.push_back(result);
    }
  }

  std::cout << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Results (compute and trajectory time is average)" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;

  real total_compute_time    = 0.0;
  real total_trajectory_time = 0.0;
  real total_convergence     = 0.0;
  for (int i = 0; i < tasks.size(); i++) {
    real compute_time    = 0.0;
    real trajectory_time = 0.0;
    real convergence     = 0.0;
    for (int j = 0; j < n_tests; j++) {
      if (res[i * n_tests + j].success && !res[i * n_tests + j].success_false) {
        compute_time += res[i * n_tests + j].compute_time;
        trajectory_time += res[i * n_tests + j].x.back();
        convergence++;
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
  std::cout << std::endl
            << "Printing trajectory..." << std::endl;
  Matrix trajectory(0, manip.n_joints);
  for (int i = 0; i < tasks.size(); i++) {
    for (int j = 0; j < n_tests; j++) {
      if (res[i * n_tests + j].success && !res[i * n_tests + j].success_false) {
        int points_more = (int) std::ceil(res[i * n_tests + j].x.back() * 1000.0) + 1;

        Bspline bspline_val_more(res[i * n_tests + j].opt->bspline.n_ctrl, points_more, res[i * n_tests + j].opt->bspline.p, res[i * n_tests + j].opt->manip.n_joints); // todo: this is expensive
        bspline_val_more.compute_trajectory(res[i * n_tests + j].x, tasks[i]);

        trajectory = append_rows(trajectory, transpose(bspline_val_more.traj.pos));
        break;
      }
    }
  }
  print_to_csv(trajectory, "../../../examples/Demo3/Statistics/trajectory_kitchen_Link6.csv");
  std::cout << "Trajectory printed." << std::endl;

  std::cout << "Printing results..." << std::endl;
  print_to_json(res, n_tests, "../../../examples/Demo3/Statistics/results_kitchen_Link6.json");
  std::cout << std::endl;
  std::cout << "Results printed." << std::endl;

  return 0;
}
