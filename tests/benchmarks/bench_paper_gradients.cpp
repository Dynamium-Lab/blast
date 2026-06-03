//
// Created by nikos on 2025-09-05.
//
#include <blast>
// #include "../tests/test_helper/test_helper.hpp"

#define BLAST_TRACE_LEVEL 0
// #define BLAST_USE_NATIVE_SQP

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include <fcntl.h> // _O_WRONLY
#include <io.h>    // _open, _dup, _dup2, _close

using namespace blast;

struct IOSilencer {
  int             saved_stdout_fd = -1;
  int             null_fd         = -1;
  std::streambuf* old_cout        = nullptr;
  std::streambuf* old_cerr        = nullptr;
  std::ofstream   null_stream;

  IOSilencer() {
    // Flush buffers
    std::cout.flush();
    std::cerr.flush();
    fflush(stdout);

    // Redirect std::cout and std::cerr to /dev/null
    null_stream.open("NUL"); // Windows null device
    old_cout = std::cout.rdbuf();
    old_cerr = std::cerr.rdbuf();
    std::cout.rdbuf(null_stream.rdbuf());
    std::cerr.rdbuf(null_stream.rdbuf());

    // Redirect printf (C stdout)
    saved_stdout_fd = _dup(_fileno(stdout));
    null_fd         = _open("NUL", _O_WRONLY);
    _dup2(null_fd, _fileno(stdout));
  }

  ~IOSilencer() {
    std::cout.rdbuf(old_cout);
    std::cerr.rdbuf(old_cerr);
    fflush(stdout);

    _dup2(saved_stdout_fd, _fileno(stdout));
    _close(null_fd);
    _close(saved_stdout_fd);
  }
};

struct Config {
  int  n_ctrl   = 0;
  int  n_points = 0;
  int  task_idx = 0;
  int  n_optim  = 0;
  Task task;
};

// We create a config list that will be used for all benchmarks
inline void fill_config_list(std::array<Config, 13>& config_list) {

  std::array<std::tuple<int, int>, 1> bspline_config_list = {std::make_tuple(16, 110)}; // These give rounded number of n_points_per_segment
  std::vector<Task>                   tasks               = get_Link6_demo1_tasks();    // 13 tasks

  int config_idx = 0;
  int task_id    = 0;

  for (auto& task: tasks) {
    for (auto& bspline_config: bspline_config_list) {
      auto [n_ctrl, n_points]          = bspline_config;
      config_list[config_idx].n_ctrl   = n_ctrl;
      config_list[config_idx].n_points = n_points;
      config_list[config_idx].task_idx = task_id;
      config_list[config_idx].n_optim  = 1;
      config_list[config_idx++].task   = task;
    }
    task_id++;
  }
}

TEST_CASE("Acceleration tests", "[ICRA]") {
  std::array<Config, 13> config_list;
  fill_config_list(config_list);

  real eps = 0.01;

  for (int config_id = 0; config_id < config_list.size(); config_id++) {
    auto n_ctrl   = config_list[config_id].n_ctrl;
    auto n_points = config_list[config_id].n_points;
    auto task_idx = config_list[config_id].task_idx;
    auto n_optim  = config_list[config_id].n_optim;
    auto task     = config_list[config_id].task;
    // Optimization opt(get_generic_gen3(), get_gen3_task());
    Optimization opt(get_generic_Link6(), task);

    opt.world = get_lab_world();

    Bspline bspline(n_ctrl, n_points, 5, 6);
    opt.bspline = bspline;

    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tool_speed          = true;
    opt.constraints.self_collisions     = true;
    opt.constraints.external_collisions = true;

    opt.constraints.n_collision_skip = 1;

    opt.max_tries         = 1;
    opt.success_tolerance = 0.01;

    opt.guess.type      = Guess::custom;
    opt.guess.initial_x = guess_random((opt.bspline), opt.task);

    cout << "Config ID:                  " << config_id << endl;
    cout << "Task id:                    " << task_idx << endl;
    cout << "n_ctrl:                     " << n_ctrl << endl;
    cout << "n_points:                   " << n_points << endl;

    Result result(&opt);
    Result result_spdup1(&opt);
    Result result_spdup2(&opt);
    Result result_spdup3(&opt);
    for (int iter = 0; iter < n_optim; iter++) {
      {
        {
          IOSilencer _;
          result = optimize(&opt, OptimizationMethod::baseline);
        }

        {
          IOSilencer _;
          result_spdup1 = optimize(&opt, OptimizationMethod::with_analytical_pva);
        }

        {
          IOSilencer _;
          result_spdup2 = optimize_acc2(&opt, OptimizationMethod::with_analytical_dynamics);
        }

        {
          IOSilencer _;
          result_spdup3 = optimize(&opt, OptimizationMethod::with_segments);
        }
        CHECK(result_spdup1.success == result.success);
        CHECK(result_spdup2.success == result.success);
        CHECK(result_spdup3.success == result.success);
        CHECK(result_spdup1.success_false == result.success_false);
        CHECK(result_spdup2.success_false == result.success_false);
        CHECK(result_spdup3.success_false == result.success_false);
        CHECK(std::abs(result.x.back() - result_spdup1.x.back()) < 0.01);
        CHECK(std::abs(result.x.back() - result_spdup2.x.back()) < 0.01);
        CHECK(std::abs(result.x.back() - result_spdup3.x.back()) < 0.01);
      }
    }
  }
};
