#define BLAST_USE_NATIVE_SQP

#include <blast>
#include "test_helper/test_helper.hpp"

#include <iomanip>

using blast::get_generic_link6_opt;
using blast::Optimization;
using blast::OptimizationMethod;
using blast::optimize;
using blast::real;
using blast::Result;

struct MethodStats {
  real success_rate         = 0;
  real mean_compute_time_ms = 0;
  real mean_traj_time_s     = 0;
  real mean_max_constraint  = 0;
};

inline MethodStats run_trials(OptimizationMethod method, int n_trials) {
  using namespace blast;
  int  n_success            = 0;
  real total_compute_time   = 0;
  real total_traj_time      = 0;
  real total_max_constraint = 0;


  // Optimization opt(get_generic_gen3(), get_gen3_task());
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 110, 5, 6);
  opt.bspline = bspline;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = true;
  opt.constraints.self_collisions     = true;
  opt.constraints.external_collisions = true;

  opt.constraints.n_collision_constraints = 1;
  opt.max_tries                           = 1;
  opt.success_tolerance                   = 0.01;

  opt.guess.type = Guess::custom;

  auto   t1 = get_tick_us();
  Result result(&opt);

  for (int i = 0; i < n_trials; i++) {
    opt.guess.initial_x = guess_random(opt.bspline, opt.task);

    result = optimize(&opt, method);

    if (result.success) {
      n_success++;
      total_traj_time += result.x.back();
      total_compute_time += result.compute_time;
      total_max_constraint = std::max(total_max_constraint, result.max_constraint_value);
    }
  }

  MethodStats s;
  s.success_rate         = 100.0 * (real) n_success / n_trials;
  s.mean_compute_time_ms = total_compute_time / n_trials;
  s.mean_traj_time_s     = n_success > 0 ? total_traj_time / n_success : 0;
  s.mean_max_constraint  = total_max_constraint / n_trials;
  return s;
}

int main() {
  constexpr int n_trials = 10;

  struct Entry {
    OptimizationMethod method;
    const char*        name;
  };

  constexpr std::array<Entry, 4> methods = {{
          {OptimizationMethod::with_segments, "with_segments"},
          {OptimizationMethod::baseline, "baseline"},
          {OptimizationMethod::with_analytical_pva, "with_analytical_pva"},
          {OptimizationMethod::with_analytical_dynamics, "with_analytical_dynamics"},
  }};

  std::cout << "\n--- Optimization Method Benchmark (" << n_trials << " trials each) ---\n";
  std::cout << std::left
            << std::setw(32) << "Method"
            << std::setw(14) << "Success"
            << std::setw(20) << "Compute time (ms)"
            << std::setw(16) << "Traj time (s)"
            << "Max constraint\n";
  std::cout << std::string(92, '-') << "\n";

  for (const auto& e: methods) {
    auto s = run_trials(e.method, n_trials);
    std::cout << std::left
              << std::setw(32) << e.name
              << std::setw(14) << s.success_rate
              << std::setw(20) << s.mean_compute_time_ms
              << std::setw(16) << s.mean_traj_time_s
              << s.mean_max_constraint << std::endl;
  }

  return 0;
}
