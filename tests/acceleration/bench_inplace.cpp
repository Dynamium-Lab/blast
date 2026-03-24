#include <blast>
#include "test_helper/test_helper.hpp"

#define BLAST_TRACE_LEVEL 3

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"

#include <fcntl.h> // _O_WRONLY
#include <io.h>    // _open, _dup, _dup2, _close

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

namespace blast {

TEST_CASE("compute_constraints vs constraints_with_segments", "[acceleration]") {
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 110, 5, 6);
  opt.bspline = bspline;

  auto x_len    = opt.bspline.x_len(opt.task);
  auto n_joints = opt.manip.n_joints;

  opt.guess.type = Guess::random;


  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = false;
  opt.constraints.self_collisions     = false;
  opt.constraints.external_collisions = false;

  opt.max_tries         = 1;
  opt.success_tolerance = 0.01;

  real eps = 1e-5;

  initialize_optimization(&opt);
  n_con(&opt);

  auto n_constraints = opt.constraints.n_constraints;

  auto x = init_guess(&opt);
  bspline.compute_trajectory(x, opt.task);

  BENCHMARK_ADVANCED("Original nlopt_constraints function")(Catch::Benchmark::Chronometer meter) {
    Array result(n_constraints);
    Array grad(x_len * n_constraints);
    meter.measure([&] {
      nlopt_constraints(n_constraints, result.data, x_len, x.data, grad.data, &opt);
      return result;
    });
  };

  BENCHMARK_ADVANCED("New nlopt_constraints_with_segments function")(Catch::Benchmark::Chronometer meter) {
    auto  n_constraints_with_segments = (n_joints * 4) * ((int) opt.bspline.n_ctrl - (int) opt.bspline.degree);
    Array result(n_constraints_with_segments);
    Array grad(x_len * n_constraints_with_segments);
    meter.measure([&] {
      nlopt_constraints_with_segments(n_constraints_with_segments, result.data, x_len, x.data, grad.data, &opt);
      return result;
    });
  };
};

TEST_CASE("Optimization with compute_constraints vs constraints_and_gradients_with_segments", "[acceleration]") {
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 110, 5, 6);
  opt.bspline = bspline;

  auto x_len    = opt.bspline.x_len(opt.task);
  auto n_joints = opt.manip.n_joints;

  // opt.guess.type = Guess::random;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = false;
  opt.constraints.self_collisions     = false;
  opt.constraints.external_collisions = false;

  opt.max_tries         = 1;
  opt.success_tolerance = 0.01;

  // initialize_optimization(&opt);
  // n_con(&opt);
  // note: use the same initial guess
  opt.guess.type = Guess::custom;
  opt.guess.initial_x   = Array(x_len, 2.0);
  // note: this only works with 16 n_ctrl
  // opt.guess.initial_x   = Array{0.152355, -0.597958, -0.756348, 0.558375, -0.946557, -0.10268, 0.901938, -0.145201, -0.436917, -0.827993,
  //                      -0.210473, -0.0664338, 0.370776, -0.29955, -0.139142, 0.757138, -0.157328, 0.909626, -0.214084, 0.745555,
  //                      -0.0416009, 0.375275, 0.126153, 0.402337, 0.723458, -0.070362, -0.444889, 0.784517, 0.655188, 0.216073,
  //                      0.889643, -0.969162, -0.0279569, 0.496376, -0.891977, 0.0247388, -0.332795, -0.239096, -0.418835, 0.922935,
  //                      0.124817, 0.34831, -0.00905602, -0.157159, 0.492551, 0.0325292, -0.299891, 0.828649, -0.909505, -0.294613,
  //                      -0.768063, -0.489629, 0.0918699, -0.603901, -0.391773, -0.919239, 0.472415, 0.591886, -0.494069, 0.0488416,
  //                      4.81184};

  Array x_orig = opt.guess.initial_x;

  Result result(&opt);
  Result result_with_segments(&opt);

  // result_with_segments = optimize_with_segments(&opt);

  // cout << "Compute time:        " << result_with_segments.compute_time << endl;
  // cout << "Function evaluations:" << result_with_segments.num_eval << endl;
  // cout << "Trajectory time:     " << result_with_segments.x[result_with_segments.x.size - 1] << endl;
  // cout << "Success:             " << result_with_segments.success << endl;
  // cout << "False success:       " << result_with_segments.success_false << endl;
  // cout << "Number of tries:     " << result_with_segments.num_tries << endl;
  // cout << "Stopping criteria:   " << result_with_segments.nlopt_exit_criteria << endl;

  // CHECK(result_with_segments.x.size > 0);

  // nlopt_result sqp(blast::Array& x, blast::Optimization& opt, nlopt_stopping* stop);
  BENCHMARK_ADVANCED("OG")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      IOSilencer _;
      result = optimize(&opt);
      return result;
    });
  };

  // opt.guess.type = Guess::custom;
  opt.guess.initial_x = x_orig;

  BENCHMARK_ADVANCED("New OG")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      IOSilencer _;
      result_with_segments = optimize_with_segments(&opt);
      return result;
    });
  };
  cout << endl;
  cout << "=============== Optimization Parameters ===============" << endl;
  cout << " n_ctrl : " << opt.bspline.n_ctrl << endl;
  cout << " n_points : " << opt.bspline.n_points << endl;
  cout << endl;
  cout << "================== Result of last OG ==================" << endl;
  cout << "Compute time:        " << result.compute_time << endl;
  cout << "Function evaluations:" << result.num_eval << endl;
  cout << "Trajectory time:     " << result.x[result.x.size - 1] << endl;
  cout << "Success:             " << result.success << endl;
  cout << "False success:       " << result.success_false << endl;
  cout << "Number of tries:     " << result.num_tries << endl;
  cout << "Stopping criteria:   " << result.nlopt_exit_criteria << endl;

  cout << endl;
  cout << "============= Result of last with_segment =============" << endl;
  cout << "Compute time:        " << result_with_segments.compute_time << endl;
  cout << "Function evaluations:" << result_with_segments.num_eval << endl;
  cout << "Trajectory time:     " << result_with_segments.x[result_with_segments.x.size - 1] << endl;
  cout << "Success:             " << result_with_segments.success << endl;
  cout << "False success:       " << result_with_segments.success_false << endl;
  cout << "Number of tries:     " << result_with_segments.num_tries << endl;
  cout << "Stopping criteria:   " << result_with_segments.nlopt_exit_criteria << endl;


  CHECK(std::abs(result.x.back() - result_with_segments.x.back()) < 1e-5);
};

} // namespace blast
