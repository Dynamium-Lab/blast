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

TEST_CASE("Benchmark accelerations", "[accelerations]") {

  ManipulatorTempData manip_data;

  // Optimization opt(get_generic_gen3(), get_gen3_task());
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 100, 5, 6);
  opt.bspline = bspline;

  opt.guess.type = Guess::custom;
  opt.guess.initial_x   = Array(opt.bspline.x_len(opt.task), 2.0);
  // opt.guess.n_random_shots = 100;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = false;
  opt.constraints.self_collisions     = false;
  opt.constraints.external_collisions = false;

  opt.max_tries         = 1;
  opt.success_tolerance = 0.01;

  Result result(&opt);
  Result result_dev(&opt);
  Result result_dev_new(&opt);

  BENCHMARK_ADVANCED("optimize")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      IOSilencer _;
      result = optimize(&opt, OptimizationMethod::baseline);
      return result;
    });
  };

  BENCHMARK_ADVANCED("optimize with_analytical_pva")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      IOSilencer _;
      result_dev = optimize(&opt, OptimizationMethod::with_analytical_pva);
      return result_dev;
    });
  };

  // BENCHMARK_ADVANCED("optimize with_analytical_dynamics")(Catch::Benchmark::Chronometer meter) {
  //   meter.measure([&] {
  //     IOSilencer _;
  //     result_dev_new = optimize(&opt, OptimizationMethod::with_analytical_dynamics);
  //     return result_dev_new;
  //   });
  // };

  CHECK(std::abs(result.x.back() - result_dev.x.back()) < 1e-5);
  // CHECK(std::abs(result.x.back() - result_dev_new.x.back()) < 1e-5);
}

} // namespace blast
