
#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include <blast>
#include "test_helper/test_helper.hpp"

namespace blast {

TEST_CASE("Test manipulator benchmark", "[blast]") {
  using namespace blast;

  // --- Test characteristics ---

  std::vector<real> gen_compute_time;
  std::vector<real> gen_trajectory_time;

  Array x0(37);
  fill_random(x0, PI);
  x0.back() = std::abs(x0.back()) * 5 + 0.1;

  Array gen_x(37);
  // BENCHMARK("Test manipulator benchmark") {
  auto opt_gen = get_generic_link6_opt();
  // opt_gen.guess.change_to_custom(x0);
  auto t1         = get_tick_us();
  auto result_gen = blast::optimize(&opt_gen, blast::OptimizationMethod::baseline);
  auto t2         = get_tick_us();
  // gen_compute_time.push_back(result_gen.compute_time);
  gen_compute_time.push_back((t2 - t1) / 1000);
  // gen_trajectory_time.push_back(result_gen.x[result_gen.x.size-1]);
  gen_trajectory_time.push_back(result_gen.x[result_gen.x.size - 1]);
  gen_x = result_gen.x;
  // };

  std::cout << "Success gen: " << result_gen.success << " / " << result_gen.success_false << std::endl;
  std::cout << "Compute time gen: " << mean(gen_compute_time) << " ms" << std::endl;
  std::cout << "Traj time gen: " << mean(gen_trajectory_time) << std::endl;

  std::vector<real> hc_compute_time;
  std::vector<real> hc_trajectory_time;
  Array             hc_x(37);
  // BENCHMARK("Hard-coded Link6 manipulator benchmark") {
  auto opt_hc = get_hardcoded_link6_opt();
  // opt_gen.guess.change_to_custom(x0);
  t1             = get_tick_us();
  auto result_hc = blast::optimize(&opt_hc, blast::OptimizationMethod::baseline);
  t2             = get_tick_us();
  // hc_compute_time.push_back(result_hc.compute_time);
  hc_compute_time.push_back((t2 - t1) / 1000);
  // hc_trajectory_time.push_back(result_hc.x[result_hc.x.size-1]);
  hc_trajectory_time.push_back(result_hc.x[result_hc.x.size - 1]);
  hc_x = result_hc.x;
  // };

  std::cout << "Success hc: " << result_hc.success << " / " << result_hc.success_false << std::endl;
  std::cout << "Compute time hc: " << mean(hc_compute_time) << std::endl;
  std::cout << "Traj time hc: " << mean(hc_trajectory_time) << std::endl;

  std::cout << "Resulting optimization vector : " << std::endl;
  print(gen_x);
  print(hc_x);
}
} // namespace blast
