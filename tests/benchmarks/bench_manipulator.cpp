#include <blast>
#include "test_helper/test_helper.hpp"

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"


namespace blast {

TEST_CASE("Manipulator benchmarks", "[manipulator]") {

  BENCHMARK_ADVANCED("Dynamics on a Gen3")(Catch::Benchmark::Chronometer meter) {
    constexpr int n_joints  = 7;
    constexpr int n_points  = 1000;

    auto manip = get_generic_gen3();

    // Pre-generate random joint positions, velocities, and accelerations (not timed)
    Trajectory traj(n_points, n_joints);
    fill_random(traj.pos, PI);
    fill_random(traj.vel, 2.0);
    fill_random(traj.acc, 5.0);

    // Optional warm-up (not timed)
    for (int i = 0; i < 10; ++i) {
      auto q   = traj.pos.col(i);
      auto qd  = traj.vel.col(i);
      auto qdd = traj.acc.col(i);

      forward_kinematics(manip, q);
      dynamics(manip, qd, qdd);
    }

    meter.measure([&] {
      // Accumulate a checksum to keep computations observable
      real checksum = 0.0;
      for (int i = 0; i < n_points; ++i) {
        auto q   = traj.pos.col(i);
        auto qd  = traj.vel.col(i);
        auto qdd = traj.acc.col(i);

        forward_kinematics(manip, q);
        dynamics(manip, qd, qdd);
        for (int j = 0; j < n_joints; ++j) {
          checksum += manip._efforts[j];
        }
      }
      return checksum;
    });
  };
}

} // namespace blast
