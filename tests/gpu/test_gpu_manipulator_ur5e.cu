// GPU UR5e manipulator correctness test.
//
// cuUR5eMultiTraj checks position, velocity, and acceleration (not torque
// -- UR5e has a real, finite acceleration limit, unlike Gen3 which has
// none, so acceleration is the natural third constraint here and needs no
// Newton-Euler recursion at all). This test cross-checks its constraint
// formula against the Manipulator object's own position_max/velocity_max/
// acceleration_max, to catch a transcription error in
// cuUR5eMultiTraj::_init()'s hardcoded limits.

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "blast/gpu/manipulator.cuh"

TEST_CASE("GpuCpuUR5eManipCorrectness", "[Manipulator][UR5e]") {
  using namespace blast;

  Manipulator ur5e = make_UR5e();
  REQUIRE(ur5e.n_joints == 6);

  cuUR5eMultiTraj gpu_manip(1, 1, 0);

  const int n_trials = 30;
  for (int trial = 0; trial < n_trials; trial++) {
    real pos[6], vel[6], acc[6];
    for (int j = 0; j < 6; j++) {
      // Random values spanning well past each limit (both directions), so
      // this exercises both the "inside limit" (negative con) and "outside
      // limit" (positive con) branches, not just one.
      pos[j] = 1.5 * (real) ur5e.position_max[j] * (2.0 * random_real() - 1.0);
      vel[j] = 1.5 * (real) ur5e.velocity_max[j] * (2.0 * random_real() - 1.0);
      acc[j] = 1.5 * (real) ur5e.acceleration_max[j] * (2.0 * random_real() - 1.0);
    }

    real con[18];
    gpu_manip.compute_constraints_point(pos, vel, acc, con);

    for (int j = 0; j < 6; j++) {
      CAPTURE(trial, j, pos[j], vel[j], acc[j]);

      real expected_pos_con = (std::abs(pos[j]) - (real) ur5e.position_max[j]) / (real) ur5e.position_max[j];
      REQUIRE(std::abs(con[j] - expected_pos_con) < 1e-6);

      real expected_vel_con = (std::abs(vel[j]) - (real) ur5e.velocity_max[j]) / (real) ur5e.velocity_max[j];
      REQUIRE(std::abs(con[6 + j] - expected_vel_con) < 1e-6);

      real expected_acc_con = (std::abs(acc[j]) - (real) ur5e.acceleration_max[j]) / (real) ur5e.acceleration_max[j];
      REQUIRE(std::abs(con[12 + j] - expected_acc_con) < 1e-6);
    }
  }
}
