#include <blast>
#include "test_helper/test_helper.hpp"

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"


namespace blast {

TEST_CASE("Manipulator benchmarks", "[manipulator]") {

  BENCHMARK_ADVANCED("Dynamics on a Gen3")(Catch::Benchmark::Chronometer meter) {

    const u32 n_points = 256;
    const u32 n_joints = 7;
    const u32 p = 5;
    const u32 ncontrol = 24;
    auto manip = get_generic_gen3();
    Bspline bspline(ncontrol, n_points, p, n_joints);
    Matrix efforts(n_joints, n_points);

    // random task
    real amp = 10;
    Matrix task(n_joints, 6);
    for (u32 i = 0; i < task.rows; i++)
    {
      for (u32 j = 0; j < task.cols; j++)
      {
        task(i, j) = amp * get_random();
      }
    }


    // random optimization vector
    Array x(bspline.x_len(task));
    for (u32 i = 0; i < x.size; i++)
      x[i] = amp * get_random();
    x[x.size - 1] = 3.f;
    bspline.compute_trajectory(x, task);

    meter.measure([&] {
      for (int i = 0; i < n_points; ++i) {
        auto q   = bspline.traj.pos.col(i);
        auto qd  = bspline.traj.vel.col(i);
        auto qdd = bspline.traj.acc.col(i);

        forward_kinematics(manip, q);
        dynamics(manip, qd, qdd);
      }
      return manip._efforts;
    });
  };
}

} // namespace blast
