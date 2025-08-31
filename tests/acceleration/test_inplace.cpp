#include <blast>
#include "test_helper/test_helper.hpp"

#define BLAST_TRACE_LEVEL 3

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"

namespace blast {


TEST_CASE("Replace INF_REAL with huge value", "[acceleration]") {
  int n_iterations = 10000;

  // both pmax and pmin are +/-INF_REAL
  real pmax = INF_REAL;
  real pmin = -INF_REAL;

  pmax = 1e300;
  pmin = -1e300;

  real c          = (pmax + pmin) / 2;
  real b          = (pmax - pmin) / 2;
  real one_over_b = 1 / b;

  Array p(n_iterations);
  p = fill_random(p, 360);

  for (int i = 0; i < n_iterations; i++) {
    real q_prime    = p[i] - c;
    real constraint = (std::abs(q_prime) - b) / b;
    real gradient   = sign(q_prime) * 1 / b;

    CHECK(constraint == -1.0);
    CHECK(std::abs(gradient) <= 1e-300);
  }
};

TEST_CASE("Analytic dv/dT & da/dT tests", "[acceleration]") {
  int n_iterations = 1000;

  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 100, 5, 6);
  opt.bspline = bspline;

  opt.guess.type = Guess::custom;
  opt.guess.x0   = Array(opt.bspline.x_len(opt.task), 2.0);
  // opt.guess.n_shot = 100;

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

  auto xlen = opt.bspline.x_len(opt.task);

  for (int iter = 0; iter < n_iterations; iter++) {

    initialize_optimization(&opt);
    n_con(&opt);

    auto x = init_guess(&opt);
    bspline.compute_trajectory(x, opt.task);

    Array result(opt.constraints.n_constraints);
    compute_constraints(result.data, x, &opt);

    Array x_plus(xlen);
    x_plus = x;
    Array grad(opt.constraints.n_constraints);
    Array grad_analytic(opt.constraints.n_constraints);

    // last point T
    u32 j = xlen - 1;
    // memcpy(x_plus.data, x, xlen * sizeof(real)); todo: check if necessary
    x_plus[j] += eps;

    // compute constraints
    Array r_plus_T(opt.constraints.n_constraints);
    compute_constraints(r_plus_T.data, x_plus, &opt);

    for (u32 i = 0; i < opt.constraints.n_constraints; i++) {
      grad[i] = (r_plus_T[i] - result[i]) / eps;
    }

    int fill_idx = 24;
    for (int point = 1; point < opt.bspline.n_points; point++) {
      // Position is all zeros
      fill_idx += opt.manip.n_joints;
      // Velocity
      for (int joint = 0; joint < opt.manip.n_joints; joint++) {
        grad_analytic[fill_idx] = -(result[fill_idx] + 1) / x[xlen - 1];
        CHECK(std::abs(grad[fill_idx] - grad_analytic[fill_idx]) < 1e-4);
        fill_idx++;
      }
      // Acceleration
      for (int joint = 0; joint < opt.manip.n_joints; joint++) {
        grad_analytic[fill_idx] = -2.0 * (result[fill_idx] + 1) / x[xlen - 1];
        CHECK(std::abs(grad[fill_idx] - grad_analytic[fill_idx]) < 1e-4);
        fill_idx++;
      }

      for (int joint = 0; joint < opt.manip.n_joints; joint++) {
        fill_idx++;
      }
    }
  }
};

} // namespace blast
