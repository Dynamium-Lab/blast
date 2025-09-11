//
// Created by nikos on 2025-09-09.
//
#include <blast>

using namespace blast;

int main() {
  int p      = 5;
  int points = 1000;
  int joints = 1;
  int n_ctrl = 10;

  Bspline spline(n_ctrl, points, p, joints);

  print_to_csv(transpose(spline.basis_p), "../../../../examples/bsplines/basis_function_p.csv");
  print_to_csv(transpose(spline.basis_v), "../../../../examples/bsplines/basis_function_v.csv");
  print_to_csv(transpose(spline.basis_a), "../../../../examples/bsplines/basis_function_a.csv");

  Matrix task(1, 6);
  task(0, 0) = 1;
  task(0, 3) = -1;
  Array x(spline.x_len(task));
  x        = fill_random(x, 1);
  x.back() = 0.5;
  spline.compute_trajectory(x, task);
  print_to_csv(transpose(spline.traj.pos), "../../../../examples/bsplines/bspline_p.csv");
  print_to_csv(transpose(spline.traj.vel), "../../../../examples/bsplines/bspline_v.csv");
  print_to_csv(transpose(spline.traj.acc), "../../../../examples/bsplines/bspline_a.csv");
}
