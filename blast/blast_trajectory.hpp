#pragma once

#include <blast>

namespace blast {

struct Trajectory;
struct Bspline;
struct Task;

struct Trajectory {
  Matrix pos; // n_joints x n_points
  Matrix vel; // n_joints x n_points
  Matrix acc; // n_joints x n_points
  Array  t;   // n_points

  Trajectory(u32 n_points, u32 n_joints) :
      pos(n_joints, n_points),
      vel(n_joints, n_points),
      acc(n_joints, n_points),
      t(n_points) {}

  Trajectory() = default;
};
Trajectory compute_5order_trajectory(real T, Matrix& task);

struct Bspline {
  Trajectory traj;
  Matrix     control; // n_ctrl x n_joints
  Matrix     basis_p; // n_ctrl x n_points
  Matrix     basis_v; // n_ctrl x n_points
  Matrix     basis_a; // n_ctrl x n_points
  std::vector<Matrix> basis; // n_ctrl x n_points, storing for derivatives up to desired "d"
  u32        n_joints;
  u32        n_points;
  u32        n_ctrl;
  u32        degree;

  std::vector<u32> lower_bounds; // nctrl
  std::vector<u32> upper_bounds; // nctrl

  inline blast_fn Bspline() = delete;

  inline blast_fn Bspline(u32 n_control, u32 n_points, u32 degree, u32 n_joints);

  inline blast_fn explicit Bspline(u32 n_joints) :
      Bspline(12, 100, 5, n_joints) {}

  inline blast_fn Bspline(const Bspline& other) = default;

  inline blast_fn Bspline(Bspline&& other) = default;

  inline blast_fn Bspline& operator=(const Bspline& other) = default;

  inline blast_fn Bspline& operator=(Bspline&& other) = default;

  // Compute a trajectory from the given optimization vector
  //  - note: fastest when 'n_control' is a multiple of 4 (SIMD)
  inline blast_fn void compute_trajectory(const Array& x, const Matrix& task);
  inline blast_fn u32  x_len(const Matrix& task) const;

  inline blast_fn void compute_basis();
  inline blast_fn void compute_basis_derivative(int d);
  inline blast_fn void compute_basis_open();
  inline blast_fn void compute_control(const Array& x, const Matrix& task, real* dst) const;
};

} // namespace blast

#include "trajectory/bspline.hpp"
#include "trajectory/polynomial_5order.hpp"
