#pragma once

#include <blast>

namespace blast {


// Rotation a fraction t of the way from r_from to r_to along the shortest arc:
// r_from * exp(t * log(r_from^T r_to)), by Rodrigues' formula. A linear blend of the two matrices
// is not a rotation (it shrinks and shears the box), which made the box's bounding volume
// disagree with its distance and the broadphase miss the closest pair.
inline host_fn Mat3 interpolate_rotation(const Mat3& r_from, const Mat3& r_to, real t) {
  const Mat3 relative  = transpose(r_from) * r_to;
  const real cos_angle = clamp((relative(0, 0) + relative(1, 1) + relative(2, 2) - 1) / 2, -1, 1);
  const real angle     = std::acos(cos_angle);
  Vec3       axis      = {relative(2, 1) - relative(1, 2), relative(0, 2) - relative(2, 0), relative(1, 0) - relative(0, 1)};
  const real axis_norm = norm(axis); // 2 sin(angle)
  if (axis_norm < 1e-9) // no rotation between the keyframes (a half turn, axis ambiguous, is not supported)
    return r_from;
  axis = axis / axis_norm;
  Mat3 skew;
  skew(0, 1) = -axis.z, skew(0, 2) = axis.y;
  skew(1, 0) = axis.z, skew(1, 2) = -axis.x;
  skew(2, 0) = -axis.y, skew(2, 1) = axis.x;
  const real step = t * angle;
  return r_from * (eye() + std::sin(step) * skew + (1 - std::cos(step)) * (skew * skew));
}

inline host_fn Box DynamicBox::lookup(real time) const {
  Assert(trajectory.size() == n_points);

  real fraction = (time - start_time) / (end_time - start_time);
  fraction      = clamp(fraction, 0, 1);

  const real increment = fraction * (n_points - 1);

  const int  increment_low = (int) floor(increment);
  const int  increment_up  = fraction == 1 ? increment_low : increment_low + 1;
  const real inc_rest      = increment - (real) increment_low;

  const Vec3  c_low = trajectory[increment_low].center;
  const Mat3& R_low = trajectory[increment_low].rotation;
  const Vec3  e_low = trajectory[increment_low].extents;

  const Vec3  c_up = trajectory[increment_up].center;
  const Mat3& R_up = trajectory[increment_up].rotation;
  const Vec3  e_up = trajectory[increment_up].extents;

  Box box;
  box.center   = inc_rest * c_up + (1 - inc_rest) * c_low;
  box.rotation = interpolate_rotation(R_low, R_up, inc_rest);
  box.extents  = inc_rest * e_up + (1 - inc_rest) * e_low;

  return box;
}


} // namespace blast
