#pragma once

#include <blast>

namespace blast {

inline host_fn Capsule DynamicCapsule::lookup(const blast::real time) const {
  Assert(trajectory.size() == n_points);

  blast::real fraction = end_time == start_time ? 1 : (time - start_time) / (end_time - start_time);
  fraction             = blast::clamp(fraction, 0, 1);

  blast::real increment = fraction * (n_points - 1);

  int         increment_low = (int) floor(increment);
  int         increment_up  = fraction == 1 ? increment_low : increment_low + 1;
  blast::real inc_rest      = increment - (blast::real) increment_low;

  Vec3        p1_low = trajectory[increment_low].p1;
  Vec3        p2_low = trajectory[increment_low].p2;
  blast::real r_low  = trajectory[increment_low].radius;

  Vec3        p1_up = trajectory[increment_up].p1;
  Vec3        p2_up = trajectory[increment_up].p2;
  blast::real r_up  = trajectory[increment_up].radius;

  blast::Capsule capsule;
  capsule.p1     = inc_rest * p1_up + (1 - inc_rest) * p1_low;
  capsule.p2     = inc_rest * p2_up + (1 - inc_rest) * p2_low;
  capsule.radius = inc_rest * r_up + (1 - inc_rest) * r_low;
  ;

  return capsule;
}

} // namespace blast
