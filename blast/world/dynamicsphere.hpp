#pragma once

#include <blast>

namespace blast {

inline host_fn Sphere DynamicSphere::lookup(const blast::real time) const {
    Assert(trajectory.size() == n_points);

    blast::real fraction = (time - start_time)/(end_time - start_time);
    fraction = blast::clamp(fraction, 0, 1);

    blast::real increment = fraction * (n_points-1);

    int increment_low = (int)floor(increment);
    int increment_up = fraction == 1 ? increment_low : increment_low + 1;
    blast::real inc_rest = increment - (blast::real)increment_low;

    Vec3 c_low = trajectory[increment_low].center;
    Vec3 c_up = trajectory[increment_up].center;

    blast::real r_low = trajectory[increment_low].radius;
    blast::real r_up = trajectory[increment_up].radius;

    blast::Sphere sphere;
    sphere.center = inc_rest*c_up + (1-inc_rest)*c_low;
    sphere.radius = inc_rest*r_up + (1-inc_rest)*r_low;

    return sphere;
}


} // namespace blast