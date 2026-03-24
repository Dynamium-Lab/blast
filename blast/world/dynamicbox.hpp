#pragma once

#include <blast>

namespace blast {


inline host_fn Box DynamicBox::lookup(real time) const {
    Assert(trajectory.size() == n_points);

    real fraction = (time - start_time)/(end_time - start_time);
    fraction = clamp(fraction, 0, 1);

    const real increment = fraction * (n_points-1);

    const int increment_low = (int)floor(increment);
    const int increment_up = fraction == 1 ? increment_low : increment_low + 1;
    const real inc_rest = increment - (real)increment_low;

    const Vec3 c_low = trajectory[increment_low].center;
    const Mat3& R_low = trajectory[increment_low].rotation;
    const Vec3 e_low = trajectory[increment_low].extents;

    const Vec3 c_up = trajectory[increment_up].center;
    const Mat3& R_up = trajectory[increment_up].rotation;
    const Vec3 e_up = trajectory[increment_up].extents;

    Box box;
    box.center = inc_rest*c_up + (1-inc_rest)*c_low;
    box.rotation = inc_rest*R_up + (1-inc_rest)*R_low; // todo: this is wrong!
    box.extents = inc_rest*e_up + (1-inc_rest)*e_low;

    return box;
}


}

