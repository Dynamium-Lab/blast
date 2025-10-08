#pragma once

#include <blast>

namespace blast {


inline host_fn Box DynamicBox::lookup(real time) const {
    Assert(trajectory.size() == n_pts);

    real fraction = (time - t0)/(tf - t0);
    fraction = clamp(fraction, 0, 1);

    const real increment = fraction * (n_pts-1);

    const int increment_low = (int)floor(increment);
    const int increment_up = fraction == 1 ? increment_low : increment_low + 1;
    const real inc_rest = increment - (real)increment_low;

    const Vec3 c_low = trajectory[increment_low].c;
    const Mat3& R_low = trajectory[increment_low].R;
    const Vec3 e_low = trajectory[increment_low].e;

    const Vec3 c_up = trajectory[increment_up].c;
    const Mat3& R_up = trajectory[increment_up].R;
    const Vec3 e_up = trajectory[increment_up].e;

    Box box;
    box.c = inc_rest*c_up + (1-inc_rest)*c_low;
    box.R = inc_rest*R_up + (1-inc_rest)*R_low; // todo: this is wrong!
    box.e = inc_rest*e_up + (1-inc_rest)*e_low;

    return box;
}


}

