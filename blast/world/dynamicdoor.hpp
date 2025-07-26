#pragma once

#include <blast>

namespace blast {


inline blast_fn Box DynamicDoor::lookup(const real t) const {
    // progression from 0 to 1
    real progression = t < t0 ? 0 : (t > tf ? 1 : (t-t0) / (tf-t0));
    Box result;
    result.e = e;
    real current_angle = start_angle*(1-progression) + end_angle*progression;
    result.R = rpy2rotation(current_angle*axis);
    result.c = hinge + result.R*static_c_from_hinge;
    return result;
}

} // namespace blast