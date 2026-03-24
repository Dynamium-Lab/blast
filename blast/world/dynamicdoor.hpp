#pragma once

#include <blast>

namespace blast {


inline blast_fn Box DynamicDoor::lookup(const real time) const {
    // progression from 0 to 1
    real progression = time < start_time ? 0 : (time > end_time ? 1 : (time-start_time) / (end_time-start_time));
    Box result;
    result.extents = extents;
    real current_angle = start_angle*(1-progression) + end_angle*progression;
    result.rotation = rpy2rotation(current_angle*axis);
    result.center = hinge + result.rotation*static_c_from_hinge;
    return result;
}

} // namespace blast