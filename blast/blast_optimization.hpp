#pragma once
#include "blast.hpp"


namespace blast {



#ifdef __NVCC__
__global__ void test_kernel(cuPvaBspline pva);
#endif




#ifdef __NVCC__
__global__ void test_kernel(cuPvaBspline pva) {
    const u32 point = blockIdx.x * blockDim.x + threadIdx.x;

    blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;
    pva.compute_trajectory(point);

    const u32 pva_offset = point * 7;
    const u32 constraints_offset = point * 21;

    const auto pos = pva.device_pos + pva_offset;
    const auto vel = pva.device_vel + pva_offset;
    const auto acc = pva.device_acc + pva_offset;

    const auto con = manip->device_constraints + constraints_offset;
    manip->compute_constraints(pos, vel, acc, con);
}


#endif
} // namespace blast
