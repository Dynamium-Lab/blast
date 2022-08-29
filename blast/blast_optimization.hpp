#pragma once
#include "blast.hpp"


namespace blast {



#ifdef __NVCC__

// compute the trajectory and the constraints (slower, but access to trajectory)
__global__ void pva_constraints_kernel(cuPvaBspline pva);

// compute the constraints, but don't store the trajectory (faster, but no access to trajectory)
__global__ void constraints_no_pva_kernel(cuPvaBspline pva);
#endif




#ifdef __NVCC__


__global__ void pva_constraints_kernel(cuPvaBspline pva) {
    const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
    pva.compute_trajectory(point);

    blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;
    const u32 pva_offset = point * 7;
    const u32 constraints_offset = point * 21;
    manip->compute_constraints(
        pva.device_pos+pva_offset,
        pva.device_vel+pva_offset,
        pva.device_acc+pva_offset,
        manip->device_constraints + constraints_offset
    );
}

__global__ void constraints_no_pva_kernel(cuPvaBspline pva) {
    const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
    const u32 nctrl = pva.ncontrol;
    const real one_over_T = pva.one_over_T;
    const real one_over_T2 = pva.one_over_T2;

    if (point == 0)
        pva.device_t[point] = pva.dt * point;
    __syncthreads();

    real pos[7];
    real vel[7];
    real acc[7];
    const auto bp = &pva.device_basis_p[point*nctrl];
    const auto bv = &pva.device_basis_v[point*nctrl];
    const auto ba = &pva.device_basis_a[point*nctrl];
    for (int joint = 0; joint < 7; joint++) {
        const auto c = pva.device_control + joint*nctrl;
        real position = 0;
        real velocity = 0;
        real acceleration = 0;
        for (int i = 0; i < nctrl; i++) {
            position     += c[i] * bp[i];
            velocity     += c[i] * bv[i];
            acceleration += c[i] * ba[i];
        }
        pos[joint] = position;
        vel[joint] = velocity     * one_over_T;
        acc[joint] = acceleration * one_over_T2;
    }

    blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;
    const u32 constraints_offset = point * 21;
    manip->compute_constraints(pos, vel, acc, manip->device_constraints + constraints_offset);
}





#endif
} // namespace blast
