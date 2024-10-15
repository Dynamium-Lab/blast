#pragma once
#include "gpu/manipulator.hpp"
#include "gpu/trajectory.hpp"


//------ KERNEL FUNCTIONS ------------------------------------------------------------------------------------

__global__ void pva_constraints_kernel(cuBspline pva) {
    const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
    pva.compute_trajectory(point);
    blast::cuGen3MultiTraj* manip = (blast::cuGen3MultiTraj*)blast::manip_broadcast_arena;
    const u32 pva_offset = point * 7;
    const u32 constraints_offset = point * 21;
    manip->compute_constraints_point(
        pva.device_pos+pva_offset,
        pva.device_vel+pva_offset,
        pva.device_acc+pva_offset,
        manip->dev_constraints + constraints_offset
    );
}

__global__ void constraints_no_pva_kernel(cuBspline pva) {
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
    blast::cuGen3MultiTraj* manip = (blast::cuGen3MultiTraj*)blast::manip_broadcast_arena;
    const u32 constraints_offset = point * 21;
    manip->compute_constraints_point(pos, vel, acc, manip->dev_constraints + constraints_offset);
}

// kernel that uses shared memory to speed up constraint computation
// todo: WIP not completed
__global__ void constraints_shared_kernel(cuBspline pva) {
    const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
    const u32 nctrl = pva.ncontrol;

    // __shared__ real bp[];
    // __shared__ real bv[];
    // __shared__ real ba[];
    // __shared__ real control[12*7];
    // if (threadIdx.x == 0)
    //     cudaMemcpyAsync(control, pva.device_control, sizeof(real) * 12 * 7, cudaMemcpyDeviceToDevice);

    // __syncthreads();

    // cudaDeviceSynchronize();

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
        vel[joint] = velocity     * pva.one_over_T;
        acc[joint] = acceleration * pva.one_over_T2;
    }

    blast::cuGen3MultiTraj* manip = (blast::cuGen3MultiTraj*)blast::manip_broadcast_arena;
    const u32 constraints_offset = point * 21;
    manip->compute_constraints_point(pos, vel, acc, manip->dev_constraints + constraints_offset);
}

