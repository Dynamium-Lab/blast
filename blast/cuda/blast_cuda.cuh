#pragma once

#include <cuda_runtime.h>

#include "cuda/blast_cu_math.cuh"
#include "cuda/blast_cu_utils.cuh"
#include "cuda/blast_cu_manipulator.cuh"
#include "cuda/blast_cu_trajectory.cuh"

namespace blast {

__global__ void test_kernel(cuPvaBspline pva);
__global__ void no_object_kernel(
    real dt, real one_over_T, real one_over_T2, unsigned joints, unsigned ncontrol,
    real* basis_p, real* basis_v, real* basis_a, real* control, real* device_pos, real* device_vel, real* device_acc, real* device_t);




__global__
void test_kernel(cuPvaBspline pva) {
    const u32 point = blockIdx.x * blockDim.x + threadIdx.x;
    if (point > 255)
        printf("Error in kernel: point %d is not valid\n", point);

    blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;
    if (point == 0) {
        printf("Manip->I[5].data[5] = %f\n", manip->I[5].data[5]);
    }

    pva.compute_trajectory(point);

    const u32 pva_offset = point * 7;
    const u32 constraints_offset = point * 21;
    if (constraints_offset > 255*21)
        printf("Error in kernel: constraints offset %d is not valid\n", constraints_offset);

    const auto pos = pva.device_pos + pva_offset;
    const auto vel = pva.device_vel + pva_offset;
    const auto acc = pva.device_acc + pva_offset;

    const auto con = manip->device_constraints + constraints_offset;
    manip->compute_constraints(pos, vel, acc, con);



    // const unsigned idx = threadIdx.x;
    // const unsigned dof = pva.joints;
    // const unsigned nctrl = pva.ncontrol;
    // const auto bp = pva.device_basis_p + threadIdx.x*nctrl;
    // const auto bv = pva.device_basis_v + threadIdx.x*nctrl;
    // const auto ba = pva.device_basis_a + threadIdx.x*nctrl;

    // const double one_over_T = pva.one_over_T;
    // const double one_over_T2 = pva.one_over_T2;
    // double *control = pva.device_control;
    // for (u32 joint = 0; joint < dof; joint++) {
    //     double position = 0;
    //     double velocity = 0;
    //     double acceleration = 0;
    //     for (int i = 0; i < nctrl; i++) {
    //         position += control[dof*i + joint] * bp[i];
    //         velocity += control[dof*i + joint] * bv[i];
    //         acceleration += control[dof*i + joint] * ba[i];
    //     }
    //     pva.device_pos[idx*dof + joint] = position;
    //     pva.device_vel[idx*dof + joint] = velocity * one_over_T;
    //     pva.device_acc[idx*dof + joint] = acceleration * one_over_T2;
    // }

}

__global__
void no_object_kernel(
    real dt, real one_over_T, real one_over_T2, unsigned joints, unsigned ncontrol,
    real* basis_p, real* basis_v, real* basis_a, real* control, real* device_pos, real* device_vel, real* device_acc, real* device_t
) {
    const u32 point = threadIdx.x;
    auto bp = basis_p + point*ncontrol;
    auto bv = basis_v + point*ncontrol;
    auto ba = basis_a + point*ncontrol;
    for (int joint = 0; joint < joints; joint++) {
        auto c = control + joint*ncontrol;
        real position = 0;
        real velocity = 0;
        real acceleration = 0;
        for (int i = 0; i < ncontrol; i++) {
            position += c[i*joints + joint] * bp[i];
            velocity += c[i*joints + joint] * bv[i];
            acceleration += c[i*joints + joint] * ba[i];
        }
        device_pos[point*joints + joint] = position;
        device_vel[point*joints + joint] = velocity * one_over_T;
        device_acc[point*joints + joint] = acceleration * one_over_T2;
    }
    device_t[point] = dt * point;
}

}
