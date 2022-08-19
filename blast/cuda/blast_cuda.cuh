#pragma once

#include <cuda_runtime.h>

#include "cuda/blast_cu_math.cuh"
#include "cuda/blast_cu_utils.cuh"
#include "cuda/blast_cu_manipulator.cuh"
#include "cuda/blast_cu_trajectory.cuh"

namespace blast {

__global__ void test_kernal(cuPvaBspline pva);




__global__ void test_kernal(cuPvaBspline pva) {
    // blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;

    if (threadIdx.x == 0) {
        // printf("Device control points:\n");
        // cuPrint(pva.device_control, 8, 7);

        printf("Device dt:%f, 1/T: %f, 1/T^2: %f\n", pva.dt, pva.one_over_T, pva.one_over_T2);
    }

    __syncthreads();

    pva.compute_trajectory(threadIdx.x);

    __syncthreads();
}

}
