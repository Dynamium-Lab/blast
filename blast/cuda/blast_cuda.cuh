#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#ifdef __INTELLISENSE__
void __syncthreads();
#endif

#include "cuda/blast_cu_math.cuh"
#include "cuda/blast_cu_utils.cuh"
#include "cuda/blast_cu_manipulator.cuh"
#include "cuda/blast_cu_trajectory.cuh"


__global__ void test_kernal();




__global__ void test_kernal() {
    blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;
    if (blockIdx.x == 0) {
        printf("%f\n", manip->pmax[0]);
        printf("%f\n", manip->pmax[1]);
        printf("%f\n", manip->pmax[2]);
        printf("%f\n", manip->pmax[3]);
        printf("%f\n", manip->pmax[4]);
        printf("%f\n", manip->pmax[5]);
        printf("%f\n", manip->pmax[6]);
    }
    __syncthreads();
}
