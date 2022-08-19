#pragma once

#include <cuda_runtime.h>

#include "cuda/blast_cu_math.cuh"
#include "cuda/blast_cu_utils.cuh"
#include "cuda/blast_cu_manipulator.cuh"
#include "cuda/blast_cu_trajectory.cuh"

namespace blast {

__global__ void test_kernal(cuPvaBspline pva);




__global__ void test_kernal(cuPvaBspline pva) {
    blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;
    pva.compute_trajectory(threadIdx.x);
}

}
