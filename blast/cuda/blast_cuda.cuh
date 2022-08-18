#pragma once

#include <cuda_runtime.h>

#include "cuda/blast_cu_math.cuh"
#include "cuda/blast_cu_utils.cuh"
#include "cuda/blast_cu_manipulator.cuh"
#include "cuda/blast_cu_trajectory.cuh"


__global__ void test_kernal();




__global__ void test_kernal() {
    blast::cuGen3_7DOF* manip = (blast::cuGen3_7DOF*)blast::manip_broadcast_arena;
}
