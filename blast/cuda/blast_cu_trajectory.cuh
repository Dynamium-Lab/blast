#pragma once

#include "cuda/blast_cuda.cuh"


namespace blast {


__device__
struct cuPvaBspline {
    unsigned dof;
    unsigned npoints;
    unsigned p;
    unsigned ncontrol;


    // compute the basis functions
    __host__ void init();

    // map to GPU using cudaMemcpy
    __host__ void map_to_gpu(double* arena);

};


} // namespace blast
