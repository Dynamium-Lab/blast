#pragma once

#include "cuda/blast_cuda.cuh"
#include "blast_math.hpp"

#include <cstdio>

#define cuda_check(expr) \
{ \
    auto code = (expr); \
    if (code != cudaSuccess) \
        fprintf(stderr, "Cuda error: %s. In file: %s (%d)", cudaGetErrorString(code), __FILE__,__LINE__);\
}

namespace blast {

__host__
inline void print_device_properties() {
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", prop.name);
        printf("  Multiprocessors: %d\n", prop.multiProcessorCount);
        printf("  Threads per multiprocessor: %d\n", prop.maxThreadsPerMultiProcessor);
        printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
        printf("  Clock Rate (KHz): %d\n", prop.clockRate);
        printf("  Concurrent Kernels: %d\n", prop.concurrentKernels);
        printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GB/s): %f\n\n", 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    }
}



__device__
inline void cuPrint(double* data, unsigned rows, unsigned cols) {

    for (u32 i = 0; i < rows; i++) {
        printf("[");
        for (u32 j = 0; j < cols-1; j++)
            printf("%0.4f, ", data[i + rows*j]);
        printf("%0.4f]\n", data[i + rows*(cols-1)]);
    }
}

} // namespace blast
