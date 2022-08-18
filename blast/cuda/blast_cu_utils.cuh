#pragma once

#include "cuda/blast_cuda.cuh"

#include <cstdio>

namespace blast {

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

} // namespace blast
