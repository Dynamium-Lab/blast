#pragma once

#if defined(_MSC_VER)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#undef NOMINMAX
#undef WIN32_LEAN_AND_MEAN
#else
#include <time.h>
#endif
#include <cmath>
#include <string>
#include "blast_math.hpp"

namespace blast {

blast_fn void print(Vec3 v) {
    printf("[% 0.4f, % 0.4f, % 0.4f]\n", v.x, v.y, v.z);
}

blast_fn void print(Mat3 m) {
    printf("\n[% 0.4f, % 0.4f, % 0.4f]\n[% 0.4f, % 0.4f, % 0.4f]\n[% 0.4f, % 0.4f, % 0.4f]\n",
           m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2), m(2, 0), m(2, 1), m(2, 2));
}

blast_fn void print(const Array& a) {
    using namespace std;
    if(a.size == 0)
        return;
    printf("[");
    for (u32 i = 0; i < a.size-1; i++)
        printf("% 0.4f, ", a[i]);
    printf("% 0.4f]\n", a[a.size-1]);
}

blast_fn void print(const Matrix& m) {
    if(m.size == 0)
        return;
    for (u32 i = 0; i < m.rows; i++) {
        printf("[");
        for (u32 j = 0; j < m.cols-1; j++)
            printf("% 0.4f, ", m(i, j));
        printf("% 0.4f]\n", m(i, m.cols-1));
    }
}

blast_fn void print(const MatrixXd &m) {
    if(m.size() == 0)
        return;
    for (u32 i = 0; i < m.rows(); i++) {
        printf("[");
        for (u32 j = 0; j < m.cols()-1; j++)
            printf("% 0.4f, ", m(i, j));
        printf("% 0.4f]\n", m(i, m.cols()-1));
    }
}

blast_fn void print(double* data, unsigned size) {
    printf("[");
    for (u32 i = 0; i < size-1; i++)
        printf("% 0.4f, ", data[i]);
    printf("% 0.4f]\n", data[size-1]);
}

blast_fn void print(float* data, unsigned size) {
    printf("[");
    for (u32 i = 0; i < size-1; i++)
        printf("% 0.4f, ", data[i]);
    printf("% 0.4f]\n", data[size-1]);
}

blast_fn void print(double* data, unsigned rows, unsigned cols) {
    for (u32 i = 0; i < rows; i++) {
        printf("[");
        for (u32 j = 0; j < cols-1; j++)
            printf("% 0.4f, ", data[i + rows*j]);
        printf("% 0.4f]\n", data[i + rows*(cols-1)]);
    }
}

blast_fn void print(float* data, unsigned rows, unsigned cols) {
    for (u32 i = 0; i < rows; i++) {
        printf("[");
        for (u32 j = 0; j < cols-1; j++)
            printf("% 0.4f, ", data[i + rows*j]);
        printf("% 0.4f]\n", data[i + rows*(cols-1)]);
    }
}

// get the time
host_fn int64_t get_tick_us() {
#if defined(_MSC_VER)
    LARGE_INTEGER start, frequency;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start);
    return (start.QuadPart * 1000000) / frequency.QuadPart;
#else
    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);
    return (start.tv_sec * 1000000LLU) + (start.tv_nsec / 1000);
#endif
}

// note: CUDA stuff, only enabled if compiling for Nvidia GPUs
#ifdef __NVCC__

// print the properties of all connected GPUs
host_fn void print_device_properties() {
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: ............................... %d\n", i);
        printf("  Device name: ............................... %s\n", prop.name);
        printf("  Multiprocessors: ........................... %d\n", prop.multiProcessorCount);
        printf("  Threads per multiprocessor: ................ %d\n", prop.maxThreadsPerMultiProcessor);
        printf("  Threads per block: ......................... %d\n", prop.maxThreadsPerBlock);
        printf("  Async Engine Count: ........................ %d\n", prop.asyncEngineCount);
        printf("  Registers per block: ....................... %d\n", prop.regsPerBlock);
        printf("  Registers per multiprocessor: .............. %d\n", prop.regsPerMultiprocessor);
        printf("  Shared memory per block (KB): .............. %f\n", prop.sharedMemPerBlock/1024.0);
        printf("  Shared memory per multiprocessor (KB): ..... %f\n", prop.sharedMemPerMultiprocessor/1024.0);
        printf("  Concurrent Kernels: ........................ %d\n", prop.concurrentKernels);
        printf("  Clock Rate (KHz): .......................... %d\n", prop.clockRate);
        printf("  Memory Clock Rate (KHz): ................... %d\n", prop.memoryClockRate);
        printf("  Memory Bus Width (bits): ................... %d\n", prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GB/s): .............. %f\n", 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
        printf("  Compute capabilities : ..................... %d,%d\n", prop.major, prop.minor);
    }
}
#endif

} // namespace blast
