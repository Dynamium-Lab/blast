#pragma once

#include <blast>
#include "blast_trajectory.hpp"

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
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include "json.hpp"

namespace blast {

inline blast_fn void print(Vec3 v) {
    printf("[% 0.4f, % 0.4f, % 0.4f]\n", v.x, v.y, v.z);
}

inline blast_fn void print(Mat3 m) {
    printf("\n[% 0.4f, % 0.4f, % 0.4f]\n[% 0.4f, % 0.4f, % 0.4f]\n[% 0.4f, % 0.4f, % 0.4f]\n",
           m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2), m(2, 0), m(2, 1), m(2, 2));
}

inline blast_fn void print(const Array& a) {
    using namespace std;
    if(a.size == 0)
        return;
    printf("[");
    for (u32 i = 0; i < a.size-1; i++)
        printf("% 0.4f, ", a[i]);
    printf("% 0.4f]\n", a[a.size-1]);
}

inline blast_fn void print(const Matrix& m) {
    if(m.size == 0)
        return;
    for (u32 i = 0; i < m.rows; i++) {
        printf("[");
        for (u32 j = 0; j < m.cols-1; j++)
            printf("% 0.4f, ", m(i, j));
        printf("% 0.4f]\n", m(i, m.cols-1));
    }
}

inline blast_fn void print(const double* data, unsigned size) {
    printf("[");
    for (u32 i = 0; i < size-1; i++)
        printf("% 0.4f, ", data[i]);
    printf("% 0.4f]\n", data[size-1]);
}

inline blast_fn void print(const float* data, unsigned size) {
    printf("[");
    for (u32 i = 0; i < size-1; i++)
        printf("% 0.4f, ", data[i]);
    printf("% 0.4f]\n", data[size-1]);
}

inline blast_fn void print(const double* data, unsigned rows, unsigned cols) {
    for (u32 i = 0; i < rows; i++) {
        printf("[");
        for (u32 j = 0; j < cols-1; j++)
            printf("% 0.4f, ", data[i + rows*j]);
        printf("% 0.4f]\n", data[i + rows*(cols-1)]);
    }
}

inline blast_fn void print(const float* data, unsigned rows, unsigned cols) {
    for (u32 i = 0; i < rows; i++) {
        printf("[");
        for (u32 j = 0; j < cols-1; j++)
            printf("% 0.4f, ", data[i + rows*j]);
        printf("% 0.4f]\n", data[i + rows*(cols-1)]);
    }
}

inline host_fn void print_to_csv(const Matrix& m, const std::string &filename) {
    std::ofstream file;
    file.open(filename);
    for (u32 i = 0; i < m.rows; i++) {
        for (u32 j = 0; j < m.cols-1; j++)
            file << m(i, j) << ",";
        file << m(i, m.cols-1) << std::endl;
    }
    file.close();
}

inline host_fn void print_to_csv(const Trajectory& traj, const std::string &filename) {
    std::ofstream file;
    file.open(filename);

    // print header
    {
        file << "t";
        for (u32 i = 0; i < traj.pos.rows; i++)
            file << ",p[" << i << "]";
        for (u32 i = 0; i < traj.vel.rows; i++)
            file << ",v[" << i << "]";
        for (u32 i = 0; i < traj.acc.rows; i++)
            file << ",a[" << i << "]";
        file << std::endl;
    }

    // print data
    for (u32 i = 0; i < traj.t.size; i++) {
        file << traj.t[i];
        for (u32 j = 0; j < traj.pos.rows; j++)
            file << "," << traj.pos(j, i);
        for (u32 j = 0; j < traj.vel.rows; j++)
            file << "," << traj.vel(j, i);
        for (u32 j = 0; j < traj.acc.rows; j++)
            file << "," << traj.acc(j, i);
        file << std::endl;
    }
    file.close();
}

inline host_fn int64_t get_tick_us() {
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

inline host_fn blast::Matrix read_csv_matrix(const std::string& filename, const char* csv_sep = ",") {
    std::ifstream file(filename);
    std::cout << "Reading from file: " << filename << std::endl;
    if (!file.is_open()) {
        Assert(false);
        std::cerr << "ERROR: File is unavailable" << std::endl;
        return blast::Matrix(0, 0);
    }

    std::string line;
    u32 num_rows = 0, num_cols = 0;

    std::getline(file, line); //discard first line
    // Determine number of rows and columns
    if (std::getline(file, line)) {
        num_cols = std::count(line.begin(), line.end(), *csv_sep) + 1;
        num_rows++;
    }
    while (std::getline(file, line)) {
        num_rows++;
    }

    // Rewind file to read again
    file.clear();
    file.seekg(0);

    // Create matrix with detected dimensions
    blast::Matrix pos(num_rows, num_cols);

    std::getline(file, line); //discard first line
    // Read the data into the matrix
    for (u32 i = 0; i < num_rows; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        for (u32 j = 0; j < num_cols; j++) {
            std::string value;
            if (std::getline(iss, value, *csv_sep)) {
                pos(i, j) = std::stod(value);
            }
        }
    }
    return pos;
}

inline host_fn blast::Matrix read_csv_matrix_no_header(const std::string& filename, const char* csv_sep = ",") {
    std::ifstream file(filename);
    std::cout << "Reading from file: " << filename << std::endl;
    if (!file.is_open()) {
        Assert(false);
        std::cerr << "ERROR: File is unavailable" << std::endl;
        return blast::Matrix(0, 0);
    }

    std::string line;
    u32 num_rows = 0, num_cols = 0;

    // Determine number of rows and columns
    if (std::getline(file, line)) {
        num_cols = std::count(line.begin(), line.end(), *csv_sep) + 1;
        num_rows++;
    }
    while (std::getline(file, line)) {
        num_rows++;
    }

    // Rewind file to read again
    file.clear();
    file.seekg(0);

    // Create matrix with detected dimensions
    blast::Matrix pos(num_rows, num_cols);

    // Read the data into the matrix
    for (u32 i = 0; i < num_rows; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        for (u32 j = 0; j < num_cols; j++) {
            std::string value;
            if (std::getline(iss, value, *csv_sep)) {
                pos(i, j) = std::stod(value);
            }
        }
    }

    return pos;
}

// note: CUDA stuff, only enabled if compiling for Nvidia GPUs
#ifdef __NVCC__
// print the properties of all connected GPUs
inline host_fn void print_device_properties() {
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