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
#include "blast.h"
#include <fstream>

namespace blast {

void print(Vec3 v);

void print(Mat3 m);

void print(const Array& a);

void print(const Matrix& m);

void print(double* data, unsigned size);

void print(float* data, unsigned size);

void print(double* data, unsigned rows, unsigned cols);

void print(float* data, unsigned rows, unsigned cols);

void print_to_csv(const Matrix& m, const std::string &filename);

// get the time
int64_t get_tick_us();

// note: CUDA stuff, only enabled if compiling for Nvidia GPUs
#ifdef __NVCC__

// print the properties of all connected GPUs
void print_device_properties();
#endif

} // namespace blast