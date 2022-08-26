#pragma once
#include <stdint.h>
#include <iostream>

//--- pre-processor directives to hand SIMD, CUDA, Etc

// use doubles by default
#ifndef BLAST_USE_DOUBLES
#define BLAST_USE_DOUBLES 1
#endif

// determine if the compilation enables building for nvidia GPUs
#ifdef __NVCC__
#define blast_fn __host__ __device__ inline
#define host_fn inline
#define dev_fn __device__ inline
#else
#pragma message("WARN, Not compiling with 'nvcc', the Nvidia GPU will not be used to speed up computation")
#define blast_fn inline
#define host_fn inline
#endif


//--- Other headers ---
#include "blast_error.hpp"
#include "blast_math.hpp"
#include "blast_trajectory.hpp"
#include "blast_manipulator.hpp"
#include "blast_optimization.hpp"
