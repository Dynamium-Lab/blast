#pragma once

#include <vector>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <immintrin.h>

// determine if the compilation enables building for nvidia GPUs
#ifdef __NVCC__
#define blast_fn    __host__ __device__
#define host_fn     __host__
#define device_fn   __device__
#else
#define blast_fn
#define host_fn
#define device_fn
#endif

namespace blast {

using u8    = uint8_t;
using u16   = uint16_t;
using u32   = uint32_t;
using u64   = uint64_t;
using i8    = int8_t;
using i16   = int16_t;
using i32   = int32_t;
using i64   = int64_t;

} // namespace blast

#include "blast_error.hpp"
#include "blast_math.hpp"

#include "blast_world.h"
#include "blast_trajectory.h"
#include "blast_manipulator.h"
#include "blast_optimization.h"
