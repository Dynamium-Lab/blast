#pragma once

#include <vector>
#include <cstdint>
#include <immintrin.h>
#include <cmath>

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

#include "blast_profile.hpp"

#include "blast_math.h"
#include "blast_world.h"
#include "blast_trajectory.h"
#include "blast_manipulator.h"
#include "blast_optimization.h"