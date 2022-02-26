#pragma once
#include <stdint.h>

//--- Basic types used everywhere ---
namespace blast {

using u8    = uint8_t;
using u16   = uint16_t;
using u32   = uint32_t;
using i8    = int8_t;
using i16   = int16_t;
using i32   = int32_t;
using real  = double; // for testing float vs double speed and precision

} // namespace blast

// todo: trigger breakpoint or something else
#ifdef _DEBUG
#define Assert(expr) \
    if (!(expr)){\
        fprintf(stderr, "Assertion failed in function: %s. File: %s(%d).\n", __FUNCTION__, __FILE__, __LINE__); \
        exit(-1); \
    }
#else
#define Assert(expr) do{} while(0)
#endif

//--- Other headers ---
#include "blast_math.hpp"
#include "blast_trajectory.hpp"
#include "blast_manipulator.hpp"
