#pragma once
#include <stdint.h>

//--- Basic types used everywhere ---
namespace blast {

using u32 = uint32_t;
using i32 = int32_t;
using real = float; // for testing float vs double speed and precision

} // namespace blast

#ifdef _DEBUG
#define Assert(expr) \
    if (!(expr)){\
        fprintf(stderr, "Assertion failed in function: %s: %s. File: %s(%d).\n", __FUNCTION__, __LINE__, __FILE__, __LINE__); \
        exit(-1); \
    }
#else
#define Assert(expr) do{} while(0)
#endif

//--- Other headers ---
#include "blast_math.hpp"
#include "blast_trajectory.hpp"
