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
#include <random>


#include <cmath>
#include <stdint.h>

namespace blast {
inline real get_random() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 e2(rd());
    static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
    return dis(e2);
}

inline int64_t get_tick() {
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

} // namespace blast

