#pragma once
#include "blast.h"
#include <random>

namespace blast {

inline blast_fn real sign(real v) {
    return v > 0 ? (real)1: v == (real)0 ? 0: (real)-1;
}

inline blast_fn real wrap2pi(real r) {
    while (r < -(real)3.1415)
        r += 2*(real)3.1415;
    while (r > (real)3.1415)
        r-= 2*(real)3.1415;
    return r;
}

inline blast_fn float wrap_to_180(float r) {
    while (r < -180)
        r += 360;
    while (r > 180)
        r-= 360;
    return r;
}

inline blast_fn real deg2rad(real r) {
    return r * (real)3.1415/180;
}

inline blast_fn real rad2deg(real r) {
    return r * 180/(real)3.1415;
}

inline blast_fn Array rad2deg(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * 180/(real)3.1415;
    return r;
}

inline blast_fn Array deg2rad(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * (real)3.1415/180;
    return r;
}

inline blast_fn real clamp(real val, real mini, real maxi) {
    return val > maxi ? maxi: val < mini ? mini: val;
}

inline blast_fn real& clamp_inplace(real& val, real mini, real maxi) {
    val = val > maxi ? maxi: val < mini ? mini: val;
    return val;
}

#ifndef __CUDA_ARCH__
inline host_fn real get_random() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 e2(rd());
    static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
    return dis(e2);
}
#endif

} // namespace blast