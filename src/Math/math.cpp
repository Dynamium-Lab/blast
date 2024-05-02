#include "blast.h"
#include <random>

namespace blast {

real sign(real v) {
    return v > 0 ? (real)1: v == (real)0 ? 0: (real)-1;
}

real wrap2pi(real r) {
    while (r < -(real)3.1415)
        r += 2*(real)3.1415;
    while (r > (real)3.1415)
        r-= 2*(real)3.1415;
    return r;
}

float wrap_to_180(float r) {
    while (r < -180)
        r += 360;
    while (r > 180)
        r-= 360;
    return r;
}

real deg2rad(real r) {
    return r * (real)3.1415/180;
}

real rad2deg(real r) {
    return r * 180/(real)3.1415;
}

Array rad2deg(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * 180/(real)3.1415;
    return r;
}

Array deg2rad(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * (real)3.1415/180;
    return r;
}

real clamp(real val, real mini, real maxi) {
    return std::min(maxi, std::max(val, mini));
}

real& clamp_inplace(real& val, real mini, real maxi) {
    val = clamp(val, mini, maxi);
    return val;
}

real get_random() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 e2(rd());
    static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
    return dis(e2);
}


} // namespace blast