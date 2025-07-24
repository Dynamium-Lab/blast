#pragma once
#include <blast>
#include <random>

namespace blast {


// Returns a rotation matrix from roll, pitch, yaw angles
inline Mat3 rpy2rotation(Vec3 rpy) {
    // Pre-compute cosines and sines.
    double cx = cos(rpy.x);
    double sx = sin(rpy.x);
    double cy = cos(rpy.y);
    double sy = sin(rpy.y);
    double cz = cos(rpy.z);
    double sz = sin(rpy.z);

    // Compute the rotation matrix components.
    Mat3 R;
    R(0, 0) = cy * cz;
    R(0, 1) = -cy * sz;
    R(0, 2) = sy;

    R(1, 0) = sx * sy * cz + cx * sz;
    R(1, 1) = -sx * sy * sz + cx * cz;
    R(1, 2) = -sx * cy;

    R(2, 0) = -cx * sy * cz + sx * sz;
    R(2, 1) = cx * sy * sz + sx * cz;
    R(2, 2) = cx * cy;

    return R;
}

// Returns the sign of a real number (zero if it is zero).
inline blast_fn real sign(real v) {
    return v > 0 ? (real)1: v == (real)0 ? 0: (real)-1;
}

// Returns the sign of a real number (zero will return 1).
inline blast_fn real sign_no_zero(real v) {
    return v >= 0 ? (real)1: (real)-1;
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