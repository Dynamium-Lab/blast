#pragma once
#include <random>

namespace blast {


// https://stackoverflow.com/a/49943540
inline host_fn double simd_hadd(__m256d v) {
  __m128d vlow  = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1);
  vlow          = _mm_add_pd(vlow, vhigh);

  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  return _mm_cvtsd_f64(_mm_add_sd(vlow, high64));
}

// https://stackoverflow.com/a/35270026
inline host_fn float simd_hadd(__m128 v) {
  __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
  __m128 sums = _mm_add_ps(v, shuf);
  shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
  sums        = _mm_add_ss(sums, shuf);
  return _mm_cvtss_f32(sums);
}

inline host_fn float simd_hadd(__m256 v) {
  __m128 vlow  = _mm256_castps256_ps128(v);
  __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
  vlow         = _mm_add_ps(vlow, vhigh);     // add the low 128
  return simd_hadd(vlow);
}


// Returns a rotation matrix from roll, pitch, yaw angles
inline Mat3 rpy2rotation(Vec3 rpy) {
  // Pre-compute cosines and sines.
  const double cx = cos(rpy.x);
  const double sx = sin(rpy.x);
  const double cy = cos(rpy.y);
  const double sy = sin(rpy.y);
  const double cz = cos(rpy.z);
  const double sz = sin(rpy.z);

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

inline Vec3 rotation2rpy(Mat3 rotation) {
  u32  condition = std::abs(rotation(0, 2)) != 1 ? 1 : rotation(0, 2) == 1 ? 2
                                                                           : 3; // allow for gimble lock
  real rx        = condition == 1 ? atan2(-rotation(1, 2), rotation(2, 2))
                                  : 0.0;


  real ry = condition == 1   ? asin(rotation(0, 2))
            : condition == 2 ? PI / 2
                             : -PI / 2;


  real rz = condition == 1   ? atan2(-rotation(0, 1), rotation(0, 0))
            : condition == 2 ? atan2(rotation(1, 0), rotation(1, 1))
                             : atan2(-rotation(1, 0), -rotation(1, 1));

  return {rx, ry, rz};
}

// Returns the sign of a real number (zero if it is zero).
inline blast_fn real sign(real v) {
  return v > 0 ? 1 : v == 0 ? 0
                            : -1;
}

// Returns the sign of a real number (zero will return 1).
inline blast_fn real sign_no_zero(real v) {
  return v >= 0 ? 1 : -1;
}

inline blast_fn real wrap2pi(real r) {
  while (r < -PI)
    r += 2 * PI;
  while (r > PI)
    r -= 2 * PI;
  return r;
}

// return a copy of a with all members wrapped to $\pi
inline blast_fn Array wrap2pi(Array a) {
  for (u32 i = 0; i < a.size; i++) {
    a[i] = wrap2pi(a[i]);
  }
  return a;
}

inline blast_fn real wrap_to_180(real r) {
  while (r < -180.0)
    r += 360.0;
  while (r > 180.0)
    r -= 360.0;
  return r;
}

inline blast_fn float wrap_to_180(float r) {
  while (r < -180)
    r += 360;
  while (r > 180)
    r -= 360;
  return r;
}

inline blast_fn real deg2rad(real r) {
  return r * PI / 180;
}

inline blast_fn real rad2deg(real r) {
  return r * 180 / PI;
}

inline blast_fn Array rad2deg(const Array& a) {
  Array r(a.size);
  for (u32 i = 0; i < a.size; i++)
    r[i] = a[i] * 180 / PI;
  return r;
}

inline blast_fn Array deg2rad(const Array& a) {
  Array r(a.size);
  for (u32 i = 0; i < a.size; i++)
    r[i] = a[i] * PI / 180;
  return r;
}

inline blast_fn real clamp(real val, real mini, real maxi) {
  return val > maxi ? maxi : val < mini ? mini
                                        : val;
}

inline blast_fn real& clamp_inplace(real& val, real mini, real maxi) {
  val = val > maxi ? maxi : val < mini ? mini
                                       : val;
  return val;
}

#ifndef __CUDA_ARCH__
inline host_fn real random_real() {
  static thread_local std::random_device                          rd;
  static thread_local std::mt19937                                e2(rd());
  static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
  return dis(e2);
}

inline host_fn u32 random_int(u32 min, u32 max) {
  // thread_local: the generator carries mutable state, so a shared static
  // would be a data race when guesses are generated from Taskflow workers.
  // Mirrors random_real() above; each thread seeds its own generator once.
  static thread_local std::random_device rd;
  static thread_local std::mt19937       gen(rd());

  // Distribution is cheap and stateless across calls — keep it local so the
  // [min, max] range can vary per call.
  std::uniform_int_distribution<u32> dist(min, max);

  return dist(gen);
}
#endif

} // namespace blast
