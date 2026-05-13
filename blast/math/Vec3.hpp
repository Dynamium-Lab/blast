#pragma once

#include <blast>

namespace blast {

//------ Vec3 ---------------------

// NOTE: Unsafe operator[]. Use only when i ∈ [0,2]. No runtime check in release builds.
inline blast_fn real& Vec3::operator[](int i) {
  Assert(i < 3);
  return *(&x + i); // note: Don't try this at home
}

// NOTE: Unsafe operator[]. Use only when i ∈ [0,2]. No runtime check in release builds.
inline blast_fn real Vec3::operator[](int i) const {
  Assert(i < 3);
  return *(&x + i); // note: Don't try this at home
}

inline blast_fn Vec3::Vec3(real x, real y, real z) :
    x(x),
    y(y),
    z(z),
    _pad(0) {
}

inline blast_fn bool Vec3::operator==(const Vec3& a) const {
  return is_close(*this, a);
}

inline blast_fn Vec3 cross(Vec3 a, Vec3 b) {
  Vec3 r;
  r.x = a.y * b.z - a.z * b.y;
  r.y = a.z * b.x - a.x * b.z;
  r.z = a.x * b.y - a.y * b.x;
  return r;
}

inline blast_fn real dot(Vec3 a, Vec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline blast_fn real norm(Vec3 a) {
  return (real) sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

inline blast_fn bool is_close(const Vec3& a, const Vec3& b, real eps) {
  return (std::abs(a.x - b.x) < eps) && (std::abs(a.y - b.y) < eps) && (std::abs(a.z - b.z) < eps);
}

inline blast_fn bool is_small(const Vec3& a, real eps) {
  return (std::abs(a.x) < eps) && (std::abs(a.y) < eps) && (std::abs(a.z) < eps);
}

inline blast_fn Vec3& zero(Vec3& a) {
  a.x = 0;
  a.y = 0;
  a.z = 0;
  return a;
}

inline blast_fn Vec3& constant(Vec3& a, real val) {
  a.x = val;
  a.y = val;
  a.z = val;
  return a;
}

inline blast_fn Vec3 operator-(Vec3 a) {
  Vec3 r;
  r.x = -a.x;
  r.y = -a.y;
  r.z = -a.z;
  return r;
}

inline blast_fn Vec3 operator-(Vec3 a, Vec3 b) {
  return Vec3{
          a.x - b.x,
          a.y - b.y,
          a.z - b.z};
}

inline blast_fn Vec3 operator+(Vec3 a, Vec3 b) {
  return Vec3{
          a.x + b.x,
          a.y + b.y,
          a.z + b.z};
}

inline blast_fn Vec3 operator*(real a, Vec3 b) {
  return Vec3{
          a * b.x,
          a * b.y,
          a * b.z};
}

inline blast_fn Vec3 operator*(Vec3 a, real b) {
  return Vec3{
          a.x * b,
          a.y * b,
          a.z * b};
}

inline blast_fn Vec3 operator/(Vec3 a, real b) {
  return Vec3{
          a.x / b,
          a.y / b,
          a.z / b};
}

inline blast_fn Vec3& operator+=(Vec3& v1, const Vec3& v2) {
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
  return v1;
}

inline blast_fn Vec3& operator-=(Vec3& v1, const Vec3& v2) {
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
  return v1;
}

inline blast_fn Vec3& operator*=(Vec3& v, real a) {
  v.x *= a;
  v.y *= a;
  v.z *= a;
  return v;
}

inline blast_fn Vec3& operator*=(Vec3& v, real a) {
  v.x *= a;
  v.y *= a;
  v.z *= a;
  return v;
}

inline blast_fn Vec3& operator/=(Vec3& v, real a) {
  v.x /= a;
  v.y /= a;
  v.z /= a;
  return v;
}
} // namespace blast
