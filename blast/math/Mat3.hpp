#pragma once

namespace blast {


inline blast_fn Mat3::Mat3(const Mat3& m) {
  memcpy(data, m.data, 9 * sizeof(real));
}

inline blast_fn Mat3::Mat3(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3) {
  data[0] = x1;
  data[1] = y1;
  data[2] = z1;
  data[3] = x2;
  data[4] = y2;
  data[5] = z2;
  data[6] = x3;
  data[7] = y3;
  data[8] = z3;
}

inline blast_fn real& Mat3::operator()(u32 row, u32 col) {
  return data[3 * col + row];
}

inline real Mat3::operator()(u32 row, u32 col) const {
  return data[3 * col + row];
}

inline blast_fn real& Mat3::operator[](u32 i) {
  Assert(i < 9);
  return data[i];
}

inline blast_fn real Mat3::operator[](u32 i) const {
  Assert(i < 9);
  return data[i];
}

inline blast_fn Mat3& Mat3::zero() {
  memset(data, 0, 9 * sizeof(real));
  return *this;
}

inline blast_fn Vec3 Mat3::col_copy(int c) const {
  Assert(c < 3);
  return Vec3(data[c * 3 + 0], data[c * 3 + 1], data[c * 3 + 2]);
}

inline blast_fn Array Mat3::col(int c) {
  Assert((c >= 0 && c < 3));
  Array result(data + 3 * c, 3);
  Assert(result.is_alias);
  return result;
}

inline blast_fn Mat3 operator*(Mat3 m, Mat3 rhs) {
  Mat3 r;
  r.data[0] = m.data[0] * rhs.data[0] + m.data[3] * rhs.data[1] + m.data[6] * rhs.data[2];
  r.data[1] = m.data[1] * rhs.data[0] + m.data[4] * rhs.data[1] + m.data[7] * rhs.data[2];
  r.data[2] = m.data[2] * rhs.data[0] + m.data[5] * rhs.data[1] + m.data[8] * rhs.data[2];
  r.data[3] = m.data[0] * rhs.data[3] + m.data[3] * rhs.data[4] + m.data[6] * rhs.data[5];
  r.data[4] = m.data[1] * rhs.data[3] + m.data[4] * rhs.data[4] + m.data[7] * rhs.data[5];
  r.data[5] = m.data[2] * rhs.data[3] + m.data[5] * rhs.data[4] + m.data[8] * rhs.data[5];
  r.data[6] = m.data[0] * rhs.data[6] + m.data[3] * rhs.data[7] + m.data[6] * rhs.data[8];
  r.data[7] = m.data[1] * rhs.data[6] + m.data[4] * rhs.data[7] + m.data[7] * rhs.data[8];
  r.data[8] = m.data[2] * rhs.data[6] + m.data[5] * rhs.data[7] + m.data[8] * rhs.data[8];
  return r;
}

inline blast_fn Mat3 operator*(real x, const Mat3& m) {
  Mat3 r;
  r.data[0] = x * m.data[0];
  r.data[1] = x * m.data[1];
  r.data[2] = x * m.data[2];
  r.data[3] = x * m.data[3];
  r.data[4] = x * m.data[4];
  r.data[5] = x * m.data[5];
  r.data[6] = x * m.data[6];
  r.data[7] = x * m.data[7];
  r.data[8] = x * m.data[8];
  return r;
}

inline blast_fn Mat3& operator*=(Mat3& lhs, Mat3 rhs) {
  lhs.data[0] = lhs.data[0] * rhs.data[0] + lhs.data[3] * rhs.data[1] + lhs.data[6] * rhs.data[2];
  lhs.data[1] = lhs.data[1] * rhs.data[0] + lhs.data[4] * rhs.data[1] + lhs.data[7] * rhs.data[2];
  lhs.data[2] = lhs.data[2] * rhs.data[0] + lhs.data[5] * rhs.data[1] + lhs.data[8] * rhs.data[2];
  lhs.data[3] = lhs.data[0] * rhs.data[3] + lhs.data[3] * rhs.data[4] + lhs.data[6] * rhs.data[5];
  lhs.data[4] = lhs.data[1] * rhs.data[3] + lhs.data[4] * rhs.data[4] + lhs.data[7] * rhs.data[5];
  lhs.data[5] = lhs.data[2] * rhs.data[3] + lhs.data[5] * rhs.data[4] + lhs.data[8] * rhs.data[5];
  lhs.data[6] = lhs.data[0] * rhs.data[6] + lhs.data[3] * rhs.data[7] + lhs.data[6] * rhs.data[8];
  lhs.data[7] = lhs.data[1] * rhs.data[6] + lhs.data[4] * rhs.data[7] + lhs.data[7] * rhs.data[8];
  lhs.data[8] = lhs.data[2] * rhs.data[6] + lhs.data[5] * rhs.data[7] + lhs.data[8] * rhs.data[8];
  return lhs;
}

inline blast_fn Mat3 operator+(Mat3 lhs, Mat3 rhs) {
  Mat3 r;
  for (u32 i = 0; i < 9; i++)
    r.data[i] = lhs.data[i] + rhs.data[i];
  return r;
}

inline blast_fn Mat3& operator+=(Mat3& lhs, Mat3 rhs) {
  for (u32 i = 0; i < 9; i++)
    lhs.data[i] += rhs.data[i];
  return lhs;
}

inline blast_fn Vec3 operator*(Mat3 m, Vec3 v) {
  Vec3 r;
  r.x = m.data[0] * v.x + m.data[3] * v.y + m.data[6] * v.z;
  r.y = m.data[1] * v.x + m.data[4] * v.y + m.data[7] * v.z;
  r.z = m.data[2] * v.x + m.data[5] * v.y + m.data[8] * v.z;
  return r;
}

inline blast_fn Mat3 operator-(Mat3 m1, Mat3 m2) {
  Mat3 r;
  for (u32 i = 0; i < 9; i++)
    r.data[i] = m1.data[i] - m2.data[i];
  return r;
}

inline blast_fn Mat3 transpose(Mat3 m) {
  Mat3 result{
          m.data[0],
          m.data[3],
          m.data[6],
          m.data[1],
          m.data[4],
          m.data[7],
          m.data[2],
          m.data[5],
          m.data[8]};
  return result;
}

inline blast_fn Mat3 eye() {
  return Mat3{
          1, 0, 0,
          0, 1, 0,
          0, 0, 1};
}

inline blast_fn Mat3& transpose_inplace(Mat3& m) {
  real tmp = m[1];
  m[1]     = m[3];
  m[3]     = tmp;
  tmp      = m[2];
  m[2]     = m[6];
  m[6]     = tmp;
  tmp      = m[5];
  m[5]     = m[7];
  m[7]     = tmp;
  return m;
}

inline blast_fn Mat3& zero(Mat3& m) {
  m.zero();
  return m;
}

inline blast_fn Mat3& constant(Mat3& m, real val) {
  for (u32 i = 0; i < 9; i++)
    m.data[i] = val;
  return m;
}

inline blast_fn bool operator==(const Mat3& a, const Mat3& b) {
  return is_close(a, b);
}

inline blast_fn bool is_close(const Mat3& a, const Mat3& b, real eps) {
  for (int i = 0; i < 9; i++)
    if (std::abs(a.data[i] - b.data[i]) > eps)
      return false;
  return true;
}

inline blast_fn bool is_small(const Mat3& m, real eps) {
  for (int i = 0; i < 9; i++)
    if (std::abs(m.data[i]) > eps)
      return false;
  return true;
}


} // namespace blast
