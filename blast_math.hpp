#pragma once

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace blast {

static const real inf = std::numeric_limits<real>::infinity();
static const real nan = std::numeric_limits<real>::quiet_NaN();

// 3D vector
struct Vec3 {
    real x = 0;
    real y = 0;
    real z = 0;
};
void zero(Vec3&);
Vec3 cross(Vec3 a, Vec3 b);
real dot(Vec3 a, Vec3 b);
Vec3 operator-(Vec3 a, Vec3 b);
Vec3 operator+(Vec3 a, Vec3 b);
Vec3 operator*(real a, Vec3 b);
Vec3 operator*(Vec3 a, real b);



// 3x3 matrix
struct Mat3 {
    real data[9] = {0};

    real& operator()(u32 row, u32 col);
};
void zero(Mat3&);
Mat3 transpose(Mat3&);
Vec3 operator*(Mat3& m, Vec3 v);
Mat3 operator*(Mat3& m, Mat3 rhs);



// Array of real numbers
struct Array {
    real * data = nullptr;
    u32 size = 0;

    Array(u32 size);
    real& operator[](u32 i);
};
void zero(Array&);
void free(Array&);



// Matrix of real numbers
struct Matrix {
    real * data = nullptr;
    u32 size = 0;
    u32 rows = 0;
    u32 cols = 0;

    Matrix(u32 r, u32 c);
    real& operator()(u32 row, u32 col);
};
void free(Matrix&);



//--- Utility and debug functions ---
void print(Vec3);
void print(Mat3);
void print(Array&);
void print(Matrix&);








//-----------------------------------------------------------------------
inline Mat3 transpose(Mat3& m) {
    Mat3 result = {
        m(0, 0),
        m(0, 1),
        m(0, 2),
        m(1, 0),
        m(1, 1),
        m(1, 2),
        m(2, 0),
        m(2, 1),
        m(2, 2)
    };
    return result;
}

inline Vec3 cross(Vec3 a, Vec3 b) {
    Vec3 r;
    r.x = a.y*b.z - a.z*b.y;
    r.y = a.z*b.x - a.x*b.z;
    r.z = a.x*b.y - a.y*b.x;
    return r;
}

inline real dot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Vec3 operator-(Vec3 a, Vec3 b) {
    return Vec3{
        a.x - b.x,
        a.y - b.y,
        a.z - b.z
    };
}

inline Vec3 operator+(Vec3 a, Vec3 b) {
    return Vec3{
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    };
}

inline Vec3 operator*(real a, Vec3 b) {
    return Vec3{
        a * b.x,
        a * b.y,
        a * b.z
    };
}

inline Vec3 operator*(Vec3 a, real b) {
    return Vec3{
        a.x * b,
        a.y * b,
        a.z * b
    };
}

inline real& Mat3::operator()(u32 row, u32 col) {
    return data[3*col + row];
}

inline void zero(Vec3& v) {
    v.x = v.y = v.z = 0;
}

inline void zero(Mat3& m) {
    memset(m.data, 0, 9*sizeof(real));
}

inline void zero(Array& a) {
    memset(a.data, 0, a.size*sizeof(real));
}

inline void zero(Matrix& m) {
    memset(m.data, 0, m.size*sizeof(real));
}

inline Vec3 operator*(Mat3& m, Vec3 v) {
    Vec3 r;
    r.x = m.data[0]*v.x + m.data[3]*v.y + m.data[6]*v.z;
    r.y = m.data[1]*v.x + m.data[4]*v.y + m.data[7]*v.z;
    r.z = m.data[2]*v.x + m.data[5]*v.y + m.data[8]*v.z;
    return r;
}

inline Mat3 operator*(Mat3& m, Mat3 rhs) {
    Mat3 r;
    r.data[0] = m.data[0]*rhs.data[0] + m.data[3]*rhs.data[1] + m.data[6]*rhs.data[2];
    r.data[1] = m.data[1]*rhs.data[0] + m.data[4]*rhs.data[1] + m.data[7]*rhs.data[2];
    r.data[2] = m.data[2]*rhs.data[0] + m.data[5]*rhs.data[1] + m.data[8]*rhs.data[2];
    r.data[3] = m.data[0]*rhs.data[3] + m.data[3]*rhs.data[4] + m.data[6]*rhs.data[5];
    r.data[4] = m.data[1]*rhs.data[3] + m.data[4]*rhs.data[4] + m.data[7]*rhs.data[5];
    r.data[5] = m.data[2]*rhs.data[3] + m.data[5]*rhs.data[4] + m.data[8]*rhs.data[5];
    r.data[6] = m.data[0]*rhs.data[6] + m.data[3]*rhs.data[7] + m.data[6]*rhs.data[8];
    r.data[7] = m.data[1]*rhs.data[6] + m.data[4]*rhs.data[7] + m.data[7]*rhs.data[8];
    r.data[8] = m.data[2]*rhs.data[6] + m.data[5]*rhs.data[7] + m.data[8]*rhs.data[8];
    return r;
}

inline Array::Array(u32 size) {
    this->size = size;
    if (size)
        data = (real*)calloc(size, sizeof(real));
}

inline real& Array::operator[](u32 i) {
    Assert(i < size);
    return data[i];
}

inline Matrix::Matrix(u32 r, u32 c) {
    size = r * c;
    rows = r;
    cols = c;
    if (size)
        data = (real*)calloc(size, sizeof(real));
}

inline real& Matrix::operator()(u32 row, u32 col) {
    Assert(row < this->rows && col < this->cols);
    return this->data[row + this->rows*col];
}

inline void free(Array& v) {
    std::free(v.data);
    v.size = 0;
    v.data = nullptr;
}

inline void free(Matrix& m) {
    std::free(m.data);
    m.size = 0;
    m.data = nullptr;
}

inline void print(Vec3 v) {
    printf("[%f, %f, %f]\n", v.x, v.y, v.z);
}

inline void print(Mat3 m) {
    printf("\n[%f, %f, %f]\n[%f, %f, %f]\n[%f, %f, %f]\n",
           m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2), m(2, 0), m(2, 1), m(2, 2));
}

inline void print(Array& a) {
    printf("[");
    for (u32 i = 0; i < a.size-1; i++)
        printf("%0.4f, ", a[i]);
    printf("%0.4f]\n", a[a.size-1]);
}

inline void print(Matrix& m) {
    printf("\n");
    for (u32 i = 0; i < m.rows; i++) {
        printf("[");
        for (u32 j = 0; j < m.cols-1; j++) {
            printf("%0.4f, ", m(i, j));
        }
        printf("%0.4f]\n", m(i, m.cols-1));
    }
}

} // namespace blast
