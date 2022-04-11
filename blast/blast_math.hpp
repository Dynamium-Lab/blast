#pragma once
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <utility>
#include <algorithm>

#include "blast.hpp"
#include "blast_simd.hpp"

namespace blast {

static const real inf = std::numeric_limits<real>::infinity();
static const real nan = std::numeric_limits<real>::quiet_NaN();
static const real pi  = 3.1415;

// 3D vector
struct Vec3 {
    real x = 0;
    real y = 0;
    real z = 0;
    real _pad = 0;

    Vec3() = default;
    Vec3(real x, real y, real z);
};
void zero(Vec3&);
Vec3 cross(Vec3 a, Vec3 b);
real dot(Vec3 a, Vec3 b);
Vec3 operator-(Vec3 a, Vec3 b);
Vec3 operator+(Vec3 a, Vec3 b);
Vec3 operator*(real a, Vec3 b);
Vec3 operator*(Vec3 a, real b);
Vec3& operator*=(Vec3&v, real a);

// 3x3 matrix
struct Mat3 {
    real data[9] = {0};

    Mat3() = default;
    Mat3(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3);
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

    Array() {}
    Array(u32 new_size);    // normal constructor
    Array(const Array&);    // copy constructor
    Array(Array&&);         // move constructor
    Array& operator=(const Array&); // copy assignment
    Array& operator=(Array&&);      // move assignment
    Array& operator=(const std::initializer_list<real>& other);
    ~Array();

    real& operator[](u32 i);
    real operator[](u32 i) const;
    Array operator-();

    void resize(u32 new_size);
    void zero();
    real& back();
    real back() const;
};
void zero(Array&);

// Matrix of real numbers
struct Matrix {
    real * data = nullptr;
    u32 size = 0;
    u32 rows = 0;
    u32 cols = 0;

    Matrix() {}
    Matrix(u32 r, u32 c);   // normal constructor
    Matrix(const Matrix&);  // copy constructor
    Matrix(Matrix&&);       // move constructor
    Matrix& operator=(const Matrix&); // copy assignment
    Matrix& operator=(Matrix&&);      // move assignment
    ~Matrix();

    real& operator()(u32 row, u32 col);
    real operator()(u32 row, u32 col) const;
    void zero();
};
void zero(Matrix&);

//--- Utility and debug functions ---
void print(Vec3);
void print(Mat3);
void print(Array&);
void print(Matrix&);










//-----------------------------------------------------------------------
inline Vec3::Vec3(real x, real y, real z)
    : x(x), y(y), z(z), _pad(0) {

}

inline Mat3::Mat3(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3) {
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

inline Mat3 transpose(Mat3& m) {
    Mat3 result {
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

inline Vec3& operator*=(Vec3&v, real a) {
    v.x += a;
    v.y += a;
    v.z += a;
    return v;
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
    a.zero();
}

inline void zero(Matrix& m) {
    m.zero();
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

// inline Mat3 operator*(Mat3& A, Mat3 B) {
//     Mat3 r;
//     __m256d A1 = _mm256_loadu_pd(A.data);
//     __m256d A2 = _mm256_loadu_pd(A.data+3);
//     __m256d A3 = _mm256_loadu_pd(A.data+6);
//     auto t1 = _mm256_add_pd(
//                   _mm256_add_pd(
//                       _mm256_mul_pd(A1, _mm256_set1_pd(B(0, 0))),
//                       _mm256_mul_pd(A2, _mm256_set1_pd(B(1, 0)))),
//                   _mm256_mul_pd(A3, _mm256_set1_pd(B(2, 0))));
//     auto t2 = _mm256_add_pd(
//                   _mm256_add_pd(
//                       _mm256_mul_pd(A1, _mm256_set1_pd(B(0, 1))),
//                       _mm256_mul_pd(A2, _mm256_set1_pd(B(1, 1)))),
//                   _mm256_mul_pd(A3, _mm256_set1_pd(B(2, 1))));
//     auto t3 = _mm256_add_pd(
//                   _mm256_add_pd(
//                       _mm256_mul_pd(A1, _mm256_set1_pd(B(0, 2))),
//                       _mm256_mul_pd(A2, _mm256_set1_pd(B(1, 2)))),
//                   _mm256_mul_pd(A3, _mm256_set1_pd(B(2, 2))));
//     _mm256_storeu_pd(r.data, t1);
//     _mm256_storeu_pd(r.data+3, t2);
//     _mm256_storeu_pd(r.data+6, t3);
//     return r;
// }

inline Array::Array(u32 size) : size(size) {
    if (size)
        data = (real*)calloc(size, sizeof(real));
}

inline Array::Array(const Array& a) : size(a.size) {
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, a.data, size*sizeof(real));
    }
}

inline Array::Array(Array&& a) : data(a.data), size(a.size) {
    a.data = nullptr;
    a.size = 0;
}

inline Array& Array::operator=(const Array& a) {
    if (this != &a) {
        if (!data) {
            data = (real*)malloc(a.size * sizeof(real));
        }
        else if (a.size >= size) {
            std::free(data);
            data = (real*)malloc(a.size * sizeof(real));
        }
        Assert(data);

        size = a.size;
        std::memcpy(data, a.data, size*sizeof(real));
    }
    return *this;
}

inline Array& Array::operator=(Array&& a) {
    if (this != &a) {
        if (data)
            std::free(data);
        data = a.data;
        size = a.size;
        a.data = nullptr;
        a.size = 0;
    }
    return *this;
}

inline Array& Array::operator=(const std::initializer_list<real>& other) {
    Assert(other.size() <= size);
    std::copy_n(other.begin(), size, data);
    return *this;
}

inline Array Array::operator-() {
    Array result(size);
    for (u32 i = 0; i < size; i++) {
        result[i] = -data[i];
    }
    return std::move(result);
}

inline Array::~Array() {
    if (data)
        std::free(data);
}

inline real& Array::operator[](u32 i) {
    Assert(i < size);
    return data[i];
}

inline real Array::operator[](u32 i) const {
    Assert(i < size);
    return data[i];
}

inline void Array::resize(u32 new_size) {
    data = (real*)realloc(data, new_size*sizeof(real));
}

inline void Array::zero() {
    if(data)
        for (u32 i = 0; i<size; i++)
            data[i] = 0;
}

inline real& Array::back() {
    return data[size-1];
}

inline real Array::back() const {
    return data[size-1];
}

inline Matrix::Matrix(u32 r, u32 c) {
    size = r * c;
    rows = r;
    cols = c;
    if (size)
        data = (real*)calloc(size, sizeof(real));
}

inline Matrix::Matrix(const Matrix& m) : size(m.size), cols(m.cols), rows(m.rows) {
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, m.data, size*sizeof(real));
    }
}

inline Matrix::Matrix(Matrix&& m) : data(m.data), size(m.size), cols(m.cols), rows(m.rows) {
    m.data = nullptr;
    m.size = 0;
    m.rows = 0;
    m.cols = 0;
}

inline Matrix& Matrix::operator=(const Matrix& m) {
    if (this != &m) {
        if (!data) {
            data = (real*)malloc(m.size * sizeof(real));
        }
        else if (m.size >= size) {
            std::free(data);
            data = (real*)malloc(m.size * sizeof(real));
        }
        Assert(data);

        size = m.size;
        cols = m.cols;
        rows = m.rows;
        std::memcpy(data, m.data, size*sizeof(real));
    }
    return *this;
}

inline Matrix& Matrix::operator=(Matrix&& m) {
    if (this != &m) {
        if (data)
            std::free(data);
        data = m.data;
        size = m.size;
        rows = m.rows;
        cols = m.cols;
        m.data = nullptr;
        m.size = 0;
        m.rows = 0;
        m.cols = 0;
    }
    return *this;
}

inline Matrix::~Matrix() {
    if (data)
        std::free(data);
}

inline void Matrix::zero() {
    memset(data, 0, size*sizeof(real));
}

inline real& Matrix::operator()(u32 row, u32 col) {
    Assert(row < this->rows && col < this->cols);
    return this->data[row + this->rows*col];
}

inline real Matrix::operator()(u32 row, u32 col) const {
    Assert(row < this->rows && col < this->cols);
    return this->data[row + this->rows*col];
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
