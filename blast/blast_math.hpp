#pragma once
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include <random>

#include "blast.hpp"
#include "blast_simd.hpp"
#include "blast_error.hpp"
#ifdef __CUDA_ARCH__
#include <math_constants.h>
#endif

#include "mipp/mipp.h"

namespace blast {

// uses doubles by default unless BLAST_USE_DOUBLES is set to 0
#if BLAST_USE_DOUBLES
using real  = double;
#define BLAST_SIZEOF_REAL 8
#else
using real  = float;
#define BLAST_SIZEOF_REAL 4
#endif
static_assert(sizeof(real) == BLAST_SIZEOF_REAL);

// various type aliases
using u8    = uint8_t;
using u16   = uint16_t;
using u32   = uint32_t;
using i8    = int8_t;
using i16   = int16_t;
using i32   = int32_t;
using i64   = int64_t;
using svector = std::vector<real>;

// forward declarations
struct Vec3;
struct Mat3;
struct Mat4;
struct Array;
struct Matrix;

// useful constants
#ifdef __CUDA_ARCH__
#if BLAST_USE_DOUBLES
#define NAN_REAL CUDART_NAN
#define INF_REAL CUDART_INF
#else
#define NAN_REAL CUDART_NAN_F
#define INF_REAL CUDART_INF_F
#endif
#else
#define NAN_REAL std::numeric_limits<real>::quiet_NaN()
#define INF_REAL std::numeric_limits<real>::infinity()
#endif
#define PI (real)3.1415

// 3D vector
struct Vec3 {
    real x = 0;
    real y = 0;
    real z = 0;
    real _pad = 0;

    Vec3() = default;

    blast_fn real& operator[](u32 i);
    blast_fn Vec3(real x, real y, real z);
};

// 3x3 matrix
struct Mat3 {
    real data[9];

    Mat3() = default;
    blast_fn Mat3(const Mat3& m);
    blast_fn Mat3(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3);

    blast_fn real& operator()(u32 row, u32 col);
    blast_fn real& operator[](u32 i);
    blast_fn real operator[](u32 i) const ;

    blast_fn Vec3 col(u32 c);
};

// 4x4 matrix
struct alignas(32) Mat4 {
    real data[16];

    // default
    Mat4() = default;

    // copy constructor
    blast_fn Mat4(const Mat4&);

    // create a diagonal matrix with the value
    blast_fn Mat4(real v);

    // copy list into the matrix
    //  - note: the rest are NOT initialized if the list does not have 16 elements
    blast_fn Mat4(const std::initializer_list<real>&);

    // set every value to 0
    blast_fn void zero();

    // copy assignment
    blast_fn Mat4& operator=(const Mat4&);

    // copy list into the matrix
    //  - note: the rest are NOT initialized if the list does not have 16 elements
    blast_fn Mat4& operator=(const std::initializer_list<real>&);

    // access the given element
    blast_fn real& operator()(u32 row, u32 col);

    // access the given element (value)
    blast_fn real operator()(u32 row, u32 col) const;

    // access the element
    //  - note: does not check out of bounds
    blast_fn real& operator[](u32 i);

    // access value of the element
    //  - note: does not check out of bounds
    blast_fn real operator[](u32 i) const;

    // unary minus
    blast_fn Mat4 operator-();

    // matrix multiply inplace
    blast_fn Mat4& operator*=(Mat4&);

    // multiply each element by the value inplace
    blast_fn Mat4& operator*=(real);

    // matrix multiply inplace
    blast_fn Mat4& operator+=(Mat4&);

    // matrix multiply inplace
    blast_fn Mat4& operator-=(Mat4&);

    // check if all values of another array are the same
    blast_fn bool operator==(const Mat4&);

    // return a 4D array pointing to the column
    //  - note: the resulting array is an alias
    blast_fn Array col(u32 c);

    // return a 4D array copying the given colum
    //  - note: Copies data because we are const
    //  - note: resulting Array is NOT an alias
    blast_fn Array col(u32 c) const;

    // // interpret as a Matrix
    // //  note: the resulting Matrix is an alias
    // operator Matrix&();

    // // interpret as a 16D Array
    // //  note: the resulting Array is an alias
    // operator Array&();

};

// Array of real numbers
struct Array {
    real* data = nullptr; // pointer to the heap allocated memory
    u32   size = 0; // size in number of elements
    bool  is_alias = false; // if true, the Array does not own the data and will not free it. note: should be as short-lived as possible

    // default constructor
    Array() = default;

    // construct and allocate memory for 'n' elements
    blast_fn Array(u32 n);

    // construct and fill with value
    blast_fn Array(u32 n, real value);

    // copy constructor
    blast_fn Array(const Array&);

    // move constructor
    blast_fn Array(Array&&);

    // create an Array from pre-existing data
    //  - note: becomes an alias
    blast_fn Array(real*, u32 n);

    // create an Array from a const std::vector<real
    //  - note: copies data
    host_fn Array(const svector&);

    // free memory if not an alias
    blast_fn ~Array();

    // copy assignment
    blast_fn Array& operator=(const Array&);

    // move assignment
    blast_fn Array& operator=(Array&&);

    // copy data from a list
    //  - note: must be the same size
    blast_fn Array& operator=(const std::initializer_list<real>& other);

    // access the element
    //  - note: does not check out of bounds
    blast_fn real& operator[](u32 i);

    // access value of the element
    //  - note: does not check out of bounds
    blast_fn real operator[](u32 i) const;

    // unary minus
    blast_fn Array operator-();

    // multiply with another vector (in place)
    blast_fn Array& operator*=(real);

    // check if all values of another array are very close the this one (1e-06)
    blast_fn bool operator==(Array&);

    // map the array to the data of the given std::vector<real>
    //  - note: becomes an alias
    host_fn Array& alias(std::vector<real>&);

    // map the array to the data of the given matrix
    //  - note: interpret all the data of the matrix as one long array, becomes an alias)
    blast_fn Array& alias(Matrix& m);

    // map the Array to the given data
    //  - note: becomes an alias
    blast_fn Array& alias(real*, u32);
    blast_fn Array& alias(const real*, u32);

    // resize the array
    //  - note: old pointers (aliases) to this data may be invalidated
    //  - note: fails if the array is an alias
    blast_fn void resize(u32 new_size);

    // access the last element
    blast_fn real& back();

    // access the value of the last element
    blast_fn real back() const;
};

// Matrix of real numbers
struct Matrix {
    real* data = nullptr; // pointer to the heap allocated memory
    u32 size = 0; // size in number of elements
    u32 rows = 0; // number of rows in the matrix
    u32 cols = 0; // number of columns in the matrix
    bool is_alias = false; // if true, the Matrix does not own the data and will not free it.
    // note: Alias matrices should be as short-lived as possible

    // default constructor
    Matrix() = default;

    // construct and allocate memory for 'r x c' elements
    blast_fn Matrix(u32 r, u32 c);

    // copy constructor
    blast_fn Matrix(const Matrix&);

    // move constructor
    blast_fn Matrix(Matrix&&);

    // create a Matrix from pre-existing data
    //  - note: becomes an alias
    blast_fn Matrix(real*, u32 r, u32 c);

    // create an matrix of a single columns from an array
    //  - note: creates a copy
    blast_fn Matrix(const Array&);

    // free memory if not an alias
    blast_fn ~Matrix();

    // resize the matrix
    host_fn void resize(u32 r, u32 c);

    // copy assignment
    blast_fn Matrix& operator=(const Matrix&);

    // move assignment
    blast_fn Matrix& operator=(Matrix&&);

    // map the Matrix to the data of the given Array
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    blast_fn Matrix& alias(Array&);

    // map the Matrix to the data of the given std::vector<real>
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    host_fn Matrix& alias(std::vector<real>&);

    // map the Matrix to the given data
    //  - note: becomes an alias
    blast_fn Matrix& alias(real*, u32 r, u32 c);

    // access element
    //  - note: fails if the matrix is not big enough
    blast_fn real& operator()(u32 row, u32 col);

    // access value of element
    //  - note: fails if the matrix is not big enough
    blast_fn real operator()(u32 row, u32 col) const;

    // get the address of an element
    blast_fn real* address(u32 row, u32 col) const;

    // return an array accessing the given colum
    //  - note: new Array is aliasing our data
    blast_fn Array col(u32 c) const;
};

//------ MISC ---------------------

// wrap to -pi to pi
blast_fn real wrap2pi(real r) {
    while (r < -(real)3.1415)
        r += 2*(real)3.1415;
    while (r > (real)3.1415)
        r-= 2*(real)3.1415;
    return r;
}

// wrap to -180 to 180
blast_fn float wrap_to_180(float r) {
    while (r < -180)
        r += 360;
    while (r > 180)
        r-= 360;
    return r;
}

// convert degrees to radians
blast_fn real deg2rad(real r) {
    return r * (real)3.1415/180;
}

// convert radians to degrees
blast_fn real rad2deg(real r) {
    return r * 180/(real)3.1415;
}

// return an array with the values converted from radians to degrees
blast_fn Array rad2deg(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * 180/(real)3.1415;
    return r;
}

// return an array with the values converted from degrees to radians
blast_fn Array deg2rad(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * (real)3.1415/180;
    return r;
}

// fill the vector with zeros
blast_fn void zero(Vec3& v) {
    v.x = v.y = v.z = 0;
}

// fill the matrix with zeros
blast_fn void zero(Mat3& m) {
    memset(m.data, 0, 9*sizeof(real));
}

// fill the array with zeros
blast_fn void zero(Array& a) {
    if(a.data)
        memset(a.data, 0, a.size*sizeof(real));
}

// fill the matrix with zeros
blast_fn void zero(Matrix& m) {
    if(m.data)
        memset(m.data, 0, m.size*sizeof(real));
}

// fill the vector with a constant value
blast_fn void constant(Vec3& v, real val) {
    v.x = val;
    v.y = val;
    v.z = val;
}

// fill the matrix with a constant value
blast_fn void constant(Mat3& m, real val) {
    for (u32 i = 0; i < 9; i++)
        m.data[i] = val;
}

// fill the matrix with a constant value
blast_fn void constant(Mat4& m, real val) {
    for (u32 i = 0; i < 16; i++)
        m.data[i] = val;
}

// fill the array with a constant value
blast_fn void constant(Array& a, real val) {
    for (u32 i = 0; i < a.size; i++)
        a[i] = val;
}

// fill the matrix with a constant value
blast_fn void constant(Matrix& m, real val) {
    for (u32 i = 0; i < m.size; i++)
        m.data[i] = val;
}

// return the minimum value of the array (sign matters)
blast_fn real array_min(const Array& a) {
    real result = INF_REAL;
    for(u32 i = 0; i < a.size; i++)
        result = a[i] < result ? a[i] : result;
    return result;
}

// return the argument of the minimum value of the array
blast_fn unsigned array_min_id(const Array& a) {
    real result = INF_REAL;
    real result_old;
    unsigned i_min;
    for(u32 i = 0; i < a.size; i++) {
        result_old = result;
        result = a[i] < result ? a[i] : result;
        i_min = a[i] < result_old ? i : i_min;
    }
    return i_min;
}

// return the maximum value of the array (sign matters)
blast_fn real array_max(const Array& a) {
    real result = -INF_REAL;
    for(u32 i = 0; i < a.size; i++)
        result = a[i] > result ? a[i] : result;
    return result;
}

// return the minimum value of the matrix (sign matters)
blast_fn real matrix_min(const Matrix& m) {
    real result = INF_REAL;
    for(u32 i = 0; i < m.size; i++)
        result = m.data[i] < result ? m.data[i] : result;
    return result;
}

// return the argument of the minimum value of the matrix
blast_fn unsigned matrix_min_id(const Matrix& m) {
    real result = INF_REAL;
    real result_old;
    unsigned i_min;
    for(u32 i = 0; i < m.size; i++) {
        result_old = result;
        result = m.data[i] < result ? m.data[i] : result;
        i_min = m.data[i] < result_old ? i : i_min;
    }
    return i_min;
}

// return the maximum value of the matrix (sign matters)
blast_fn real matrix_max(const Matrix& m) {
    real result = -INF_REAL;
    for(u32 i = 0; i < m.size; i++)
        result = m.data[i] > result ? m.data[i] : result;
    return result;
}

// return true if the difference of all elements of each array are close to zero
blast_fn bool is_close(const Array& a1, const Array& a2, real eps = 1e-05) {
    Assert(a1.size == a2.size);
    for (u32 i =0; i < a1.size; i++)
        if(a1[i] - a2[i] > eps || a1[i] - a2[i] < -eps)
            return false;
    return true;
}

// return true if all elements are close to zero
blast_fn bool is_small(const Array& a, real eps = 1e-05) {
    for (u32 i =0; i < a.size; i++)
        if(abs(a[i]) > eps)
            return false;
    return true;
}

// return 1 if v > 0, 0 if v == 0, -1 if v < 0
blast_fn real sign(real v) {
    return v > 0 ? 1: v == 0 ? 0: -1;
}

// return a random number between -1 and 1
host_fn real get_random() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 e2(rd());
    static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
    return dis(e2);
}

// fill the given array with random values between -A and A
host_fn void fill_random(Array& v, real A) {
    for (int i = 0; i < (int)v.size; i++)
        v[i] = A * get_random();
    // todo: performance
}

// fill the given matrix with random values between -A and A
host_fn void fill_random(Mat3& m, real A) {
    for (int i = 0; i < 9; i++)
        m[i] = A * get_random();
}

// create a new array with random values between -A and A
host_fn Array random_array(u32 n, real A) {
    Array result(n);
    for (int i = 0; i < (int)n; i++)
        result[i] = A * get_random();
    return result;
}

// return the value clamped to [mini, maxi]
blast_fn real clamp(real val, real mini, real maxi) {
    real r = val < mini ? mini : val;
    r = r > maxi ? maxi : r;
    return r;
}

// return a new array with each value clamped to [mini, maxi]
blast_fn Array clamp(const Array& a, real mini, real maxi) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = clamp(a[i], mini, maxi);
    return r;
}

// return a new array with each value clamped to [lb, ub] (each element can have different lb and ub)
blast_fn Array clamp(const Array& a, const Array& lb, const Array& ub) {
    Assert(a.size == lb.size && a.size == ub.size);
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = clamp(a[i], lb[i], ub[i]);
    return r;
}

// modify the value to clamp to [mini, maxi]
blast_fn void clamp_inplace(real& val, real mini, real maxi) {
    val = clamp(val, mini, maxi);
}

// modify the array to clamp each value to [mini, maxi]
blast_fn void clamp_inplace(Array& a, real mini, real maxi) {
    for (u32 i = 0; i < a.size; i++)
        clamp_inplace(a[i], mini, maxi);
}

// clamp inplace with lower bound and upper bound
blast_fn void clamp_inplace(Array& a, const Array& lb, const Array& ub) {
    Assert(a.size == lb.size && a.size == ub.size);
    for (u32 i = 0; i < a.size; i++)
        a[i] = clamp(a[i], lb[i], ub[i]);
}

// return the sum of the array
blast_fn real sum(const Array& a) {
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += a[i];
    return result;
}

// return the mean of the array
blast_fn real mean(const Array& a) {
    return sum(a) / a.size;
}

// return the euclidean norm (2-norm)
blast_fn real norm(const Array& a) {
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += a[i] * a[i];
    return sqrt(result);
}

// return the sum of the square of each element
blast_fn real norm_sqr(const Array& a) {
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += a[i] * a[i];
    return result;
}

// return the infinity norm
blast_fn real norm_inf(const Array& a) {
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result = abs(a[i]) > result ? abs(a[i]) : result;
    return result;
}

// return the 1-norm
blast_fn real norm1(const Array& a) {
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += ::abs(a[i]);
    return result;
}

// return an array of the same size with each element being the absolute value of the corresponding element
blast_fn Array array_abs(const Array& a) {
    Array result(a.size);
    for (u32 i = 0; i < a.size; i++)
        result[i] = ::abs(a[i]);
    return result;
}

// return an array of the same size with each element being the square of the corresponding element
blast_fn Array array_abs2(const Array& a) {
    Array result(a.size);
    for (u32 i = 0; i < a.size; i++)
        result[i] = a[i] * a[i];
    return result;
}

// modify the array in place with each element being the absolute value of the corresponding element
blast_fn void array_abs_inplace(Array& a) {
    for (u32 i = 0; i < a.size; i++)
        a[i] = ::abs(a[i]);
}

// modify the array in place with each element being the square of the corresponding element
blast_fn void array_abs2_inplace(Array& a) {
    for (u32 i = 0; i < a.size; i++)
        a[i] = a[i] * a[i];
}

// fill the destination pointer with a matrix of the same size as m with each element (i,j) being a[i] - m(i,j)
blast_fn void minus_insert(const Array& a, const Matrix& m, real* dst) {
    Assert(m.rows == a.size);
    auto m_data = m.data;
    for (u32 c = 0; c < m.cols; c++) {
        for (u32 r = 0; r < a.size; r++) {
            *dst = a[r] - *m_data;
            m_data++;
            dst++;
        }
    }
}

// fill the destination pointer with a matrix of the same size as m with each element (i,j) being m(i,j) - a[i]
blast_fn void minus_insert(const Matrix& m, const Array& a, real* dst) {
    Assert(m.rows == a.size);
    auto m_data = m.data;
    for (u32 c = 0; c < m.cols; c++) {
        for (u32 r = 0; r < a.size; r++) {
            *dst = *m_data - a[r];
            m_data++;
            dst++;
        }
    }
}

//------ Vec3 ---------------------
blast_fn real& Vec3::operator[](u32 i) {
    Assert(i < 3);
    return *(&x + i); //note: Don't try this at home
}

blast_fn Vec3::Vec3(real x, real y, real z)
    : x(x), y(y), z(z), _pad(0) {
}

blast_fn Vec3 cross(Vec3 a, Vec3 b) {
    Vec3 r;
    r.x = a.y*b.z - a.z*b.y;
    r.y = a.z*b.x - a.x*b.z;
    r.z = a.x*b.y - a.y*b.x;
    return r;
}

blast_fn real dot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

blast_fn real norm(Vec3 a) {
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

blast_fn Vec3 operator-(const Vec3& a) {
    Vec3 r;
    r.x = -a.x;
    r.y = -a.y;
    r.z = -a.z;
    return r;
}

blast_fn Vec3 operator-(Vec3 a, Vec3 b) {
    return Vec3{
        a.x - b.x,
        a.y - b.y,
        a.z - b.z
    };
}

blast_fn Vec3 operator+(Vec3 a, Vec3 b) {
    return Vec3{
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    };
}

blast_fn Vec3 operator*(real a, Vec3 b) {
    return Vec3{
        a * b.x,
        a * b.y,
        a * b.z
    };
}

blast_fn Vec3 operator*(Vec3 a, real b) {
    return Vec3{
        a.x * b,
        a.y * b,
        a.z * b
    };
}

blast_fn Vec3& operator+=(Vec3& v1, const Vec3& v2) {
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

blast_fn Vec3& operator*=(Vec3&v, real a) {
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}

//------ Mat3 ---------------------

blast_fn Mat3::Mat3(const Mat3& m) {
    memcpy(data, m.data, 9*sizeof(real));
}

blast_fn Mat3::Mat3(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3) {
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

blast_fn real& Mat3::operator()(u32 row, u32 col) {
    return data[3*col + row];
}

blast_fn Mat3 operator*(const Mat3& m, const Mat3 rhs) {
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

blast_fn Mat3& operator*=(Mat3& lhs, const Mat3& rhs) {
    lhs = lhs * rhs;
    return lhs;
}

blast_fn Mat3 operator+(const Mat3& lhs, const Mat3& rhs) {
    Mat3 r;
    for (u32 i = 0; i < 9; i++)
        r.data[i] = lhs.data[i] + rhs.data[i];
    return r;
}

blast_fn Mat3& operator+=(Mat3& lhs, const Mat3& rhs) {
    lhs = lhs + rhs;
    return lhs;
}

blast_fn Mat3 transpose(const Mat3& m) {
    Mat3 result {
        m.data[0],
        m.data[3],
        m.data[6],
        m.data[1],
        m.data[4],
        m.data[7],
        m.data[2],
        m.data[5],
        m.data[8]
    };
    return result;
}

blast_fn real& Mat3::operator[](u32 i) {
    Assert(i < 9);
    return data[i];
}

blast_fn real Mat3::operator[](u32 i) const {
    Assert(i < 9);
    return data[i];
}

blast_fn Vec3 Mat3::col(u32 c) {
    Assert(c < 3);
    return Vec3(data[c*3+0], data[c*3+1], data[c*3+2]);
}

//------ Mat4 ---------------------

blast_fn Mat4::Mat4(const Mat4& m) {
    memcpy(data, m.data, 16*sizeof(real));
}

blast_fn Mat4::Mat4(real v) {
    zero();
    data[0] = v;
    data[4] = v;
    data[8] = v;
    data[12] = v;
}

blast_fn void Mat4::zero() {
    memset(data, 0, 16*sizeof(real));
}
// todo: other zeros in misc ?

blast_fn Mat4::Mat4(const std::initializer_list<real>& l) {
    u32 i = 0;
    memset(data, 0, sizeof(*this));
    for (i = 0; i < l.size() && i < 16; i++)
        data[i] = l.begin()[i];
}

blast_fn Mat4& Mat4::operator=(const Mat4& m) {
    memcpy(data, m.data, 16*sizeof(real));
    return *this;
}

blast_fn Mat4& Mat4::operator=(const std::initializer_list<real>& l) {
    u32 i;
    for (i = 0; i < l.size() && i < 16; i++)
        data[i] = l.begin()[i];
    for (; i < 16; i++)
        data[i]=(real)0.0;
    return *this;
}

blast_fn real& Mat4::operator()(u32 row, u32 col) {
    Assert(row<4 && col<4);
    return data[4*col + row];
}

blast_fn real Mat4::operator()(u32 row, u32 col) const {
    Assert(row<4 && col<4);
    return data[4*col + row];
}

blast_fn real& Mat4::operator[](u32 i) {
    Assert(i < 16);
    return data[i];
}

blast_fn real Mat4::operator[](u32 i) const {
    Assert(i < 16);
    return data[i];
}

blast_fn Mat4 Mat4::operator-() {
    Mat4 r;
    for (u32 i = 0; i < 16; i++)
        r.data[i] = -data[i];
    return r;
}

blast_fn Mat4& Mat4::operator*=(Mat4& rhs) {
    // naive implementation (no SIMD)
    Mat4 tmp = *this;
    for (u32 i = 0; i < 4; i++) {
        for (u32 j = 0; j < 4; j++) {
            (*this)(i, j) = tmp(i, 0)*rhs(0, j) + tmp(i, 1)*rhs(1, j) + tmp(i, 2)*rhs(2, j) + tmp(i, 3)*rhs(3, j);
        }
    }
    return *this;

}

blast_fn Mat4& Mat4::operator*=(real v) {
    for (int i = 0; i < 16; i++)
        data[i] *= v;
    return *this;
}

blast_fn Mat4& Mat4::operator+=(Mat4& rhs) {
    for (u32 i = 0; i < 16; i++)
        data[i] += rhs[i];
    return *this;
}

blast_fn Mat4& Mat4::operator-=(Mat4& rhs) {
    for (u32 i = 0; i < 16; i++)
        data[i] -= rhs[i];
    return *this;
}

blast_fn bool Mat4::operator==(const Mat4& rhs) {
    for (u32 i = 0; i < 16; i++)
        if (data[i] != rhs.data[i])
            return false;
    return true;
}

blast_fn Mat4 operator*(const Mat4& lhs, const Mat4& rhs) {
    Mat4 result;
    // naive implementation (no SIMD)
    //  note: order of for loops is important!
    for (u32 j = 0; j < 4; j++) {
        for (u32 i = 0; i < 4; i++) {
            result(i, j) = lhs(i, 0)*rhs(0, j) + lhs(i, 1)*rhs(1, j) + lhs(i, 2)*rhs(2, j) + lhs(i, 3)*rhs(3, j);
        }
    }
    return result;
}

blast_fn Array Mat4::col(u32 c) {
    Assert(c < 4);
    Array result;
    result.data = data + 4*c;
    result.size = 4;
    result.is_alias = true;
    return result;
}

blast_fn Array Mat4::col(u32 c) const {
    Assert(c < 4);
    Array result(4);
    memcpy(result.data, data + c*4, 4*sizeof(real));
    return result;
}

//------ Array ---------------------

blast_fn Array::Array(u32 new_size) : size(new_size) {
    if (size) {
        data = (real*)malloc(size * sizeof(real));
        memset(data, 0, size * sizeof(real));
    }
}

blast_fn Array::Array(u32 n, real value) : size(n) {
    if (n) {
        data = (real*)malloc(size * sizeof(real));
        for (u32 i = 0; i < n; i++)
            data[i] = value;
    }
}

blast_fn Array::Array(const Array& a) : size(a.size) {
    if (size) {
        data = (real*)malloc(size * sizeof(real));
        memcpy(data, a.data, size*sizeof(real));
    }
}

blast_fn Array::Array(Array&& a) : data(a.data), size(a.size), is_alias(a.is_alias) {
    a.is_alias = true;
}

blast_fn Array::Array(real* d, u32 n) {
    data = d;
    size = n;
    is_alias = true;
}

host_fn Array::Array(const svector& v) {
    size = (u32)v.size();
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, v.data(), size*sizeof(real));
    }
}

blast_fn Array::~Array() {
    if (!is_alias && data) {
        free(data);
    }
}

blast_fn Array& Array::operator=(const Array& a) {
    if (this != &a) {
        if (data && !is_alias)
            free(data);
        size = a.size;
        if (size) {
            data = (real*)malloc(a.size * sizeof(real));
            memcpy(data, a.data, size*sizeof(real));
        }
        is_alias = false;
    }
    return *this;
}

blast_fn Array& Array::operator=(Array&& a) {
    if (this != &a) {
        if (data && !is_alias)
            std::free(data);
        data = a.data;
        size = a.size;
        is_alias = a.is_alias;
        a.data = nullptr;
        a.size = 0;
    }
    return *this;
}

blast_fn Array& Array::operator=(const std::initializer_list<real>& other) {
    Assert(other.size() <= size);
    memcpy(data, other.begin(), other.size() * sizeof(real));
    return *this;
}

blast_fn Array Array::operator-() {
    Array result(size);
    for (u32 i = 0; i < size; i++)
        result[i] = -data[i];
    return std::move(result);
}

blast_fn bool Array::operator==(Array& a) {
    Assert(size == a.size);
    return is_close(*this, a);
}

blast_fn Array& Array::operator*=(real n) {
    for (u32 i = 0; i < size; i++)
        data[i] *= n;
    return *this;
}

blast_fn real& Array::operator[](u32 i) {
    Assert(i < size);
    return data[i];
}

blast_fn real Array::operator[](u32 i) const {
    Assert(i < size);
    return data[i];
}

host_fn Array& Array::alias(svector& v) {
    if (data && !is_alias)
        free(data);
    data = v.data();
    size = (u32)v.size();
    is_alias = true;
    return *this;
}

blast_fn Array& Array::alias(Matrix& m) {
    if (data && !is_alias)
        free(data);
    data = m.data;
    size = m.size;
    is_alias = true;
    return *this;
}

blast_fn Array& Array::alias(real* p, u32 n) {
    Assert(p);
    if (data && !is_alias)
        free(data);
    data = p;
    size = n;
    is_alias = true;
    return *this;
}

blast_fn Array& Array::alias(const real* p, u32 n) {
    Assert(p);
    Assert(data == nullptr);

    data = const_cast<real*>(p);
    size = n;
    is_alias = true;
    return *this;
}

blast_fn void Array::resize(u32 new_size) {
    Assert(!is_alias);
    if (new_size > size) {
        data = (real*)realloc(data, new_size*sizeof(real));
        memset(data + size, 0, (new_size-size)*sizeof(real));
    }
    size = new_size;
}

blast_fn real& Array::back() {
    Assert(size);
    return data[size-1];
}

blast_fn real Array::back() const {
    Assert(size);
    return data[size-1];
}

blast_fn Vec3 operator*(const Mat3& m, const Vec3 v) {
    Vec3 r;
    r.x = m.data[0]*v.x + m.data[3]*v.y + m.data[6]*v.z;
    r.y = m.data[1]*v.x + m.data[4]*v.y + m.data[7]*v.z;
    r.z = m.data[2]*v.x + m.data[5]*v.y + m.data[8]*v.z;
    return r;
}

blast_fn Array operator-(const Array& v1, const Array& v2) {
    Array r = v1;
    for (u32 i = 0; i < v1.size; i++)
        r[i] -= v2[i];
    return r;
}

blast_fn Array operator+(const Array& v1, const Array& v2) {
    Assert(v1.size == v2.size);
    Array r = v1;
    for (u32 i = 0; i < v1.size; i++)
        r[i] += v2[i];
    return r;
}

blast_fn Array operator*(const Array& a, real b) {
    Array r(a);
    r *= b;
    return r;
}

blast_fn Array operator*(real b, const Array& a) {
    Array r(a);
    r *= b;
    return r;
}

// Compute the dot product of the given arrays
//  - note: fastest when the number of elements are a factor of 4 or even 8
blast_fn real dot(Array& a, Array& b) {
    Assert(a.size == b.size);
    real r = 0;
#if defined(__CUDA_ARCH__)
    for (int i = 0; i < (int)a.size; i++)
        r += a[i] * b[i];
#else
    int i = 0;

    // SIMD what you can
    mipp::Reg<real> ra; // wide register for part of the 'a' array
    mipp::Reg<real> rb; // wide register for part of the 'b' array
    mipp::Reg<real> accum = 0.0;
    auto vecLoopSize = (a.size / mipp::N<real>()) * mipp::N<real>();
    for (; i < vecLoopSize; i += mipp::N<real>()) {
        ra.load(&a[i]);
        rb.load(&b[i]);
        accum = mipp::fmadd(ra, rb, accum);
    }
    r = accum.hadd();
    for (; i < (int)a.size; i++)
        r += a[i] * b[i];
#endif
    return r;
}

// Compute the sine and the cosine of every element
//  - note: fastest when the number of elements are a factor of 4 (or even 8 if real is float)
blast_fn void sincos(const Array& angles, Array& sines, Array& cosines) {
    Assert(angles.size == sines.size && angles.size == cosines.size);

#if defined(__CUDA_ARCH__)
#if BLAST_USE_DOUBLES
    for (int i = 0; i < angles.size; i++)
        ::sincos(angles[i], &sines[i], &cosines[i]);
#else
    for (int i = 0; i < angles.size; i++)
        ::sincosf(angles[i], &sines[i], &cosines[i]);
#endif
#elif defined(__GNUC__)
    int i = 0;
    for (; i < (int)angles.size; i++) {
        sines[i] = sin(angles[i]);
        cosines[i] = cos(angles[i]);
    }
#else
    int i = 0;
    // simd what you can
    mipp::Reg<real> a;
    // mipp::Reg<real> c;
    // mipp::Reg<real> s;
    auto vecLoopSize = (angles.size / mipp::N<real>()) * mipp::N<real>();
    for (; i < vecLoopSize; i += mipp::N<real>()) {
        a.load(angles.data + i);
        auto r = mipp::sincos(a);
        r[0].store(sines.data + i);
        r[1].store(cosines.data + i);
    }
    // serialize the rest
    for (; i < (int)angles.size; i++) {
        sines[i] = sin(angles[i]);
        cosines[i] = cos(angles[i]);
    }
#endif
}

//------ Matrix ---------------------

blast_fn Matrix::Matrix(u32 r, u32 c) {
    size = r * c;
    rows = r;
    cols = c;
    if (size) {
        data = (real*)malloc(size * sizeof(real));
        memset(data, 0, size * sizeof(real));
    }
}

blast_fn Matrix::Matrix(const Matrix& m) : size(m.size), cols(m.cols), rows(m.rows), is_alias(false) {
    if (size) {
        data = (real*)malloc(size * sizeof(real));
        memcpy(data, m.data, size*sizeof(real));
    }
}

blast_fn Matrix::Matrix(Matrix&& m) : data(m.data), size(m.size), cols(m.cols), rows(m.rows), is_alias(m.is_alias) {
    m.data = nullptr;
    m.size = 0;
    m.rows = 0;
    m.cols = 0;
}

blast_fn Matrix::Matrix(real* d, u32 r, u32 c) {
    data = d;
    rows = r;
    cols = c;
    size = r*c;
    is_alias = true;
}

blast_fn Matrix::Matrix(const Array& v) : size(v.size), cols(1), rows(v.size), is_alias(false) {
    if (size) {
        data = (real*)malloc(size * sizeof(real));
        memcpy(data, v.data, size*sizeof(real));
    }
}

blast_fn Matrix::~Matrix() {
    if (!is_alias && data)
        free(data);
}

host_fn void Matrix::resize(u32 r, u32 c) {
    Assert(!is_alias);

    if (data == nullptr)
        data = (real*)malloc(r*c * sizeof(real));
    else if (size < r*c)
        data = (real*)realloc(data, r*c * sizeof(real));
    size = r*c;
    rows = r;
    cols = c;
    Assert(data);
}

blast_fn Matrix& Matrix::operator=(const Matrix& m) {
    if (this != &m) {
        if (m.size == 0) {
            size = 0;
            rows = 0;
            cols = 0;
            if (data)
                free(data);
            data = nullptr;
            return *this;
        }
        else if (data == nullptr || is_alias) {
            data = (real*)malloc(m.size * sizeof(real));
        }
        else {
            free(data);
            data = (real*)malloc(m.size * sizeof(real));
        }
        Assert(data);
        size = m.size;
        cols = m.cols;
        rows = m.rows;
        is_alias = false;
        memcpy(data, m.data, size*sizeof(real));
    }
    return *this;
}

blast_fn Matrix& Matrix::operator=(Matrix&& m) {
    if (this != &m) {
        if (data && !is_alias)
            free(data);
        data = m.data;
        size = m.size;
        rows = m.rows;
        cols = m.cols;
        is_alias = m.is_alias;
        m.is_alias = true; // m is about to be destroid but we don't want it to free the memory
    }
    return *this;
}

blast_fn Matrix& Matrix::alias(Array& a) {
    Assert(a.data);
    if (data && !is_alias)
        free(data);
    size = a.size;
    rows = size;
    cols = 1;
    data = a.data;
    is_alias = true;
    return *this;
}

host_fn Matrix& Matrix::alias(std::vector<real>& v) {
    Assert(!v.empty());
    if (data && !is_alias)
        free(data);
    size = (u32)v.size();
    rows = size;
    cols = 1;
    data = v.data();
    is_alias = true;
    return *this;
}

blast_fn Matrix& Matrix::alias(real* p, u32 r, u32 c) {
    Assert(p);
    if (data && !is_alias)
        free(data);
    size = r*c;
    rows = r;
    cols = c;
    data = p;
    is_alias = true;
    return *this;
}

blast_fn real& Matrix::operator()(u32 row, u32 col) {
    Assert(row < rows && col < cols);
    return data[row + rows*col];
}

blast_fn real Matrix::operator()(u32 row, u32 col) const {
    Assert(row < rows && col < cols);
    return data[row + rows*col];
}

blast_fn real* Matrix::address(u32 row, u32 col) const {
    Assert(row < rows && col < cols);
    return &data[row + rows*col];
}

blast_fn Array Matrix::col(u32 c) const {
    Assert(c < this->cols);
    Array result;
    result.data = data + rows*c;
    result.size = rows;
    result.is_alias = true;
    return result;
}

blast_fn Matrix operator+(const Matrix& m1, const Matrix& m2) {
    Assert(m1.cols == m2.cols && m1.rows == m2.rows);
    Matrix result(m1.rows, m1.cols);
    for (u32 i = 0; i < m1.size; i++)
        result.data[i] = m1.data[i] + m2.data[i];

    return result;
}

blast_fn Matrix operator-(const Matrix& m1, const Matrix& m2) {
    Assert(m1.cols == m2.cols && m1.rows == m2.rows);
    Matrix result(m1.rows, m1.cols);
    for (u32 i = 0; i < m1.size; i++)
        result.data[i] = m1.data[i] - m2.data[i];
    return result;
}

blast_fn Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
    Assert(lhs.cols == rhs.rows);
    Matrix result(lhs.rows, rhs.cols);
    for (u32 i = 0; i < lhs.rows; i++) {
        for (u32 j = 0; j < rhs.cols; j++) {
            real sum = 0;
            for (u32 k = 0; k < lhs.cols; k++) {
                sum += lhs(i, k) * rhs(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

blast_fn Array operator*(const Matrix& m, const Array& v) {
    Assert(m.cols == v.size);
    Array result(m.rows);
    for (u32 i = 0; i < m.rows; i++) {
        real sum = 0;
        for (u32 j = 0; j < m.cols; j++) {
            sum += m(i, j) * v[j];
        }
        result[i] = sum;
    }
    return result;
}

blast_fn Matrix operator*(const Matrix& m, real& r) {
    Matrix result(m.rows, m.cols);
    for (u32 i = 0; i < m.size; i++)
        result.data[i] = m.data[i]*r;
    return result;
}

blast_fn Matrix pw_mult(const Matrix& m1, const Matrix& m2) {
    Matrix result(m1.rows, m1.cols);
    Assert(m1.size == m2.size);
    for (u32 i = 0; i < m1.size; i++)
        result.data[i] = m1.data[i] * m2.data[i];

    return result;
}

blast_fn Matrix transpose(const Matrix& m) {
    Matrix result(m.cols, m.rows);
    for (u32 i = 0; i < m.rows; ++i) {
        for (u32 j = 0; j < m.cols; ++j) {
            result(j, i) = m(i, j);
        }
    }
    return result;
}

blast_fn real determinant(const Matrix & m) {
    int n = m.rows;
    Matrix LU = m;

    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            LU(i, k) /= LU(k, k);
            for (int j = k + 1; j < n; j++) {
                LU(i, j) -= LU(i, k) * LU(k, j);
            }
        }
    }

    real det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= LU(i, i);
    }

    return det;
}

blast_fn Matrix LU_decomp(const Matrix& m) {
    int n = m.rows;
    Matrix LU = m;

    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            LU(i, k) /= LU(k, k);
            for (int j = k + 1; j < n; j++) {
                LU(i, j) -= LU(i, k) * LU(k, j);
            }
        }
    }
    return LU;
}

blast_fn Array solveLU(const Matrix& LU, const Array& b) {
    int n = LU.rows;
    Array x(n, 0.0);

    for (int i = 0; i < n; i++) {
        real sum = 0.0;
        for (int j = 0; j < i; j++)
            sum += LU(i, j) * x[j];

        x[i] = b[i] - sum;
    }

    // Solve Ux = y (backsubstitution)
    for (int i = n - 1; i >= 0; i--) {
        real sum = 0.0;
        for (int j = i + 1; j < n; j++)
            sum += LU(i, j) * x[j];

        x[i] = (x[i] - sum) / LU(i, i);
    }

    return x;
}

blast_fn Matrix inverse(const Matrix& m) {
    Assert(m.rows == m.cols); // square matrix

    int n = m.rows;
    Matrix LU = LU_decomp(m);
    Matrix inv(n, n);

    for (int i = 0; i < n; i++) {
        // unit vector in one dimension
        Array unit(n);
        unit[i] = 1.0;

        // Solve for the i-th column of the inverse
        Array column = solveLU(LU, unit);

        // insert the column into the inverse matrix
        inv.col(i) = column;
    }
    return inv;
}

blast_fn Matrix pinv(const Matrix& m) {
    Matrix res(m.cols, m.rows);
    if (m.cols >= m.rows) {
        res = transpose(m) * inverse(m * transpose(m));
    }
    else {
        res = inverse(transpose(m) * m) * transpose(m);
    }
    return res;
}

blast_fn Matrix eye(const u32 s) {
    Matrix result(s, s);
    for (u32 i = 0; i < s; ++i) {
        for (u32 j = 0; j < s; ++j) {
            result(j, i) = i == j ? 1 : 0;
        }
    }
    return result;
}

//------ Collision ---------------------

blast_fn real clamped_root(real slope, real h0, real h1) {
//note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    real r;
    if (h0 < 0) {
        if (h1 > 0) {
            r = -h0 / slope;
            if (r > 1)
                r = 0.5;
        }
        else
            r = 1;
    }
    else
        r = 0;
    return r;
}

blast_fn void compute_intersection(real* sValue, i32* classify, real b, real f00, real f10, i32* edge, real end[][2]) {
// note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    real const zero = 0;
    real const half = (real)0.5;
    real const one = 1;
    if (classify[0] < 0) {
        edge[0] = 0;
        end[0][0] = zero;
        end[0][1] = f00 / b;
        if (end[0][1] < zero || end[0][1] > one)
            end[0][1] = half;
        if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
    }
    else if (classify[0] == 0) {
        edge[0] = 2;
        end[0][0] = sValue[0];
        end[0][1] = zero;
        if (classify[1] < 0) {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
        else if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
    }
    else {
        edge[0] = 1;
        end[0][0] = one;
        end[0][1] = f10 / b;
        if (end[0][1] < zero || end[0][1] > one)
            end[0][1] = half;
        if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
    }
}

blast_fn void compute_minimum_parameters(i32* edge, real end[][2], real b, real c, real e, real g00, real g10, real g01, real g11, real* parameter) {
// note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    real const zero = 0;
    real const one = 1;
    real const delta = end[1][1] - end[0][1];
    real h0 = delta * (-b * end[0][0] + c * end[0][1] - e);
    if (h0 >= zero) {
        if (edge[0] == 0) {
            parameter[0] = zero;
            parameter[1] = clamped_root(c, g00, g01);
        }
        else if (edge[0] == 1) {
            parameter[0] = one;
            parameter[1] = clamped_root(c, g10, g11);
        }
        else {
            parameter[0] = end[0][0];
            parameter[1] = end[0][1];
        }
    }
    else {
        real h1 = delta * (-b * end[1][0] + c * end[1][1] - e);
        if (h1 <= zero) {
            if (edge[1] == 0) {
                parameter[0] = zero;
                parameter[1] = clamped_root(c, g00, g01);
            }
            else if (edge[1] == 1) {
                parameter[0] = one;
                parameter[1] = clamped_root(c, g10, g11);
            }
            else {
                parameter[0] = end[1][0];
                parameter[1] = end[1][1];
            }
        }
        else {
            real z = clamp( h0/(h0 - h1), 0, 1 );
            real omz = one - z;
            parameter[0] = omz * end[0][0] + z * end[1][0];
            parameter[1] = omz * end[0][1] + z * end[1][1];
        }
    }
}

blast_fn real two_segment_distance_sqr(Vec3 P0, Vec3 P1, Vec3 Q0, Vec3 Q1) {
// note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    auto const P1mP0 = P1 - P0;
    auto const Q1mQ0 = Q1 - Q0;
    auto const P0mQ0 = P0 - Q0;
    real a = dot(P1mP0, P1mP0);
    real b = dot(P1mP0, Q1mQ0);
    real c = dot(Q1mQ0, Q1mQ0);
    real d = dot(P1mP0, P0mQ0);
    real e = dot(Q1mQ0, P0mQ0);
    real f00 = d;
    real f10 = f00 + a;
    real f01 = f00 - b;
    real f11 = f10 - b;
    real g00 = -e;
    real g10 = g00 - b;
    real g01 = g00 + c;
    real g11 = g10 + c;
    real parameter[2] = {0, 0};
    if (a > 0 && c > 0) {
        real sValue[2] = {
            clamped_root(a, f00, f10),
            clamped_root(a, f01, f11)
        };
        i32 classify[2] = {0, 0};
        for (size_t i = 0; i < 2; ++i) {
            if (sValue[i] <= 0)
                classify[i] = -1;
            else if (sValue[i] >= 1)
                classify[i] = 1;
            else
                classify[i] = 0;
        }
        if (classify[0] == -1 && classify[1] == -1) {
            parameter[0] = 0;
            parameter[1] = clamped_root(c, g00, g01);
        }
        else if (classify[0] == +1 && classify[1] == +1) {
            parameter[0] = 1;
            parameter[1] = clamped_root(c, g10, g11);
        }
        else {
            i32 edge[2] = { 0, 0 };
            real end[2][2];
            compute_intersection(sValue, classify, b, f00, f10, edge, end);
            compute_minimum_parameters(edge, end, b, c, e, g00, g10, g01, g11, parameter);
        }
    }
    else    {
        if (a > 0) {
            parameter[0] = clamped_root(a, f00, f10);
            parameter[1] = 0;
        }
        else     if (c > 0) {
            parameter[0] = 0;
            parameter[1] = clamped_root(c, g00, g01);
        }
        else {
            parameter[0] = 0;
            parameter[1] = 0;
        }
    }
    Vec3 closest0 = P0 + parameter[0]*P1mP0;
    Vec3 closest1 = Q0 + parameter[1]*Q1mQ0;
    Vec3 diff = closest0 - closest1;

    // auto result = sqrt(dot(diff, diff));
    return dot(diff, diff);
}

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
TEST_CASE("Arrays", "[Math]") {
    using namespace blast;
    SECTION("Basic operations") {
        Array a(8);
        a = { 1, 1, 1, 1, 1, 1, 1, 1 };
        Array b(8);
        b = { 1, 1, 1, 1, 1, 1, 2, 2 };
        auto r = dot(a, b);
        REQUIRE(r == (real)10.0);
    }

    SECTION("Correct aliasing functionnality") {
        Array a({1, 2, 3, 4, 5});
        {
            Array b(a);
            Array c(a.data, 4);
            b[1] = 12;
            REQUIRE(float(a[1]) == 2.0f);
            c[1] = 14;
            REQUIRE(float(a[1]) == 14.0f);
            REQUIRE_FALSE(a.is_alias);
            REQUIRE_FALSE(b.is_alias);
            REQUIRE(c.is_alias);
        }
    }

    SECTION("Correct Vec3 addressing") {
        Vec3 a = {1, 2, 3};
        REQUIRE(a.x == a[0]);
        REQUIRE(a.y == a[1]);
        REQUIRE(a.z == a[2]);
        REQUIRE(a.x == 1);
        REQUIRE(a.y == 2);
        REQUIRE(a.z == 3);
    }
}

TEST_CASE("Mat3", "[Math]") {
    using namespace blast;

    // mult
    {
        Mat3 m1;  fill_random(m1, 10);
        Mat3 m2;  fill_random(m2, 5);

        Mat3 m3;
        for (u32 i = 0; i < 3; i++) {
            for (u32 j = 0; j < 3; j++) {
                m3(i, j) = m1(i, 0)*m2(0, j) + m1(i, 1)*m2(1, j) + m1(i, 2)*m2(2, j);
            }
        }

        Mat3 m4 = m1 * m2;
        m1 *= m2;

        // m1, m3 and m4 should now all be the same
        for (u32 i = 0; i < 9; i++) {
            REQUIRE((float)m1[i] == (float)m3[i]);
            REQUIRE((float)m1[i] == (float)m4[i]);
        }
    }

    // add
    {
        Mat3 m1;
        Mat3 m2;
        for (u32 i = 0; i < 9; i++) {
            m1[i] = 10*get_random();
            m2[i] = 5*get_random();
        }

        Mat3 m3;
        for (u32 i = 0; i < 3; i++) {
            for (u32 j = 0; j < 3; j++) {
                m3(i, j) = m1(i, j) + m2(i, j);
            }
        }

        Mat3 m4 = m1 + m2;
        m1 += m2;

        // m1, m3 and m4 should now all be the same
        for (u32 i = 0; i < 9; i++) {
            REQUIRE((float)m1[i] == (float)m3[i]);
            REQUIRE((float)m1[i] == (float)m4[i]);
        }
    }

    // col
    {
        Mat3 m;
        for (u32 i = 0; i < 9; i++) {
            m[i] = 10*get_random();
        }

        Vec3 v[3];
        v[0] = m.col(0);
        v[1] = m.col(1);
        v[2] = m.col(2);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                REQUIRE(m(i, j) == v[j][i]);
            }
        }
    }
}

TEST_CASE("Mat4Operations", "[Math]") {
    using namespace blast;

    Mat4 m1;
    Mat4 m2;
    for (u32 i = 0; i < 16; i++) {
        m1[i] = 10*get_random();
        m2[i] = 5*get_random();
    }

    Mat4 m3;
    for (u32 i = 0; i < 4; i++) {
        for (u32 j = 0; j < 4; j++) {
            m3(i, j) = m1(i, 0)*m2(0, j) + m1(i, 1)*m2(1, j) + m1(i, 2)*m2(2, j) + m1(i, 3)*m2(3, j);
        }
    }

    Mat4 m4 = m1 * m2;
    m1 *= m2;

    // m1, m3 and m4 should now all be the same
    for (u32 i = 0; i < 16; i++) {
        REQUIRE((float)m1[i] == (float)m3[i]);
        REQUIRE((float)m1[i] == (float)m4[i]);
    }
}

TEST_CASE("TwoSegmentDist", "[Math]") {
    auto dist_sqr_test1 = blast::two_segment_distance_sqr({1, 1, 1}, {2, 2, 2}, {1, 0, 0}, {2, 0, 0});
    auto dist_sqr_test2 = blast::two_segment_distance_sqr({1, 1, 1}, {2, 2, 2}, {1, 0, 0}, {2, 3, 2});
    auto dist_sqr_test3 = blast::two_segment_distance_sqr({1.5, 3, 2}, {7, 0, 4}, {8.2, 0, 5}, {2, 2, 0});
    auto dist_sqr_test4 = blast::two_segment_distance_sqr({1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4});

    REQUIRE((float)dist_sqr_test1 == 2.f);
    REQUIRE((float)dist_sqr_test2 == 0.16666666666f);
    REQUIRE((float)dist_sqr_test3 == 0.07709516434f);
    REQUIRE((float)dist_sqr_test4 == 3.f);
}

TEST_CASE("MatrixOperations", "[Math]") {
    using namespace blast;

    SECTION("transpose") {
        Matrix m(3, 2);
        for (int i = 0; i < m.size; i++)
            m.data[i] = i;
        auto mT = transpose(m);

        REQUIRE(mT(0, 0) == m(0, 0));
        REQUIRE(mT(0, 1) == m(1, 0));
        REQUIRE(mT(0, 2) == m(2, 0));
        REQUIRE(mT(1, 0) == m(0, 1));
        REQUIRE(mT(1, 1) == m(1, 1));
        REQUIRE(mT(1, 2) == m(2, 1));
    }

    SECTION("Addition & substraction & real multiplication") {
        Matrix m1(3, 3);
        Matrix m2(3, 3);
        real r = 2;
        for (int i = 0; i < m1.size; i++) {
            m1.data[i] = i;
            m2.data[i] = i;
        }
        auto m3 = m1 + m2;
        auto m4 = m1 - m2;

        REQUIRE(m3.cols == 3);
        REQUIRE(m3.rows == 3);
        REQUIRE(m3(0, 0) == 0);
        REQUIRE(m3(0, 1) == 6);
        REQUIRE(m3(0, 2) == 12);
        REQUIRE(m3(1, 0) == 2);
        REQUIRE(m3(1, 1) == 8);
        REQUIRE(m3(1, 2) == 14);
        REQUIRE(m3(2, 0) == 4);
        REQUIRE(m3(2, 1) == 10);
        REQUIRE(m3(2, 2) == 16);

        REQUIRE(m4.cols == 3);
        REQUIRE(m4.rows == 3);
        REQUIRE(m4(0, 0) == 0);
        REQUIRE(m4(0, 1) == 0);
        REQUIRE(m4(0, 2) == 0);
        REQUIRE(m4(1, 0) == 0);
        REQUIRE(m4(1, 1) == 0);
        REQUIRE(m4(1, 2) == 0);
        REQUIRE(m4(2, 0) == 0);
        REQUIRE(m4(2, 1) == 0);
        REQUIRE(m4(2, 2) == 0);

        auto m5 = m3*r;
        REQUIRE(m5.cols == 3);
        REQUIRE(m5.rows == 3);
        REQUIRE(m5(0, 0) == 0);
        REQUIRE(m5(0, 1) == 12);
        REQUIRE(m5(0, 2) == 24);
        REQUIRE(m5(1, 0) == 4);
        REQUIRE(m5(1, 1) == 16);
        REQUIRE(m5(1, 2) == 28);
        REQUIRE(m5(2, 0) == 8);
        REQUIRE(m5(2, 1) == 20);
        REQUIRE(m5(2, 2) == 32);
    }

    SECTION("Multiplication") {
        Matrix m2(2, 3);
        Matrix m(3, 2);
        for (int i = 0; i < m.size; i++)
            m.data[i] = i;
        for (int i = 0; i < m2.size; i++)
            m2.data[i] = i+1;

        auto m3 = m * m2;
        REQUIRE(m3.cols == 3);
        REQUIRE(m3.rows == 3);
        REQUIRE(m3(0, 0) == 6);
        REQUIRE(m3(0, 1) == 12);
        REQUIRE(m3(0, 2) == 18);
        REQUIRE(m3(1, 0) == 9);
        REQUIRE(m3(1, 1) == 19);
        REQUIRE(m3(1, 2) == 29);
        REQUIRE(m3(2, 0) == 12);
        REQUIRE(m3(2, 1) == 26);
        REQUIRE(m3(2, 2) == 40);
    }

    SECTION("Piecewies multiplication") {
        Matrix m1(3, 3);
        Matrix m2(3, 3);
        for (int i = 0; i < m1.size; i++) {
            m1.data[i] = i;
            m2.data[i] = i;
        }

        auto m3 = pw_mult(m1, m2);
        REQUIRE(m3.cols == 3);
        REQUIRE(m3.rows == 3);
        REQUIRE(m3(0, 0) == 0);
        REQUIRE(m3(0, 1) == 9);
        REQUIRE(m3(0, 2) == 36);
        REQUIRE(m3(1, 0) == 1);
        REQUIRE(m3(1, 1) == 16);
        REQUIRE(m3(1, 2) == 49);
        REQUIRE(m3(2, 0) == 4);
        REQUIRE(m3(2, 1) == 25);
        REQUIRE(m3(2, 2) == 64);
    }
    SECTION("Identity matrix") {
        Matrix identity;
        identity = eye(4);

        REQUIRE(identity(0, 0) == 1);
        REQUIRE(identity(0, 1) == 0);
        REQUIRE(identity(0, 2) == 0);
        REQUIRE(identity(0, 3) == 0);
        REQUIRE(identity(1, 0) == 0);
        REQUIRE(identity(1, 1) == 1);
        REQUIRE(identity(1, 2) == 0);
        REQUIRE(identity(1, 3) == 0);
        REQUIRE(identity(2, 0) == 0);
        REQUIRE(identity(2, 1) == 0);
        REQUIRE(identity(2, 2) == 1);
        REQUIRE(identity(2, 3) == 0);
        REQUIRE(identity(3, 0) == 0);
        REQUIRE(identity(3, 1) == 0);
        REQUIRE(identity(3, 2) == 0);
        REQUIRE(identity(3, 3) == 1);
    }
}

TEST_CASE("Min/Max") {
    using namespace blast;

    SECTION("min/max array") {
        Array a(4);
        a = {1, 2, 3, 0};

        auto a_min = array_min(a);
        REQUIRE(a_min == 0);
        auto i_a_min = array_min_id(a);
        REQUIRE(i_a_min == 3);
        auto a_max = array_max(a);
        REQUIRE(a_max == 3);
    }

    SECTION("min/max matrix") {
        Matrix m(2, 2);
        m(0, 0) = 1;
        m(0, 1) = 2;
        m(1, 0) = 3;
        m(1, 1) = 0;

        auto m_min = matrix_min(m);
        REQUIRE(m_min == 0);
        auto i_m_min = matrix_min_id(m);
        REQUIRE(i_m_min == 3);
        auto m_max = matrix_max(m);
        REQUIRE(m_max == 3);
    }
}
#endif
