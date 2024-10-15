#pragma once

#include <immintrin.h>
#include <limits>
#include <vector>

#ifdef __CUDA_ARCH__
#include <cuda_runtime.h>
#include <math_constants.h>
#endif

namespace blast {

const u32 ALIGN = 64;

// Types declared in this file
struct Vec3;
struct Mat3;
struct Mat4;
struct Array;
struct Matrix;


// uses doubles by default unless BLAST_USE_DOUBLES is set to 0
#ifndef BLAST_USE_DOUBLES
#define BLAST_USE_DOUBLES 1
#endif
#if BLAST_USE_DOUBLES
using real  = double;
#define BLAST_SIZEOF_REAL 8
#else
using real  = float;
#define BLAST_SIZEOF_REAL 4
#endif
static_assert(sizeof(real) == BLAST_SIZEOF_REAL);


// Constants
#ifdef __CUDA_ARCH__
#if BLAST_USE_DOUBLES
const real NAN_REAL = CUDART_NAN;
const real INF_REAL = CUDART_INF;
#else
const real NAN_REAL = CUDART_NAN_F;
const real INF_REAL = CUDART_INF_F;
#endif
#else
const real NAN_REAL = std::numeric_limits<real>::quiet_NaN();
const real INF_REAL = std::numeric_limits<real>::infinity();
#endif
const real PI = (real)3.1415;


// Type aliases
using svector = std::vector<real>;


// 3D vector
struct Vec3 {
    real x = 0;
    real y = 0;
    real z = 0;
    real _pad = 0;
    Vec3() = default;
    blast_fn Vec3(real x, real y, real z);
    blast_fn real& operator[](int i);
    blast_fn bool operator==(const Vec3&) const;
};
blast_fn Vec3  cross(Vec3, Vec3);
blast_fn Vec3& zero(Vec3&);
blast_fn Vec3& constant(Vec3&, real value);
blast_fn real  dot(Vec3, Vec3);
blast_fn real  norm(Vec3);
blast_fn Vec3  operator-(Vec3);
blast_fn Vec3  operator-(Vec3, Vec3);
blast_fn Vec3  operator+(Vec3, Vec3);
blast_fn Vec3  operator*(real, Vec3);
blast_fn Vec3  operator*(Vec3, real);
blast_fn Vec3  operator/(Vec3, real);
blast_fn Vec3& operator+=(Vec3&, const Vec3&);
blast_fn Vec3& operator*=(Vec3&, real);


// 3x3 matrix
struct Mat3 {
    real data[9] = {0};
    Mat3() = default;
    blast_fn Mat3(const Mat3& m);
    blast_fn Mat3(real x1, real y1, real z1,
                  real x2, real y2, real z2,
                  real x3, real y3, real z3);
    blast_fn real& operator()(u32 row, u32 col);
    blast_fn real  operator()(u32 row, u32 col) const;
    blast_fn real& operator[](u32 i);
    blast_fn real  operator[](u32 i) const ;
    blast_fn Mat3& zero();
    blast_fn Array col(int c);
    blast_fn Vec3  col_copy(int c) const;
};
blast_fn Mat3& zero(Mat3&);
blast_fn Mat3& transpose_inplace(Mat3& m);
blast_fn Mat3  transpose(const Mat3& m);
blast_fn Mat3  operator*(const Mat3& m, const Mat3 rhs);
blast_fn Mat3 operator*(const real x, const Mat3& m);
blast_fn Mat3& operator*=(Mat3& lhs, const Mat3& rhs);
blast_fn Mat3  operator+(const Mat3& lhs, const Mat3& rhs);
blast_fn Mat3& operator+=(Mat3& lhs, const Mat3& rhs);
blast_fn Vec3  operator*(const Mat3& m, const Vec3 v);


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

    // initializer constructor
    //  - note: not available on CUDA
    blast_fn Array(const std::initializer_list<real>&);

    // create an Array from pre-existing data
    //  - note: becomes an alias
    blast_fn Array(real*, u32 n);

    // create an Array from a const std::vector<real
    //  - note: copies data
    //  - note: not available on CUDA
    host_fn Array(const std::vector<real>&);

    // free memory if not an alias
    blast_fn ~Array();

    // copy assignment
    blast_fn Array& operator=(const Array&);

    // move assignment
    blast_fn Array& operator=(Array&&);

    // copy data from a list
    //  - note: must be the same size
    host_fn Array& operator=(const std::initializer_list<real>& other);

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

    // divide with another vector (in place)
    blast_fn Array& operator/=(real n);

    // check if all values of another array are very close the this one (1e-06)
    blast_fn bool operator==(Array&);

    // map the array to the data of the given std::vector<real>
    //  - note: becomes an alias
    //  - note: not available on CUDA
    host_fn Array& alias(svector&);

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
    //  - note: not available on CUDA
    blast_fn void resize(u32 new_size);

    // access the last element
    blast_fn real& back();

    // access the value of the last element
    blast_fn real back() const;
};
blast_fn Array operator-(const Array& v1, const Array& v2);
blast_fn Array operator+(const Array& v1, const Array& v2);
blast_fn Array operator/(const Array& a, real b);
blast_fn Array operator*(const Array& a, real b);
blast_fn Array operator*(real b, const Array& a);
blast_fn real  dot(const Array& a, const Array& b);



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
    //  - note: not available on CUDA
    blast_fn void resize(u32 r, u32 c);

    // copy assignment
    blast_fn Matrix& operator=(const Matrix&);

    // move assignment
    blast_fn Matrix& operator=(Matrix&&);

    blast_fn Matrix operator-();

    blast_fn bool operator==(const Matrix&) const;

    blast_fn bool operator!=(const Matrix&) const;

    // map the Matrix to the data of the given Array
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    blast_fn Matrix& alias(Array&);

    // map the Matrix to the data of the given std::vector<real>
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    //  - note: not available on CUDA
    host_fn Matrix& alias(svector&);

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
blast_fn Matrix operator+(const Matrix& m1, const Matrix& m2);
blast_fn Matrix operator-(const Matrix& m1, const Matrix& m2);
blast_fn Matrix operator*(const Matrix& lhs, const Matrix& rhs);
blast_fn Array operator*(const Matrix& m, const Array& v);
blast_fn Matrix operator*(real& r, const Matrix& m);
blast_fn Matrix operator/(const Matrix& m, real& r);
blast_fn Matrix pw_mult(const Matrix& m1, const Matrix& m2);
blast_fn Matrix eye(const u32 s);
blast_fn Matrix transpose(const Matrix& m);
blast_fn Matrix pinv_svd(const Matrix& m);

// Conversion functions
blast_fn real  wrap2pi(real);
blast_fn real  wrap_to_180(real);
blast_fn real  deg2rad(real);
blast_fn real  rad2deg(real);
blast_fn Array deg2rad(const Array&);
blast_fn Array rad2deg(const Array&);

// fill container with zeros (does not resize)
blast_fn Vec3&   zero(Vec3&);
blast_fn Mat3&   zero(Mat3&);
blast_fn Mat4&   zero(Mat4&);
blast_fn Array&  zero(Array&);
blast_fn Matrix& zero(Matrix&);

// fill container with given value (does not resize)
blast_fn Vec3&   constant(Vec3&,    real val);
blast_fn Mat3&   constant(Mat3&,    real val);
blast_fn Mat4&   constant(Mat4&,    real val);
blast_fn Array&  constant(Array&,   real val);
blast_fn Matrix& constant(Matrix&,  real val);

// min/max functions
blast_fn real min(const Array&);
blast_fn real min(const Matrix&);
blast_fn u32  argmin(const Array&);
blast_fn u32  argmin(const Matrix&);
blast_fn real max(const Array&);
blast_fn real max(const Matrix&);
blast_fn u32  argmax(const Array&);
blast_fn u32  argmax(const Matrix&);

// Check if all values of two containers are (within eps)
blast_fn bool is_close(const Vec3&,   const Vec3&,   real eps = 1e-05);
blast_fn bool is_close(const Array&,  const Array&,  real eps = 1e-05);
blast_fn bool is_close(const Matrix&, const Matrix&, real eps = 1e-05);
blast_fn bool is_close(const Mat3&,   const Mat3&,   real eps = 1e-05);
blast_fn bool is_close(const Mat4&,   const Mat4&,   real eps = 1e-05);

// Check if all values of the container are small (within eps)
blast_fn bool is_small(const Array&,  real eps = 1e-05);
blast_fn bool is_small(const Matrix&, real eps = 1e-05);
blast_fn bool is_small(const Vec3&,   real eps = 1e-05);
blast_fn bool is_small(const Mat3&,   real eps = 1e-05);
blast_fn bool is_small(const Mat4&,   real eps = 1e-05);

// Random
blast_fn real    get_random();
blast_fn Array&  fill_random(Array&,     real amplitude);
blast_fn Matrix& fill_random(Matrix&,    real amplitude);
blast_fn Array   random_array(u32 size,  real amplitude);

// Clamp
blast_fn real    clamp(real val, real mini, real maxi);
blast_fn Array   clamp(const Array& a, real mini, real maxi);
blast_fn Array   clamp(const Array& a, const Array& lb, const Array& ub);
blast_fn real&   clamp_inplace(real& val, real mini, real maxi);
blast_fn Array&  clamp_inplace(Array& a,  real mini, real maxi);
blast_fn Array&  clamp_inplace(Array& a, const Array& lb, const Array& ub);

// Array operations
blast_fn real    sum(const Array&);
blast_fn real    mean(const Array&);
blast_fn real    norm(const Array&);
blast_fn real    norm_sqr(const Array&);
blast_fn real    norm_inf(const Array&);
blast_fn real    norm_1(const Array&);
blast_fn Array   abs(const Array&);
blast_fn Array&  abs_inplace(Array&);
blast_fn Array   sqr(const Array&); // Square each element
blast_fn Array&  sqr_inplace(Array&); // Square each element (modify Array)

// Misc
blast_fn real sign(real v);
blast_fn void sincos(const Array& angles, Array& sines, Array& cosines);
blast_fn Matrix pinv(const Matrix& m);

// https://stackoverflow.com/a/49943540
host_fn inline double simd_hadd(__m256d v) {
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);
    vlow  = _mm_add_pd(vlow, vhigh);

    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));
}

// https://stackoverflow.com/a/35270026
host_fn inline float simd_hadd(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

host_fn inline float simd_hadd(__m256 v) {
    __m128 vlow  = _mm256_castps256_ps128(v);
    __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
    vlow  = _mm_add_ps(vlow, vhigh);     // add the low 128
    return simd_hadd(vlow);
}

} // namespace blast

#if defined(__CUDA_ARCH__)
#define Malloc(a, s) malloc(s)
#define Free free
#elif defined(_MSC_VER)
#define Malloc(a, s) _aligned_malloc(s, a)
#define Free _aligned_free
#else
#define Malloc aligned_alloc
#define Free free
#endif

#include "math/Array.hpp"
#include "math/Matrix.hpp"
#include "math/Vec3.hpp"
#include "math/Mat3.hpp"
#include "math/misc.hpp"
