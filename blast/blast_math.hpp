#pragma once

#include <immintrin.h>
#include <limits>
#include <vector>

#ifdef __CUDA_ARCH__
#include <cuda_runtime.h>
#include <math_constants.h>
#endif
namespace blast {
// Types declared in this file
struct Vec3;
struct Mat3;
struct Array;
struct Matrix;


// Blast uses doubles by default unless BLAST_USE_DOUBLES is set to 0
#ifndef BLAST_USE_DOUBLES
#define BLAST_USE_DOUBLES 1
#endif
#if BLAST_USE_DOUBLES
using real = double;
#define BLAST_SIZEOF_REAL 8
#else
using real = float;
#define BLAST_SIZEOF_REAL 4
#endif
static_assert(sizeof(real) == BLAST_SIZEOF_REAL);


// Constants
constexpr u32 ALIGN = 64;
constexpr real BLAST_EPSILON = 0.00001;
constexpr real PI = 3.141592653589793;
#ifdef __CUDA_ARCH__
#if BLAST_USE_DOUBLES
const real NAN_REAL = CUDART_NAN;
const real INF_REAL = CUDART_INF;
#else
const real NAN_REAL = CUDART_NAN_F;
const real INF_REAL = CUDART_INF_F;
#endif
#else
constexpr real NAN_REAL = std::numeric_limits<real>::quiet_NaN();
constexpr real INF_REAL = std::numeric_limits<real>::infinity();
#endif


// 3D vector
struct Vec3 {
  real x = 0;
  real y = 0;
  real z = 0;
  real _pad = 0;
  Vec3() = default;
  inline blast_fn Vec3(real x, real y, real z);
  inline blast_fn real &operator[](int i);
  inline blast_fn real operator[](int i) const;
  inline blast_fn bool operator==(const Vec3 &) const;
};

inline blast_fn Vec3 operator-(Vec3);
inline blast_fn Vec3 operator-(Vec3, Vec3);
inline blast_fn Vec3 operator+(Vec3, Vec3);
inline blast_fn Vec3 operator*(real, Vec3);
inline blast_fn Vec3 operator*(Vec3, real);
inline blast_fn Vec3 operator/(Vec3, real);
inline blast_fn Vec3 &operator+=(Vec3 &, const Vec3 &);
inline blast_fn Vec3 &operator*=(Vec3 &, real);
inline blast_fn Vec3 cross(Vec3, Vec3);
inline blast_fn Vec3 &zero(Vec3 &);
inline blast_fn real dot(Vec3, Vec3);
inline blast_fn real norm(Vec3);
inline blast_fn bool is_small(const Vec3 &, real eps = 1e-05);
inline blast_fn bool is_close(const Vec3 &, const Vec3 &, real eps = 1e-05);
inline blast_fn Vec3 &constant(Vec3 &, real val);


// 3x3 matrix
struct Mat3 {
  real data[9] = {0};
  Mat3() = default;
  inline blast_fn Mat3(const Mat3 &m);
  inline blast_fn Mat3(real x1, real y1, real z1,
                       real x2, real y2, real z2,
                       real x3, real y3, real z3);
  inline blast_fn real &operator()(u32 row, u32 col);
  inline blast_fn real operator()(u32 row, u32 col) const;
  inline blast_fn real &operator[](u32 i);
  inline blast_fn real operator[](u32 i) const;
  inline blast_fn Mat3 &zero();
  inline blast_fn Array col(int c);
  inline blast_fn Vec3 col_copy(int c) const;
};

inline blast_fn Mat3 operator*(const Mat3 &m, Mat3 rhs);
inline blast_fn Mat3 operator*(real x, const Mat3 &m);
inline blast_fn Vec3 operator*(const Mat3 &m, Vec3 v);
inline blast_fn Mat3 &operator*=(Mat3 &lhs, const Mat3 &rhs);
inline blast_fn Mat3 operator+(const Mat3 &lhs, const Mat3 &rhs);
inline blast_fn Mat3 &operator+=(Mat3 &lhs, const Mat3 &rhs);
inline blast_fn Mat3 operator-(const Mat3 &m, const Mat3 &m2);
inline blast_fn Mat3 &zero(Mat3 &);
inline blast_fn Mat3 &transpose_inplace(Mat3 &m);
inline blast_fn Mat3 transpose(const Mat3 &m);
inline blast_fn Mat3 eye();
inline blast_fn Mat3 &constant(Mat3 &, real val);
inline blast_fn bool is_close(const Mat3 &, const Mat3 &, real eps = 1e-05);
inline blast_fn bool is_small(const Mat3 &, real eps = 1e-05);
inline blast_fn Mat3 rpy2rotation(Vec3 rpy);


// Array of real numbers
struct Array {
  real *data = nullptr; // pointer to the heap allocated memory
  u32 size = 0; // size in number of elements
  bool is_alias = false; // if true, the Array does not own the data and will not free it. note: should be as short-lived as possible

  // default constructor
  Array() = default;

  // construct and allocate memory for 'n' elements
  inline blast_fn explicit Array(u32 new_size);

  // construct and fill with value
  inline blast_fn Array(u32 n, real value);

  // copy constructor
  inline blast_fn Array(const Array &);

  // move constructor
  inline blast_fn Array(Array &&) noexcept;

  // initializer constructor
  //  - note: not available on CUDA
  host_fn Array(const std::initializer_list<real> &);

  // create an Array from pre-existing data
  //  - note: becomes an alias
  inline blast_fn Array(real *, u32 n);

  // create an Array from a const std::vector<real
  //  - note: copies data
  //  - note: not available on CUDA
  host_fn explicit Array(const std::vector<real> &);

  // free memory if not an alias
  inline blast_fn ~Array();

  // copy assignment
  inline blast_fn Array &operator=(const Array &);

  // move assignment
  inline blast_fn Array &operator=(Array &&) noexcept;

  // copy data from a list
  //  - note: must be the same size
  inline blast_fn Array &operator=(const std::initializer_list<real> &other);

  // access the element
  //  - note: does not check out of bounds
  inline blast_fn real &operator[](u32 i);

  // access value of the element
  //  - note: does not check out of bounds
  inline blast_fn real operator[](u32 i) const;

  // unary minus
  inline blast_fn Array operator-() const;

  // multiply with another vector (in place)
  inline blast_fn Array &operator*=(real);

  // divide with another vector (in place)
  inline blast_fn Array &operator/=(real n);

  // check if all values of another array are very close to this one (1e-06)
  inline blast_fn bool operator==(const Array &) const;

  // map the array to the data of the given std::vector<real>
  //  - note: becomes an alias
  //  - note: not available on CUDA
  inline blast_fn Array &alias(std::vector<real> &);

  // map the array to the data of the given matrix
  //  - note: interpret all the data of the matrix as one long array, becomes an alias
  inline blast_fn Array &alias(Matrix &m);

  // map the Array to the given data
  //  - note: becomes an alias
  inline blast_fn Array &alias(real *, u32);
  inline blast_fn Array &alias(const real *, u32);

  // resize the array
  //  - note: old pointers (aliases) to this data may be invalidated
  //  - note: fails if the array is an alias
  //  - note: not available on CUDA
  inline blast_fn void resize(u32 new_size);

  // access the last element
  inline blast_fn real &back();

  // access the value of the last element
  inline blast_fn real back() const;
};

inline blast_fn Array random_array(u32 n, real A);
inline blast_fn Array operator-(const Array &v1, const Array &v2);
inline blast_fn Array operator+(const Array &v1, const Array &v2);
inline blast_fn Array operator/(const Array &a, real b);
inline blast_fn Array operator*(const Array &a, real b);
inline blast_fn Array operator*(real b, const Array &a);
inline blast_fn real dot(const Array &a, const Array &b);
inline blast_fn Array deg2rad(const Array &);
inline blast_fn Array rad2deg(const Array &);
inline blast_fn Array &constant(Array &, real val);
inline blast_fn Array &fill_random(Array &, real A);
inline blast_fn Array &zero(Array &);
inline blast_fn real min(const Array &);
inline blast_fn real max(const Array &);
inline blast_fn u32 argmin(const Array &);
inline blast_fn u32 argmax(const Array &);
inline blast_fn bool is_close(const Array &, const Array &, real eps = 1e-05);
inline blast_fn bool is_small(const Array &, real eps = 1e-05);
inline blast_fn real sum(const Array &);
inline blast_fn real mean(const Array &);
inline blast_fn real norm(const Array &);
inline blast_fn real norm_sqr(const Array &);
inline blast_fn real norm_inf(const Array &);
inline blast_fn real norm_1(const Array &);
inline blast_fn Array abs(const Array &);
inline blast_fn Array &abs_inplace(Array &);
inline blast_fn Array sqr(const Array &); // Square each element
inline blast_fn Array &sqr_inplace(Array &); // Square each element (modify Array)
inline blast_fn void sincos(const Array &angles, Array &sines, Array &cosines);
inline blast_fn Array clamp(const Array &a, real mini, real maxi);
inline blast_fn Array clamp(const Array &a, const Array &lb, const Array &ub);
inline blast_fn Array &clamp_inplace(Array &a, real mini, real maxi);
inline blast_fn Array &clamp_inplace(Array &a, const Array &lb, const Array &ub);


// Matrix of real numbers
struct Matrix {
  real *data = nullptr; // pointer to the heap allocated memory
  u32 size = 0; // size in number of elements
  u32 rows = 0; // number of rows in the matrix
  u32 cols = 0; // number of columns in the matrix
  bool is_alias = false; // if true, the Matrix does not own the data and will not free it.
  // note: Alias matrices should be as short-lived as possible

  // default constructor
  Matrix() = default;

  // construct and allocate memory for 'r x c' elements
  inline blast_fn Matrix(u32 r, u32 c);

  // copy constructor
  inline blast_fn Matrix(const Matrix &);

  // move constructor
  inline blast_fn Matrix(Matrix &&) noexcept;

  // create a Matrix from pre-existing data
  //  - note: becomes an alias
  inline blast_fn Matrix(real *, u32 r, u32 c);

  // create a matrix of a single columns from an array
  //  - note: creates a copy
  inline blast_fn explicit Matrix(const Array &);

  // free memory if not an alias
  inline blast_fn ~Matrix();

  // resize the matrix
  //  - note: not available on CUDA
  inline blast_fn void resize(u32 r, u32 c);

  // copy assignment
  inline blast_fn Matrix &operator=(const Matrix &);

  // move assignment
  inline blast_fn Matrix &operator=(Matrix &&) noexcept;

  inline blast_fn Matrix operator-();

  inline blast_fn bool operator==(const Matrix &) const;

  inline blast_fn bool operator!=(const Matrix &) const;

  // map the Matrix to the data of the given Array
  //  - note: interpret as an n x 1 matrix
  //  - note: becomes an alias
  inline blast_fn Matrix &alias(Array &);

  // map the Matrix to the data of the given std::vector<real>
  //  - note: interpret as an n x 1 matrix
  //  - note: becomes an alias
  //  - note: not available on CUDA
  inline blast_fn Matrix &alias(std::vector<real> &);

  // map the Matrix to the given data
  //  - note: becomes an alias
  inline blast_fn Matrix &alias(real *, u32 r, u32 c);

  // access element
  //  - note: fails if the matrix is not big enough
  inline blast_fn real &operator()(u32 row, u32 col);

  // access value of element
  //  - note: fails if the matrix is not big enough
  inline blast_fn real operator()(u32 row, u32 col) const;

  // get the address of an element
  inline blast_fn real *address(u32 row, u32 col) const;

  // return an array accessing the given colum
  //  - note: new Array is aliasing our data
  inline blast_fn Array col(u32 c) const;
};

inline blast_fn Matrix operator-(const Matrix &m1, const Matrix &m2);
inline blast_fn Matrix operator*(const Matrix &lhs, const Matrix &rhs);
inline blast_fn Array operator*(const Matrix &m, const Array &v);
inline blast_fn Matrix operator*(real &r, const Matrix &m);
inline blast_fn Matrix operator/(const Matrix &m, real &r);
inline blast_fn Matrix operator+(const Matrix &m1, const Matrix &m2);
inline blast_fn Matrix pw_mult(const Matrix &m1, const Matrix &m2);
inline blast_fn Matrix pinv(const Matrix &m);
inline blast_fn Matrix eye(int s);
inline blast_fn Matrix transpose(const Matrix &m);
inline blast_fn Matrix &zero(Matrix &);
inline blast_fn Matrix &constant(Matrix &, real val);
inline blast_fn bool is_close(const Matrix &, const Matrix &, real eps = 1e-05);
inline blast_fn bool is_small(const Matrix &, real eps = 1e-05);
inline blast_fn real min(const Matrix &);
inline blast_fn real max(const Matrix &);
inline blast_fn std::tuple<u32, u32> argmin(const Matrix &); // return row, column of smallest element
inline blast_fn std::tuple<u32, u32> argmax(const Matrix &); // return row, column of largest element
inline blast_fn Matrix &fill_random(Matrix &, real A);


// Functions on real floating point values
inline blast_fn real wrap2pi(real);
inline blast_fn real wrap_to_180(real);
inline blast_fn real deg2rad(real);
inline blast_fn real rad2deg(real);
inline blast_fn real get_random();
inline blast_fn real clamp(real val, real mini, real maxi);
inline blast_fn real &clamp_inplace(real &val, real mini, real maxi);
inline blast_fn real sign(real v);


// SIMD
inline host_fn double simd_hadd(__m256d v);
inline host_fn float simd_hadd(__m128 v);
inline host_fn float simd_hadd(__m256 v);

} // namespace blast


#if defined(_MSC_VER)
#define Malloc(a, s) _aligned_malloc(s, a)
#define Free _aligned_free
#elif defined(__NVCC__)
#define Malloc(a, s) malloc(s)
#define Free free
#else
#define Malloc aligned_alloc
#define Free free
#endif

#include "math/Array.hpp"
#include "math/Mat3.hpp"
#include "math/Matrix.hpp"
#include "math/Vec3.hpp"
#include "math/misc.hpp"
