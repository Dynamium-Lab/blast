#pragma once

#include <limits>
#include <vector>

namespace blast {


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
    Vec3(real x, real y, real z);
    real& operator[](int i);
    bool operator==(const Vec3&) const;
};
Vec3  cross(Vec3, Vec3);
Vec3& zero(Vec3&);
Vec3& constant(Vec3&, real value);
real  dot(Vec3, Vec3);
real  norm(Vec3);
Vec3  operator-(const Vec3&);
Vec3  operator-(Vec3, Vec3);
Vec3  operator+(Vec3, Vec3);
Vec3  operator*(real, Vec3);
Vec3  operator*(Vec3, real);
Vec3  operator/(Vec3&, real);
Vec3& operator+=(Vec3&, const Vec3&);
Vec3& operator*=(Vec3&, real);


// 3x3 matrix
struct Mat3 {
    real data[9] = {0};
    Mat3() = default;
    Mat3(const Mat3& m);
    Mat3(real x1, real y1, real z1,
         real x2, real y2, real z2,
         real x3, real y3, real z3);
    real& operator()(u32 row, u32 col);
    real& operator[](u32 i);
    real  operator[](u32 i) const ;
    Mat3& zero();
    Array col(int c);
    Vec3  col_copy(int c) const;
};
Mat3& zero(Mat3&);
Mat3& transpose_inplace(Mat3& m);
Mat3  transpose(const Mat3& m);
Mat3  operator*(const Mat3& m, const Mat3 rhs);
Mat3 operator*(const real x, const Mat3& m);
Mat3& operator*=(Mat3& lhs, const Mat3& rhs);
Mat3  operator+(const Mat3& lhs, const Mat3& rhs);
Mat3& operator+=(Mat3& lhs, const Mat3& rhs);
Vec3  operator*(const Mat3& m, const Vec3 v);


// 4x4 matrix
struct Mat4 {
    real data[16] = {0};

    // default
    Mat4() = default;

    // copy constructor
    Mat4(const Mat4&);

    // create a diagonal matrix with the value
    Mat4(real v);

    // copy list into the matrix
    //  - note: the rest are NOT initialized if the list does not have 16 elements
    Mat4(const std::initializer_list<real>&);

    // set every value to 0
    void zero();

    // copy assignment
    Mat4& operator=(const Mat4&);

    // copy list into the matrix
    //  - note: the rest are NOT initialized if the list does not have 16 elements
    Mat4& operator=(const std::initializer_list<real>&);

    // access the given element
    real& operator()(u32 row, u32 col);

    // access the given element (value)
    real operator()(u32 row, u32 col) const;

    // access the element
    //  - note: does not check out of bounds
    real& operator[](u32 i);

    // access value of the element
    //  - note: does not check out of bounds
    real operator[](u32 i) const;

    // unary minus
    Mat4 operator-();

    // matrix multiply inplace
    Mat4& operator*=(Mat4&);

    // multiply each element by the value inplace
    Mat4& operator*=(real);

    // matrix multiply inplace
    Mat4& operator+=(Mat4&);

    // matrix multiply inplace
    Mat4& operator-=(Mat4&);

    // check if all values of another array are the same
    bool operator==(const Mat4&);

    // return a 4D array pointing to the column
    //  - note: the resulting array is an alias
    Array col(u32 c);

    // return a 4D array copying the given colum
    //  - note: Copies data because we are const
    //  - note: resulting Array is NOT an alias
    Array col(u32 c) const;

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
    Array(u32 n);

    // construct and fill with value
    Array(u32 n, real value);

    // copy constructor
    Array(const Array&);

    // move constructor
    Array(Array&&);

    // initializer constructor
    //  - note: not available on CUDA
    Array(const std::initializer_list<real>&);

    // create an Array from pre-existing data
    //  - note: becomes an alias
    Array(real*, u32 n);

    // create an Array from a const std::vector<real
    //  - note: copies data
    //  - note: not available on CUDA
    Array(const svector&);

    // free memory if not an alias
    ~Array();

    // copy assignment
    Array& operator=(const Array&);

    // move assignment
    Array& operator=(Array&&);

    // copy data from a list
    //  - note: must be the same size
    Array& operator=(const std::initializer_list<real>& other);

    // access the element
    //  - note: does not check out of bounds
    real& operator[](u32 i);

    // access value of the element
    //  - note: does not check out of bounds
    real operator[](u32 i) const;

    // unary minus
    Array operator-();

    // multiply with another vector (in place)
    Array& operator*=(real);

    // divide with another vector (in place)
    Array& operator/=(real n);

    // check if all values of another array are very close the this one (1e-06)
    bool operator==(Array&);

    // map the array to the data of the given std::vector<real>
    //  - note: becomes an alias
    //  - note: not available on CUDA
    Array& alias(svector&);

    // map the array to the data of the given matrix
    //  - note: interpret all the data of the matrix as one long array, becomes an alias)
    Array& alias(Matrix& m);

    // map the Array to the given data
    //  - note: becomes an alias
    Array& alias(real*, u32);
    Array& alias(const real*, u32);

    // resize the array
    //  - note: old pointers (aliases) to this data may be invalidated
    //  - note: fails if the array is an alias
    //  - note: not available on CUDA
    void resize(u32 new_size);

    // access the last element
    real& back();

    // access the value of the last element
    real back() const;
};
Array operator-(const Array& v1, const Array& v2);
Array operator+(const Array& v1, const Array& v2);
Array operator/(const Array& a, real b);
Array operator*(const Array& a, real b);
Array operator*(real b, const Array& a);
real  dot(const Array& a, const Array& b);



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
    Matrix(u32 r, u32 c);

    // copy constructor
    Matrix(const Matrix&);

    // move constructor
    Matrix(Matrix&&);

    // create a Matrix from pre-existing data
    //  - note: becomes an alias
    Matrix(real*, u32 r, u32 c);

    // create an matrix of a single columns from an array
    //  - note: creates a copy
    Matrix(const Array&);

    // free memory if not an alias
    ~Matrix();

    // resize the matrix
    //  - note: not available on CUDA
    void resize(u32 r, u32 c);

    // copy assignment
    Matrix& operator=(const Matrix&);

    // move assignment
    Matrix& operator=(Matrix&&);

    Matrix operator-();

    bool operator==(const Matrix&) const;

    bool operator!=(const Matrix&) const;

    // map the Matrix to the data of the given Array
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    Matrix& alias(Array&);

    // map the Matrix to the data of the given std::vector<real>
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    //  - note: not available on CUDA
    Matrix& alias(svector&);

    // map the Matrix to the given data
    //  - note: becomes an alias
    Matrix& alias(real*, u32 r, u32 c);

    // access element
    //  - note: fails if the matrix is not big enough
    real& operator()(u32 row, u32 col);

    // access value of element
    //  - note: fails if the matrix is not big enough
    real operator()(u32 row, u32 col) const;

    // get the address of an element
    real* address(u32 row, u32 col) const;

    // return an array accessing the given colum
    //  - note: new Array is aliasing our data
    Array col(u32 c) const;
};
Matrix operator+(const Matrix& m1, const Matrix& m2);
Matrix operator-(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& lhs, const Matrix& rhs);
Array operator*(const Matrix& m, const Array& v);
Matrix operator*(real& r, const Matrix& m);
Matrix operator/(const Matrix& m, real& r);
Matrix pw_mult(const Matrix& m1, const Matrix& m2);
Matrix eye(const u32 s);
Matrix transpose(const Matrix& m);
Matrix pinv_svd(const Matrix& m);

// Conversion functions
real  wrap2pi(real);
real  wrap_to_180(real);
real  deg2rad(real);
real  rad2deg(real);
Array deg2rad(const Array&);
Array rad2deg(const Array&);

// fill container with zeros (does not resize)
Vec3&   zero(Vec3&);
Mat3&   zero(Mat3&);
Mat4&   zero(Mat4&);
Array&  zero(Array&);
Matrix& zero(Matrix&);

// fill container with given value (does not resize)
Vec3&   constant(Vec3&,    real val);
Mat3&   constant(Mat3&,    real val);
Mat4&   constant(Mat4&,    real val);
Array&  constant(Array&,   real val);
Matrix& constant(Matrix&,  real val);

// min/max functions
real min(const Array&);
real min(const Matrix&);
u32  argmin(const Array&);
u32  argmin(const Matrix&);
real max(const Array&);
real max(const Matrix&);
u32  argmax(const Array&);
u32  argmax(const Matrix&);

// Check if all values of two containers are (within eps)
bool is_close(const Vec3&,   const Vec3&,   real eps = 1e-05);
bool is_close(const Array&,  const Array&,  real eps = 1e-05);
bool is_close(const Matrix&, const Matrix&, real eps = 1e-05);
bool is_close(const Mat3&,   const Mat3&,   real eps = 1e-05);
bool is_close(const Mat4&,   const Mat4&,   real eps = 1e-05);

// Check if all values of the container are small (within eps)
bool is_small(const Array&,  real eps = 1e-05);
bool is_small(const Matrix&, real eps = 1e-05);
bool is_small(const Vec3&,   real eps = 1e-05);
bool is_small(const Mat3&,   real eps = 1e-05);
bool is_small(const Mat4&,   real eps = 1e-05);

// Random
real    get_random();
Array&  fill_random(Array&,     real amplitude);
Matrix& fill_random(Matrix&,    real amplitude);
Array   random_array(u32 size,  real amplitude);

// Clamp
real    clamp(real val, real mini, real maxi);
Array   clamp(const Array& a, real mini, real maxi);
Array   clamp(const Array& a, const Array& lb, const Array& ub);
real&   clamp_inplace(real& val, real mini, real maxi);
Array&  clamp_inplace(Array& a,  real mini, real maxi);
Array&  clamp_inplace(Array& a, const Array& lb, const Array& ub);

// Array operations
real    sum(const Array&);
real    mean(const Array&);
real    norm(const Array&);
real    norm_sqr(const Array&);
real    norm_inf(const Array&);
real    norm_1(const Array&);
Array   abs(const Array&);
Array&  abs_inplace(Array&);
Array   sqr(const Array&); // Square each element
Array&  sqr_inplace(Array&); // Square each element (modify Array)

// Misc
real sign(real v);
void sincos(const Array& angles, Array& sines, Array& cosines);
Matrix pinv(const Matrix& m);

} // namespace blast