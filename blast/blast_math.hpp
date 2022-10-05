#pragma once
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>

#include "blast_simd.hpp"
#include "blast_error.hpp"

namespace blast {

using u8    = uint8_t;
using u16   = uint16_t;
using u32   = uint32_t;
using i8    = int8_t;
using i16   = int16_t;
using i32   = int32_t;
using real  = double; // for testing float vs double speed and precision
#define BLAST_SIZEOF_REAL 8
static_assert(sizeof(real) == BLAST_SIZEOF_REAL);

using svector = std::vector<real>;

// forward declarations
struct Vec3;
struct Mat3;
struct Mat4;
struct Array;
struct Matrix;

// useful constants
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

// 3x3 matrix
struct Mat3 {
    real data[9];

    Mat3() = default;
    Mat3(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3);

    real& operator()(u32 row, u32 col);
    real& operator[](u32 i);
    real operator[](u32 i) const ;
};

// 4x4 matrix
struct alignas(32) Mat4 {
    union {
        real    data[16];
        __m256d ymm[4];
    };

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
    inline void zero();

    // copy assignment
    Mat4& operator=(const Mat4&);

    // copy list into the matrix
    //  - note: the rest are NOT initialized if the list does not have 16 elements
    Mat4& operator=(const std::initializer_list<real>&);

    // copy packed doubles list into the matrix
    //  - note: the rest are 0 if the list does not have 4 elements
    Mat4& operator=(const std::initializer_list<__m256d>&);

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

    // copy constructor
    Array(const Array&);

    // move constructor
    Array(Array&&);

    // create an Array from pre-existing data
    //  - note: becomes an alias
    Array(real*, u32 n);

    // create an Array from a const std::vector<real
    //  - note: copies data
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

    // check if all values of another array are very close the this one (1e-06)
    bool operator==(Array&);

    // map the array to the data of the given std::vector<real>
    //  - note: becomes an alias
    Array& alias(std::vector<real>&);

    // map the array to the data of the given matrix
    //  - note: interpret all the data of the matrix as one long array, becomes an alias)
    Array& alias(Matrix& m);

    // map the Array to the given data
    //  - note: becomes an alias
    Array& alias(real*, u32);

    // resize the array
    //  - note: old pointers (aliases) to this data may be invalidated
    //  - note: fails if the array is an alias
    void resize(u32 new_size);

    // access the last element
    real& back();

    // access the value of the last element
    real back() const;
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

    // copy assignment
    Matrix& operator=(const Matrix&);

    // move assignment
    Matrix& operator=(Matrix&&);

    // map the Matrix to the data of the given Array
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    Matrix& alias(Array&);

    // map the Matrix to the data of the given std::vector<real>
    //  - note: interpret as a n x 1 matrix
    //  - note: becomes an alias
    Matrix& alias(std::vector<real>&);

    // map the Matrix to the given data
    //  - note: becomes an alias
    Matrix& alias(real*, u32 r, u32 c);

    // access element
    //  - note: fails if the matrix is not big enough
    real& operator()(u32 row, u32 col);

    // access value of element
    //  - note: fails if the matrix is not big enough
    real operator()(u32 row, u32 col) const;

    // return an array accessing the given colum
    //  - note: new Array is aliasing our data
    Array col(u32 c) const;
};



//------ MISC ---------------------

inline real wrap2pi(real r) {
    while (r < -pi)
        r += 2*pi;
    while (r > pi)
        r-= 2*pi;
    return r;
}

inline float wrap_to_180(float r) {
    while (r < -180)
        r += 360;
    while (r > 180)
        r-= 360;
    return r;
}

inline real deg2rad(real r) {
    return r * pi/180;
}

inline real rad2deg(real r) {
    return r * 180/pi;
}

inline Array rad2deg(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * 180/pi;
    return r;
}

inline Array deg2rad(const Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * pi/180;
    return r;
}

inline void zero(Vec3& v) {
    v.x = v.y = v.z = 0;
}

inline void zero(Mat3& m) {
    memset(m.data, 0, 9*sizeof(real));
}

inline void zero(Array& a) {
    if(a.data)
        memset(a.data, 0, a.size*sizeof(real));
}

inline void zero(Matrix& m) {
    if(m.data)
        memset(m.data, 0, m.size*sizeof(real));
}

inline void constant(Vec3& v, real val) {
    v.x = val;
    v.y = val;
    v.z = val;
}

inline void constant(Mat3& m, real val) {
    for (u32 i = 0; i < 9; i++)
        m.data[i] = val;
}

inline void constant(Mat4& m, real val) {
    for (u32 i = 0; i < 16; i++)
        m.data[i] = val;
}

inline void constant(Array& a, real val) {
    for (u32 i = 0; i < a.size; i++)
        a[i] = val;
}

inline void constant(Matrix& m, real val) {
    for (u32 i = 0; i < m.size; i++)
        m.data[i] = val;
}

inline void minus_insert(const Array& a, const Matrix& m, real* dst) {
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

inline void minus_insert(const Matrix& m, const Array& a, real* dst) {
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

inline real array_min(const Array& a) {
    real result = inf;
    for(u32 i = 0; i < a.size; i++)
        result = a[i] < result ? a[i] : result;
    return result;
}

inline real array_max(const Array& a) {
    real result = -inf;
    for(u32 i = 0; i < a.size; i++)
        result = a[i] > result ? a[i] : result;
    return result;
}

inline bool close(const Array& a1, const Array& a2, real eps = 1e-05) {
    Assert(a1.size == a2.size);
    for (u32 i =0; i < a1.size; i++)
        if(a1[i] - a2[i] > eps || a1[i] - a2[i] < -eps)
            return false;
    return true;
}


//------ Vec3 ---------------------

inline Vec3::Vec3(real x, real y, real z)
    : x(x), y(y), z(z), _pad(0) {

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

inline Vec3& operator+=(Vec3& v1, Vec3& v2) {
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

inline Vec3& operator*=(Vec3&v, real a) {
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}


//------ Mat3 ---------------------

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

inline real& Mat3::operator()(u32 row, u32 col) {
    return data[3*col + row];
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

inline Mat3& operator*=(Mat3& lhs, Mat3& rhs) {
    lhs = lhs * rhs;
    return lhs;
}

inline Mat3 operator+(Mat3& lhs, Mat3& rhs) {
    Mat3 r;
    for (u32 i = 0; i < 9; i++)
        r.data[i] = lhs.data[i] + rhs.data[i];
    return r;
}

inline Mat3& operator+=(Mat3& lhs, Mat3& rhs) {
    lhs = lhs + rhs;
    return lhs;
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

inline real& Mat3::operator[](u32 i) {
    Assert(i < 9);
    return data[i];
}

inline real Mat3::operator[](u32 i) const {
    Assert(i < 9);
    return data[i];
}


//------ Mat4 ---------------------

inline Mat4::Mat4(const Mat4& m) {
    memcpy(data, m.data, 16*sizeof(real));
}

inline Mat4::Mat4(real v) {
    zero();
    data[0] = v;
    data[4] = v;
    data[8] = v;
    data[12] = v;
}

inline void Mat4::zero() {
    memset(data, 0, 16*sizeof(real));
}

inline Mat4::Mat4(const std::initializer_list<real>& l) {
    u32 i = 0;
    memset(data, 0, sizeof(*this));
    for (i = 0; i < l.size() && i < 16; i++)
        data[i] = l.begin()[i];
}

inline Mat4& Mat4::operator=(const Mat4& m) {
    memcpy(data, m.data, 16*sizeof(real));
    return *this;
}

inline Mat4& Mat4::operator=(const std::initializer_list<real>& l) {
    u32 i;
    for (i = 0; i < l.size() && i < 16; i++)
        data[i] = l.begin()[i];
    for (; i < 16; i++)
        data[i]=(real)0.0;
    return *this;
}

inline Mat4& Mat4::operator=(const std::initializer_list<__m256d>& l) {
    u32 i;
    for (i = 0; i < l.size() && i < 4; i++)
        ymm[i] = l.begin()[i];
    for (; i < 4; i++)
        ymm[i] = _mm256_setzero_pd();
    return *this;
}

inline real& Mat4::operator()(u32 row, u32 col) {
    Assert(row<4 && col<4);
    return data[4*col + row];
}

inline real Mat4::operator()(u32 row, u32 col) const {
    Assert(row<4 && col<4);
    return data[4*col + row];
}

inline real& Mat4::operator[](u32 i) {
    Assert(i < 16);
    return data[i];
}

inline real Mat4::operator[](u32 i) const {
    Assert(i < 16);
    return data[i];
}

inline Mat4 Mat4::operator-() {
    Mat4 r;
    for (u32 i = 0; i < 16; i++)
        r.data[i] = -data[i];
    return r;
}

inline Mat4& Mat4::operator*=(Mat4& rhs) {
#if 0
    // naive implementation (no SIMD)
    Mat4 tmp = *this;
    for (u32 i = 0; i < 4; i++) {
        for (u32 j = 0; j < 4; j++) {
            (*this)(i, j) = tmp(i, 0)*rhs(0, j) + tmp(i, 1)*rhs(1, j) + tmp(i, 2)*rhs(2, j) + tmp(i, 3)*rhs(3, j);
        }
    }
#elif BLAST_SIZEOF_REAL == 8
    const auto a0 = ymm[0];
    const auto a1 = ymm[1];
    const auto a2 = ymm[2];
    const auto a3 = ymm[3];

    for (u32 i = 0; i < 4; i++) {
        const auto _b0 = _mm256_set1_pd(rhs(0, i));
        const auto _b1 = _mm256_set1_pd(rhs(1, i));
        const auto _b2 = _mm256_set1_pd(rhs(2, i));
        const auto _b3 = _mm256_set1_pd(rhs(3, i));
        ymm[i] = _mm256_fmadd_pd(a0, _b0,
                                 _mm256_fmadd_pd(a1, _b1,
                                         _mm256_fmadd_pd(a2, _b2,
                                                 _mm256_mul_pd(a3, _b3))));
    }

#else
#error Mat4 multiplication is not (yet) implemented for floats
    // todo: implement
#endif
    return *this;
}

inline Mat4& Mat4::operator*=(real v) {
#if BLAST_SIZEOF_REAL == 8
    const __m256d _v = _mm256_set1_pd(v);
    for (u32 i = 0; i < 4; i++)
        ymm[i] =_mm256_mul_pd(_v, ymm[i]);
#else
    const __m256 _v = _mm256_set1_ps(v);
    _mm256_store_ps(data, _mm256_mul_ps(_v, _mm256_load_ps(data)));
    _mm256_store_ps(data+8, _mm256_mul_ps(_v, _mm256_load_ps(data+8)));
#endif
    return *this;
}

inline Mat4& Mat4::operator+=(Mat4& rhs) {
#if 0
    for (u32 i = 0; i < 16; i++)
        data[i] += rhs[i];

#elif BLAST_SIZEOF_REAL == 8
    for (u32 i = 0; i < 4; i++)
        ymm[i] = _mm256_add_pd(ymm[i], rhs.ymm[i]);

#else
#error not yet implemented
#endif
    return *this;
}

inline Mat4& Mat4::operator-=(Mat4& rhs) {
#if 0
    for (u32 i = 0; i < 16; i++)
        data[i] -= rhs[i];

#elif BLAST_SIZEOF_REAL == 8
    for (u32 i = 0; i < 4; i++)
        ymm[i] = _mm256_sub_pd(ymm[i], rhs.ymm[i]);

#else
#error not yet implemented
#endif
    return *this;
}

inline bool Mat4::operator==(const Mat4& rhs) {
    for (u32 i = 0; i < 16; i++)
        if (data[i] != rhs.data[i])
            return false;
    return true;
}

inline Mat4 operator*(const Mat4& lhs, const Mat4& rhs) {
    Mat4 result;
#if 0
    // naive implementation (no SIMD)
    //  note: order of for loops is important!
    for (u32 j = 0; j < 4; j++) {
        for (u32 i = 0; i < 4; i++) {
            result(i, j) = lhs(i, 0)*rhs(0, j) + lhs(i, 1)*rhs(1, j) + lhs(i, 2)*rhs(2, j) + lhs(i, 3)*rhs(3, j);
        }
    }
#elif BLAST_SIZEOF_REAL == 8
// note: this implementation is approx 20% faster than above serialized code
    const auto a0 = lhs.ymm[0];
    const auto a1 = lhs.ymm[1];
    const auto a2 = lhs.ymm[2];
    const auto a3 = lhs.ymm[3];

    for (u32 i = 0; i < 4; i++) {
        const auto _b0 = _mm256_set1_pd(rhs(0, i));
        const auto _b1 = _mm256_set1_pd(rhs(1, i));
        const auto _b2 = _mm256_set1_pd(rhs(2, i));
        const auto _b3 = _mm256_set1_pd(rhs(3, i));
        result.ymm[i] = _mm256_fmadd_pd(a0, _b0,
                                        _mm256_fmadd_pd(a1, _b1,
                                                _mm256_fmadd_pd(a2, _b2,
                                                        _mm256_mul_pd(a3, _b3))));
    }

#else
#error Mat4 multiplication is not (yet) implemented for floats
    // todo: implement
#endif
    return result;
}

inline Array Mat4::col(u32 c) {
    Assert(c < 4);
    Array result;
    result.data = data + 4*c;
    result.size = 4;
    result.is_alias = true;
    return result;
}

inline Array Mat4::col(u32 c) const {
    Assert(c < 4);
    Array result(4);
    memcpy(result.data, data + c*4, 4*sizeof(real));
    return result;
}

// inline Mat4::operator Matrix&() {
//     Matrix r;
//     return r; // todo: implement
// }

// inline Mat4::operator Array&() {
//     Array r;
//     return r; // todo: implement
// }



//------ Array ---------------------

inline Array::Array(u32 new_size) : size(new_size) {
    if (new_size)
        data = (real*)calloc(new_size, sizeof(real));
}

inline Array::Array(const Array& a) : size(a.size) {
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, a.data, size*sizeof(real));
    }
}

inline Array::Array(Array&& a) : data(a.data), size(a.size), is_alias(a.is_alias) {
    a.data = nullptr;
    a.size = 0;
}

inline Array::Array(real* d, u32 n) {
    data = d;
    size = n;
    is_alias = true;
}

inline Array::Array(const svector& v) {
    size = (u32)v.size();
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, v.data(), size*sizeof(real));
    }
}

inline Array::~Array() {
    if (!is_alias && data)
        free(data);
}

inline Array& Array::operator=(const Array& a) {
    if (this != &a) {
        if (data && !is_alias)
            free(data);
        size = a.size;
        if (size) {
            data = (real*)malloc(a.size * sizeof(real));
            std::copy_n(a.data, size, data);
        }
        is_alias = false;
    }
    return *this;
}

inline Array& Array::operator=(Array&& a) {
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

inline Array& Array::operator=(const std::initializer_list<real>& other) {
    Assert(other.size() <= size);
    std::copy_n(other.begin(), size, data);
    return *this;
}

inline Array Array::operator-() {
    Array result(size);
    for (u32 i = 0; i < size; i++)
        result[i] = -data[i];
    return std::move(result);
}

inline bool Array::operator==(Array& a) {
    Assert(size == a.size);
    return close(*this, a);
}

inline Array& Array::operator*=(real n) {
    for (u32 i = 0; i < size; i++)
        data[i] *= n;
    return *this;
}

inline real& Array::operator[](u32 i) {
    Assert(i < size);
    return data[i];
}

inline real Array::operator[](u32 i) const {
    Assert(i < size);
    return data[i];
}

inline Array& Array::alias(svector& v) {
    if (data && !is_alias)
        free(data);
    data = v.data();
    size = (u32)v.size();
    is_alias = true;
    return *this;
}

inline Array& Array::alias(Matrix& m) {
    if (data && !is_alias)
        free(data);
    data = m.data;
    size = m.size;
    is_alias = true;
    return *this;
}

inline Array& Array::alias(real* p, u32 n) {
    Assert(p);
    if (data && !is_alias)
        free(data);
    data = p;
    size = n;
    is_alias = true;
    return *this;
}

inline void Array::resize(u32 new_size) {
    Assert(!is_alias);
    if (new_size > size) {
        data = (real*)realloc(data, new_size*sizeof(real));
        memset(data + size, 0, (new_size-size)*sizeof(real));
    }
    size = new_size;
}

inline real& Array::back() {
    Assert(size);
    return data[size-1];
}

inline real Array::back() const {
    Assert(size);
    return data[size-1];
}

inline Vec3 operator*(Mat3& m, Vec3 v) {
    Vec3 r;
    r.x = m.data[0]*v.x + m.data[3]*v.y + m.data[6]*v.z;
    r.y = m.data[1]*v.x + m.data[4]*v.y + m.data[7]*v.z;
    r.z = m.data[2]*v.x + m.data[5]*v.y + m.data[8]*v.z;
    return r;
}

inline Array operator-(Array& v1, Array& v2) {
    Array r = v1;
    for (u32 i = 0; i < v1.size; i++)
        r[i] -= v2[i];
    return r;
}

inline Array operator*(Array& a, real b) {
    Array r(a);
    r *= b;
    return r;
}

inline Array operator*(real b, Array& a) {
    Array r(a);
    r *= b;
    return r;
}

// Compute the dot product of the given arrays
//  - note: fastest when the number of elements are a factor of 4 or even 8
inline real dot(Array& a, Array& b) {
    Assert(a.size == b.size);

    real r = (real)0.0;
    int i = 0;

    // SIMD what you can
#if BLAST_SIZEOF_REAL == 8
    auto accum = _mm256_setzero_pd();
    for (; i < (int)a.size-3; i += 4) {
        const auto a_v = _mm256_loadu_pd(&a[i]);
        const auto b_v = _mm256_loadu_pd(&b[i]);
        accum = _mm256_fmadd_pd(a_v, b_v, accum);
    }
    r = simd_hadd(accum);
#elif BLAST_SIZEOF_REAL == 4
    auto accum = _mm256_setzero_ps();
    for (; i < (int)a.size-7; i += 8) {
        const auto a_v = _mm256_loadu_ps(&a[i]);
        const auto b_v = _mm256_loadu_ps(&b[i]);
        accum = _mm256_fmadd_ps(a_v, b_v, accum);
    }
    r = simd_hadd(accum);
#endif
    // serialize the rest
    for (; i < (int)a.size; i++)
        r += a[i] * b[i];

    return r;
}

// Compute the sine and the cosine of every element
//  - note: fastest when the number of elements are a factor of 4 (or even 8 if real is float)
//  - note: doing this manually in your function is still faster by about 10%-20%
inline void sincos(const Array& angles, Array& sines, Array& cosines) {
    Assert(angles.size == sines.size && angles.size == cosines.size);
    // SIMD what we can
    int i = 0;
#if BLAST_SIZEOF_REAL == 4
    for (; i < angles.size-7; i += 8) {
        __m256 s_tmp;
        __m256 c_tmp;
        __m256 angle_v = _mm256_load_ps(&angles.data[i]);
        s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
        _mm256_storeu_ps(&sines.data[i], s_tmp);
        _mm256_storeu_ps(&cosines.data[i], c_tmp);
    }
    for (; i < angles.size-3; i += 4) {
        __m128 s_tmp;
        __m128 c_tmp;
        __m128 angle_v = _mm_load_ps(&angles.data[i]);
        s_tmp = _mm_sincos_ps(&c_tmp, angle_v);
        _mm_storeu_ps(&sines.data[i], s_tmp);
        _mm_storeu_ps(&cosines.data[i], c_tmp);
    }
#elif BLAST_SIZEOF_REAL == 8
    for (; i < (int)angles.size-3; i += 4) {
        __m256d s_tmp;
        __m256d c_tmp;
        __m256d angle_v = _mm256_load_pd(angles.data + i);
        s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
        _mm256_storeu_pd(&sines.data[i], s_tmp);
        _mm256_storeu_pd(&cosines.data[i], c_tmp);
    }
    for (; i < (int)angles.size-1; i += 2) {
        __m128d s_tmp;
        __m128d c_tmp;
        __m128d angle_v = _mm_load_pd(&angles.data[i]);
        s_tmp = _mm_sincos_pd(&c_tmp, angle_v);
        _mm_storeu_pd(&sines.data[i], s_tmp);
        _mm_storeu_pd(&cosines.data[i], c_tmp);
    }
#endif
    // serialize the rest
    for (; i < (int)angles.size; i++) {
        sines[i] = sin(angles[i]);
        cosines[i] = cos(angles[i]);
    }
}

//------ Matrix ---------------------

inline Matrix::Matrix(u32 r, u32 c) {
    size = r * c;
    rows = r;
    cols = c;
    if (size)
        data = (real*)calloc(size, sizeof(real));
}

inline Matrix::Matrix(const Matrix& m) : size(m.size), cols(m.cols), rows(m.rows), is_alias(false) {
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, m.data, size*sizeof(real));
    }
}

inline Matrix::Matrix(Matrix&& m) : data(m.data), size(m.size), cols(m.cols), rows(m.rows), is_alias(m.is_alias) {
    m.data = nullptr;
    m.size = 0;
    m.rows = 0;
    m.cols = 0;
}

inline Matrix::Matrix(real* d, u32 r, u32 c) {
    data = d;
    rows = r;
    cols = c;
    size = r*c;
    is_alias = true;
}

inline Matrix::Matrix(const Array& v) : size(v.size), cols(1), rows(v.size), is_alias(false) {
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, v.data, size*sizeof(real));
    }
}

inline Matrix::~Matrix() {
    if (!is_alias && data)
        free(data);
}

inline Matrix& Matrix::operator=(const Matrix& m) {
    Assert(m.data);
    if (this != &m) {
        if (data == nullptr || is_alias) {
            data = (real*)malloc(m.size * sizeof(real));
        }
        else {
            data = (real*)realloc(data, m.size * sizeof(real));
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

inline Matrix& Matrix::operator=(Matrix&& m) {
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

inline Matrix& Matrix::alias(Array& a) {
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

inline Matrix& Matrix::alias(std::vector<real>& v) {
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

inline Matrix& Matrix::alias(real* p, u32 r, u32 c) {
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

inline real& Matrix::operator()(u32 row, u32 col) {
    Assert(row < rows && col < cols);
    return data[row + rows*col];
}

inline real Matrix::operator()(u32 row, u32 col) const {
    Assert(row < rows && col < cols);
    return data[row + rows*col];
}

inline Array Matrix::col(u32 c) const {
    Assert(c < this->cols);
    Array result;
    result.data = data + rows*c;
    result.size = rows;
    result.is_alias = true;
    return result;
}



//------ Collision ---------------------

static real clamped_root(real slope, real h0, real h1) {
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

static void compute_intersection(real* sValue, i32* classify, real b, real f00, real f10, i32* edge, real end[][2]) {
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

static void compute_minimum_parameters(i32* edge, real end[][2], real b, real c, real e, real g00, real g10, real g01, real g11, real* parameter) {
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
            real z = std::min(std::max(h0 / (h0 - h1), zero), one);
            real omz = one - z;
            parameter[0] = omz * end[0][0] + z * end[1][0];
            parameter[1] = omz * end[0][1] + z * end[1][1];
        }
    }
}

inline real two_segment_distance_sqr(Vec3 P0, Vec3 P1, Vec3 Q0, Vec3 Q1) { // todo: Fix funcion to work with gradients
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
