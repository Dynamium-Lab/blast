#pragma once
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>

#include "blast.hpp"
#include "blast_simd.hpp"

namespace blast {

using svector = std::vector<real>;

// forward declarations
struct Vec3;
struct Mat3;
struct Array;
struct ArrayAlias;
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
    Vec3& operator+=(Vec3&v);
};

// 3x3 matrix
struct Mat3 {
    real data[9] = {0};

    Mat3() = default;
    Mat3(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3);
    real& operator()(u32 row, u32 col);
    Mat3& operator*=(Mat3& rhs);
};

// Array of real numbers
struct Array {
    real * data = nullptr;
    u32 size = 0;

    Array() = default;
    Array(u32 new_size);            // normal constructor
    Array(const Array&);            // copy constructor
    Array(Array&&);                 // move constructor
    Array(const Matrix& m);         // create an array from a matrix (note:flattens and copies)
    Array(const svector&);
    virtual ~Array();
    Array&  operator=(const Array&); // copy assignment
    Array&  operator=(Array&&);      // move assignment
    Array&  operator=(const std::initializer_list<real>& other);
    real&   operator[](u32 i);
    real    operator[](u32 i) const;
    Array   operator-();
    Array&  operator*=(real);
    bool    operator==(Array&);
    void    resize(u32 new_size);
    real&   back();
    real    back() const;

    // interop with std::vector<real>
    operator        svector&();
    operator const  svector&() const;
};

// Array that does not own it's internal memory (does not free upon destruction)
struct ArrayAlias : Array {
    ArrayAlias() = default;
    ArrayAlias(Matrix& m);
    ArrayAlias(Array&);
    ArrayAlias(ArrayAlias&&);
    ArrayAlias(svector&);
    virtual ~ArrayAlias() override;

    const ArrayAlias& operator=(const Array&);
    ArrayAlias& operator=(Array&);

    //--- deleted functions ---
    ArrayAlias(u32) = delete;
    ArrayAlias(const Array&) = delete;
    ArrayAlias(const Matrix&) = delete;
    ArrayAlias(const svector&) = delete;
    ArrayAlias& operator=(ArrayAlias&&) = delete;
    ArrayAlias& operator=(const std::initializer_list<real>& other) = delete;
    void        resize(u32 new_size) = delete;
};

// Matrix of real numbers
struct Matrix {
    real * data = nullptr;
    u32 size = 0;
    u32 rows = 0;
    u32 cols = 0;

    Matrix() = default;
    Matrix(u32 r, u32 c);   // normal constructor
    Matrix(const Matrix&);  // copy constructor
    Matrix(Matrix&&);       // move constructor
    Matrix(const Array&);   // construct a matrix by copying an array into the first collumn (copy)
    ~Matrix();
    Matrix& operator=(const Matrix&); // copy assignment
    Matrix& operator=(Matrix&&);      // move assignment
    real&   operator()(u32 row, u32 col);
    real    operator()(u32 row, u32 col) const;
    ArrayAlias col(u32 c);
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

inline Array rad2deg(Array& a) {
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = a[i] * 180/pi;
    return r;
}

inline Array deg2rad(Array& a) {
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

inline void print(Vec3 v) {
    printf("[%f, %f, %f]\n", v.x, v.y, v.z);
}

inline void print(Mat3 m) {
    printf("\n[%f, %f, %f]\n[%f, %f, %f]\n[%f, %f, %f]\n",
           m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2), m(2, 0), m(2, 1), m(2, 2));
}

inline void print(Array& a) {
    if(a.size == 0)
        return;
    printf("[");
    for (u32 i = 0; i < a.size-1; i++)
        printf("%0.4f, ", a[i]);
    printf("%0.4f]\n", a[a.size-1]);
}

inline void print(Matrix& m) {
    if(m.size == 0)
        return;
    printf("\n");
    for (u32 i = 0; i < m.rows; i++) {
        printf("[");
        for (u32 j = 0; j < m.cols-1; j++)
            printf("%0.4f, ", m(i, j));
        printf("%0.4f]\n", m(i, m.cols-1));
    }
}

inline void minus_insert(Array& a, Matrix& m, real* dst) {
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

inline void minus_insert(Matrix& m, Array& a, real* dst) {
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

inline real array_min(Array& a) {
    real result = inf;
    for(u32 i = 0; i < a.size; i++)
        result = a[i] < result ? a[i] : result;
    return result;
}

inline real array_max(Array& a) {
    real result = -inf;
    for(u32 i = 0; i < a.size; i++)
        result = a[i] > result ? a[i] : result;
    return result;
}

inline bool close(Array& a1, Array& a2, real eps = 1e-05) {
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

inline Vec3& Vec3::operator+=(Vec3&v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
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

inline Mat3& Mat3::operator*=(Mat3& rhs) {
    data[0] = data[0]*rhs.data[0] + data[3]*rhs.data[1] + data[6]*rhs.data[2];
    data[1] = data[1]*rhs.data[0] + data[4]*rhs.data[1] + data[7]*rhs.data[2];
    data[2] = data[2]*rhs.data[0] + data[5]*rhs.data[1] + data[8]*rhs.data[2];
    data[3] = data[0]*rhs.data[3] + data[3]*rhs.data[4] + data[6]*rhs.data[5];
    data[4] = data[1]*rhs.data[3] + data[4]*rhs.data[4] + data[7]*rhs.data[5];
    data[5] = data[2]*rhs.data[3] + data[5]*rhs.data[4] + data[8]*rhs.data[5];
    data[6] = data[0]*rhs.data[6] + data[3]*rhs.data[7] + data[6]*rhs.data[8];
    data[7] = data[1]*rhs.data[6] + data[4]*rhs.data[7] + data[7]*rhs.data[8];
    data[8] = data[2]*rhs.data[6] + data[5]*rhs.data[7] + data[8]*rhs.data[8];
    return *this;
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

inline Array::Array(Array&& a) : data(a.data), size(a.size) {
    a.data = nullptr;
    a.size = 0;
}

inline Array::Array(const Matrix& m) : size(m.size) {
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, m.data, size*sizeof(real));
    }
}

inline Array::Array(const svector& v) {
    Assert(!v.empty());
    size = (u32)v.size();
    data = (real*)malloc(size*sizeof(real));
    std::copy(v.begin(), v.end(), data);
}

inline Array::~Array() {
    if (data)
        free(data);
}

inline Array& Array::operator=(const Array& a) {
    if (this != &a) {
        if (!data)
            data = (real*)malloc(a.size * sizeof(real));
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
    for (u32 i = 0; i < size; i++)
        result[i] = -data[i];
    return std::move(result);
}

inline bool  Array::operator==(Array& a) {
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

inline void Array::resize(u32 new_size) {
    data = (real*)realloc(data, new_size*sizeof(real));
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

// todo: should we copy???
inline Array::operator svector&() {
    svector r(size);
    memcpy(r.data(), data, size*sizeof(real));
    return r;
}

inline Array::operator const svector&() const {
    svector r(size);
    memcpy(r.data(), data, size*sizeof(real));
    return r;
}


//------ ArrayAlias ---------------------

inline ArrayAlias::ArrayAlias(Matrix& m) {
    data = m.data;
    size = m.size;
}

inline ArrayAlias::ArrayAlias(Array& a) {
    data = a.data;
    size = a.size;
}

inline ArrayAlias::ArrayAlias(ArrayAlias&& a) {
    data = a.data;
    size = a.size;
}

inline ArrayAlias::ArrayAlias(svector& v) {
    Assert(!v.empty());
    size = (u32)v.size();
    data = v.data();
}

inline ArrayAlias::~ArrayAlias() {}

inline const ArrayAlias& ArrayAlias::operator=(const Array& a) {
    data = a.data;
    size = a.size;
    return *this;
}

inline ArrayAlias& ArrayAlias::operator=(Array& a) {
    data = a.data;
    size = a.size;
    return *this;
}


//------ Matrix ---------------------

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

inline Matrix::~Matrix() {
    if (data)
        std::free(data);
}

inline Matrix& Matrix::operator=(const Matrix& m) {
    if (this != &m) {
        if (!data)
            data = (real*)malloc(m.size * sizeof(real));
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

inline Matrix::Matrix(const Array& v) : size(v.size), cols(1), rows(v.size) {
    if (size) {
        data = (real*)calloc(size, sizeof(real));
        memcpy(data, v.data, size*sizeof(real));
    }
}

inline real& Matrix::operator()(u32 row, u32 col) {
    Assert(row < this->rows && col < this->cols);
    return this->data[row + this->rows*col];
}

inline real Matrix::operator()(u32 row, u32 col) const {
    Assert(row < this->rows && col < this->cols);
    return this->data[row + this->rows*col];
}

inline ArrayAlias Matrix::col(u32 c) {
    Assert(c < this->cols);
    ArrayAlias result;
    result.data = data + rows*c;
    result.size = rows;
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

inline real two_segment_distance_sqr(Vec3 P0, Vec3 P1, Vec3 Q0, Vec3 Q1) {
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

            compute_minimum_parameters(edge, end, b, c, e, g00, g10,
                                       g01, g11, parameter);
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
    return dot(diff, diff);
}

} // namespace blast
