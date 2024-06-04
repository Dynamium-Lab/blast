#include "blast.h"
#include "blast_error.h"
#include "mipp.h"
#include <cstdlib>

namespace blast {

#if defined(_MSC_VER)
#define Malloc(a, s) _aligned_malloc(s, a)
#define Free _aligned_free
#else
#define Malloc aligned_alloc
#define Free free
#endif
const u32 ALIGN = 64;

// Member functions
Array::Array(u32 new_size) : size(new_size) {
    if (size) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        memset(data, 0, size * sizeof(real));
    }
}

Array::Array(u32 n, real value) : size(n) {
    if (n) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        for (u32 i = 0; i < n; i++)
            data[i] = value;
    }
}

Array::Array(const Array& a) : size(a.size) {
    if (size) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        memcpy(data, a.data, size*sizeof(real));
    }
}

Array::Array(Array&& a) : data(a.data), size(a.size), is_alias(a.is_alias) {
    a.is_alias = true;
}

Array::Array(const std::initializer_list<real>& list) : size((u32)list.size()) {
    if (size) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        memcpy(data, list.begin(), list.size() * sizeof(real));
    }
}

Array::Array(real* d, u32 n) {
    data = d;
    size = n;
    is_alias = true;
}

Array::Array(const svector& v) {
    size = (u32)v.size();
    if (size) {
        data = (real*)Malloc(ALIGN, size * sizeof(real));
        memcpy(data, v.data(), size*sizeof(real));
    }
}

Array::~Array() {
    if (!is_alias && data) {
        Free(data);
    }
}

Array& Array::operator=(const Array& a) {
    if (this != &a) {
        if (data && !is_alias)
            Free(data);
        size = a.size;
        if (size) {
            data = (real*)Malloc(ALIGN, a.size * sizeof(real));
            memcpy(data, a.data, size*sizeof(real));
        }
        is_alias = false;
    }
    return *this;
}

Array& Array::operator=(Array&& a) {
    if (this != &a) {
        if (data && !is_alias)
            Free(data);
        data = a.data;
        size = a.size;
        is_alias = a.is_alias;
        a.data = nullptr;
        a.size = 0;
    }
    return *this;
}

Array& Array::operator=(const std::initializer_list<real>& other) {
    Assert(other.size() <= size);
    memcpy(data, other.begin(), other.size() * sizeof(real));
    return *this;
}

Array Array::operator-() {
    Array result(size);
    for (u32 i = 0; i < size; i++)
        result[i] = -data[i];
    return std::move(result);
}

bool Array::operator==(Array& a) {
    Assert(size == a.size);
    return is_close(*this, a);
}

Array& Array::operator*=(real n) {
    for (u32 i = 0; i < size; i++)
        data[i] *= n;
    return *this;
}

Array& Array::operator/=(real n) {
    for (u32 i = 0; i < size; i++)
        data[i] /= n;
    return *this;
}

real& Array::operator[](u32 i) {
    Assert(i < size);
    return data[i];
}

real Array::operator[](u32 i) const {
    Assert(i < size);
    return data[i];
}

Array& Array::alias(svector& v) {
    if (data && !is_alias)
        Free(data);
    data = v.data();
    size = (u32)v.size();
    is_alias = true;
    return *this;
}

Array& Array::alias(Matrix& m) {
    if (data && !is_alias)
        Free(data);
    data = m.data;
    size = m.size;
    is_alias = true;
    return *this;
}

Array& Array::alias(real* p, u32 n) {
    Assert(p);
    if (data && !is_alias)
        Free(data);
    data = p;
    size = n;
    is_alias = true;
    return *this;
}

Array& Array::alias(const real* p, u32 n) {
    Assert(p);
    Assert(data == nullptr);

    data = const_cast<real*>(p);
    size = n;
    is_alias = true;
    return *this;
}

void Array::resize(u32 new_size) {
    Assert(!is_alias);
    if (new_size > size) {
        real* tmp = (real*)Malloc(ALIGN, new_size * sizeof(real));
        memset(tmp, 0, new_size * sizeof(real));
        memcpy(tmp, data, size * sizeof(real));
        Free(data);
        data = tmp;
    }
    size = new_size;
}

real& Array::back() {
    Assert(size);
    return data[size-1];
}

real Array::back() const {
    Assert(size);
    return data[size-1];
}


// Other functions

Array operator-(const Array& v1, const Array& v2) {
    Array r = v1;
    for (u32 i = 0; i < v1.size; i++)
        r[i] -= v2[i];
    return r;
}

Array operator+(const Array& v1, const Array& v2) {
    Assert(v1.size == v2.size);
    Array r = v1;
    for (u32 i = 0; i < v1.size; i++)
        r[i] += v2[i];
    return r;
}

Array operator/(const Array& a, real b) {
    Array r(a);
    r /= b;
    return r;
}

Array operator*(const Array& a, real b) {
    Array r(a);
    r *= b;
    return r;
}

Array operator*(real b, const Array& a) {
    Array r(a);
    r *= b;
    return r;
}

real sum(const Array& a) {
    Assert(a.size > 0);
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += a[i];
    return result;
}

real mean(const Array& a) {
    Assert(a.size > 0);
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += a[i];
    result/=a.size;
    return result;
}

real norm(const Array& a) {
    Assert(a.size > 0);
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += a[i] * a[i];
    return sqrt(result);
}

real norm_sqr(const Array& a) {
    Assert(a.size > 0);
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += a[i] * a[i];
    return result;
}

real norm_1(const Array& a) {
    Assert(a.size > 0);
    real result = 0;
    for (u32 i = 0; i < a.size; i++)
        result += std::abs(a[i]);
    return result;
}

real norm_inf(const Array& a) {
    Assert(a.size > 0);
    return blast::max(blast::abs(a));
}

real min(const Array& a) {
    Assert(a.size > 0);
    real result = INF_REAL;
    for (u32 i = 0; i < a.size; i++)
        result = std::min(result, a[i]);
    return result;
}

real max(const Array& a) {
    Assert(a.size > 0);
    real result = -INF_REAL;
    for (u32 i = 0; i < a.size; i++)
        result = std::max(result, a[i]);
    return result;
}

u32 argmin(const Array& a) {
    Assert(a.size > 0);
    u32 idx = 0;
    real min_tmp = INF_REAL;
    for (u32 i = 0; i < a.size; i++)
        if (a[i] < min_tmp) {
            idx = i;
            min_tmp = a[i];
        }
    return idx;
}

u32 argmax(const Array& a) {
    Assert(a.size > 0);
    u32 idx = 0;
    real max_tmp = -INF_REAL;
    for (u32 i = 0; i < a.size; i++)
        if (a[i] > max_tmp) {
            idx = i;
            max_tmp = a[i];
        }
    return idx;
}

Array abs(const Array&a) {
    Array result(a.size);
    for (u32 i = 0; i < a.size; i++)
        result[i] = std::abs(a[i]);
    return result;
}

Array& abs_inplace(Array& a) {
    for (u32 i = 0; i < a.size; i++)
        a[i] = std::abs(a[i]);
    return a;
}

Array sqr(const Array& a) {
    Array result = a;
    for (u32 i = 0; i < a.size; i++)
        result[i] *= result[i];
    return result;
}

Array& sqr_inplace(Array& a) {
    for (u32 i = 0; i < a.size; i++)
        a[i] *= a[i];
    return a;
}







real dot(const Array& a, const Array& b) {
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
    for (; i < (int)vecLoopSize; i += mipp::N<real>()) {
        ra.load(&a.data[i]);
        rb.load(&b.data[i]);
        accum = mipp::fmadd(ra, rb, accum);
    }
    r = accum.hadd();
    for (; i < (int)a.size; i++)
        r += a[i] * b[i];
#endif
    return r;
}

void sincos(const Array& angles, Array& sines, Array& cosines) {
    Assert(sines.size >= angles.size && cosines.size >= angles.size);

#if defined(__CUDA_ARCH__)
#if BLAST_USE_DOUBLES
    for (int i = 0; i < angles.size; i++)
        ::sincos(angles[i], &sines[i], &cosines[i]);
#else
    for (int i = 0; i < angles.size; i++)
        ::sincosf(angles[i], &sines[i], &cosines[i]);
#endif
#elif defined(__INTEL_COMPILER)
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
#else
    for (int i = 0; i < (int)angles.size; i++) {
        sines[i] = sin(angles[i]);
        cosines[i] = cos(angles[i]);
    }
#endif
}

bool is_close(const Array& a1, const Array& a2, real eps) {
    Assert(a1.size == a2.size);
    for (u32 i =0; i < a1.size; i++)
        if(a1[i] - a2[i] > eps || a1[i] - a2[i] < -eps)
            return false;
    return true;
}

Array& zero(Array& a) {
    if(a.data)
        memset(a.data, 0, a.size*sizeof(real));
    return a;
}

Array& constant(Array& a, real val) {
    for (u32 i = 0; i < a.size; i++)
        a[i] = val;
    return a;
}

bool is_small(const Array& a, real eps) {
    for (u32 i =0; i < a.size; i++)
        if(std::abs(a[i]) > eps)
            return false;
    return true;
}

Array clamp(const Array& a, const Array& lb, const Array& ub) {
    Assert(a.size == lb.size && a.size == ub.size);
    Array r(a.size);
    for (u32 i = 0; i < a.size; i++)
        r[i] = clamp(a[i], lb[i], ub[i]);
    return r;
}

Array clamp(const Array& a, real mini, real maxi) {
    Array r = a;
    for (u32 i = 0; i < a.size; i++)
        clamp_inplace(r[i], mini, maxi);
    return r;
}

Array& clamp_inplace(Array& a,  real mini, real maxi) {
    for (u32 i = 0; i < a.size; i++)
        clamp_inplace(a[i], mini, maxi);
    return a;
}

Array& clamp_inplace(Array& a, const Array& lb, const Array& ub) {
    for (u32 i = 0; i < a.size; i++)
        clamp_inplace(a[i], lb[i], ub[i]);
    return a;
}

Array random_array(u32 n, real A) {
    Array result(n);
    for (int i = 0; i < (int)n; i++)
        result[i] = A * get_random();
    return result;
}

Array& fill_random(Array& a, real A) {
    for (int i = 0; i < (int)a.size; i++)
        a[i] = A * get_random();
    return a;
}



} // namespace blast