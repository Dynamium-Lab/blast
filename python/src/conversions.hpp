#pragma once
// Type-conversion helpers between blast C++ types and numpy arrays.
// All functions are inline — include this header in every bind_*.cpp.

#include <blast>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;
using blast::real;

// ---------------------------------------------------------------------------
// 1-D Array  (blast::Array <-> numpy 1-D float array)
// ---------------------------------------------------------------------------

/// Zero-copy view into a blast::Array (no explicit owner).
/// Use with nb::rv_policy::reference_internal on the property getter:
/// nanobind will keep the parent Python object alive as the ndarray's owner.
inline nb::ndarray<nb::numpy, real>
array_to_view(blast::Array* arr) {
    size_t shape = (size_t)arr->size;
    return nb::ndarray<nb::numpy, real>(arr->data, 1, &shape);
}

/// Copy a blast::Array into a freshly allocated numpy array.
/// Uses blast's copy constructor so the capsule owns the entire Array object.
inline nb::ndarray<nb::numpy, real>
array_to_copy(const blast::Array& arr) {
    auto* copy = new blast::Array(arr);
    nb::capsule owner(copy, [](void* p) noexcept { delete static_cast<blast::Array*>(p); });
    size_t shape = (size_t)copy->size;
    return nb::ndarray<nb::numpy, real>(copy->data, 1, &shape, owner);
}

/// Copy a C-contiguous numpy 1-D float array into a blast::Array.
inline blast::Array
array_from_numpy(nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> arr) {
    blast::Array result((blast::u32)arr.size());
    std::memcpy(result.data, arr.data(), arr.size() * sizeof(real));
    return result;
}

// ---------------------------------------------------------------------------
// 2-D Matrix  (blast::Matrix <-> numpy 2-D float array, Fortran-order)
// ---------------------------------------------------------------------------
// blast::Matrix is column-major: element (r,c) lives at data[r + rows*c].
// We expose this as a Fortran-contiguous numpy array of shape (rows, cols)
// so that Python indexing arr[r, c] correctly maps to the same element.

/// Zero-copy Fortran-order view into a blast::Matrix (no explicit owner).
/// Use with nb::rv_policy::reference_internal on the property getter:
/// nanobind will keep the parent Python object alive as the ndarray's owner.
inline nb::ndarray<nb::numpy, real>
matrix_to_view(blast::Matrix* m) {
    size_t  shape[2]   = {(size_t)m->rows, (size_t)m->cols};
    int64_t strides[2] = {1, (int64_t)m->rows};  // element strides (Fortran/column-major)
    return nb::ndarray<nb::numpy, real>(m->data, 2, shape, {}, strides);
}

/// Copy a blast::Matrix into a freshly allocated Fortran-order numpy array.
/// Uses blast's copy constructor so the capsule owns the entire Matrix object.
inline nb::ndarray<nb::numpy, real>
matrix_to_copy(const blast::Matrix& m) {
    auto* copy = new blast::Matrix(m);
    nb::capsule owner(copy, [](void* p) noexcept { delete static_cast<blast::Matrix*>(p); });
    size_t  shape[2]   = {(size_t)copy->rows, (size_t)copy->cols};
    int64_t strides[2] = {1, (int64_t)copy->rows};  // element strides (Fortran/column-major)
    return nb::ndarray<nb::numpy, real>(copy->data, 2, shape, owner, strides);
}

// ---------------------------------------------------------------------------
// Vec3 / Mat3  (always copied — they are small value types)
// ---------------------------------------------------------------------------

/// Copy a Vec3 into a (3,) numpy array.
inline nb::ndarray<nb::numpy, real>
vec3_to_numpy(const blast::Vec3& v) {
    real* buf = new real[3]{v.x, v.y, v.z};
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<real*>(p); });
    size_t shape = 3;
    return nb::ndarray<nb::numpy, real>(buf, 1, &shape, owner);
}

/// Build a Vec3 from a length-3 numpy array.
inline blast::Vec3
vec3_from_numpy(nb::ndarray<nb::numpy, real, nb::ndim<1>> arr) {
    if (arr.size() != 3)
        throw nb::value_error("Vec3 requires an array of length 3");
    const real* d = arr.data();
    return blast::Vec3{d[0], d[1], d[2]};
}

/// Copy a Mat3 (column-major) into a C-order (3,3) numpy array.
inline nb::ndarray<nb::numpy, real>
mat3_to_numpy(const blast::Mat3& m) {
    real* buf = new real[9];
    // mat3 is column-major: m(r,c) = data[3*c+r].
    // Output is C-order: buf[r*3+c] = m(r,c).
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            buf[r * 3 + c] = m(r, c);
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<real*>(p); });
    size_t shape[2] = {3, 3};
    return nb::ndarray<nb::numpy, real>(buf, 2, shape, owner);
}

/// Build a Mat3 from a (3,3) numpy array.
inline blast::Mat3
mat3_from_numpy(nb::ndarray<nb::numpy, real, nb::ndim<2>> arr) {
    if (arr.shape(0) != 3 || arr.shape(1) != 3)
        throw nb::value_error("Mat3 requires a (3,3) array");
    blast::Mat3 m;
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            m(r, c) = arr(r, c);
    return m;
}

// ---------------------------------------------------------------------------
// std::array<real, N> helpers
// ---------------------------------------------------------------------------

template<size_t N>
inline nb::ndarray<nb::numpy, real>
real_array_to_numpy(const std::array<real, N>& a, size_t n) {
    real* buf = new real[n];
    std::memcpy(buf, a.data(), n * sizeof(real));
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<real*>(p); });
    return nb::ndarray<nb::numpy, real>(buf, 1, &n, owner);
}

template<size_t N>
inline void
real_array_from_numpy(std::array<real, N>& a,
                      nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> arr) {
    if (arr.size() > N)
        throw nb::value_error("Array length exceeds MAX_JOINTS");
    std::memcpy(a.data(), arr.data(), arr.size() * sizeof(real));
}
