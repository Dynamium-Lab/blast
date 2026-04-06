#include <blast>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <sstream>
#include "conversions.hpp"

namespace nb = nanobind;
using namespace blast;

void bind_math(nb::module_& m) {
    // ------------------------------------------------------------------
    // Vec3
    // ------------------------------------------------------------------
    nb::class_<Vec3>(m, "Vec3")
        .def(nb::init<>())
        .def(nb::init<real, real, real>(), nb::arg("x"), nb::arg("y"), nb::arg("z"))
        .def("__init__", [](Vec3* self, nb::ndarray<nb::numpy, real, nb::ndim<1>> arr) {
            if (arr.size() != 3)
                throw nb::value_error("Vec3 requires an array of length 3");
            const real* d = arr.data();
            new (self) Vec3{d[0], d[1], d[2]};
        }, nb::arg("array"))
        .def_rw("x", &Vec3::x)
        .def_rw("y", &Vec3::y)
        .def_rw("z", &Vec3::z)
        .def("__getitem__", [](const Vec3& v, int i) -> real {
            if (i < 0) i += 3;
            if (i < 0 || i >= 3) throw nb::index_error("Vec3 index out of range");
            return v[i];
        })
        .def("__setitem__", [](Vec3& v, int i, real val) {
            if (i < 0) i += 3;
            if (i < 0 || i >= 3) throw nb::index_error("Vec3 index out of range");
            v[i] = val;
        })
        .def("__len__", [](const Vec3&) { return 3; })
        .def("to_numpy", [](const Vec3& v) { return vec3_to_numpy(v); },
             "Return a copy as a (3,) numpy array.")
        .def_static("from_numpy",
             [](nb::ndarray<nb::numpy, real, nb::ndim<1>> arr) {
                 return vec3_from_numpy(arr);
             },
             nb::arg("array"), "Construct a Vec3 from a (3,) numpy array.")
        .def("__repr__", [](const Vec3& v) {
            std::ostringstream s;
            s << "Vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
            return s.str();
        });

    // Allow numpy (3,) arrays to be passed where a Vec3 is expected.
    nb::implicitly_convertible<nb::ndarray<nb::numpy, real, nb::ndim<1>>, Vec3>();

    // ------------------------------------------------------------------
    // Mat3
    // ------------------------------------------------------------------
    nb::class_<Mat3>(m, "Mat3")
        .def(nb::init<>())
        // Nine-value constructor (column-by-column, matching C++ Mat3 storage):
        //   Mat3(c0r0, c0r1, c0r2, c1r0, c1r1, c1r2, c2r0, c2r1, c2r2)
        .def("__init__",
             [](Mat3* self, real a0, real a1, real a2,
                            real a3, real a4, real a5,
                            real a6, real a7, real a8) {
                 new (self) Mat3{a0, a1, a2, a3, a4, a5, a6, a7, a8};
             },
             nb::arg("a0"), nb::arg("a1"), nb::arg("a2"),
             nb::arg("a3"), nb::arg("a4"), nb::arg("a5"),
             nb::arg("a6"), nb::arg("a7"), nb::arg("a8"),
             "Construct from 9 values stored column-by-column:\n"
             "  Mat3(col0_row0, col0_row1, col0_row2,\n"
             "       col1_row0, col1_row1, col1_row2,\n"
             "       col2_row0, col2_row1, col2_row2)")
        .def("__init__", [](Mat3* self, nb::ndarray<nb::numpy, real, nb::ndim<2>> arr) {
            new (self) Mat3(mat3_from_numpy(arr));
        }, nb::arg("array"), "Construct from a (3,3) numpy array (row-major convention).")
        .def("__call__", [](const Mat3& m, int r, int c) -> real {
            if (r < 0 || r >= 3 || c < 0 || c >= 3)
                throw nb::index_error("Mat3 index out of range");
            return m(r, c);
        }, nb::arg("row"), nb::arg("col"))
        .def("to_numpy", [](const Mat3& m) { return mat3_to_numpy(m); },
             "Return a copy as a C-order (3,3) numpy array.")
        .def_static("from_numpy",
             [](nb::ndarray<nb::numpy, real, nb::ndim<2>> arr) {
                 return mat3_from_numpy(arr);
             },
             nb::arg("array"), "Construct a Mat3 from a (3,3) numpy array.")
        .def("__repr__", [](const Mat3& m) {
            std::ostringstream s;
            s << "Mat3([[" << m(0,0) << ", " << m(0,1) << ", " << m(0,2) << "],\n"
              << "      [" << m(1,0) << ", " << m(1,1) << ", " << m(1,2) << "],\n"
              << "      [" << m(2,0) << ", " << m(2,1) << ", " << m(2,2) << "]])";
            return s.str();
        });

    nb::implicitly_convertible<nb::ndarray<nb::numpy, real, nb::ndim<2>>, Mat3>();

    // ------------------------------------------------------------------
    // Module-level constants
    // ------------------------------------------------------------------
    m.attr("MAX_JOINTS")   = blast::MAX_JOINTS;
    m.attr("MAX_CAPSULES") = blast::MAX_CAPSULES;
}
