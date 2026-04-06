#include <blast>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include "conversions.hpp"

namespace nb = nanobind;
using namespace blast;

void bind_trajectory(nb::module_& m) {
    // ------------------------------------------------------------------
    // Trajectory
    // ------------------------------------------------------------------
    // Trajectory.pos/vel/acc are blast::Matrix (column-major, shape n_joints x n_points).
    // We return them as Fortran-order numpy views so arr[joint, time_idx] is natural.
    //
    // nb::rv_policy::reference_internal (the default for property getters) is
    // exactly right here: nanobind will keep the Trajectory Python object alive
    // as the owner of the returned ndarray view, preventing the C++ data from
    // being freed while the view is in use.
    nb::class_<Trajectory>(m, "Trajectory")
        .def(nb::init<u32, u32>(), nb::arg("n_points"), nb::arg("n_joints"))
        .def(nb::init<>())
        .def_prop_ro("pos",
            [](const Trajectory& t) { return matrix_to_copy(t.pos); },
            nb::rv_policy::automatic,
            "Joint positions, shape (n_joints, n_points), Fortran-contiguous float64 array.")
        .def_prop_ro("vel",
            [](const Trajectory& t) { return matrix_to_copy(t.vel); },
            nb::rv_policy::automatic,
            "Joint velocities, shape (n_joints, n_points), Fortran-contiguous float64 array.")
        .def_prop_ro("acc",
            [](const Trajectory& t) { return matrix_to_copy(t.acc); },
            nb::rv_policy::automatic,
            "Joint accelerations, shape (n_joints, n_points), Fortran-contiguous float64 array.")
        .def_prop_ro("t",
            [](const Trajectory& t) { return array_to_copy(t.t); },
            nb::rv_policy::automatic,
            "Time samples, shape (n_points,), float64 array.")
        .def("__repr__", [](const Trajectory& t) {
            return "<Trajectory n_joints=" + std::to_string(t.pos.rows) +
                   " n_points=" + std::to_string(t.pos.cols) + ">";
        });

    // ------------------------------------------------------------------
    // Bspline
    // ------------------------------------------------------------------
    nb::class_<Bspline>(m, "Bspline")
        .def(nb::init<u32>(), nb::arg("n_joints"),
             "Construct a B-spline with default parameters (12 ctrl pts, 100 eval pts, degree 5).")
        .def(nb::init<u32, u32, u32, u32>(),
             nb::arg("n_control"), nb::arg("n_points"), nb::arg("degree"), nb::arg("n_joints"),
             "Construct a B-spline with explicit parameters.\n\n"
             "Parameters\n"
             "----------\n"
             "n_control : int\n"
             "    Number of control points. Multiples of 4 are fastest (SIMD).\n"
             "n_points  : int\n"
             "    Number of evaluation points (trajectory resolution).\n"
             "degree    : int\n"
             "    B-spline degree (5 = quintic, C4-continuous).\n"
             "n_joints  : int\n"
             "    Number of robot joints.")
        .def_ro("n_joints",  &Bspline::n_joints,  "Number of joints.")
        .def_ro("n_points",  &Bspline::n_points,  "Number of evaluation points.")
        .def_ro("n_ctrl",    &Bspline::n_ctrl,    "Number of control points.")
        .def_ro("degree",    &Bspline::degree,    "Spline degree.")
        .def_prop_ro("traj",
            [](Bspline& b) -> Trajectory& { return b.traj; },
            nb::rv_policy::reference_internal,
            "Trajectory computed from the last optimize() call.")
        .def("__repr__", [](const Bspline& b) {
            return "<Bspline n_ctrl=" + std::to_string(b.n_ctrl) +
                   " n_points=" + std::to_string(b.n_points) +
                   " degree=" + std::to_string(b.degree) +
                   " n_joints=" + std::to_string(b.n_joints) + ">";
        });
}
