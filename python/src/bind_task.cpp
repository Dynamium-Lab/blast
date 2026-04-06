#include <blast>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include "conversions.hpp"

namespace nb = nanobind;
using namespace blast;

void bind_task(nb::module_& m) {
    // ------------------------------------------------------------------
    // TaskBoundary
    // ------------------------------------------------------------------
    nb::class_<TaskBoundary>(m, "TaskBoundary")
        .def(nb::init<>())
        // Getters return copies (capsule-owned) → use rv_policy::move so
        // nanobind does not also try to apply reference_internal (double-owner).
        .def_prop_rw("position",
            [](const TaskBoundary& b) { return array_to_copy(b.position); },
            [](TaskBoundary& b, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) {
                b.position = array_from_numpy(a);
            },
            nb::rv_policy::move,
            "Joint positions (rad).")
        .def_prop_rw("velocity",
            [](const TaskBoundary& b) { return array_to_copy(b.velocity); },
            [](TaskBoundary& b, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) {
                b.velocity = array_from_numpy(a);
            },
            nb::rv_policy::move,
            "Joint velocities (rad/s).")
        .def_prop_rw("acceleration",
            [](const TaskBoundary& b) { return array_to_copy(b.acceleration); },
            [](TaskBoundary& b, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) {
                b.acceleration = array_from_numpy(a);
            },
            nb::rv_policy::move,
            "Joint accelerations (rad/s²).");

    // ------------------------------------------------------------------
    // Task
    // ------------------------------------------------------------------
    nb::class_<Task>(m, "Task")
        .def(nb::init<>())
        .def(nb::init<u32>(), nb::arg("n_joints"),
             "Construct a zero-initialised task for n_joints joints.")
        .def_static("stop_to_stop",
            [](nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> start,
               nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> goal) {
                return Task::stop_to_stop(array_from_numpy(start),
                                          array_from_numpy(goal));
            },
            nb::arg("start_pos"), nb::arg("goal_pos"),
            "Create a point-to-point task that starts and ends at rest.\n\n"
            "Parameters\n"
            "----------\n"
            "start_pos : array-like, shape (n_joints,)\n"
            "    Starting joint angles (rad).\n"
            "goal_pos : array-like, shape (n_joints,)\n"
            "    Goal joint angles (rad).")
        .def_rw("start", &Task::start, "Start boundary condition (TaskBoundary).")
        .def_rw("goal",  &Task::goal,  "Goal boundary condition (TaskBoundary).")
        .def("__repr__", [](const Task& t) {
            return "<Task n_joints=" + std::to_string(t.start.position.size) + ">";
        });
}
