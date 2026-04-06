#include <blast>
#include <nanobind/nanobind.h>

namespace nb = nanobind;

// Forward declarations — each bind_*.cpp defines one of these.
void bind_math(nb::module_& m);
void bind_world(nb::module_& m);
void bind_manipulator(nb::module_& m);
void bind_task(nb::module_& m);
void bind_trajectory(nb::module_& m);
void bind_optimization(nb::module_& m);

NB_MODULE(_blast_core, m) {
    m.doc() = "Blast: fast trajectory optimization for robot manipulators.";

    // Registration order matters: later modules may use types from earlier ones.
    bind_math(m);          // Vec3, Mat3, MAX_JOINTS, MAX_CAPSULES
    bind_world(m);         // Box, Sphere, Capsule, World, Dynamic*
    bind_manipulator(m);   // ManipulatorLimits, ManipulatorKinematics, ..., Manipulator
    bind_task(m);          // TaskBoundary, Task
    bind_trajectory(m);    // Trajectory, Bspline
    bind_optimization(m);  // Guess, ConstraintSelection, Optimization, Result, optimize()

    // Expose the scalar dtype so Python code can always allocate compatible arrays:
    //   np.array([1.0, 2.0], dtype=blast.REAL_DTYPE)
    nb::object numpy = nb::module_::import_("numpy");
    m.attr("REAL_DTYPE")    = numpy.attr(BLAST_USE_DOUBLES ? "float64" : "float32");
    m.attr("REAL_IS_DOUBLE") = (bool)(BLAST_USE_DOUBLES);
}
