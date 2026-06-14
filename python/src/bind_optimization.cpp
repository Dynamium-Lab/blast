#include <blast>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include "conversions.hpp"

namespace nb = nanobind;
using namespace blast;

void bind_optimization(nb::module_& m) {
  // ------------------------------------------------------------------
  // GuessType enum
  // ------------------------------------------------------------------
  nb::enum_<Guess::GuessType>(m, "GuessType")
          .value("custom", Guess::custom,
                 "Use a user-provided initial solution vector.")
          .value("random", Guess::random,
                 "Use a randomly generated initial guess.")
          .value("from_list", Guess::from_list,
                 "Try each column of the candidates matrix as an initial guess.")
          .value("shotgun", Guess::shotgun,
                 "Try n_random_shots random initial guesses and keep the best.")
          .export_values();

  // ------------------------------------------------------------------
  // OptimizationMethod enum
  // ------------------------------------------------------------------
  nb::enum_<OptimizationMethod>(m, "OptimizationMethod")
          .value("with_segments", OptimizationMethod::with_segments,
                 "Segment-based constraints (default, fastest).")
          .value("baseline", OptimizationMethod::baseline,
                 "Point-based constraints with finite-difference gradients.")
          .value("with_analytical_pva", OptimizationMethod::with_analytical_pva,
                 "Point-based with analytical gradients for pos/vel/acc.")
          .value("with_analytical_dynamics", OptimizationMethod::with_analytical_dynamics,
                 "Point-based with analytical gradients for PVA + torque dynamics.")
          .export_values();

  // ------------------------------------------------------------------
  // ConstraintSelection
  // ------------------------------------------------------------------
  nb::class_<ConstraintSelection>(m, "ConstraintSelection")
          .def(nb::init<>())
          .def_rw("position", &ConstraintSelection::position,
                  "Enforce joint position limits.")
          .def_rw("velocity", &ConstraintSelection::velocity,
                  "Enforce joint velocity limits.")
          .def_rw("acceleration", &ConstraintSelection::acceleration,
                  "Enforce joint acceleration limits.")
          .def_rw("torque", &ConstraintSelection::torque,
                  "Enforce joint torque limits (requires dynamics).")
          .def_rw("tool_speed", &ConstraintSelection::tool_speed,
                  "Enforce maximum tool-center-point speed.")
          .def_rw("self_collisions", &ConstraintSelection::self_collisions,
                  "Enforce self-collision avoidance (requires capsule model).")
          .def_rw("external_collisions", &ConstraintSelection::external_collisions,
                  "Enforce collision avoidance with world obstacles.")
          .def_rw("n_collision_constraints", &ConstraintSelection::n_collision_constraints)
          .def_rw("n_collision_skip", &ConstraintSelection::n_collision_skip)
          .def("__repr__", [](const ConstraintSelection& c) {
            std::string s = "<ConstraintSelection";
            if (c.position)
              s += " pos";
            if (c.velocity)
              s += " vel";
            if (c.acceleration)
              s += " acc";
            if (c.torque)
              s += " torque";
            if (c.tool_speed)
              s += " tool_speed";
            if (c.self_collisions)
              s += " self_col";
            if (c.external_collisions)
              s += " ext_col";
            return s + ">";
          });

  // ------------------------------------------------------------------
  // Objective
  // ------------------------------------------------------------------
  nb::class_<Objective>(m, "Objective")
          .def(nb::init<>())
          .def_rw("time_weight", &Objective::time_weight,
                  "Weight for the trajectory time objective (default 1.0).");

  // ------------------------------------------------------------------
  // Guess
  // ------------------------------------------------------------------
  nb::class_<Guess>(m, "Guess")
          .def(nb::init<>())
          .def(nb::init<u32>(), nb::arg("n_random_shots"),
               "Create a random-shots guess with the specified number of attempts.")
          .def_rw("type", &Guess::type,
                  "Strategy for generating the initial solution (GuessType).")
          .def_rw("n_random_shots", &Guess::n_random_shots,
                  "Number of random initial guesses to try (GuessType.random).")
          .def_rw("parameter", &Guess::parameter,
                  "Scalar parameter (interpretation depends on GuessType).")
          .def_prop_rw("initial_x", [](Guess& g) {
                if (g.initial_x.size == 0) return nb::ndarray<nb::numpy, real>{};
                return array_to_view(&g.initial_x); }, [](Guess& g, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { g.initial_x = array_from_numpy(a); }, "Initial solution vector for GuessType.custom (float64 view).")
          .def("__repr__", [](const Guess& g) { return "<Guess type=" + std::to_string((int) g.type) +
                                                       " n_random_shots=" + std::to_string(g.n_random_shots) + ">"; });

  // ------------------------------------------------------------------
  // Optimization
  // ------------------------------------------------------------------
  nb::class_<Optimization>(m, "Optimization")
          .def(nb::init<const Manipulator&, const Task&>(),
               nb::arg("manip"), nb::arg("task"),
               "Construct an optimization problem. Default constraints: pos/vel/acc/tool_speed.")
          .def(nb::init<const Manipulator&, const Task&, const Bspline&>(),
               nb::arg("manip"), nb::arg("task"), nb::arg("bspline"),
               "Construct with an explicit B-spline parametrisation.")
          .def_rw("bspline", &Optimization::bspline)
          .def_rw("guess", &Optimization::guess)
          .def_rw("constraints", &Optimization::constraints)
          .def_rw("objective", &Optimization::objective)
          .def_rw("world", &Optimization::world)
          .def_rw("success_tolerance", &Optimization::success_tolerance,
                  "Fraction of constraint violation still considered a success (default 0.0).")
          .def_rw("max_tries", &Optimization::max_tries,
                  "Maximum number of optimization retries (default 1).")
          .def_rw("max_eval", &Optimization::max_eval,
                  "Maximum NLopt function evaluations per try (default 1000).")
          .def_rw("max_time", &Optimization::max_time,
                  "Maximum wall time per NLopt call in seconds (default 30).")
          .def_rw("trajectory_start_time", &Optimization::trajectory_start_time,
                  "Trajectory start time for dynamic obstacle synchronisation.")
          .def("set_task", &Optimization::set_task, nb::arg("task"))
          .def("set_world", &Optimization::set_world, nb::arg("world"))
          .def("x_len", &Optimization::x_len,
               "Total number of decision variables (B-spline control points + time).")
          .def("__repr__", [](const Optimization& o) {
            return "<Optimization n_joints=" + std::to_string(o.manip.n_joints) +
                   " x_len=" + std::to_string(o.x_len()) + ">";
          });

  // ------------------------------------------------------------------
  // Result
  // ------------------------------------------------------------------
  // Note: Result::opt (raw Optimization*) is intentionally not exposed —
  // it is marked with "// todo: fix this" in the C++ source and can
  // dangle easily from Python.
  nb::class_<Result>(m, "Result")
          .def_ro("success", &Result::success,
                  "True if the optimizer found a feasible trajectory.")
          .def_ro("compute_time", &Result::compute_time,
                  "Wall time for the optimization (milliseconds).")
          .def_ro("num_eval", &Result::num_eval,
                  "Number of NLopt function evaluations.")
          .def_ro("num_tries", &Result::num_tries,
                  "Number of optimization retries performed.")
          .def_ro("max_constraint_value", &Result::max_constraint_value,
                  "Maximum constraint violation in the final solution.")
          .def_ro("max_constraint_idx", &Result::max_constraint_idx,
                  "Index of the most-violated constraint.")
          .def_prop_ro("x", [](Result& r) { return array_to_view(&r.x); }, "Optimized B-spline control vector (float64 view).")
          .def_prop_ro("x0", [](Result& r) { return array_to_view(&r.x0); }, "Initial guess vector that was used (float64 view).")
          .def_prop_ro("trajectory", [](Result& self) -> Trajectory& { return self.trajectory; }, nb::rv_policy::reference_internal, "Computed Trajectory (pos/vel/acc/t). The views it returns keep "
                                                                                                                                     "the Result alive.")
          .def("__repr__", [](const Result& r) { return "<Result success=" + std::string(r.success ? "True" : "False") +
                                                        " compute_time=" + std::to_string(r.compute_time) + "ms>"; });

  // ------------------------------------------------------------------
  // optimize()  — main entry point
  // ------------------------------------------------------------------
  m.def("optimize", [](Optimization* opt, u32 output_steps_ms) -> Result {
            nb::gil_scoped_release release;
            return blast::optimize(opt, output_steps_ms); }, nb::arg("opt"), nb::arg("output_steps_ms") = 1u, "Run trajectory optimization.\n\n"
                                                                                                                               "Parameters\n"
                                                                                                                               "----------\n"
                                                                                                                               "opt : Optimization\n"
                                                                                                                               "    Fully configured optimization problem. Do NOT access opt from\n"
                                                                                                                               "    another thread while this call is in progress.\n"
                                                                                                                               "method : OptimizationMethod\n"
                                                                                                                               "    Gradient/constraint evaluation strategy.\n"
                                                                                                                               "output_steps_ms : int\n"
                                                                                                                               "    Print progress every N milliseconds (0 = silent).\n\n"
                                                                                                                               "Returns\n"
                                                                                                                               "-------\n"
                                                                                                                               "Result\n"
                                                                                                                               "    Optimization result including trajectory and diagnostics.\n\n"
                                                                                                                               "Notes\n"
                                                                                                                               "-----\n"
                                                                                                                               "The GIL is released for the duration of the NLopt call, allowing\n"
                                                                                                                               "other Python threads to run concurrently.");
}
