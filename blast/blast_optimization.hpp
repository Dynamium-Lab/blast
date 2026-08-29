#pragma once

#include <blast>
#include <cstring>
#include <functional>

// A float build cannot use external NLopt (its C API is double-only); it requires
// the in-house SQP solver. Mirror the CMake FATAL_ERROR for header-only consumers
// who set the macros directly instead of going through Blast's CMake options.
#if !BLAST_USE_DOUBLES && !defined(BLAST_USE_NATIVE_SQP)
#error "Blast float build (BLAST_USE_DOUBLES=0) requires BLAST_USE_NATIVE_SQP: external NLopt is double-only."
#endif

namespace blast {
struct Optimization;
struct Objective;
struct ConstraintSelection;
struct Guess;

using ConstraintFunctionVector = std::vector<void (*)(real* result, Optimization* opt)>;
using ConstraintFunction       = void (*)(real*, Optimization*);


using ObjectiveFunction = real (*)(Optimization* opt);

struct Objective {
  real time_weight = 1.0;
  // real energy_weight = 0; // not supported atm
  // real jerk_weight = 0; // not supported atm
  // real obstacle_avoidance_weight = 0; // not supported atm

  // Setup function to return the NLopt-compatible objective function
  std::vector<real>                    k_extra_objectives{};
  std::vector<real (*)(Optimization*)> extra_objectives{};

  inline void add_custom_objective(ObjectiveFunction function, real k);
};

struct Guess {
  enum GuessType : u32 {
    custom,
    random,
    shotgun,
    from_list
  };

  GuessType type           = random;
  real      parameter      = 0;
  u32       n_random_shots = 100;
  Array     initial_x;
  Matrix    candidates;

  Guess() = default;

  // Constructor for Guess::custom, initializing initial_x
  explicit Guess(Array x_0) :
      type(Guess::custom),
      initial_x(std::move(x_0)) {
  }

  // Constructor for Guess::from_list, initializing candidates
  explicit Guess(Matrix m) :
      type(Guess::from_list),
      candidates(std::move(m)) {
  }

  // Constructor for Guess::shotgun, initializing n_random_shots
  explicit Guess(u32 shots) :
      type(Guess::shotgun),
      n_random_shots(shots) {
  }
};

enum ConstraintType : unsigned {
  position     = 1 << 0,
  velocity     = 1 << 1,
  acceleration = 1 << 2,
  collision    = 1 << 3,
};

enum class OptimizationMethod : u32 {
  with_segments,            // segment-based constraints (default)
  baseline,                 // point-based constraints, finite-difference gradients
  with_analytical_pva,      // point-based, analytical gradients for position/velocity/acceleration
  with_analytical_dynamics, // point-based, analytical gradients for PVA + torque dynamics
};

struct ConstraintSelection {
  bool position            = false;
  bool velocity            = false;
  bool acceleration        = false;
  bool torque              = false;
  bool jerk                = false; // not supported atm
  bool tool_speed          = false;
  bool self_collisions     = false;
  bool external_collisions = false;

  int n_collision_constraints   = 5; // todo: remove because of new paradigm
  int n_collision_skip          = 2; // todo: remove because of new paradigm
  int n_constraints             = 0;
  int n_constraints_per_segment = 0;

  bool               collect_x_each_iteration = false; // set to true to record optimization vector for all iterations
  std::vector<Array> x_list;

  ConstraintFunctionVector extra_constraints   = {};
  std::vector<u32>         n_extra_constraints = {};
  void                     add_constraint(ConstraintFunction, int n_con);
};

struct Optimization {
  OptimizationMethod  method = OptimizationMethod::with_segments;
  Manipulator         manip;
  Bspline             bspline;
  Guess               guess;
  ConstraintSelection constraints;
  Objective           objective;
  Matrix              task;
  World               world;
  real                trajectory_start_time = 0.0;
  real                success_tolerance     = 0.01; // constraint violation after optimization that is still considered a success
  bool                tighten_for_tolerance = true; // ABLATION SWITCH. false makes
                                                    // tighten_for_success_tolerance() capture its snapshot and
                                                    // return without modifying anything, so the solver sees the
                                                    // TRUE limits and geometry and the gate is applied to
                                                    // untightened constraints. Isolates the tightening from
                                                    // everything else in a build -- the only way to attribute a
                                                    // throughput change to it rather than to a neighbouring
                                                    // commit. Not a production setting: with it off, an accepted
                                                    // solve may violate the real limits by up to
                                                    // success_tolerance.
  real collision_buffer = 0.001;                    // HOW CLOSE, in metres, the trajectory may come to real
                                                    // contact. This is the collision analogue of
                                                    // success_tolerance, which is dimensionless (a fraction of a
                                                    // limit) and so cannot express a distance: comparing it
                                                    // against a raw metre value is what made the default bar 1 cm
                                                    // of penetration. tighten_for_success_tolerance() adds
                                                    // collision_buffer/2 to every capsule and every obstacle, so
                                                    // a pairwise check carries exactly this much margin, and
                                                    // derives the unit change below that makes success_tolerance
                                                    // land on it. Set 0 (or negative) to fall back to
                                                    // success_tolerance itself, which is the pre-buffer
                                                    // behaviour exactly.
  real collision_scale = 1.0;                       // DERIVED -- do not set. tighten_for_success_tolerance() writes
                                                    // success_tolerance/collision_buffer here and
                                                    // restore_from_tolerance() puts it back to 1, so outside the
                                                    // tightening window the collision constraint is raw metres.
                                                    // constraints.hpp multiplies the collision rows by it.
  int  max_tries = 1;                               // Maximum number of tries in the optimization loop.
  int  max_eval  = 1000;                            // Maximum number of function evaluations for a single NLopt call.
  real max_time  = 30.0;                            // Maximum time (seconds) for a single NLopt call.

  void* custom_data;

  Array lb;
  Array ub;

  Optimization() = delete;

  Optimization(const Manipulator& new_manip, const Task& new_task);

  Optimization(const Manipulator& new_manip, const Task& new_task, const Bspline& new_bspline);

  int x_len() const;

  void set_manip(Manipulator new_manip);
  void set_bspline(Bspline new_bspline);
  void set_guess(Guess new_guess);
  void set_constraints(ConstraintSelection new_constraints);
  void set_objective(Objective new_objective);
  void set_task(const Task& new_task);
  void set_world(World new_world);
};

inline void constraints_and_gradients_with_segments(const Array& x, Optimization& opt, Array& constraints,
                                                    Matrix& grad);
// inline void compute_constraints_with_segments(const Array& x, Optimization& opt, Array& constraints);
inline void nlopt_constraints_with_segments(unsigned m, real* result, unsigned x_len, const real* x, real* grad,
                                            void* f_data);

inline void compute_constraints(real* result, const Array& x, Optimization* opt);
inline void nlopt_constraints(unsigned m, real* result, unsigned x_len, const real* x, real* grad,
                              void* f_data);
inline real compute_objective(Array& current_x, Optimization* opt);
inline bool validate_task(Optimization* opt);

inline void nlopt_constraints_with_analytical_pva(unsigned m, real* result, unsigned xlen, const real* x, real* grad, void* f_data);
inline void nlopt_constraints_with_analytical_dynamics(unsigned m, real* result, unsigned xlen, const real* x, real* grad, void* f_data);


// inline real   bound_constraint(const real& value, const real& value_min, const real& value_max);
// inline Matrix get_J_tool(const Optimization* opt);
} // namespace blast

#include "optimization/constraints.hpp"
#include "optimization/initial_guess.hpp"
#include "optimization/objective.hpp"
#include "optimization/optimization.hpp"
