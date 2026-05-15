#pragma once

#include <blast>
#include <cstring>
#include <functional>

namespace blast {
struct Optimization;
struct Objective;
struct ConstraintSelection;
struct Guess;

using ConstraintFunctionVector = std::vector<void (*)(double* result, Optimization* opt)>;
using ConstraintFunction       = void (*)(double*, Optimization*);


using ObjectiveFunction = double (*)(Optimization* opt);

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
    from_list,
    rrt_connect
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

  // Constructor for Guess::random, initializing n_random_shots
  explicit Guess(u32 shots) :
      type(Guess::random),
      n_random_shots(shots) {
  }

  // Constructor for Guess::rrt_connect, initializing parameter
  explicit Guess(real param) :
      type(Guess::rrt_connect),
      parameter(param) {
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
  Manipulator         manip;
  Bspline             bspline;
  Guess               guess;
  ConstraintSelection constraints;
  Objective           objective;
  Matrix              task;
  World               world;
  real                trajectory_start_time = 0.0;
  real                success_tolerance     = 0.0;  // % of constraint violation after optimization that is still considered a success
  int                 max_tries             = 1;    // Maximum number of tries in the optimization loop.
  int                 max_eval              = 1000; // Maximum number of function evaluations for a single NLopt call.
  real                max_time              = 30.0; // Maximum time (seconds) for a single NLopt call.

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
inline void nlopt_constraints_with_segments(unsigned m, double* result, unsigned x_len, const double* x, double* grad,
                                            void* f_data);

inline void   compute_constraints(double* result, const Array& x, Optimization* opt);
inline void   nlopt_constraints(unsigned m, double* result, unsigned x_len, const double* x, double* grad,
                                void* f_data);
inline double compute_objective(Array& current_x, Optimization* opt);
inline bool   validate_task(Optimization* opt);

inline void nlopt_constraints_with_analytical_pva(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data);
inline void nlopt_constraints_with_analytical_dynamics(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data);


// inline real   bound_constraint(const real& value, const real& value_min, const real& value_max);
// inline Matrix get_J_tool(const Optimization* opt);
// inline bool   validate_task(Optimization* opt);
} // namespace blast

#include "optimization/constraints.hpp"
#include "optimization/initial_guess.hpp"
#include "optimization/objective.hpp"
#include "optimization/optimization.hpp"
