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
  real K_time = 1.0;
  // real K_energy = 0; // not supported atm
  // real K_jerk = 0; // not supported atm
  // real K_obstacle_avoidance = 0; // not supported atm

  // Setup function to return the NLopt-compatible objective function
  std::vector<real>                    k_extra_objectives;
  std::vector<real (*)(Optimization*)> extra_objectives;

  inline void add_custom_objective(ObjectiveFunction function, real k);
};

struct Guess {
  enum GuessType : u32 {
    custom,
    random,
    from_list,
    rrt_connect
  };
  GuessType type      = random;
  real      parameter = 0;
  u32       n_shot    = 100;
  Array     x0;
  Matrix    candidates;

  Guess() = default;

  // Constructor for Guess::custom, initializing x0
  explicit Guess(Array x_0) :
      type(Guess::custom),
      x0(std::move(x_0)) {}

  // Constructor for Guess::from_list, initializing candidates
  explicit Guess(Matrix m) :
      type(Guess::from_list),
      candidates(std::move(m)) {}

  // Constructor for Guess::random, initializing n_shot
  explicit Guess(u32 shots) :
      type(Guess::random),
      n_shot(shots) {}

  // Constructor for Guess::rrt_connect, initializing parameter
  explicit Guess(real param) :
      type(Guess::rrt_connect),
      parameter(param) {}
};

enum ConstraintType : unsigned {
  position     = 1 << 0,
  velocity     = 1 << 1,
  acceleration = 1 << 2,
};

struct ConstraintSelection {
  bool position            = false;
  bool velocity            = false;
  bool acceleration        = false;
  bool torque              = false;
  bool jerk                = false; // not supported atm
  bool tcp_speed           = false;
  bool self_collisions     = false;
  bool external_collisions = false;

  u32 n_collision_constraints = 5;
  u32 n_collision_skip        = 2;
  u32 n_constraints           = 0;

  // Added more info for testing
  bool                show_info = false;
  std::vector<Matrix> grad_list;
  std::vector<Array>  constr_list;
  std::vector<Array>  x_list;

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
  real                success_tolerance     = 0.0; // % of constraint violation after optimization that is still considered a success
  int                 max_restart           = 1;

  void* custom_data;

  Optimization() = delete;

  Optimization(const Manipulator& new_manip, const Matrix& new_task);

  Optimization(const Manipulator& new_manip, const Matrix& new_task, const Bspline& new_bspline);

  void set_manip(Manipulator new_manip);
  void set_bspline(Bspline new_bspline);
  void set_guess(Guess new_guess);
  void set_constraints(ConstraintSelection new_constraints);
  void set_objective(Objective new_objective);
  void set_task(Matrix new_task);
  void set_world(World new_world);
};

inline void   compute_constraints(double* result, const Array& x, Optimization* opt);
inline void   nlopt_constraints(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data);
inline double compute_objective(Array& current_x, Optimization* opt);
inline bool   validate_task(Optimization* opt);

// inline real   bound_constraint(const real& value, const real& value_min, const real& value_max);
// inline Matrix get_J_tool(const Optimization* opt);
// inline bool   validate_task(Optimization* opt);


} // namespace blast

#include "optimization/constraints.hpp"
#include "optimization/initial_guess.hpp"
#include "optimization/objective.hpp"
#include "optimization/optimization.hpp"
