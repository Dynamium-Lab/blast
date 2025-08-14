#pragma once

#include <blast>
#include <utility>
#include "nlopt.h"
#include "optimization/constraints.hpp"
#include "optimization/initial_guess.hpp"
#include "optimization/objective.hpp"

namespace blast {

struct Result {
  bool         success       = false;
  bool         success_false = false;
  real         compute_time;
  Optimization opt;
  Array        x;
  Array        x0;
  nlopt_result nlopt_exit_criteria = NLOPT_SUCCESS;
  int          num_eval            = 0;
  int          num_restart         = 0;
  Trajectory   trajectory;

  Result() = delete;

  explicit Result(Optimization optim) :
      compute_time(0.0),
      opt(std::move(optim)) {}
};

inline Optimization::Optimization(const Manipulator& new_manip, const Matrix& new_task) :
    manip(new_manip),
    bspline(12, 256, 5, new_manip.joints),
    task(new_task),
    custom_data(nullptr) {
  // Default values
  guess.type   = Guess::random;
  guess.n_shot = 100;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.tcp_speed    = true;

  objective.K_time = 1.0;
}

inline Optimization::Optimization(const Manipulator& new_manip, const Matrix& new_task, const Bspline& new_bspline) :
    manip(new_manip),
    bspline(new_bspline),
    task(new_task),
    custom_data(nullptr) {
  guess.type   = Guess::random;
  guess.n_shot = 100;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.tcp_speed    = true;

  objective.K_time = 1.0;
}

inline void Optimization::set_manip(Manipulator new_manip) {
  manip = std::move(new_manip);
}

inline void Optimization::set_bspline(Bspline new_bspline) {
  bspline = std::move(new_bspline);
}

inline void Optimization::set_guess(Guess new_guess) {
  guess = std::move(new_guess);
}

inline void Optimization::set_constraints(ConstraintSelection new_constraints) {
  constraints = std::move(new_constraints);
}

inline void Optimization::set_objective(Objective new_objective) {
  objective = std::move(new_objective);
}

inline void Optimization::set_task(Matrix new_task) {
  task = std::move(new_task);
}

inline void Optimization::set_world(World new_world) {
  world = std::move(new_world);
}

inline void configure_internal_data(Optimization* opt) {
  opt->manip._efforts.resize(opt->manip.joints);
  opt->manip._rotations.resize(opt->manip.joints);
  opt->manip._rotations_mult.resize(opt->manip.joints);
  opt->manip._p_j.resize(opt->manip.joints + 1);
  if (opt->manip._n_caps != 0)
    opt->manip._capsule_list.resize(opt->manip._n_caps);
  // if constexpr (std::is_same_v<T_manip, GenericManipulator>) {
  //     opt->manip._n_caps = opt->manip._collision_model.size();
  //     opt->manip._n_internal_collisions = 0;
  //     for (int i = 0; i < opt->manip._collision_matrix.cols; i++) {
  //         for (int j = i+1; j < opt->manip._n_caps; j++) {
  //             opt->manip._n_internal_collisions += opt->manip._collision_matrix(j, i) == 0 ? 0 : 1;
  //         }
  //         opt->manip._n_internal_collisions += opt->manip._collision_base[i] == 0 ? 0 : 1;
  //     }
  // }
}

inline void ncon(Optimization* opt) {
  const auto n_points            = opt->bspline.n_points;
  const auto n_joints            = opt->manip.joints;
  const auto n_constraints_basic = n_points * n_joints;
  opt->constraints.n_constraints = 0;
  if (opt->constraints.position)
    opt->constraints.n_constraints += n_constraints_basic;
  if (opt->constraints.velocity)
    opt->constraints.n_constraints += n_constraints_basic;
  if (opt->constraints.acceleration)
    opt->constraints.n_constraints += n_constraints_basic;
  if (opt->constraints.torque)
    opt->constraints.n_constraints += n_constraints_basic;

  if (opt->constraints.tcp_speed)
    opt->constraints.n_constraints += n_points;

  if (opt->constraints.external_collisions)
    opt->constraints.n_constraints += opt->constraints.n_collision_constraints;
  if (opt->constraints.self_collisions) {
    // if constexpr (std::is_same_v<T_manip, GenericManipulator>) {
    //     // auto col_matrix = opt->manip._collision_matrix;
    //     for (int i = 0; i < opt->manip._collision_matrix.cols; i++) {
    //         for (int j = i+1; j < opt->manip._n_caps; j++) {
    //             opt->manip._n_internal_collisions += opt->manip._collision_matrix(j, i) == 0 ? 0 : 1;
    //         }
    //         opt->manip._n_internal_collisions += opt->manip._collision_base[i] == 0 ? 0 : 1;
    //     }
    // }
    // else {
    opt->constraints.n_constraints += opt->manip._n_internal_collisions * n_points;
    // }
  }

  for (auto& n_con: opt->constraints.n_extra_constraints)
    opt->constraints.n_constraints += n_con;
}

/*template <typename T_manip>
inline Optimization setup_optimization(const OptimConfig* conf) {
    Assert(conf);

    Optimization result;

    // Must be given
    Assert(conf->manip);
    Assert(conf->task);
    result.task = conf->task;
    result.custom_data = conf->custom_data;

    auto manip = conf->manip;

    // Initialize bspline according to number of points
    result.bspline = new Bspline(conf->n_ctrl, conf->n_points, conf->p, manip.joints);

    // Configure manipulator according to number of points
    result.manip = manip;
    configure_internal_data(&result);

    // Default guess if none is given
    if (conf->guess) {
        result.guess = conf->guess;
    }
    else {
        result.guess = new Guess((u32)100); // Guess::random,
    }

    // Empty world if none is given
    if (conf->world) {
        result.world = conf->world;
    }
    else {
        result.world = new World();
    }

    // Time objective if none is given
    if (conf->objective) {
        result.objective = conf->objective;
    }
    else {
        result.objective = new Objective();
        result.objective.K_time = 1.0;
    }

    // pva + torque constraints if none is given
    if (conf->constraints) {
        result.constraints = conf->constraints;
    }
    else {
        result.constraints = new Constraints();
        result.constraints.position = true;
        result.constraints.velocity = true;
        result.constraints.acceleration = true;
        result.constraints.torque = true;
    }

    // Set number of constraints (extra constraints are taken into account as they are added)
    ncon(result);

    return result;
}*/

inline void initialize_optimization(Optimization* opt) {
  // Constraints
  if (!opt->constraints.external_collisions) {
    opt->constraints.n_collision_constraints = 0;
  }

  // Task
  auto task = opt->task;
  for (int i = 0; i < opt->manip.joints; i++) {
    while (std::abs(task(i, 0) - task(i, 3)) > PI) { // PI since the highest angle between two points on a circle is halfway through, so PI
      // Try updating start value
      real tmp_value = task(i, 0);
      if (task(i, 0) < task(i, 3)) {
        tmp_value += 2 * PI;
      } else {
        tmp_value -= 2 * PI;
      }
      if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
        task(i, 0) = tmp_value;
      } else {
        // Try updating end value
        tmp_value = task(i, 3);
        if (task(i, 3) < task(i, 0)) {
          tmp_value += 2 * PI;
        } else {
          tmp_value -= 2 * PI;
        }
        if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
          task(i, 0) = tmp_value;
        } else {
          break; // Nothing to be done
        }
      }
    }
  }
  opt->task = task;
}

inline Result optimize(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto T1            = get_tick_us();
  Result result(*opt);

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  ncon(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, 0.001);
  Array x_tol(n, 0.000001);


  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);

  bool is_valid_more = false;
  bool is_valid      = false;

  auto start_guess   = opt->guess;
  int  restart_count = 0; // declared here intentionally
  int  max_restart   = 5;
  for (; restart_count <= max_restart; restart_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
    auto x         = init_guess(opt);
    result.x0 = x;

    // launch optimization
    double f;
    result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
    result.num_eval            = nlopt_get_numevals(o);

    // check validity
    Bspline bspline_val(opt->bspline.n_ctrl, opt->bspline.n_points, opt->bspline.p, opt->manip.joints);
    bspline_val.compute_trajectory(x, opt->task);
    auto opt_val(*opt);
    opt_val.set_bspline(bspline_val);
    ncon(&opt_val);
    Array constraints_points(opt_val.constraints.n_constraints);
    compute_constraints(constraints_points.data, x, &opt_val);
    auto max_c = max(constraints_points);
    is_valid   = max_c < 0.01;

    ncon(opt);
    Array constraints_points_2(opt->constraints.n_constraints);
    compute_constraints(constraints_points_2.data, x, opt);


    u64 steps_ms    = static_cast<u64>(std::ceil(x.back() * 1e3 / output_steps_ms));
    x.back()        = static_cast<blast::real>(std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
    int points_more = static_cast<int>(steps_ms + 1);

    Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.joints);
    bspline_val_more.compute_trajectory(x, opt->task);
    auto opt_val_more(*opt);
    opt_val_more.set_bspline(bspline_val_more);
    ncon(&opt_val_more);
    Array constraints_more_points(opt_val_more.constraints.n_constraints);
    compute_constraints(constraints_more_points.data, x, &opt_val_more);
    is_valid_more = max(constraints_more_points) < 0.01;

    result.x = x;

    if (is_valid && is_valid_more) {
      result.trajectory = bspline_val_more.traj;
      break;
    } else if (opt->guess.type != Guess::random && restart_count == 0) {
      opt->guess.type = Guess::random;
    }
  }

  opt->guess = start_guess;

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = *opt;
  result.num_restart   = restart_count;

  nlopt_destroy(o);

  return result;
}


} // namespace blast
