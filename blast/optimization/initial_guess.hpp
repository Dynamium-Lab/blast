#pragma once

#include <blast>

namespace blast {

// Create n pseudo-random arrays of d dimensions form : d-by-n
// (https://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/#GeneralizingGoldenRatio)
inline blast::Matrix lds_golden_ratio(const u32 n, const u32 d) {
  blast::Matrix result(d, n);

  blast::real g = 1.6180339887498948482; // Recommanded by wikipedia article

  // precompute factors a
  std::vector<blast::real> a;
  a.reserve(d);
  a.push_back(1.0 / g);
  for (u32 i = 1; i < d; i++)
    a.push_back(a.at(i - 1) * a.at(0));

  for (u32 i = 0; i < n; i++) {
    for (u32 j = 0; j < d; j++) {
      result(j, i) = fmod(0.5 + a[j] * (i + 1), 1);
    }
  }

  return result;
}

Array guess_random(Bspline& bspline, Matrix& task) {
  Array x(bspline.x_len(task));
  fill_random(x, 1);
  x.back() = std::abs(x.back()) * 2 + 0.1;
  return x;
}

inline Matrix get_random_guesses(const u32 n, const u32 d) {
  Array  x(d);
  Matrix result(d, n);
  for (int i = 0; i < n; i++) {
    fill_random(x, 1);
    x.back() = std::abs(x.back()) * 5 + 0.1;
    for (int j = 0; j < d; j++) {
      result(j, i) = x[j];
    }
  }
  return result;
}

inline Array guess_shot_mean_segments(Optimization* opt) { // no collisions
  Array best_x(opt->bspline.x_len(opt->task));
  real  best_val = INF_REAL;
  for (int idx_nshot = 0; idx_nshot < opt->guess.n_random_shots; idx_nshot++) {
    auto x = guess_random((opt->bspline), opt->task);
    opt->bspline.compute_trajectory(x, opt->task);
    Array  c1(opt->constraints.n_constraints);
    Matrix gradient;
    constraints_and_gradients_with_segments(x, *opt, c1, gradient); // todo: change for segments
    real r = 0;
    for (u32 i = 0; i < c1.size; i++)
      r += std::max({c1[i], (real) 0});
    Assert(!isnan(r));
    r = r * x.back(); // todo: Evaluate time estimate impact on trajectory
    if (r < best_val) {
      best_x   = x;
      best_val = r;
    }
  }
  return best_x;
}

inline Array get_best_x_segments(Optimization* opt) {
  real  best_f = INF_REAL;
  Array best_x(opt->guess.candidates.rows);
  for (u32 i = 0; i < opt->guess.candidates.cols; i++) {
    auto current_x = opt->guess.candidates.col(i);
    auto current_f = compute_objective(current_x, opt);

    best_x = current_f < best_f ? current_x : best_x;
    best_f = current_f < best_f ? current_f : best_f;
  }

  return best_x;
}

// Deterministic guess: straight line start->goal, time from manip's vel/accel limits
// (per-joint trapezoidal/triangular, max across joints). Requires fully specified task.
inline Array guess_straight_line(Optimization* opt) {
  Array     x(opt->bspline.x_len(opt->task));
  const u32 n_free = opt->bspline.n_ctrl - 6;
  Assert(x.size == (u32) opt->manip.n_joints * n_free + 1);

  real total_time = 0.0;
  u32  k          = 0;
  for (int j = 0; j < opt->manip.n_joints; j++) {
    const real start = opt->task(j, 0);
    const real goal  = opt->task(j, 3);
    const real dist  = std::abs(goal - start);

    const real vmax = opt->manip.velocity_max[j];
    const real amax = opt->manip.acceleration_max[j];
    if (vmax > 0 && amax > 0) {
      const real d_switch = vmax * vmax / amax; // dist covered while accelerating to vmax and back down
      const real t_j      = (dist >= d_switch) ? dist / vmax + vmax / amax : 2.0 * std::sqrt(dist / amax);
      total_time          = std::max(total_time, t_j);
    }

    for (u32 i = 0; i < n_free; i++) {
      const real a = (real) (i + 1) / (real) (n_free + 1);
      x[k++]       = start + a * (goal - start);
    }
  }
  x.back() = total_time > 0 ? total_time : 2.0;
  return x;
}

inline Array init_guess_segments(Optimization* opt) {
  Array x(opt->bspline.x_len(opt->task));
  switch (opt->guess.type) {
    case Guess::random: {
      x = guess_random(opt->bspline, opt->task);
      break;
    }
    case Guess::shotgun: {
      x = guess_shot_mean_segments(opt);
      break;
    }
    case Guess::straight_line: {
      x = guess_straight_line(opt);
      break;
    }
    case Guess::custom: {
      x = opt->guess.initial_x;
      break;
    }
    case Guess::from_list: {
      x = get_best_x_segments(opt);
      break;
    }
    default:
      Assert(false);
  }
  return x;
}

// todo: Remove collisions
inline Array guess_shot_mean(Optimization* opt) { // no collisions
  Array best_x(opt->bspline.x_len(opt->task));
  real  best_val = INF_REAL;
  for (int idx_nshot = 0; idx_nshot < opt->guess.n_random_shots; idx_nshot++) {
    auto x = guess_random((opt->bspline), opt->task);
    opt->bspline.compute_trajectory(x, opt->task);
    Array c1(opt->constraints.n_constraints);
    compute_constraints(c1.data, x, opt); // todo: change for segments
    real r = 0;
    for (u32 i = 0; i < c1.size; i++)
      r += std::max({c1[i], (real) 0});
    Assert(!isnan(r));
    r = r * x.back(); // todo: Evaluate time estimate impact on trajectory
    if (r < best_val) {
      best_x   = x;
      best_val = r;
    }
  }
  return best_x;
}

inline Array get_best_x(Optimization* opt) {
  real  best_f = INF_REAL;
  Array best_x(opt->guess.candidates.rows);
  for (u32 i = 0; i < opt->guess.candidates.cols; i++) {
    auto current_x = opt->guess.candidates.col(i);
    auto current_f = compute_objective(current_x, opt);

    best_x = current_f < best_f ? current_x : best_x;
    best_f = current_f < best_f ? current_f : best_f;
  }

  return best_x;
}

inline Array init_guess(Optimization* opt) {
  Array x(opt->bspline.x_len(opt->task));
  switch (opt->guess.type) {
    case Guess::random: {
      x = guess_random(opt->bspline, opt->task);
      break;
    }
    case Guess::shotgun: {
      x = guess_shot_mean(opt);
      break;
    }
    case Guess::straight_line: {
      x = guess_straight_line(opt);
      break;
    }
    case Guess::custom: {
      x = opt->guess.initial_x;
      break;
    }
    case Guess::from_list: {
      x = get_best_x(opt);
      break;
    }
    default:
      Assert(false);
  }
  return x;
}

} // namespace blast
