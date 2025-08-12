#pragma once

#include <blast>

namespace blast {

inline double compute_objective(Array& current_x, Optimization* opt) {
  double result = 0;
  if (opt->objective.K_time > 0) {
    result += current_x.back();
  }
  for (int i = 0; i < opt->objective.extra_objectives.size(); i++) {
    result += opt->objective.k_extra_objectives[i] * opt->objective.extra_objectives[i](opt);
  }
  return result;
}

inline double objective_function(unsigned int n, const double* x, double* grad, void* data) {
  // {
  //     blast_time_block("objective_function");

  Optimization* opt = (Optimization*) data;

  Array xv;
  xv.alias(x, n);
  double result = compute_objective(xv, opt);

  if (grad) {
    Array x_plus(n);
    for (u32 j = 0; j < n; j++) {
      constexpr real eps = 1e-5;
      memcpy(x_plus.data, x, n * sizeof(real));
      x_plus[j] += eps;

      real r_plus = compute_objective(x_plus, opt);
      grad[j]     = (r_plus - result) / eps;
    }
  }
  return result;
  // }
}

inline void Objective::add_custom_objective(ObjectiveFunction function, real k) {
  k_extra_objectives.push_back(k);
  extra_objectives.push_back(function);
}


} // namespace blast
