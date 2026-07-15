#pragma once

#include <blast>

namespace blast {

inline real compute_objective(Array& current_x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  PROFILE_FUNCTION;
#endif
  real result = 0;
  if (opt->objective.time_weight > 0) {
    result += current_x.back();
  }
  for (int i = 0; i < opt->objective.extra_objectives.size(); i++) {
    result += opt->objective.k_extra_objectives[i] * opt->objective.extra_objectives[i](opt);
  }
  return result;
}

inline real objective_function(unsigned int n, const real* x, real* grad, void* data) {
#if BLAST_TRACE_LEVEL >= 1
  PROFILE_FUNCTION;
#endif

  Optimization* opt = (Optimization*) data;

  Array xv;
  xv.alias(x, n);
  real result = compute_objective(xv, opt);

  if (grad) {
    Array x_plus(n);
    for (u32 j = 0; j < n; j++) {
      constexpr real eps = BLAST_FD_STEP;
      memcpy(x_plus.data, x, n * sizeof(real));
      x_plus[j] += eps;

      real r_plus = compute_objective(x_plus, opt);
      grad[j]     = (r_plus - result) / eps;
    }
  }
  return result;
}

inline void Objective::add_custom_objective(ObjectiveFunction function, real k) {
  k_extra_objectives.push_back(k);
  extra_objectives.push_back(function);
}


} // namespace blast
