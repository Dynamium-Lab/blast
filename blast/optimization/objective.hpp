#pragma once

#include <blast>

namespace blast {

inline real compute_objective(Array& current_x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  PROFILE_FUNCTION;
#endif
  real result = 0;
  if (opt->objective.time_weight > 0) {
    result += current_x.back() * opt->objective.time_weight;
  }
  for (size_t i = 0; i < opt->objective.extra_objectives.size(); i++) {
    result += opt->objective.k_extra_objectives[i] * opt->objective.extra_objectives[i](opt, &current_x);
  }
  for (size_t i = 0; i < opt->objective.extra_objectives_with_gradient.size(); i++) {
    result += opt->objective.k_extra_objectives_with_gradient[i] * opt->objective.extra_objectives_with_gradient[i](opt, &current_x, nullptr);
  }
  return result;
}

inline real objective_function(unsigned int n, const real* x, real* grad, void* data) {
#if BLAST_TRACE_LEVEL >= 1
  PROFILE_FUNCTION;
#endif

  Optimization* opt       = (Optimization*) data;
  const auto&   objective = opt->objective;

  Array xv;
  xv.alias(x, n);

  if (!grad)
    return compute_objective(xv, opt);

  for (u32 j = 0; j < n; j++)
    grad[j] = 0;
  real result = 0;

  // duration
  if (objective.time_weight > 0) {
    result += xv.back() * objective.time_weight;
    grad[n - 1] += objective.time_weight;
  }

  // analytical gradients
  if (!objective.extra_objectives_with_gradient.empty()) {
    Array term_grad(n);
    for (size_t i = 0; i < objective.extra_objectives_with_gradient.size(); i++) {
      for (u32 j = 0; j < n; j++)
        term_grad[j] = 0;
      const real k = objective.k_extra_objectives_with_gradient[i];
      result += k * objective.extra_objectives_with_gradient[i](opt, &xv, &term_grad);
      for (u32 j = 0; j < n; j++)
        grad[j] += k * term_grad[j];
    }
  }

  // forward finite differences
  if (!objective.extra_objectives.empty()) {
    constexpr real eps = BLAST_FD_STEP;
    Array          x_plus(n);
    memcpy(x_plus.data, x, n * sizeof(real));
    for (size_t i = 0; i < objective.extra_objectives.size(); i++) {
      const real k     = objective.k_extra_objectives[i];
      const real value = k * objective.extra_objectives[i](opt, &xv);
      result += value;
      for (u32 j = 0; j < n; j++) {
        x_plus[j] += eps;
        grad[j] += (k * objective.extra_objectives[i](opt, &x_plus) - value) / eps;
        x_plus[j] = x[j];
      }
    }
  }
  return result;
}

inline void Objective::add_custom_objective(ObjectiveFunction function, real k) {
  k_extra_objectives.push_back(k);
  extra_objectives.push_back(function);
}

inline void Objective::add_custom_objective_with_gradient(ObjectiveGradientFunction function, real k) {
  Assert(function != nullptr);
  k_extra_objectives_with_gradient.push_back(k);
  extra_objectives_with_gradient.push_back(function);
}


} // namespace blast
