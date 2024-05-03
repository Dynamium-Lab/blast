#pragma once
#include "blast_optimization.hpp"

namespace blast {

struct Wolf  {
    blast::Array x;
};

blast_fn real gwo_optimize(Array& x, Optimisation& optim) {
    Assert(x.size == optim.bspline->xlen(*optim.task));
    const auto n_dimensions = x.size;
    // Initialization of GWO parameters
    const int n_wolves = 500;    // Number of wolves
    const int n_iterations = 50;    // Number of iterations
    Array x_alpha(n_dimensions);
    Array x_beta(n_dimensions);
    Array x_delta(n_dimensions);

    real alpha_fit = INF_REAL;
    real beta_fit = INF_REAL;
    real delta_fit = INF_REAL;
    std::vector<Wolf> wolf(n_wolves);
    // Initialization of wolves' positions
    for (int i = 0; i < n_wolves; i++) {
        for (int j = 0; j < (int)n_dimensions; j++) {
            wolf[i].x = random_array(n_dimensions, 1); // Random values between -1 and 1
        }
    }
    // GWO Computation
    for (int k = 0; k < n_iterations; k++) {
        real a = 2.0 - k * (2.0 / n_iterations);
        for (int i = 0; i < n_wolves; i++) {
            // Evaluating fitness and deciding which wolves are Alpha, Beta and Delta
            real fitness = penalty_obj_time(wolf[i].x, optim);
            // Alpha
            if (fitness < alpha_fit) {
                alpha_fit = fitness;
                x_alpha = wolf[i].x;
            }
            // Beta
            if (fitness > alpha_fit && fitness < beta_fit) {
                beta_fit = fitness;
                x_beta = wolf[i].x;
            }
            // Delta
            if (fitness > alpha_fit && fitness > beta_fit && fitness < delta_fit) {
                delta_fit = fitness;
                x_delta = wolf[i].x;
            }
            // Hunting
            // Approaching the prey

            // Search for the Prey (Diverging the wolves from each other)
            real a_alpha = a * ((2*(0.5 * get_random() + 0.5))-1);
            real a_beta = a * ((2*(0.5 * get_random() + 0.5))-1);
            real a_delta = a * ((2*(0.5 * get_random() + 0.5))-1);

            // Emphasising exploration
            // Effects of obstacles to approaching prey in nature
            // Works like a weight to make it harder for the wolf to reach the prey
            real c_alpha = 2* (0.5 * get_random() + 0.5);
            real c_beta = 2* (0.5 * get_random() + 0.5);
            real c_delta = 2* (0.5 * get_random() + 0.5);

            for (int j = 0; j < (int)n_dimensions; j++) {
                // Search Agent updates its position
                real d_alpha = abs(c_alpha * x_alpha[j] - wolf[i].x[j]);;
                real d_beta = abs(c_beta * x_beta[j] - wolf[i].x[j]);;
                real d_delta = abs(c_delta * x_delta[j] - wolf[i].x[j]);

                real X1 = (x_alpha[j] - a_alpha * d_alpha);
                real X2 = (x_beta[j] - a_beta * d_beta);
                real X3 = (x_delta[j] - a_delta * d_delta);

                wolf[i].x[j]= (X1 + X2 + X3) / 3;
            }
            wolf[i].x.back() = clamp(wolf[i].x.back(), 0.1, 10);
        }
    }
    real best_f = penalty_obj_time(x_alpha, optim);
    x = x_alpha;
    return best_f;
}
}