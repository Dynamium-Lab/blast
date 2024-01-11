#pragma once
#include "blast.hpp"
#include "blast_optimization.hpp"

namespace blast {

struct Wolf  {
    blast::Array x;
};

real gwo_optimize(Array& x, Optimisation& optim) {
    Assert(x.size == optim.bspline->xlen(*optim.task));
    const auto N_Dimensions = x.size;
    // Initialization of GWO parameters
    const int N_wolves = 500;    // Number of wolves    
    const int N_iterations = 50;    // Number of iterations 
    Array x_Alpha(N_Dimensions);
    Array x_Beta(N_Dimensions);
    Array x_Delta(N_Dimensions);

    real Alpha_fit = INF_REAL;
    real Beta_fit = INF_REAL;
    real Delta_fit = INF_REAL;
    std::vector<Wolf> wolf(N_wolves);
    // Initialization of wolves' positions
    for (int i = 0; i < N_wolves; i++) {
        for (int j = 0; j < N_Dimensions; j++) {
            wolf[i].x = random_array(N_Dimensions, 1); // Random values between -1 and 1
        }
    }
    // GWO Computation
    for (int k = 0; k < N_iterations; k++) {
        real a = 2.0 - k * (2.0 / N_iterations);
        for (int i = 0; i < N_wolves; i++) {
            // Evaluating fitness and deciding which wolves are Alpha, Beta and Delta
            real fitness = penalty_obj_time(wolf[i].x, optim);
            // Alpha
            if (fitness < Alpha_fit) {
                Alpha_fit = fitness; 
                x_Alpha = wolf[i].x;
            }
            // Beta
            if (fitness > Alpha_fit && fitness < Beta_fit) {
                Beta_fit = fitness; 
                x_Beta = wolf[i].x;
            }
            // Delta
            if (fitness > Alpha_fit && fitness > Beta_fit && fitness < Delta_fit) {
                Delta_fit = fitness; 
                x_Delta = wolf[i].x;
            }
            // Hunting
            // Approaching the prey 
            
            // Search for the Prey (Diverging the wolves from each other)
            real A_alpha = a * ((2*(0.5 * get_random() + 0.5))-1);
            real A_beta = a * ((2*(0.5 * get_random() + 0.5))-1);
            real A_delta = a * ((2*(0.5 * get_random() + 0.5))-1);
                
            // Emphasising exploration 
            // Effects of obstacles to approaching prey in nature
            // Works like a weight to make it harder for the wolf to reach the prey
            real C_alpha = 2* (0.5 * get_random() + 0.5);
            real C_beta = 2* (0.5 * get_random() + 0.5);
            real C_delta = 2* (0.5 * get_random() + 0.5);
                
            for (int j = 0; j < N_Dimensions; j++) {
                // Search Agent updates its position
                real D_alpha = abs(C_alpha * x_Alpha[j] - wolf[i].x[j]);;
                real D_beta = abs(C_beta * x_Beta[j] - wolf[i].x[j]);;
                real D_delta = abs(C_delta * x_Delta[j] - wolf[i].x[j]);

                real X1 = (x_Alpha[j] - A_alpha * D_alpha);
                real X2 = (x_Beta[j] - A_beta * D_beta);
                real X3 = (x_Delta[j] - A_delta * D_delta);

                wolf[i].x[j]= (X1 + X2 + X3) / 3;
            }
            wolf[i].x.back() = clamp(wolf[i].x.back(), 0.1, 10);
        }
    }
    real best_f = penalty_obj_time(x_Alpha, optim);
    x = x_Alpha;
    return best_f;
}
}