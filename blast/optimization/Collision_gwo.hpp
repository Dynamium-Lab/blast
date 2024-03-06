#pragma once
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "blast_world.hpp"
namespace blast {

struct Wolf1  {
    blast::Array x;
};

host_fn real dist_OBB_point(OBB OBB_test, Vec3 point) {
    Mat3 Rtrans = transpose(OBB_test.R);
 
    Vec3 point_OBB = Rtrans * (point - OBB_test.c);
 
    Vec3 proj = {clamp(point_OBB.x, -OBB_test.e.x, OBB_test.e.x), clamp(point_OBB.y, -OBB_test.e.y, OBB_test.e.y), clamp(point_OBB.z, -OBB_test.e.z, OBB_test.e.z)};
    Array dist_in = (abs(point.x) - OBB_test.e.x, abs(point.y) - OBB_test.e.y, abs(point.z) - OBB_test.e.z);
    real result_if_inside = array_max(dist_in);
    return result_if_inside > 0 ? norm(proj - point_OBB) : result_if_inside;
}
 
host_fn real OBJ_function(Array x, Matrix caps_list, OBB OBB) {
    real t = x[0]*(caps_list).cols;
    int t_step = t;
    int t_step_plus1 = t_step + 1;
    real s = x[1];
 
    Vec3 p1_1 = {(caps_list)(t_step, 0), (caps_list)(t_step, 1), (caps_list)(t_step, 2)};
    Vec3 p1_2 = {(caps_list)(t_step_plus1, 0), (caps_list)(t_step_plus1, 1), (caps_list)(t_step_plus1, 2)};
    Vec3 p2_1 = {(caps_list)(t_step, 3), (caps_list)(t_step, 4), (caps_list)(t_step, 5)};
    Vec3 p2_2 = {(caps_list)(t_step_plus1, 3), (caps_list)(t_step_plus1, 4), (caps_list)(t_step_plus1, 5)};
 
    Vec3 p1 = (1 - (t - t_step)) * p1_1 + (t - t_step) * p1_2;
    Vec3 p2 = (1 - (t - t_step)) * p2_1 + (t - t_step) * p2_2;
 
    Vec3 p = (1 - s) * p1 + s * p2;
 
    return dist_OBB_point(OBB, p);
}
real collision_gwo(Matrix caps_list, OBB OBB) {
    // Initialization of GWO parameters
    const int N_wolves = 500;    // Number of wolves    
    const int N_iterations = 50;    // Number of iterations 
    const auto N_Dimensions = 2;
    Array x_Alpha(N_Dimensions);
    Array x_Beta(N_Dimensions);
    Array x_Delta(N_Dimensions);

    real Alpha_fit = INF_REAL;
    real Beta_fit = INF_REAL;
    real Delta_fit = INF_REAL;
    std::vector<Wolf1> wolf(N_wolves);
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
            real fitness = OBJ_function(wolf[i].x, caps_list, OBB);
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
    real best_f = OBJ_function(x_Alpha, caps_list, OBB);
    return best_f;
}
}