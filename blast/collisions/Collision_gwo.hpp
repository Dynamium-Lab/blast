#pragma once
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "collisions/world.h"
namespace blast {

struct Wolf1  {
    blast::Array x;
};

real collision_gwo(const Matrix &caps_list, OBB OBB, const int N_wol, const int N_it) {
    // Initialization of GWO parameters
    const int N_wolves = N_wol;    // Number of wolves    
    const int N_iterations = N_it;    // Number of iterations 
    const auto N_Dimensions = 2;
    Array x_Alpha(N_Dimensions);
    Array x_Beta(N_Dimensions);
    Array x_Delta(N_Dimensions);

    real Alpha_fit = INF_REAL;
    real Beta_fit = INF_REAL;
    real Delta_fit = INF_REAL;
    std::vector<Wolf1> wolf(N_wolves);
    // Initialization of wolves' positions
    for (int i = 4; i < N_wolves; i++) {
        wolf[i].x.resize(2);
        real fraction = (N_Dimensions == 1) ? i / (N_Dimensions) : i / (N_Dimensions - 1);
        for (int j = 0; j < N_Dimensions; j++) {
            wolf[i].x[j] = fraction;
            wolf[i].x[j] = 0.5*get_random() + 0.5;
        }
    }
    wolf[0].x.resize(2);
    wolf[1].x.resize(2);
    wolf[2].x.resize(2);
    wolf[3].x.resize(2);
    wolf[0].x[0] = 0;
    wolf[0].x[1] = 0;
    wolf[1].x[0] = 0;
    wolf[1].x[1] = 1;
    wolf[2].x[0] = 1;
    wolf[2].x[1] = 0;
    wolf[3].x[0] = 1;
    wolf[3].x[1] = 1;
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
            wolf[i].x = clamp(wolf[i].x, 0, 1);
        }
    }
    real best_f = OBJ_function(x_Alpha, caps_list, OBB);
    return best_f;
}

real test_collision_gwo_OBB(const Matrix &cart_pos, objlist* world, const int N_wol, const int N_it) {
    int n_caps = cart_pos.rows/3 - 1;
    int n_points = cart_pos.cols;
    real min_dist = INF_REAL;
    real temp_dist;
    Matrix temp(6, n_points);
    for (int j = 0; j < n_caps; j++) {
        for (int i = 0; i < n_points; i++) {
            temp(0, i) = cart_pos(j*3, i); 
            temp(1, i) = cart_pos(j*3+1, i); 
            temp(2, i) = cart_pos(j*3+2, i); 
            temp(3, i) = cart_pos(j*3+3, i); 
            temp(4, i) = cart_pos(j*3+4, i); 
            temp(5, i) = cart_pos(j*3+5, i); 
        } 
        for (int i = 0; i < world->OBBlist.size(); i++) {
            temp_dist = collision_gwo(temp, world->OBBlist[i], N_wol, N_it);
            min_dist = temp_dist < min_dist ? temp_dist : min_dist;
        }
    }
    return min_dist;
}

real collision_gwo(const Matrix &caps_list, const objlist* world, const int N_wol, const int N_it) {
    // Initialization of GWO parameters
    const int N_wolves = N_wol;    // Number of wolves    
    const int N_iterations = N_it;    // Number of iterations 
    const auto N_Dimensions = 2;
    Array x_Alpha(N_Dimensions);
    Array x_Beta(N_Dimensions);
    Array x_Delta(N_Dimensions);

    real Alpha_fit = INF_REAL;
    real Beta_fit = INF_REAL;
    real Delta_fit = INF_REAL;
    std::vector<Wolf1> wolf(N_wolves);
    // Initialization of wolves' positions
    for (int i = 4; i < N_wolves; i++) {
        wolf[i].x.resize(2);
        real fraction = (N_Dimensions == 1) ? i / (N_Dimensions) : i / (N_Dimensions - 1);
        for (int j = 0; j < N_Dimensions; j++) {
            wolf[i].x[j] = fraction;
            wolf[i].x[j] = 0.5*get_random() + 0.5;
        }
    }
    wolf[0].x.resize(2);
    wolf[0].x[0] = 0;
    wolf[0].x[1] = 0;

    wolf[1].x.resize(2);
    wolf[1].x[0] = 0;
    wolf[1].x[1] = 1;

    wolf[2].x.resize(2);
    wolf[2].x[0] = 1;
    wolf[2].x[1] = 0;

    wolf[3].x.resize(2);
    wolf[3].x[0] = 1;
    wolf[3].x[1] = 1;

    // GWO Computation
    for (int k = 0; k < N_iterations; k++) {
        real a = 2.0 - k * (2.0 / N_iterations);
        for (int i = 0; i < N_wolves; i++) {
            // Evaluating fitness and deciding which wolves are Alpha, Beta and Delta
            real fitness = OBJ_function(wolf[i].x, caps_list, world);
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
            wolf[i].x = clamp(wolf[i].x, 0, 1);
        }
    }
    real best_f = OBJ_function(x_Alpha, caps_list, world);
    return best_f;
}

real test_collision_gwo_world_1caps(const Matrix &cart_pos, const objlist* world, const int N_individuals, const int N_iterations) {
    int n_caps = cart_pos.rows / 3 - 1;
    real distmin = INF_REAL;
    real current_dist;
    Matrix temp(6, cart_pos.cols);
    for (int j = 0; j < n_caps; j++) {
        for (int i = 0; i < cart_pos.cols; i++) {
            temp(0, i) = cart_pos(j * 3, i);
            temp(1, i) = cart_pos(j * 3 + 1, i);
            temp(2, i) = cart_pos(j * 3 + 2, i);
            temp(3, i) = cart_pos(j * 3 + 3, i);
            temp(4, i) = cart_pos(j * 3 + 4, i);
            temp(5, i) = cart_pos(j * 3 + 5, i);
        }
        current_dist = collision_gwo(temp, world, N_individuals, N_iterations);
        distmin = (distmin < 0) ? (current_dist > distmin ? current_dist : distmin) :
                                (current_dist < distmin ? current_dist : distmin);
    }
    return distmin;
}

real test_collision_gwo_world_full_robot(const Matrix &cart_pos, objlist* world, int N_particles, int N_iterations) {
    real distmin = collision_gwo(cart_pos, world, N_particles, N_iterations);
    return distmin;
}

} // namespace blast