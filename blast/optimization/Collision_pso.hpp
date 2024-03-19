#pragma once
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "collisions/world.h"

namespace blast {
    
host_fn real dist_OBB_point(OBB OBB_test, Vec3 point) {
    Mat3 Rtrans = transpose(OBB_test.R);
 
    Vec3 point_OBB = Rtrans * (point - OBB_test.c);
 
    Vec3 proj = {clamp(point_OBB.x, -OBB_test.e.x, OBB_test.e.x), clamp(point_OBB.y, -OBB_test.e.y, OBB_test.e.y), clamp(point_OBB.z, -OBB_test.e.z, OBB_test.e.z)};
    Array dist_in = (abs(point.x) - OBB_test.e.x, abs(point.y) - OBB_test.e.y, abs(point.z) - OBB_test.e.z);
    real result_if_inside = array_max(dist_in);
    return result_if_inside > 0 ? norm(proj - point_OBB) : result_if_inside;
}
host_fn real OBJ_function(Array& x, Matrix caps_list, OBB OBB) {
    real t = x[0]*((caps_list).cols-1);
    int t_step = t;
    int t_step_plus1 = t_step == (caps_list).cols-1 ? t_step : t_step + 1;
    real s = x[1];
 
    Vec3 p1_1 = {caps_list(0, t_step), caps_list(1,t_step), caps_list(2, t_step)};
    Vec3 p1_2 = {caps_list(0, t_step_plus1), caps_list(1, t_step_plus1), caps_list(2, t_step_plus1)};
    Vec3 p2_1 = {caps_list(3, t_step), caps_list(4, t_step), caps_list(5, t_step)};
    Vec3 p2_2 = {caps_list(3, t_step_plus1), caps_list(4, t_step_plus1), caps_list(5, t_step_plus1)};
 
    Vec3 p1 = (1 - (t - t_step)) * p1_1 + (t - t_step) * p1_2;
    Vec3 p2 = (1 - (t - t_step)) * p2_1 + (t - t_step) * p2_2;
 
    Vec3 p = (1 - s) * p1 + s * p2;
 
    return dist_OBB_point(OBB, p);
}
struct PsoParticle1 {
    Array x;          // Position of each particle i
    Array v;          // Velocity of each particle i
    Array best_x;     // Best known position of each particle i
    real  best_f;     // Best known fitness of each particle i
};

real collision_pso(Matrix caps_list, OBB OBB) {
    //Assert(x.size == optim.bspline->xlen(*optim.task));
    const auto N_Dimensions = 2;
    const int N_particles  = 400;
    const int N_iterations = 50;
    const double w_min = 0.2;          
    const double w_max = 0.9;      // Inertia Weight [0, 1]
    //const double w = 0.5;
    const double c = 2;            // Cognitive weight [0, 2]
    const double s = 2;            // Social weight [0, 2]

    Array gbest_x(N_Dimensions);
    real  gbest_f = INF_REAL;

    std::vector<PsoParticle1> particle(N_particles);
    for(int i = 0; i < N_particles; i++) {
        particle[i].x = random_array(N_Dimensions, 1);
        particle[i].x = clamp(particle[i].x, 0, 1);

        particle[i].v = random_array(N_Dimensions, 1);

        particle[i].best_x = particle[i].x; // Best position for Pi
        particle[i].best_f = OBJ_function(particle[i].x, caps_list, OBB); // Best fitness for Pi
    }

    for (int j = 0; j < N_iterations; j++) {
        auto r1 = 0.5 * get_random() + 0.5;
        auto r2 = 0.5 * get_random() + 0.5;
        real w = w_max - j * ((w_max - w_min) / N_iterations);

        for(int i = 0; i < N_particles; i++) {
            for (int k = 0; k < N_Dimensions; k++) {
                //Updating Velocity
                particle[i].v[k] = w * particle[i].v[k] + c * r1 * (particle[i].best_x[k]- particle[i].x[k]) + s* r2 * (gbest_x[k] - particle[i].x[k]);
                particle[i].v[k] = clamp(particle[i].v[k], -1, 1);

                // Updating position
                particle[i].x[k] += particle[i].v[k];
                particle[i].x[k] = clamp(particle[i].x[k], 0, 1);

  
            }
            //particle[i].x.back() = clamp(particle[i].x.back(), 0.1, 10);
            // Evaluate fitness
            real fitness = OBJ_function(particle[i].x, caps_list, OBB);
            // printf("fitness with penalty = %f \n", fitness);

            // Updating local best
            if (fitness < particle[i].best_f) {
                particle[i].best_f = fitness;
                particle[i].best_x = particle[i].x;
            }
            // Updating global best
            if (fitness < gbest_f) {
                gbest_f = fitness;
                gbest_x = particle[i].best_x;
            }
        }
    }

    // printf("global best f: %f\n", gbest_f);
    // printf("global best time: %f\n", x.back());
    // printf("Global best x: ");
    // print(x);

    return gbest_f;
}
}