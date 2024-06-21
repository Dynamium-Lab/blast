#pragma once
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "collisions/world.h"

namespace blast {
    
struct PsoParticle1 {
    Array x;          // Position of each particle i
    Array v;          // Velocity of each particle i
    Array best_x;     // Best known position of each particle i
    real  best_f;     // Best known fitness of each particle i
};

// Solves the collision distance problem with only one box at a time, returns best fitness score and best particle position
real collision_pso(const Matrix &caps_list, const box &box, const int N_particles, const int N_iterations, Array* best_x = nullptr) {
    //Assert(x.size == optim.bspline->xlen(*optim.task));
    const auto N_Dimensions = 2;
    const double w_min = 0.2;          
    const double w_max = 0.9;      // Inertia Weight [0, 1]
    //const double w = 0.5;
    const double c = 2;            // Cognitive weight [0, 2]
    const double s = 2;            // Social weight [0, 2]

    Array gbest_x(N_Dimensions);
    real  gbest_f = INF_REAL;

    // Intializing random particles
    std::vector<PsoParticle1> particle(N_particles);
    for(int i = 4; i < N_particles; i++) {
        particle[i].x.resize(2);
        real fraction = N_particles == 1 ? i / (N_particles) : i / (N_particles - 1);
        for (int j = 0; j < N_Dimensions; j++) {
            particle[i].x[j] = fraction;
        }
    }

    particle[0].x.resize(2);
    particle[0].x[0] = 0;
    particle[0].x[1] = 0;

    particle[1].x.resize(2);
    particle[1].x[0] = 0;
    particle[1].x[1] = 1;

    particle[2].x.resize(2);
    particle[2].x[0] = 1;
    particle[2].x[1] = 0;

    particle[3].x.resize(2);
    particle[3].x[0] = 1;
    particle[3].x[1] = 1;

    for (int i = 0; i < N_particles; i++) {
        particle[i].v = random_array(N_Dimensions, 1);

        particle[i].best_x = particle[i].x; // Best position for Pi
        particle[i].best_f = OBJ_function(particle[i].x, caps_list, box); // Best fitness for Pi
    }

    // Main loop
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
            real fitness = OBJ_function(particle[i].x, caps_list, box);
            // printf("fitness with penalty = %f \n", fitness);

            // Updating local best
            particle[i].best_x = (fitness < particle[i].best_f) ? particle[i].x : particle[i].best_x;
            particle[i].best_f = (fitness < particle[i].best_f) ? fitness : particle[i].best_f;
            // Updating global best
            gbest_x = (fitness < gbest_f) ? particle[i].best_x : gbest_x;
            gbest_f = (fitness < gbest_f) ? fitness : gbest_f;
        }
    }
    
    *best_x = gbest_x;
    return gbest_f;
}

// Calls collision_pso to solve collision distance problem one box at a time, one member at a time. Returns best best fitness score and best particle position
real test_collision_pso_OBB(const Matrix &cart_pos, const  objlist* world, const int N_particles, const int N_iterations, Array* best_x = nullptr, int* caps = nullptr) {
    int n_caps = cart_pos.rows/3 - 1;
    int n_points = cart_pos.cols;
    real gbest_f= INF_REAL;
    real temp_dist;
    Array x(2);
    Array* px = &x;
    Matrix temp(6, n_points);
    // for (int j = 0; j < n_caps; j++) {
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
            temp_dist = collision_pso(temp, world->OBBlist[i], N_particles, N_iterations, px);
            *caps   = (temp_dist < gbest_f) ? j : *caps;
            *best_x = (temp_dist < gbest_f) ? *px : *best_x;
            gbest_f = (temp_dist < gbest_f) ? temp_dist : gbest_f;
        }
    }
    return gbest_f;
}

// Solves the collision distance problem with the full world, returns best fitness score
real collision_pso(const Matrix &caps_list, const objlist* world, const int N_particles, const int N_iterations) {
    const auto N_Dimensions = 2;
    const double w_min = 0.2;          
    const double w_max = 0.9;      // Inertia Weight [0, 1]
    //const double w = 0.5;
    const double c = 2;            // Cognitive weight [0, 2]
    const double s = 2;            // Social weight [0, 2]

    Array gbest_x(N_Dimensions);
    real  gbest_f = INF_REAL;

    // Intializing random particles
    std::vector<PsoParticle1> particle(N_particles);
    for(int i = 4; i < N_particles; i++) {
        particle[i].x.resize(2);
        real fraction = N_particles == 1 ? i / (N_particles) : i / (N_particles - 1);
        for (int j = 0; j < N_Dimensions; j++) {
            particle[i].x[j] = fraction;
        }
    }
    particle[0].x.resize(2);
    particle[0].x[0] = 0;
    particle[0].x[1] = 0;

    particle[1].x.resize(2);
    particle[1].x[0] = 0;
    particle[1].x[1] = 1;

    particle[2].x.resize(2);
    particle[2].x[0] = 1;
    particle[2].x[1] = 0;

    particle[3].x.resize(2);
    particle[3].x[0] = 1;
    particle[3].x[1] = 1;

    for (int i = 0; i < N_particles; i++) {
        particle[i].v = random_array(N_Dimensions, 1);

        particle[i].best_x = particle[i].x; // Best position for Pi
        particle[i].best_f = OBJ_function(particle[i].x, caps_list, world); // Best fitness for Pi
    }

    // Main loop
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
            real fitness = OBJ_function(particle[i].x, caps_list, world);
            // printf("fitness with penalty = %f \n", fitness);

            // Updating local best
            particle[i].best_x = (fitness < particle[i].best_f) ? particle[i].x : particle[i].best_x;
            particle[i].best_f = (fitness < particle[i].best_f) ? fitness : particle[i].best_f;
            // Updating global best
            gbest_x = (fitness < gbest_f) ? particle[i].best_x : gbest_x;
            gbest_f = (fitness < gbest_f) ? fitness : gbest_f;
        }
    }
    return gbest_f;
}

} // namespace blast