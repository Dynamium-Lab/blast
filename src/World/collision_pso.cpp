#include "blast.h"
#include "world.h"

namespace blast {

struct PsoParticle1 {
    Array x;          // Position of each particle i
    Array v;          // Velocity of each particle i
    Array best_x;     // Best known position of each particle i
    real  best_f;     // Best known fitness of each particle i
};

// Solves the collision distance problem with only one box at a time, returns best fitness score and best particle position
real collision_pso(const Matrix &robot_cartesian_positions, const Box &box, const int n_particles, const int n_iterations, Array* best_x) {
    //Assert(x.size == optim.bspline->xlen(*optim.task));
    const auto n_dimensions = 2;
    const real w_min = (real)0.2;
    const real w_max = (real)0.9;      // Inertia Weight [0, 1]
    //const double w = 0.5;
    const real c = (real)2;            // Cognitive weight [0, 2]
    const real s = (real)2;            // Social weight [0, 2]

    Array gbest_x(n_dimensions);
    real  gbest_f = INF_REAL;

    // Intializing random particles
    std::vector<PsoParticle1> particle(n_particles);
    for(int i = 4; i < n_particles; i++) {
        particle[i].x.resize(2);
        real fraction = n_particles == 1 ? (real)i : (real)i / (n_particles - 1);
        for (int j = 0; j < n_dimensions; j++) {
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

    for (int i = 0; i < n_particles; i++) {
        particle[i].v = random_array(n_dimensions, 1);

        particle[i].best_x = particle[i].x; // Best position for Pi
        particle[i].best_f = obj_function(particle[i].x, robot_cartesian_positions, box); // Best fitness for Pi
    }

    // Main loop
    for (int j = 0; j < n_iterations; j++) {
        real r1 = (real)0.5 * get_random() + (real)0.5;
        real r2 = (real)0.5 * get_random() + (real)0.5;
        real w = w_max - j * ((w_max - w_min) / (real)n_iterations);

        for(int i = 0; i < n_particles; i++) {
            for (int k = 0; k < n_dimensions; k++) {
                //Updating Velocity
                particle[i].v[k] = w * particle[i].v[k] + c * r1 * (particle[i].best_x[k]- particle[i].x[k]) + s* r2 * (gbest_x[k] - particle[i].x[k]);
                particle[i].v[k] = clamp(particle[i].v[k], -1, 1);

                // Updating position
                particle[i].x[k] += particle[i].v[k];
                particle[i].x[k] = clamp(particle[i].x[k], 0, 1);
            }
            //particle[i].x.back() = clamp(particle[i].x.back(), 0.1, 10);
            // Evaluate fitness
            real fitness = obj_function(particle[i].x, robot_cartesian_positions, box);
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

// Solves the collision distance problem with the full world, returns best fitness score
real collision_pso(const Matrix &robot_cartesian_positions, const World* world, const int n_particles, const int n_iterations) {
    const auto n_dimensions = 2;
    const real w_min = (real)0.2;
    const real w_max = (real)0.9;      // Inertia Weight [0, 1]
    //const double w = 0.5;
    const real c = 2;            // Cognitive weight [0, 2]
    const real s = 2;            // Social weight [0, 2]

    Array gbest_x(n_dimensions);
    real  gbest_f = INF_REAL;

    // Intializing random particles
    std::vector<PsoParticle1> particle(n_particles);
    for(int i = 4; i < n_particles; i++) {
        particle[i].x.resize(2);
        real fraction = n_particles == 1 ? (real)i : (real)i / (n_particles - 1);
        for (int j = 0; j < n_dimensions; j++) {
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

    for (int i = 0; i < n_particles; i++) {
        particle[i].v = random_array(n_dimensions, 1);

        particle[i].best_x = particle[i].x; // Best position for Pi
        particle[i].best_f = obj_function(particle[i].x, robot_cartesian_positions, world); // Best fitness for Pi
    }

    // Main loop
    for (int j = 0; j < n_iterations; j++) {
        real r1 = (real)0.5 * get_random() + (real)0.5;
        real r2 = (real)0.5 * get_random() + (real)0.5;
        real w = w_max - j * ((w_max - w_min) / n_iterations);

        for(int i = 0; i < n_particles; i++) {
            for (int k = 0; k < n_dimensions; k++) {
                //Updating Velocity
                particle[i].v[k] = w * particle[i].v[k] + c * r1 * (particle[i].best_x[k]- particle[i].x[k]) + s* r2 * (gbest_x[k] - particle[i].x[k]);
                particle[i].v[k] = clamp(particle[i].v[k], -1, 1);

                // Updating position
                particle[i].x[k] += particle[i].v[k];
                particle[i].x[k] = clamp(particle[i].x[k], 0, 1);
            }
            // Evaluate fitness
            real fitness = obj_function(particle[i].x, robot_cartesian_positions, world);

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