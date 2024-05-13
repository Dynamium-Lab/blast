#include "blast.h"
#include "blast_math.h"


namespace blast {

struct PsoParticle1 {
    Array x;          // Position of each particle i
    Array v;          // Velocity of each particle i
    Array best_x;     // Best known position of each particle i
    real  best_f;     // Best known fitness of each particle i
};

// Solves the collision distance problem with only one box at a time, returns best fitness score and best particle position
real collision_pso(const Matrix &robot_cartesian_positions, const Box &box, const int n_particles, const int n_iterations, Array* best_x = nullptr) {
    //Assert(x.size == optim.bspline->xlen(*optim.task));
    const auto n_dimensions = 2;
    const double w_min = 0.2;          
    const double w_max = 0.9;      // Inertia Weight [0, 1]
    //const double w = 0.5;
    const double c = 2;            // Cognitive weight [0, 2]
    const double s = 2;            // Social weight [0, 2]

    Array gbest_x(n_dimensions);
    real  gbest_f = INF_REAL;

    // Intializing random particles
    std::vector<PsoParticle1> particle(n_particles);
    for(int i = 4; i < n_particles; i++) {
        particle[i].x.resize(2);
        real fraction = n_particles == 1 ? i / (n_particles) : i / (n_particles - 1);
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
        auto r1 = 0.5 * get_random() + 0.5;
        auto r2 = 0.5 * get_random() + 0.5;
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

// Calls collision_pso to solve collision distance problem one box at a time, one member at a time. Returns best best fitness score and best particle position
real test_collision_pso_box(const Matrix &robot_cartesian_positions, const  World* world, const int n_particles, const int n_iterations, Array* best_x = nullptr, int* robot_link = nullptr) {
    int n_caps = robot_cartesian_positions.rows/3 - 1;
    int n_points = robot_cartesian_positions.cols;
    real gbest_f= INF_REAL;
    real temp_dist;
    Array* px;
    Matrix temp(6, n_points);
    for (int j = 0; j < n_caps; j++) {
        for (int i = 0; i < n_points; i++) {
            temp(0, i) = robot_cartesian_positions(j*3, i); 
            temp(1, i) = robot_cartesian_positions(j*3+1, i); 
            temp(2, i) = robot_cartesian_positions(j*3+2, i); 
            temp(3, i) = robot_cartesian_positions(j*3+3, i); 
            temp(4, i) = robot_cartesian_positions(j*3+4, i); 
            temp(5, i) = robot_cartesian_positions(j*3+5, i); 
        }
        for (int i = 0; i < world->boxes.size(); i++) {
            temp_dist = collision_pso(temp, world->boxes[i], n_particles, n_iterations, px);
            *robot_link   = (temp_dist < gbest_f) ? j : *robot_link;
            *best_x = (temp_dist < gbest_f) ? *px : *best_x;
            gbest_f = (temp_dist < gbest_f) ? temp_dist : gbest_f;
        }
    }
    return gbest_f;
}

// Solves the collision distance problem with the full world, returns best fitness score
real collision_pso(const Matrix &robot_cartesian_positions, const World* world, const int n_particles, const int n_iterations) {
    const auto n_dimensions = 2;
    const double w_min = 0.2;          
    const double w_max = 0.9;      // Inertia Weight [0, 1]
    //const double w = 0.5;
    const double c = 2;            // Cognitive weight [0, 2]
    const double s = 2;            // Social weight [0, 2]

    Array gbest_x(n_dimensions);
    real  gbest_f = INF_REAL;

    // Intializing random particles
    std::vector<PsoParticle1> particle(n_particles);
    for(int i = 4; i < n_particles; i++) {
        particle[i].x.resize(2);
        real fraction = n_particles == 1 ? i / (n_particles) : i / (n_particles - 1);
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
        auto r1 = 0.5 * get_random() + 0.5;
        auto r2 = 0.5 * get_random() + 0.5;
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

// Calls collision_pso to solve collision distance problem with the full world, one member at a time. Returns best fitness score
real test_collision_pso_world_1caps(const Matrix &robot_cartesian_positions, const World* world, const int n_particles, const int n_iterations) {
    int n_caps = robot_cartesian_positions.rows / 3 - 1;
    real distmin = INF_REAL;
    real current_dist;
    Matrix temp(6, robot_cartesian_positions.cols);
    for (int j = 0; j < n_caps; j++) {
        for (int i = 0; i < robot_cartesian_positions.cols; i++) {
            temp(0, i) = robot_cartesian_positions(j * 3, i);
            temp(1, i) = robot_cartesian_positions(j * 3 + 1, i);
            temp(2, i) = robot_cartesian_positions(j * 3 + 2, i);
            temp(3, i) = robot_cartesian_positions(j * 3 + 3, i);
            temp(4, i) = robot_cartesian_positions(j * 3 + 4, i);
            temp(5, i) = robot_cartesian_positions(j * 3 + 5, i);
        }
        current_dist = collision_pso(temp, world, n_particles, n_iterations);
        distmin = current_dist < distmin ? current_dist : distmin;
    }
    return distmin;
}

// Calls collision_pso to solve collision distance problem with the full world, all members at once. Returns best fitness score
real test_collision_pso_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_particles, int n_iterations) {
    real distmin = collision_pso(robot_cartesian_positions, world, n_particles, n_iterations);
    return distmin;
}

} // namespace blast