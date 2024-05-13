#include "blast.h"
#include <algorithm>

namespace blast {

struct GAIndividual {
    Array x;          // Position of each individual
    real fitness;     // Fitness of the individual
};

// Solves the collision distance problem with only one box at a time, returns best fitness score
real collision_ga(const Matrix &robot_cartesian_positions, const Box &box, const int n_individuals, const int n_iterations) {
    const auto n_dimensions = 2;
    const double mutation_rate = 0.001; // Mutation rate
    const double elite_percentage = 0.2; // Percentage of elite individuals to retain
    const int elite_count = (int)(elite_percentage * n_individuals);
    Array fitness_fraction(n_individuals);

    Array gbest_x(n_dimensions);
    real gbest_f = -INF_REAL;

    // Initialize random population
    std::vector<GAIndividual> population(n_individuals);
    population.resize(n_individuals);
    for (int i = 4; i < n_individuals; ++i) {
        population[i].x.resize(2);
        for (int j = 0; j < n_dimensions; ++j) {
            population[i].x[j] = 0.5*get_random() + 0.5; // Initialize with random values between [0, 1]
        }
    }
    population[0].x.resize(2);
    population[0].x[0] = 0;
    population[0].x[1] = 0;

    population[1].x.resize(2);
    population[1].x[0] = 0;
    population[1].x[1] = 1;

    population[2].x.resize(2);
    population[2].x[0] = 1;
    population[2].x[1] = 0;

    population[3].x.resize(2);
    population[3].x[0] = 1;
    population[3].x[1] = 1;

    // Main loop
    for (int iteration = 0; iteration < n_iterations; ++iteration) {
        std::vector<GAIndividual> new_population(n_individuals);
        for (int i = 0; i < n_individuals; ++i) {
            new_population[i].x.resize(2);
        }

        // Evaluate fitness
        real fitness_total = 0;
        for (int i = 0; i < n_individuals; i++) {
            population[i].fitness = 1/obj_function(population[i].x, robot_cartesian_positions, box);
            fitness_total += population[i].fitness;
            // Update global best
            gbest_f = (population[i].fitness > gbest_f) ? population[i].fitness : gbest_f;
            gbest_x = (population[i].fitness > gbest_f) ? population[i].x : gbest_x;
        }

        for (int i = 0; i < n_individuals; i++) {
            fitness_fraction[i] = population[i].fitness / fitness_total;
        }
        for (int i = 1; i < n_individuals; i++) {
            fitness_fraction[i] = fitness_fraction[i-1] + fitness_fraction[i];
        }

        std::vector<GAIndividual> sorted_fit;
        sorted_fit = population;
        std::sort(sorted_fit.begin(), sorted_fit.end(), [](const GAIndividual& a, const GAIndividual& b) {
            return a.fitness > b.fitness;
        });

        // Keep elite
        for (int i = 0; i < elite_count/2; ++i) {
            new_population[i].x = sorted_fit[i].x; // Copy elite individuals to the next generation
            new_population[i].fitness = sorted_fit[i].fitness;
        }
        for (int i = elite_count/2; i < elite_count; ++i) {
            for (int j = 0; j < n_dimensions; j++) {
                new_population[i].x[j] = sorted_fit[i-elite_count/2].x[j] + 0.001 * get_random(); // Mutate elite individuals
            }
            new_population[i].x = clamp(new_population[i].x, 0, 1);
        }
        // Perform selection and crossover
        for (int i = elite_count; i < n_individuals; ++i) {
            // Select parents based on their fitness
            real parent_idx1 = 0.5*get_random() + 0.5; // Select parent 1
            real parent_idx2 = 0.5*get_random() + 0.5; // Select parent 2

            Array parent1(n_individuals);
            int fit_frac1;
            for (int j = 0; j < n_individuals; j++) {
                if (fitness_fraction[j] >= parent_idx1) {
                    fit_frac1 = j;
                    parent1 = population[j].x;
                    break;
                }
            }
            Array parent2(n_individuals);
            for (int j = 0; j < n_individuals; j++) {
                if (fitness_fraction[j] > parent_idx2) {
                    parent2 = j == fit_frac1 ? population[(j+1) % n_individuals].x : population[j].x;
                    break;
                }
            }

            for (int j = 0; j < n_dimensions; j++) {
                // Perform crossover
                real lambda = 0.5*get_random() + 0.5;
                new_population[i].x[j] = lambda*parent1[j] + (1-lambda)*parent2[j];
                // Perform mutation
                new_population[i].x[j] = (std::abs(get_random()) < mutation_rate) ? new_population[i].x[j] + 0.001 * get_random() : new_population[i].x[j];
            }
            new_population[i].x = clamp(new_population[i].x, 0, 1);
        }
        // Update population
        population = std::move(new_population);
    }
    return 1/gbest_f;
}

// Calls collision_pso to solve collision distance problem one box at a time, one member at a time. Returns best fitness score
real test_collision_ga_box(const Matrix &robot_cartesian_positions, const World* world, const int n_individuals, const int n_iterations) {
    u32 n_caps = robot_cartesian_positions.rows / 3 - 1;
    real distmin = INF_REAL;
    real current_dist;
    Matrix temp(6, robot_cartesian_positions.cols);
    for (u32 j = 0; j < n_caps; j++) {
        for (u32 i = 0; i < robot_cartesian_positions.cols; i++) {
            temp(0, i) = robot_cartesian_positions(j * 3, i);
            temp(1, i) = robot_cartesian_positions(j * 3 + 1, i);
            temp(2, i) = robot_cartesian_positions(j * 3 + 2, i);
            temp(3, i) = robot_cartesian_positions(j * 3 + 3, i);
            temp(4, i) = robot_cartesian_positions(j * 3 + 4, i);
            temp(5, i) = robot_cartesian_positions(j * 3 + 5, i);
        }
        for (int i = 0; i < world->boxes.size(); i++) {
            current_dist = collision_ga(temp, world->boxes[i], n_individuals, n_iterations);
            distmin = (distmin < 0) ? (current_dist > distmin ? current_dist : distmin) :
                      (current_dist < distmin ? current_dist : distmin);
        }
    }
    // std::cout << distmin << std::endl;
    return distmin;
}

// Solves the collision distance problem with the full world, returns best fitness score
real collision_ga(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations) {
    const auto n_dimensions = 2;
    const double mutation_rate = 0.001; // Mutation rate
    const double mutation_step = 0.001; // Mutation step
    const double elite_percentage = 0.2; // Percentage of elite individuals to retain
    const int elite_count = (int)(elite_percentage * n_individuals);
    Array fitness_fraction(n_individuals);

    Array gbest_x(n_dimensions);
    real gbest_f = -INF_REAL;

    // Initialize random population
    std::vector<GAIndividual> population(n_individuals);
    population.resize(n_individuals);
    for (int i = 4; i < n_individuals; ++i) {
        population[i].x.resize(2);
        for (int j = 0; j < n_dimensions; ++j) {
            population[i].x[j] = 0.5*get_random() + 0.5; // Initialize with random values between [0, 1]
        }
    }
    population[0].x.resize(2);
    population[1].x.resize(2);
    population[2].x.resize(2);
    population[3].x.resize(2);
    population[0].x[0] = 0;
    population[0].x[1] = 0;
    population[1].x[0] = 0;
    population[1].x[1] = 1;
    population[2].x[0] = 1;
    population[2].x[1] = 0;
    population[3].x[0] = 1;
    population[3].x[1] = 1;

    // Main loop
    for (int iteration = 0; iteration < n_iterations; ++iteration) {
        std::vector<GAIndividual> new_population(n_individuals);
        for (int i = 0; i < n_individuals; ++i) {
            new_population[i].x.resize(2);
        }

        // Evaluate fitness
        real fitness_total = 0;
        for (int i = 0; i < n_individuals; i++) {
            population[i].fitness = 1/obj_function(population[i].x, robot_cartesian_positions, world);
            fitness_total += population[i].fitness;
            // Update global best
            gbest_f = (population[i].fitness > gbest_f) ? population[i].fitness : gbest_f;
            gbest_x = (population[i].fitness > gbest_f) ? population[i].x : gbest_x;
        }

        for (int i = 0; i < n_individuals; i++) {
            fitness_fraction[i] = population[i].fitness / fitness_total;
        }
        for (int i = 1; i < n_individuals; i++) {
            fitness_fraction[i] = fitness_fraction[i-1] + fitness_fraction[i];
        }

        std::vector<GAIndividual> sorted_fit;
        sorted_fit = population;
        std::sort(sorted_fit.begin(), sorted_fit.end(), [](const GAIndividual& a, const GAIndividual& b) {
            return a.fitness > b.fitness;
        });

        // Keep elite
        for (int i = 0; i < elite_count/2; ++i) {
            new_population[i].x = sorted_fit[i].x; // Copy elite individuals to the next generation
            new_population[i].fitness = sorted_fit[i].fitness;
        }
        for (int i = elite_count/2; i < elite_count; ++i) {
            for (int j = 0; j < n_dimensions; j++) {
                new_population[i].x[j] = sorted_fit[i-elite_count/2].x[j] + mutation_step * get_random(); // Mutate elite individuals
            }
            new_population[i].x = clamp(new_population[i].x, 0, 1);
        }
        // Perform selection and crossover
        for (int i = elite_count; i < n_individuals; ++i) {
            // Select parents based on their fitness
            real parent_idx1 = 0.5*get_random() + 0.5; // Select parent 1
            real parent_idx2 = 0.5*get_random() + 0.5; // Select parent 2

            Array parent1(n_individuals);
            int fit_frac1;
            for (int j = 0; j < n_individuals; j++) {
                if (fitness_fraction[j] >= parent_idx1) {
                    fit_frac1 = j;
                    parent1 = population[j].x;
                    break;
                }
            }
            Array parent2(n_individuals);
            for (int j = 0; j < n_individuals; j++) {
                if (fitness_fraction[j] > parent_idx2) {
                    parent2 = j == fit_frac1 ? population[(j+1) % n_individuals].x : population[j].x;
                    break;
                }
            }

            for (int j = 0; j < n_dimensions; j++) {
                // Perform crossover
                real lambda = 0.5*get_random() + 0.5;
                new_population[i].x[j] = lambda*parent1[j] + (1-lambda)*parent2[j];
                // Perform mutation
                new_population[i].x[j] = (std::abs(get_random()) < mutation_rate) ? new_population[i].x[j] + mutation_step * get_random() : new_population[i].x[j];
            }
            new_population[i].x = clamp(new_population[i].x, 0, 1);
        }
        // Update population
        population = std::move(new_population);
    }
    return 1/gbest_f;
}

// Calls collision_ga to solve collision distance problem with the full world, one member at a time. Returns best fitness score
real test_collision_ga_world_1caps(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations) {
    u32 n_caps = robot_cartesian_positions.rows / 3 - 1;
    real distmin = INF_REAL;
    real current_dist;
    Matrix temp(6, robot_cartesian_positions.cols);
    for (u32 j = 0; j < n_caps; j++) {
        for (u32 i = 0; i < robot_cartesian_positions.cols; i++) {
            temp(0, i) = robot_cartesian_positions(j * 3, i);
            temp(1, i) = robot_cartesian_positions(j * 3 + 1, i);
            temp(2, i) = robot_cartesian_positions(j * 3 + 2, i);
            temp(3, i) = robot_cartesian_positions(j * 3 + 3, i);
            temp(4, i) = robot_cartesian_positions(j * 3 + 4, i);
            temp(5, i) = robot_cartesian_positions(j * 3 + 5, i);
        }
        current_dist = collision_ga(temp, world, n_individuals, n_iterations);
        distmin = (distmin < 0) ? (current_dist > distmin ? current_dist : distmin) :
                  (current_dist < distmin ? current_dist : distmin);
    }
    return distmin;
}

// Calls collision_ga to solve collision distance problem with the full world, all members at once. Returns best fitness score
real test_collision_ga_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations) {
    real distmin = collision_ga(robot_cartesian_positions, world, n_individuals, n_iterations);
    return distmin;
}

} // namespace blast
