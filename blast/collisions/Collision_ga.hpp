#pragma once
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "collisions/world.h"

namespace blast {

struct GAIndividual {
    Array x;          // Position of each individual
    real fitness;     // Fitness of the individual
};
// Rosenbrock function to test 
double rosenbrock(Array x) {
    const double a = 1.0;
    const double b = 100.0;
    return pow(a - x[1], 2) + b * pow(x[2] - x[1] * x[1], 2);
}
real collision_ga(Matrix caps_list, OBB OBB, int N_individuals, int N_iterations) {
    const auto N_Dimensions = 2;
    const double mutation_rate = 0.1; // Mutation rate
    const double elite_percentage = 0.1; // Percentage of elite individuals to retain

    Array gbest_x(N_Dimensions);
    real gbest_f = INF_REAL;

    // Initialize random population
    std::vector<GAIndividual> population(N_individuals);
    for (int i = 0; i < N_individuals; ++i) {
        population[i].x.resize(2);
        for (int j = 0; j < N_Dimensions; ++j) {
            population[i].x[j] = get_random(); // Initialize with random values between [0, 1]
        }
        population[i].fitness = OBJ_function(population[i].x, caps_list, OBB); 
    }
    // Main loop
    for (int iteration = 0; iteration < N_iterations; ++iteration) {
        // Sort population based on fitness
        std::sort(population.begin(), population.end(), [](const GAIndividual& a, const GAIndividual& b) {
            return a.fitness < b.fitness;
        });

        // Perform selection and crossover
        std::vector<GAIndividual> new_population(N_individuals);
        int elite_count = elite_percentage * N_individuals;
        for (int i = 0; i < elite_count; ++i) {
            new_population[i] = population[i]; // Copy elite individuals to the next generation
        }
        for (int i = elite_count; i < N_individuals; ++i) {
            // Select parents from elite individuals
            int parent_idx1 = rand() % elite_count; // Select parent 1 from elite individuals
            int parent_idx2 = rand() % elite_count; // Select parent 2 from elite individuals
            for (int j = 0; j < N_Dimensions; j++) {
                // Perform crossover
                new_population[i].x[j] = (population[parent_idx1].x[j] + population[parent_idx2].x[j]) / 2.0;
                // Perform mutation
                if ((double)rand() / RAND_MAX < mutation_rate) {
                    new_population[i].x[j] += 0.1 * (((double)rand() / RAND_MAX) - 0.5);
                }
            }
            // Evaluate fitness
            population[i].fitness = OBJ_function(population[i].x, caps_list, OBB); 
            // Update global best
            if (new_population[i].fitness < gbest_f) {
                gbest_f = new_population[i].fitness;
                gbest_x = new_population[i].x;
            }
        }
        // Update population
        population = std::move(new_population);
    }
    printf("GA finished \n");
    printf("Global best fitness = %f \n", gbest_f);    
    return gbest_f;
}

real test_collision_ga(Matrix cart_pos, objlist* world, int N_individuals, int N_iterations) {
    int n_caps = cart_pos.rows / 3 - 1;
    real min_dist = INF_REAL;
    real temp_dist;
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
        for (int i = 0; i < world->OBBlist.size(); i++) {
            temp_dist = collision_ga(temp, world->OBBlist[i], N_individuals, N_iterations);
            min_dist = std::min(temp_dist, min_dist);
        }
    }
    return min_dist;
}

}
