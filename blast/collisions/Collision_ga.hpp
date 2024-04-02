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
// double rosenbrock(Array x) {
//     const double a = 1.0;
//     const double b = 100.0;
//     return pow(a - x[1], 2) + b * pow(x[2] - x[1] * x[1], 2);
// }
real collision_ga(Matrix caps_list, OBB OBB, int N_individuals, int N_iterations) {
    const auto N_Dimensions = 2;
    const double mutation_rate = 0.001; // Mutation rate
    const double elite_percentage = 0.1; // Percentage of elite individuals to retain
    const int elite_count = elite_percentage * N_individuals;
    Array fitness_fraction(N_individuals);

    Array gbest_x(N_Dimensions);
    real gbest_f = -INF_REAL;

    // Initialize random population 
    std::vector<GAIndividual> population(N_individuals);
    population.resize(N_individuals);
    for (int i = 0; i < N_individuals; ++i) {
        population[i].x.resize(2);
        for (int j = 0; j < N_Dimensions; ++j) {
            population[i].x[j] = 0.5*get_random() + 0.5; // Initialize with random values between [0, 1]
        } 
    }

    // Main loop
    for (int iteration = 0; iteration < N_iterations; ++iteration) {
        std::vector<GAIndividual> new_population(N_individuals);
        for (int i = 0; i < N_individuals; ++i) {
            new_population[i].x.resize(2);
        }

        // Evaluate fitness
        real fitness_total = 0;
        for (int i = 0; i < N_individuals; i++) {
            population[i].fitness = 1/OBJ_function(population[i].x, caps_list, OBB);
            fitness_total += population[i].fitness;
            // Update global best
            if (population[i].fitness > gbest_f) {
                gbest_f = population[i].fitness;
                gbest_x = population[i].x;
            }
        }

        for (int i = 0; i < N_individuals; i++) {
            fitness_fraction[i] = population[i].fitness / fitness_total;
        }
        for (int i = 1; i < N_individuals; i++) {
            fitness_fraction[i] = fitness_fraction[i-1] + fitness_fraction[i];
        }

        std::vector<GAIndividual> sorted_fit;
        sorted_fit = population;
        std::sort(sorted_fit.begin(), sorted_fit.end(), [](const GAIndividual& a, const GAIndividual& b) {
            return a.fitness > b.fitness;
        });

        // Keep elite
        for (int i = 0; i < elite_count; ++i) {
            new_population[i].x = population[i].x; // Copy elite individuals to the next generation
            new_population[i].fitness = population[i].fitness;
        }
        // Perform selection and crossover
        for (int i = elite_count; i < N_individuals; ++i) {
            // Select parents based on their fitness
            real parent_idx1 = 0.5*get_random() + 0.5; // Select parent 1
            real parent_idx2 = 0.5*get_random() + 0.5; // Select parent 2

            Array parent1(N_individuals);
            parent1 = population[0].x;
            for (int j = 0; j < N_individuals; j++) {
                if (fitness_fraction[j] > parent_idx1) {
                    parent1 = population[j].x;
                    break;
                }
            }
            Array parent2(N_individuals);
            parent2 = population[0].x;
            for (int j = 0; j < N_individuals; j++) {
                if (fitness_fraction[j] > parent_idx2) {
                    parent2 = population[j].x;
                    break;
                }
            }

            for (int j = 0; j < N_Dimensions; j++) {
                // Perform crossover
                real lambda = 0.5*get_random() + 0.5;
                new_population[i].x[j] = lambda*parent1[j] + (1-lambda)*parent2[j];
                // Perform mutation
                if (abs(get_random()) < mutation_rate) {
                    new_population[i].x[j] += 0.1 * (0.5*get_random());
                }
            }
            new_population[i].x = clamp(new_population[i].x, 0, 1);
        }
        // Update population
        population = std::move(new_population);
    }  
    return 1/gbest_f;
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
