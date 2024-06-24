#pragma once

namespace blast {


real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const Box &box);
real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const World* world);
real obj_function_gradient(const Array& x, const Matrix &robot_cartesian_positions, const World* world, Array* gradient);

// PSO Functions
real collision_pso(const Matrix &robot_cartesian_positions, const Box &box, const int n_particles, const int n_iterations, Array* best_x = nullptr);
real test_collision_pso_box(const Matrix &robot_cartesian_positions, const  World* world, const int n_particles, const int n_iterations, Array* best_x = nullptr, int* robot_link = nullptr);
real test_collision_pso_world_1caps(const Matrix &robot_cartesian_positions, const World* world, const int n_particles, const int n_iterations);
real test_collision_pso_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_particles, int n_iterations);

// GWO Functions
real collision_gwo(const Matrix &robot_cartesian_positions, Box box, const int n_wolves, const int n_iterations);
real test_collision_gwo_box(const Matrix &robot_cartesian_positions, World* world, const int n_wolves, const int n_iterations);
real test_collision_gwo_world_1caps(const Matrix &robot_cartesian_positions, const World* world, const int n_wolves, const int n_iterations);
real test_collision_gwo_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_wolves, int n_iterations);

// GA Functions
real collision_ga(const Matrix &robot_cartesian_positions, const Box &box, const int n_individuals, const int n_iterations);
real test_collision_ga_box(const Matrix &robot_cartesian_positions, const World* world, const int n_individuals, const int n_iterations);
real test_collision_ga_world_1caps(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations);
real test_collision_ga_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations);


} // namespace blast