#pragma once

#include "blast.h"

namespace blast {

struct Box;
struct Capsule;
struct Sphere;
struct Cylinder;
struct World;

struct World {
    std::vector<Box>      boxes;
    std::vector<Cylinder> cylinders;
    std::vector<Sphere>   spheres;
    std::vector<Capsule>  capsules;
};

struct Box {
    Vec3 c; // Box center point
    Vec3 e; // Positive halfwidth extents of Box along each axis
    Mat3 R; // Local x-, y-, and z-axes (Rotation matrix)
};

struct Capsule {
    Vec3 p1;
    Vec3 p2;
    real r;
};

struct Sphere {
    Vec3 c;
    real r;
};

struct Cylinder {
    Vec3 p1;
    Vec3 p2;
    real r;
};


// todo: What do we do with cylinder?

void add_box(Vec3 center_point, Vec3 half_width, Mat3 rotation_matrix, World* world);
void add_sphere(Vec3 center_point, real radius, World* world);
void add_cylinder(Vec3 point1, Vec3 point2, real radius, World* world);
void add_capsule(Vec3 point1, Vec3 point2, real radius, World* world);
std::vector<real> test_collision(std::vector<Capsule>* robot_capsule_list, World* world, int n_lowest_distances);

real distance(Capsule capsule, Cylinder cylinder);
real distance(Capsule capsule, Sphere sphere);
real distance(Capsule capsule1, Capsule capsule2);
real distance(Capsule capsule, Box box);
// ======================================
//            Optimization functions
// ======================================

real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const Box &box);
real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const World* world);

// PSO Functions
real collision_pso(const Matrix &robot_cartesian_positions, const Box &box, const int n_particles, const int n_iterations, Array* best_x = nullptr);
real test_collision_pso_box(const Matrix &robot_cartesian_positions, const  World* world, const int n_particles, const int n_iterations, Array* best_x = nullptr, int* robot_link = nullptr);
real collision_pso(const Matrix &robot_cartesian_positions, const World* world, const int n_particles, const int n_iterations);
real test_collision_pso_world_1caps(const Matrix &robot_cartesian_positions, const World* world, const int n_particles, const int n_iterations);
real test_collision_pso_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_particles, int n_iterations);

// GWO Functions
real collision_gwo(const Matrix &robot_cartesian_positions, Box box, const int n_wolves, const int n_iterations);
real test_collision_gwo_box(const Matrix &robot_cartesian_positions, World* world, const int n_wolves, const int n_iterations);
real collision_gwo(const Matrix &robot_cartesian_positions, const World* world, const int n_wolves, const int n_iterations);
real test_collision_gwo_world_1caps(const Matrix &robot_cartesian_positions, const World* world, const int n_wolves, const int n_iterations);
real test_collision_gwo_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_wolves, int n_iterations);

// GA Functions
real collision_ga(const Matrix &robot_cartesian_positions, const Box &box, const int n_individuals, const int n_iterations);
real test_collision_ga_box(const Matrix &robot_cartesian_positions, const World* world, const int n_individuals, const int n_iterations);
real collision_ga(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations);
real test_collision_ga_world_1caps(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations);
real test_collision_ga_world_full_robot(const Matrix &robot_cartesian_positions, World* world, int n_individuals, int n_iterations);




} // namespace blast