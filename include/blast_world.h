#pragma once

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

// ======================================
//            Optimization functions
// ======================================

// PSO Functions
real collision_pso(const Matrix &caps_list, const box &box, const int n_particles, const int n_iterations, Array* best_x = nullptr);
real test_collision_pso_OBB(const Matrix &cart_pos, const  objlist* world, const int n_particles, const int n_iterations, Array* best_x = nullptr, int* caps = nullptr);
real collision_pso(const Matrix &caps_list, const objlist* world, const int n_particles, const int n_iterations);
real test_collision_pso_world_1caps(const Matrix &cart_pos, const objlist* world, const int n_particles, const int n_iterations);
real test_collision_pso_world_full_robot(const Matrix &cart_pos, objlist* world, int n_particles, int n_iterations);

// GWO Functions
real collision_gwo(const Matrix &caps_list, box box, const int n_wolves, const int n_iterations);
real test_collision_gwo_OBB(const Matrix &cart_pos, objlist* world, const int n_wolves, const int n_iterations);
real collision_gwo(const Matrix &caps_list, const objlist* world, const int n_wolves, const int n_iterations);
real test_collision_gwo_world_1caps(const Matrix &cart_pos, const objlist* world, const int n_wolves, const int n_iterations);
real test_collision_gwo_world_full_robot(const Matrix &cart_pos, objlist* world, int n_wolves, int n_iterations);

// GA Functions
real collision_ga(const Matrix &caps_list, const box &box, const int n_individuals, const int n_iterations);
real test_collision_ga_OBB(const Matrix &cart_pos, const objlist* world, const int n_individuals, const int n_iterations);
real collision_ga(const Matrix &caps_list, objlist* world, int n_individuals, int n_iterations);
real test_collision_ga_world_1caps(const Matrix &cart_pos, objlist* world, int n_individuals, int n_iterations);
real test_collision_ga_world_full_robot(const Matrix &cart_pos, objlist* world, int n_individuals, int n_iterations);




} // namespace blast