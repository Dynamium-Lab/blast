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
Array test_collision(std::vector<Capsule> robot_capsule_list, World* world, u32 n_lowest_distances);

real distance(Capsule capsule,  Cylinder cylinder);
real distance(Capsule capsule,  Sphere sphere);
real distance(Capsule capsule1, Capsule capsule2);
real distance(Capsule capsule,  Box box);

// ======================================
//            Optimization functions
// ======================================


real collision_pso(const Matrix &robot_cartesian_positions, const World* world, int n_particles,   int n_iterations);
real collision_gwo(const Matrix &robot_cartesian_positions, const World* world, int n_wolves,      int n_iterations);
real collision_ga(const Matrix &robot_cartesian_positions,  const World* world, int n_individuals, int n_iterations);

Vec3 get_point(const Array &x, const Matrix &robot_cartesian_positions);
real distance(const Box &box, const Vec3 &point);
real distance(const Capsule &capsule, const Vec3 &point);
real distance(const Sphere &sphere, const  Vec3 &point);


} // namespace blast