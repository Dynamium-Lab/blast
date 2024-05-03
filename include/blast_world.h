#pragma once

namespace blast {
using capslist = std::vector<Capsule>;

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

struct capslist {
    std::vector<Capsule> capsules;
};

// todo: What do we do with cylinder?

void add_box(Vec3 center_point, Vec3 half_width, Mat3 rotation_matrix, World* world);
void add_sphere(Vec3 center_point, real radius, World* world);
void add_cylinder(Vec3 point1, Vec3 point2, real radius, World* world);
void add_capsule(Vec3 point1, Vec3 point2, real radius, World* world);
std::vector<real> test_collision(capslist* robot_capsule_list, World* world, int n_lowest_distances);



} // namespace blast