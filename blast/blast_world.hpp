#pragma once

#include <blast>

namespace blast {

struct World;

struct Box;
struct Capsule;
struct Sphere;
struct Cylinder;
struct DynamicBox;
struct DynamicSphere;
struct DynamicCapsule;
struct DynamicDoor;




// todo: Add cylinder
struct World {
    std::vector<Box>            boxes;
    std::vector<Sphere>         spheres;
    std::vector<Capsule>        capsules;
    std::vector<DynamicBox>     dynamic_boxes;
    std::vector<DynamicSphere>  dynamic_spheres;
    std::vector<DynamicCapsule> dynamic_capsules;
    std::vector<DynamicDoor>    dynamic_doors;
};

struct Box {
    Vec3 c; // Box center point
    Vec3 e; // Positive halfwidth extents of Box along each axis
    Mat3 R; // Local x-, y-, and z-axes (Rotation matrix)
};

struct Capsule {
    Vec3 p1;
    Vec3 p2;
    real r = 0;
};

struct Sphere {
    Vec3 c;
    real r = 0;
};


// todo: make fast!!
struct DynamicBox {
    u32 n_pts;
    real t0;
    real tf;
    std::vector<Box> trajectory;    // Should be of size n_pts

    inline blast_fn Box lookup(real time) const;
};

struct DynamicCapsule {
    u32 n_pts;
    real t0;
    real tf;
    std::vector<Capsule> trajectory; // Should be of size n_pts

    inline blast_fn Capsule lookup(real time) const;
};

struct DynamicSphere {
    u32 n_pts;
    real t0;
    real tf;
    std::vector<Sphere> trajectory; // Should be of size n_pts

    inline blast_fn Sphere lookup(real time) const;
};

struct DynamicDoor {
    Vec3 e;
    Vec3 hinge;
    Vec3 static_c_from_hinge;
    real start_angle;
    real end_angle;
    real t0;
    real tf;
    Vec3 axis = {0, 0, 1};

    inline blast_fn Box lookup(real t) const;
};


inline blast_fn void add_box(Vec3 center_point, Vec3 half_width, Mat3 rotation_matrix, World* world);
inline blast_fn void add_sphere(Vec3 center_point, real radius, World* world);
inline blast_fn void add_capsule(Vec3 point1, Vec3 point2, real radius, World* world);
inline blast_fn Array test_collision(std::vector<Capsule> robot_capsule_list, World* world, u32 n_lowest_distances);

inline blast_fn real distance(Capsule capsule,  Sphere sphere);
inline blast_fn real distance(Capsule capsule1, Capsule capsule2);
inline blast_fn real distance(Capsule capsule,  Box box);

inline blast_fn Vec3 get_point(const Array& x, const Matrix& capsule_list);
inline blast_fn real distance(const  Box& box, const Vec3& point);
inline blast_fn real distance(const  Capsule& capsule, const Vec3& point);
inline blast_fn real distance(const  Sphere& sphere, const  Vec3& point);


inline blast_fn Box DynamicDoor::lookup(const real t) const {
    // progression from 0 to 1
    real progression = t < t0 ? 0 : (t > tf ? 1 : (t-t0) / (tf-t0));
    Box result;
    result.e = e;
    real current_angle = start_angle*(1-progression) + end_angle*progression;
    result.R = rpy2rotation(current_angle*axis);
    result.c = hinge + result.R*static_c_from_hinge;
    return result;
}

} // namespace blast


#include "world/world.hpp"
#include "world/primitive.hpp"
#include "world/CoDO.hpp"



