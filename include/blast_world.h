#pragma once

namespace blast {

struct Box;
struct Capsule;
struct Sphere;
struct Cylinder;

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
    std::vector<Capsule> caps;
};


} // namespace blast