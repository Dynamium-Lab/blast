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

struct World {
  std::vector<Box>            boxes;
  std::vector<Sphere>         spheres;
  std::vector<Capsule>        capsules;
  std::vector<DynamicBox>     dynamic_boxes;
  std::vector<DynamicSphere>  dynamic_spheres;
  std::vector<DynamicCapsule> dynamic_capsules;
  std::vector<DynamicDoor>    dynamic_doors;
  u32                         size = 0;

  host_fn void add_box(const Box& box);
  host_fn void add_box(Vec3 center_point, Vec3 half_width, Mat3 rotation_matrix);

  host_fn void add_sphere(const Sphere& sphere);
  host_fn void add_sphere(Vec3 center_point, real radius);

  host_fn void add_capsule(const Capsule& capsule);
  host_fn void add_capsule(Vec3 point1, Vec3 point2, real radius);

  host_fn void add_dynamic_box(const DynamicBox& box);
  host_fn void add_dynamic_box(const std::vector<Box>& new_boxes, u32 n_points, real start_time, real end_time);

  host_fn void add_dynamic_sphere(const DynamicSphere& sphere);
  host_fn void add_dynamic_sphere(const std::vector<Sphere>& new_spheres, u32 n_points, real start_time, real end_time);

  host_fn void add_dynamic_capsule(const DynamicCapsule& capsule);
  host_fn void add_dynamic_capsule(const std::vector<Capsule>& new_capsules, u32 n_points, real start_time, real end_time);
};

struct Box {
  Vec3 center;   // Box center point
  Vec3 extents;  // Positive halfwidth extents of Box along each axis
  Mat3 rotation; // Local x-, y-, and z-axes (Rotation matrix)
};

struct Capsule {
  Vec3 p1;
  Vec3 p2;
  real radius = 0;
};

struct Sphere {
  Vec3 center;
  real radius = 0;
};

struct DynamicBox {
  u32              n_points;
  real             start_time = 0;
  real             end_time   = 0;
  std::vector<Box> trajectory; // Should be of size n_points

  inline blast_fn Box lookup(real time) const;
};

struct DynamicCapsule {
  u32                  n_points;
  real                 start_time = 0;
  real                 end_time   = 0;
  std::vector<Capsule> trajectory; // Should be of size n_points

  inline blast_fn Capsule lookup(real time) const;
};

struct DynamicSphere {
  u32                 n_points;
  real                start_time = 0;
  real                end_time   = 0;
  std::vector<Sphere> trajectory; // Should be of size n_points

  inline blast_fn Sphere lookup(real time) const;
};

struct DynamicDoor {
  Vec3 extents;
  Vec3 hinge;
  Vec3 static_c_from_hinge;
  real start_angle = 0;
  real end_angle   = 0;
  real start_time  = 0;
  real end_time    = 0;
  Vec3 axis        = {0, 0, 1};

  inline blast_fn Box lookup(real t) const;
};

inline blast_fn Array test_collisions(const ObjMatrix<Capsule>& robot_capsules, const World* world, u32 n_lowest_distances, real start_time, real end_time);

inline blast_fn real distance(const Capsule& capsule, const Box& box);
inline blast_fn real distance(const Capsule& capsule, const Sphere& sphere);
inline blast_fn real distance(const Capsule& capsule1, const Capsule& capsule2);
inline blast_fn real distance(const Capsule& capsule, const Vec3& point);
inline blast_fn real distance(const Box& box, const Vec3& point);
inline blast_fn real distance(const Sphere& sphere, const Vec3& point);

inline blast_fn Vec3 get_point(const Array& x, const Matrix& capsule_list);


} // namespace blast


#include "world/distance.hpp"
#include "world/world.hpp"

#include "world/dynamicbox.hpp"
#include "world/dynamiccapsule.hpp"
#include "world/dynamicdoor.hpp"
#include "world/dynamicsphere.hpp"

#include "world/CoDO.hpp"
