#pragma once

#include <blast>

namespace blast {

inline host_fn void add_box(const Box &box, World *world) {
  world->boxes.push_back(box);
  world->size++;
}

inline host_fn void add_sphere(const Sphere &sphere, World *world) {
  world->spheres.push_back(sphere);
  world->size++;
}

inline host_fn void add_capsule(const Capsule &capsules, World *world) {
  world->capsules.push_back(capsules);
  world->size++;
}

inline host_fn void add_dynamic_box(const DynamicBox &box, World *world) {
  world->dynamic_boxes.push_back(box);
  world->size++;
}

inline host_fn void add_dynamic_sphere(const DynamicSphere &sphere, World *world) {
  world->dynamic_spheres.push_back(sphere);
  world->size++;
}

inline host_fn void add_dynamic_capsule(const DynamicCapsule &capsules, World *world) {
  world->dynamic_capsules.push_back(capsules);
  world->size++;
}

inline host_fn void add_box(Vec3 center_point, Vec3 half_width, Mat3 rotation_matrix, World *world) {
  Box new_box;
  new_box.c = center_point;
  new_box.e = half_width;
  new_box.R = rotation_matrix;
  add_box(new_box, world);
}

inline host_fn void add_capsule(Vec3 point1, Vec3 point2, real radius, World *world) {
  Capsule new_capsule;
  new_capsule.p1 = point1;
  new_capsule.p2 = point2;
  new_capsule.r = radius;
  add_capsule(new_capsule, world);
}

inline host_fn void add_sphere(Vec3 center_point, real radius, World *world) {
  Sphere new_sphere;
  new_sphere.c = center_point;
  new_sphere.r = radius;
  add_sphere(new_sphere, world);
}

inline host_fn void add_dynamic_box(const std::vector<blast::Box> &boxes, const u32 n_points, const blast::real start_time, const blast::real end_time, World *world) {
  DynamicBox new_dynamic_box;
  new_dynamic_box.n_pts = n_points;
  new_dynamic_box.trajectory = boxes;
  new_dynamic_box.t0 = start_time;
  new_dynamic_box.tf = end_time;
  add_dynamic_box(new_dynamic_box, world);
}

inline host_fn void add_dynamic_capsule(const std::vector<blast::Capsule> &capsules, const u32 n_points, const blast::real start_time, const blast::real end_time, World *world) {
  DynamicCapsule new_dynamic_capsule;
  new_dynamic_capsule.n_pts = n_points;
  new_dynamic_capsule.trajectory = capsules;
  new_dynamic_capsule.t0 = start_time;
  new_dynamic_capsule.tf = end_time;
  add_dynamic_capsule(new_dynamic_capsule, world);
}

inline host_fn void add_dynamic_sphere(const std::vector<blast::Sphere> &spheres, const u32 n_points, const blast::real start_time, const blast::real end_time, World *world) {
  DynamicSphere new_dynamic_sphere;
  new_dynamic_sphere.n_pts = n_points;
  new_dynamic_sphere.trajectory = spheres;
  new_dynamic_sphere.t0 = start_time;
  new_dynamic_sphere.tf = end_time;
  add_dynamic_sphere(new_dynamic_sphere, world);
}

} // namespace blast
