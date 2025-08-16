#pragma once

#include <blast>

namespace blast {

inline host_fn void World::add_box(const Box &box) {
  boxes.push_back(box);
  size++;
}
inline host_fn void World::add_box(Vec3 center_point, Vec3 half_width, Mat3 rotation_matrix) {
  add_box({center_point, half_width, rotation_matrix});
}

inline host_fn void World::add_sphere(const Sphere &sphere) {
  spheres.push_back(sphere);
  size++;
}
inline host_fn void World::add_sphere(Vec3 center_point, real radius) {
  add_sphere({center_point, radius});
}

inline host_fn void World::add_capsule(const Capsule &capsule) {
  capsules.push_back(capsule);
  size++;
}
inline host_fn void World::add_capsule(Vec3 point1, Vec3 point2, real radius) {
  add_capsule({point1, point2, radius});
}

inline host_fn void World::add_dynamic_box(const DynamicBox &box) {
  dynamic_boxes.push_back(box);
  size++;
}
inline host_fn void World::add_dynamic_box(const std::vector<Box> &new_boxes, const u32 n_points, const real start_time, const real end_time) {
  dynamic_boxes.push_back({n_points, start_time, end_time, new_boxes});
  size++;
}

inline host_fn void World::add_dynamic_sphere(const DynamicSphere &sphere) {
  dynamic_spheres.push_back(sphere);
  size++;
}
inline host_fn void World::add_dynamic_sphere(const std::vector<Sphere> &new_spheres, u32 n_points, real start_time, real end_time) {
  add_dynamic_sphere({n_points, start_time, end_time, new_spheres});
}

inline host_fn void World::add_dynamic_capsule(const DynamicCapsule &capsule) {
  dynamic_capsules.push_back(capsule);
  size++;
}
inline host_fn void World::add_dynamic_capsule(const std::vector<Capsule> &new_capsules, const u32 n_points, const blast::real start_time, const blast::real end_time) {
  add_dynamic_capsule({n_points, start_time, end_time, new_capsules});
}

} // namespace blast
