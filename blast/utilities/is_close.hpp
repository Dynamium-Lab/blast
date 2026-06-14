#pragma once
#include <blast>
#include <type_traits>

namespace blast {

template<typename T, std::size_t N>
host_fn bool is_close(const std::array<T, N>& a1, const std::array<T, N>& a2, real eps) {
  for (u32 i = 0; i < N; i++)
    if (!is_close(a1[i], a2[i], eps))
      return false;

  return true;
}

template<typename T>
host_fn bool is_close(const ObjMatrix<T>& a1, const ObjMatrix<T>& a2, real eps) {
  if (a1.size != a2.size) {
    return false;
  }
  for (u32 i = 0; i < a1.size; i++) {
    // Integral element types (e.g. ObjMatrix<u8>) compare exactly; eps is meaningless for them.
    if constexpr (std::is_integral_v<T>) {
      if (a1.data[i] != a2.data[i])
        return false;
    } else {
      if (!is_close(a1.data[i], a2.data[i], eps))
        return false;
    }
  }

  return true;
}

template<typename T>
host_fn bool is_close(const std::vector<T>& a1, const std::vector<T>& a2, real eps) {
  if (a1.size() != a2.size())
    return false;

  for (u32 i = 0; i < a1.size(); i++)
    if (!is_close(a1[i], a2[i], eps))
      return false;

  return true;
}

inline host_fn bool is_close(real r1, real r2, real eps) {
  return (r1 == INF_REAL && r2 == INF_REAL) || (r1 == -INF_REAL && r2 == -INF_REAL)
                 ? true
                 : (r1 == INF_REAL || r2 == INF_REAL || r1 == -INF_REAL || r2 == -INF_REAL
                            ? false
                            : (std::abs(r1 - r2) < eps));
}

inline host_fn bool is_close(const Box& box1, const Box& box2, real eps) {
  if (!is_close(box1.center, box2.center, eps))
    return false;
  if (!is_close(box1.extents, box2.extents, eps))
    return false;
  if (!is_close(box1.rotation, box2.rotation, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const DynamicBox& box1, const DynamicBox& box2, real eps) {
  if (box1.n_points != box2.n_points)
    return false;
  if (!is_close(box1.start_time, box2.start_time, eps))
    return false;
  if (!is_close(box1.end_time, box2.end_time, eps))
    return false;
  if (!is_close(box1.trajectory, box2.trajectory, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const Sphere& sph1, const Sphere& sph2, real eps) {
  if (!is_close(sph1.center, sph2.center, eps))
    return false;
  if (!is_close(sph1.radius, sph2.radius, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const DynamicSphere& sph1, const DynamicSphere& sph2, real eps) {
  if (sph1.n_points != sph2.n_points)
    return false;
  if (!is_close(sph1.start_time, sph2.start_time, eps))
    return false;
  if (!is_close(sph1.end_time, sph2.end_time, eps))
    return false;
  if (!is_close(sph1.trajectory, sph2.trajectory, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const Capsule& capsule1, const Capsule& capsule2, real eps) {
  if (!((is_close(capsule1.p1, capsule2.p1, eps) && is_close(capsule1.p2, capsule2.p2, eps)) ||
        (is_close(capsule1.p1, capsule2.p2, eps) && is_close(capsule1.p2, capsule2.p1, eps))))
    return false;
  if (!is_close(capsule1.radius, capsule2.radius, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const DynamicCapsule& capsule1, const DynamicCapsule& capsule2, real eps) {
  if (capsule1.n_points != capsule2.n_points)
    return false;
  if (!is_close(capsule1.start_time, capsule2.start_time, eps))
    return false;
  if (!is_close(capsule1.end_time, capsule2.end_time, eps))
    return false;
  if (!is_close(capsule1.trajectory, capsule2.trajectory, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const World& world1, const World& world2, real eps) {
  if (world1.size != world2.size)
    return false;
  if (!is_close(world1.boxes, world2.boxes, eps))
    return false;
  if (!is_close(world1.capsules, world2.capsules, eps))
    return false;
  if (!is_close(world1.dynamic_boxes, world2.dynamic_boxes, eps))
    return false;
  if (!is_close(world1.dynamic_capsules, world2.dynamic_capsules, eps))
    return false;
  if (!is_close(world1.dynamic_spheres, world2.dynamic_spheres, eps))
    return false;
  if (!is_close(world1.spheres, world2.spheres, eps))
    return false;

  return true;
}

inline host_fn bool operator==(const CollisionModelCapsule& a, const CollisionModelCapsule& b) {
  return ((a.p1 == b.p1 && a.p2 == b.p2) || (a.p1 == b.p2 && a.p2 == b.p1)) && is_close(a.radius, b.radius);
}

inline host_fn bool operator==(const Manipulator& a, const Manipulator& b) {
  if (a.n_joints != b.n_joints)
    return false;

  if (!is_close(a.tool_speed_max, b.tool_speed_max))
    return false;

  if (!(a.base_position == b.base_position))
    return false;
  if (!(a.base_rotation == b.base_rotation))
    return false;

  for (int j = 0; j < a.n_joints; j++) {
    if (!(a.joint_offsets[j] == b.joint_offsets[j]))
      return false;
    if (!(a.joint_axes[j] == b.joint_axes[j]))
      return false;
    if (!is_close(a.link_masses[j], b.link_masses[j]))
      return false;
    if (!(a.inertia_tensors[j] == b.inertia_tensors[j]))
      return false;
    if (!(a.cog_offsets[j] == b.cog_offsets[j]))
      return false;
    if (!(a.cog_from_next_joint[j] == b.cog_from_next_joint[j]))
      return false;
  }

  for (int j = 0; j < a.n_joints + 1; j++) {
    if (!(a.static_rotations[j] == b.static_rotations[j]))
      return false;
  }

  if (!(a._collision_matrix == b._collision_matrix))
    return false;
  if (!(a._collision_base == b._collision_base))
    return false;
  if (a._n_caps != b._n_caps)
    return false;
  if (a._n_internal_collisions != b._n_internal_collisions)
    return false;

  for (int c = 0; c < a._n_caps; c++) {
    if (!(a._collision_model[c] == b._collision_model[c]))
      return false;
  }
  if (!is_close(a._base_sphere, b._base_sphere))
    return false;

  return a.position_max == b.position_max && a.position_min == b.position_min && a.velocity_max == b.velocity_max && a.acceleration_max == b.acceleration_max && a.torque_max == b.torque_max;
}

inline host_fn bool is_close(const ManipulatorTempData& manip_data1, const ManipulatorTempData& manip_data2,
                             const u32 n_joints, const u32 n_caps, real eps) {
  for (int joint = 0; joint < n_joints; joint++) {
    if (!is_close(manip_data1.efforts[joint], manip_data2.efforts[joint], eps))
      return false;
    if (!is_close(manip_data1.rotations[joint], manip_data2.rotations[joint], eps))
      return false;
    if (!is_close(manip_data1.rotations_mult[joint], manip_data2.rotations_mult[joint], eps))
      return false;
  }
  for (int joint = 0; joint < n_joints + 1; joint++) {
    if (!is_close(manip_data1.p_j[joint], manip_data2.p_j[joint], eps))
      return false;
  }
  for (int capsule = 0; capsule < n_caps; capsule++) {
    if (!is_close(manip_data1.capsule_list[capsule], manip_data2.capsule_list[capsule], eps))
      return false;
  }

  return true;
}

// todo: why only test these parameters?
inline host_fn bool operator==(const Bspline& a, const Bspline& b) {
  return a.n_joints == b.n_joints && a.degree == b.degree && a.n_points == b.n_points && a.n_ctrl == b.n_ctrl;
}

inline host_fn bool operator==(const ConstraintSelection& a, const ConstraintSelection& b) {
  return a.position == b.position && a.velocity == b.velocity && a.acceleration == b.acceleration && a.torque == b.torque && a.tool_speed == b.tool_speed && a.self_collisions == b.self_collisions && a.external_collisions == b.external_collisions && a.n_collision_constraints == b.n_collision_constraints && a.n_collision_skip == b.n_collision_skip && a.n_constraints == b.n_constraints && a.extra_constraints.size() == b.extra_constraints.size() && a.n_extra_constraints == b.n_extra_constraints;
}

inline host_fn bool operator==(const Objective& a, const Objective& b) {
  return is_close(a.time_weight, b.time_weight) && is_close(a.k_extra_objectives, b.k_extra_objectives);
}

inline host_fn bool operator==(const Guess& a, const Guess& b) {
  if (a.type != b.type)
    return false;
  switch (a.type) {
    case Guess::custom:
      return a.initial_x == b.initial_x;
    case Guess::random:
      return a.n_random_shots == b.n_random_shots;
    default:
      Assert(false);
      return false;
  }
}

inline host_fn bool operator==(const Optimization& a, const Optimization& b) {
  return a.bspline == b.bspline && a.constraints == b.constraints && a.objective == b.objective && a.manip == b.manip && a.task == b.task && is_close(a.world, b.world);
}

inline host_fn bool operator==(const Result& a, const Result& b) {
  return a.success == b.success && a.success_false == b.success_false && *a.opt == *b.opt && a.x == b.x && a.x0 == b.x0 && a.nlopt_exit_criteria == b.nlopt_exit_criteria;
}


} // namespace blast
