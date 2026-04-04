#pragma once
#include <blast>
#include <random>

namespace blast {
template<typename T>
host_fn bool is_close(const T type1, const T type2, real eps = 1e-5) {
  if (type1 != type2)
    return false;
  if (!is_close(type1, type2, eps))
    return false;

  return true;
}

template<typename T, std::size_t N>
host_fn bool is_close(const std::array<T, N>& a1, const std::array<T, N>& a2, real eps = 1e-5) {
  for (u32 i = 0; i < N; i++)
    if (!is_close(a1[i], a2[i], eps))
      return false;

  return true;
}

template<typename T>
host_fn bool is_close(const ObjMatrix<T>& a1, const ObjMatrix<T>& a2, real eps = 1e-5) {
  if (a1.size != a2.size) {
    Assert(false);
    return false;
  }
  return true;
}

template<typename T>
host_fn bool is_close(const std::vector<T>& a1, const std::vector<T>& a2, real eps = 1e-5) {
  if (a1.size() != a2.size()) {
    return false;
  }

  for (u32 i = 0; i < a1.size(); i++)
    if (!is_close(a1[i], a2[i], eps))
      return false;

  return true;
}

// todo: remove?
// template<typename T>
// host_fn bool is_close(const ObjDblArray<T>& a1, const ObjDblArray<T>& a2, real eps) {
//   if (a1.rows.size != a2.rows.size) {
//     Assert(false);
//     return false;
//   }
//   for (u32 i = 0; i < a1.cols; i++) {
//     if (a1.rows[i] != a2.rows[i]) {
//       Assert(false);
//       return false;
//     }
//     for (u32 j = 0; j < a2.rows[i]; j++) {
//       if (!is_close(a1(j, i), a2(j, i), eps))
//         return false;
//     }
//   }
//
//   return true;
// }

// note: Does not use eps but necessary for consistency for usability with templates
inline blast_fn bool is_close(const u8 a1, const u8& a2, real eps = 1e-5) {
  return (a1 == a2);
}

// note: Does not use eps but necessary for consistency for usability with templates
inline blast_fn bool is_close(u32 a1, u32 a2, real eps = 1e-5) {
  return (a1 == a2);
}

inline host_fn bool is_close(real r1, real r2, real eps = 1e-5) {
  return (r1 == INF_REAL && r2 == INF_REAL) || (r1 == -INF_REAL && r2 == -INF_REAL)
                 ? true
                 : (r1 == INF_REAL || r2 == INF_REAL || r1 == -INF_REAL || r2 == -INF_REAL
                            ? false
                            : (std::abs(r1 - r2) < eps));
}

inline host_fn bool is_close(const Box& box1, const Box& box2, real eps = 1e-5) {
  if (!is_close(box1.center, box2.center, eps))
    return false;
  if (!is_close(box1.extents, box2.extents, eps))
    return false;
  if (!is_close(box1.rotation, box2.rotation, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const DynamicBox& box1, const DynamicBox& box2, real eps = 1e-5) {
  if (!is_close(box1.n_points, box2.n_points))
    return false;
  if (!is_close(box1.start_time, box2.start_time, eps))
    return false;
  if (!is_close(box1.end_time, box2.end_time, eps))
    return false;
  if (!is_close(box1.trajectory, box2.trajectory, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const Sphere& sph1, const Sphere& sph2, real eps = 1e-5) {
  if (!is_close(sph1.center, sph2.center, eps))
    return false;
  if (!is_close(sph1.radius, sph2.radius, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const DynamicSphere& sph1, const DynamicSphere& sph2, real eps = 1e-5) {
  if (!is_close(sph1.n_points, sph2.n_points))
    return false;
  if (!is_close(sph1.start_time, sph2.start_time, eps))
    return false;
  if (!is_close(sph1.end_time, sph2.end_time, eps))
    return false;
  if (!is_close(sph1.trajectory, sph2.trajectory, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const Capsule& capsule1, const Capsule& capsule2, real eps = 1e-5) {
  if (!((is_close(capsule1.p1, capsule2.p1, eps) && is_close(capsule1.p2, capsule2.p2, eps)) ||
        (is_close(capsule1.p1, capsule2.p2, eps) && is_close(capsule1.p2, capsule2.p1, eps))))
    return false;
  if (!is_close(capsule1.radius, capsule2.radius, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const DynamicCapsule& capsule1, const DynamicCapsule& capsule2, real eps = 1e-5) {
  if (!is_close(capsule1.n_points, capsule2.n_points))
    return false;
  if (!is_close(capsule1.start_time, capsule2.start_time, eps))
    return false;
  if (!is_close(capsule1.end_time, capsule2.end_time, eps))
    return false;
  if (!is_close(capsule1.trajectory, capsule2.trajectory, eps))
    return false;

  return true;
}

inline host_fn bool is_close(const World& world1, const World& world2, real eps = 1e-5) {
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
  if (!is_close(world1.size, world2.size))
    return false;
  if (!is_close(world1.spheres, world2.spheres, eps))
    return false;

  return true;
}

// todo: operator==
inline host_fn bool is_close(const CollisionModelCapsule& capsule1, const CollisionModelCapsule& capsule2,
                             real eps = 1e-5) {
  if (!((is_close(capsule1.p1, capsule2.p1, eps) && is_close(capsule1.p2, capsule2.p2, eps)) ||
        (is_close(capsule1.p1, capsule2.p2, eps) && is_close(capsule1.p2, capsule2.p1, eps))))
    return false;
  if (!is_close(capsule1.radius, capsule2.radius, eps))
    return false;
  return true;
}

// todo: operator==
inline host_fn bool is_close(const Manipulator& manip1, const Manipulator& manip2, real eps = 1e-5) {
  if (manip1.n_joints != manip2.n_joints)
    return false;

  // manipulator limits
  for (int joint = 0; joint < manip1.n_joints; joint++) {
    if (!is_close(manip1.position_max[joint], manip2.position_max[joint], eps))
      return false; // max joint position
    if (!is_close(manip1.position_min[joint], manip2.position_min[joint], eps))
      return false; // min joint position
    if (!is_close(manip1.velocity_max[joint], manip2.velocity_max[joint], eps))
      return false; // max joint velocity
    if (!is_close(manip1.acceleration_max[joint], manip2.acceleration_max[joint], eps))
      return false; // max joint acceleration
    if (!is_close(manip1.torque_max[joint], manip2.torque_max[joint], eps))
      return false; // max joint torque
  }

  if (!is_close(manip1.tcp_speed_max, manip2.tcp_speed_max, eps))
    return false; // max tcp speed

  if (!is_close(manip1.base_position, manip2.base_position, eps))
    return false;
  if (!is_close(manip1.base_rotation, manip2.base_rotation, eps))
    return false;

  // kinematic properties
  for (int joint = 0; joint < manip1.n_joints; joint++) {
    if (!is_close(manip1.joint_offsets[joint], manip2.joint_offsets[joint], eps))
      return false; // vector to next joint
    if (!is_close(manip1.joint_axes[joint], manip2.joint_axes[joint], eps))
      return false; // direction vectors of joint

    // dynamic properties
    if (!is_close(manip1.link_masses[joint], manip2.link_masses[joint], eps))
      return false; // link mass
    if (!is_close(manip1.inertia_tensors[joint], manip2.inertia_tensors[joint], eps))
      return false; // inertial tensors
    if (!is_close(manip1.cog_offsets[joint], manip2.cog_offsets[joint], eps))
      return false; // center of mass
    if (!is_close(manip1.cog_from_next_joint[joint], manip2.cog_from_next_joint[joint], eps))
      return false; // center of mass from next joint
  }
  // Q_static is size n_joint + 1
  for (int joint = 0; joint < manip1.n_joints + 1; joint++) {
    if (!is_close(manip1.static_rotations[joint], manip2.static_rotations[joint]))
      return false; // static rotation matrix
  }

  // if (!is_close(manip1._capsule_list, manip2._capsule_list, eps))
  //   return false; // put the capsules temporarily when computing the constraints
  if (!is_close(manip1._collision_matrix, manip2._collision_matrix, eps))
    return false;
  if (!is_close(manip1._collision_base, manip2._collision_base, eps))
    return false;
  if (manip1._n_caps != manip2._n_caps)
    return false;
  if (manip1._n_internal_collisions != manip2._n_internal_collisions)
    return false;

  for (int capsule = 0; capsule < manip1._n_caps; capsule++) {
    if (!is_close(manip1._collision_model[capsule], manip2._collision_model[capsule], eps))
      return false;
  }
  if (!is_close(manip1._base_sphere, manip2._base_sphere, eps))
    return false;

  // manipulator limits
  if (!is_close(manip1.position_max, manip2.position_max, eps))
    return false; // max joint position
  if (!is_close(manip1.position_min, manip2.position_min, eps))
    return false; // min joint position
  if (!is_close(manip1.velocity_max, manip2.velocity_max, eps))
    return false; // max joint velocity
  if (!is_close(manip1.acceleration_max, manip2.acceleration_max, eps))
    return false; // max joint acceleration
  if (!is_close(manip1.torque_max, manip2.torque_max, eps))
    return false; // max joint torque


  return true;
}

inline host_fn bool is_close(const ManipulatorTempData& manip_data1, const ManipulatorTempData& manip_data2,
                             const u32 n_joints, const u32 n_caps, real eps = 1e-5) {
  for (int joint = 0; joint < n_joints; joint++) {
    // internal variables
    if (!is_close(manip_data1.efforts[joint], manip_data2.efforts[joint], eps))
      return false;

    if (!is_close(manip_data1.rotations[joint], manip_data2.rotations[joint], eps))
      return false; // put the rotation matrices temporarily when computing the constraints
    if (!is_close(manip_data1.rotations_mult[joint], manip_data2.rotations_mult[joint], eps))
      return false; // put the rotation matrices multiplications temporarily when computing the constraints
  }
  for (int joint = 0; joint < n_joints + 1; joint++) {
    if (!is_close(manip_data1.p_j[joint], manip_data2.p_j[joint], eps))
      return false; // put the joint coordinates temporarily when computing the constraints
  }

  for (int capsule = 0; capsule < n_caps; capsule++) {
    if (!is_close(manip_data1.capsule_list[capsule], manip_data2.capsule_list[capsule], eps))
      return false; // put the capsules temporarily when computing the constraints
  }

  return true;
}

// todo: operator==
// todo: why only test these parameters?
inline host_fn bool is_close(const Bspline& spline1, const Bspline& spline2, real eps = 1e-5) {
  if (!is_close(spline1.n_joints, spline2.n_joints, eps))
    return false;
  if (!is_close(spline1.degree, spline2.degree, eps))
    return false;
  if (!is_close(spline1.n_points, spline2.n_points, eps))
    return false;
  if (!is_close(spline1.n_ctrl, spline2.n_ctrl, eps))
    return false;

  return true;
}

// todo: operator==
inline host_fn bool is_close(const ConstraintSelection& constraints1, const ConstraintSelection& constraints2,
                             real eps = 1e-5) {
  if (constraints1.position != constraints2.position)
    return false;
  if (constraints1.velocity != constraints2.velocity)
    return false;
  if (constraints1.acceleration != constraints2.acceleration)
    return false;
  if (constraints1.torque != constraints2.torque)
    return false;
  if (constraints1.tcp_speed != constraints2.tcp_speed)
    return false;
  if (constraints1.self_collisions != constraints2.self_collisions)
    return false;
  if (constraints1.external_collisions != constraints2.external_collisions)
    return false;
  if (!is_close((u32) constraints1.n_collision_constraints, (u32) constraints2.n_collision_constraints, eps))
    return false;
  if (!is_close((u32) constraints1.n_collision_skip, (u32) constraints2.n_collision_skip, eps))
    return false;
  if (!is_close((u32) constraints1.n_constraints, (u32) constraints2.n_constraints, eps))
    return false;
  if (!is_close((u32) constraints1.extra_constraints.size(), (u32) constraints2.extra_constraints.size(), eps))
    return false;
  if (!is_close(constraints1.n_extra_constraints, constraints2.n_extra_constraints, eps))
    return false;

  return true;
}

// todo: operator==
inline host_fn bool is_close(const Objective& objective1, const Objective& objective2, real eps = 1e-5) {
  if (!is_close(objective1.time_weight, objective2.time_weight, eps))
    return false;
  if (!is_close(objective1.k_extra_objectives, objective2.k_extra_objectives, eps))
    return false;

  return true;
}

// todo: operator==
inline host_fn bool is_close(const Guess& guess1, const Guess& guess2, real eps = 1e-5) {
  if (guess1.type != guess2.type)
    return false;
  switch (guess1.type) {
    case Guess::custom:
      if (!is_close(guess1.initial_x, guess2.initial_x))
        return false;
      break;
    case Guess::random:
      if (is_close(guess1.n_random_shots, guess2.n_random_shots))
        return false;
      break;
    case Guess::rrt_connect:
      if (!is_close(guess1.parameter, guess2.parameter))
        return false;
      break;
    default:
      Assert(false);
  }

  return true;
}

// todo: operator==
inline host_fn bool is_close(const Optimization& opt1, const Optimization& opt2, real eps = 1e-5) {
  if (!is_close(opt1.bspline, opt2.bspline, eps))
    return false;
  if (!is_close(opt1.constraints, opt2.constraints, eps))
    return false;
  if (!is_close(opt1.objective, opt2.objective, eps))
    return false;
  if (!is_close(opt1.manip, opt2.manip, eps))
    return false;
  if (!is_close(opt1.task, opt2.task, eps))
    return false;
  if (!is_close(opt1.world, opt2.world, eps))
    return false;
  // if (!is_close(opt1.guess, opt2.guess, eps)) // todo: why does this not work?
  //     return false;

  return true;
}

// todo: operator==
inline host_fn bool is_close(const nlopt_result& result1, const nlopt_result& result2, real eps = 1e-5) {
  return (result1 == result2);
}

// todo: operator==
inline host_fn bool is_close(const Result& result1, Result& result2, real eps = 1e-5) {
  if (result1.success != result2.success)
    return false;
  if (result1.success_false != result2.success_false)
    return false;
  if (!is_close(*result1.opt, *result2.opt))
    return false;
  if (!is_close(result1.x, result2.x))
    return false;
  if (!is_close(result1.x0, result2.x0))
    return false;
  if (result1.nlopt_exit_criteria != result2.nlopt_exit_criteria)
    return false;

  return true;
}

// todo: make thread safe?
inline host_fn u32 random_int(u32 min, u32 max) {
  // Create a static random number generator (seeded once)
  static std::random_device rd;
  static std::mt19937       gen(rd());

  // Create a uniform distribution for integers in the given range
  std::uniform_int_distribution<u32> dist(min, max);

  return dist(gen);
}
} // namespace blast
