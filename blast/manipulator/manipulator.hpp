#pragma once

#include <algorithm>
#include <blast>
#include <optional>
#include <utility>
// #include "nlopt.h"

namespace blast {

struct host_fn IK_opt {
  Manipulator manip;
  Array       desired_pose;

  IK_opt(Manipulator new_manip, Array new_pose) :
      manip(std::move(new_manip)),
      desired_pose(std::move(new_pose)) {}
};

inline host_fn Matrix jacobian(const Manipulator& manip, const ManipulatorTempData& temp) {
  std::vector<Vec3> r_tool(manip.n_joints);
  if (manip.has_tool) {
    r_tool[manip.n_joints - 1] = manip.tool.position + manip.tool.tool_center_position;
  } else {
    r_tool[manip.n_joints - 1] = Vec3(0, 0, 0);
  }

  for (int i = (int) manip.n_joints - 2; i >= 0; i--) {
    r_tool[i] = manip.joint_offsets[i] + temp.rotations[i + 1] * r_tool[i + 1];
  }

  for (int i = 0; i < manip.n_joints; i++) {
    r_tool[i] = temp.rotations_mult[i] * r_tool[i];
  }

  Matrix J_tool(6, manip.n_joints);
  for (int i = 0; i < manip.n_joints; i++) {
    const Vec3 e       = temp.rotations_mult[i] * manip.joint_axes[i]; // replaced e directly in function to skip copy
    const Vec3 cr_tool = cross(e, r_tool[i]);
    J_tool(3, i)       = e.x;
    J_tool(4, i)       = e.y;
    J_tool(5, i)       = e.z;
    J_tool(0, i)       = cr_tool.x;
    J_tool(1, i)       = cr_tool.y;
    J_tool(2, i)       = cr_tool.z;
  }

  return J_tool;
}

// todo: remove?
// inline host_fn void forward_kinematics(Manipulator& manip, const Array& joint_pos) {
//   manip.compute_rotation_matrices(joint_pos);
//
//   manip._rotations_mult[0] = manip.base_rotation * manip._rotations[0];
//   for (u32 j = 1; j < manip._rotations_mult.size(); j++) {
//     manip._rotations_mult[j] = manip._rotations_mult[j - 1] * manip._rotations[j];
//   }
//
//   manip._p_j[0] = manip.base_position + manip.base_rotation * manip.first_joint_position;
//   for (u32 j = 1; j < manip._p_j.size(); j++) {
//     manip._p_j[j] = manip._p_j[j - 1] + manip._rotations_mult[j - 1] * manip.joint_offsets[j - 1];
//   }
// }

inline host_fn void forward_kinematics(const Manipulator& manip, ManipulatorTempData& temp, const Array& joint_pos) {
  Assert(manip.n_joints == joint_pos.size);

  // real s[MAX_JOINTS];
  // real c[MAX_JOINTS];

  // #if 0
  // vcl::Vec4d    a, cr;
  // constexpr int vecLoopSize = MAX_JOINTS; // todo: fix for non divisible by 4
  // for (int i = 0; i < vecLoopSize; i += 4) {
  //   a.load(joint_pos.data + i);
  //   vcl::Vec4d sr = vcl::sincos(&cr, a);
  //   sr.store(s + i);
  //   cr.store(c + i);
  // }
  // #else
  //   for (int i = 0; i < n_joints; ++i) {
  //     s[i] = sin(joint_pos[i]);
  //     c[i] = cos(joint_pos[i]);
  //   }
  // #endif

  for (u32 j = 0; j < manip.n_joints; ++j) {
    auto       s      = sin(joint_pos[j]);
    auto       c      = cos(joint_pos[j]);
    const Mat3 temp_Q = {c, s, 0, -s, c, 0, 0, 0, 1};
    temp.rotations[j] = manip.static_rotations[j] * temp_Q;
  }

  temp.rotations_mult[0] = manip.base_rotation * temp.rotations[0];
  temp.p_j[0]            = manip.base_position + manip.base_rotation * manip.first_joint_position;
  for (u32 j = 1; j < manip.n_joints; j++) {
    temp.rotations_mult[j] = temp.rotations_mult[j - 1] * temp.rotations[j]; // note: add this to *= temp._rotations_mult ?
    temp.p_j[j]            = temp.p_j[j - 1] + temp.rotations_mult[j - 1] * manip.joint_offsets[j - 1];
  }
  // originally for the tool (todo: remove when agreed)
  temp.p_j[manip.n_joints] = temp.p_j[manip.n_joints - 1] + temp.rotations_mult[manip.n_joints - 1] * manip.joint_offsets[manip.n_joints - 1];

  if (manip.has_tool) {
    temp.tool_position = temp.p_j[manip.n_joints - 1] + temp.rotations_mult[manip.n_joints - 1] * manip.tool.position;
    temp.tool_rotation = temp.rotations_mult[manip.n_joints - 1] * manip.tool.rotation;
  }
  if (manip.has_payload) {
    temp.payload_position = temp.p_j[manip.n_joints - 1] + temp.rotations_mult[manip.n_joints - 1] * manip.payload.position;
    temp.payload_rotation = temp.rotations_mult[manip.n_joints - 1] * manip.payload.rotation;
  }
}

// todo: add tool and payload
inline host_fn void dynamics(const Manipulator& manip, ManipulatorTempData& temp, const Array& vel, const Array& acc) {
  Assert(vel.size == acc.size);
  Assert(vel.size == manip.n_joints);
  const auto joints = manip.n_joints;

  Vec3 w[MAX_JOINTS];
  Vec3 wd[MAX_JOINTS];
  Vec3 cdd[MAX_JOINTS];
  Mat3 Qt[MAX_JOINTS];
  Vec3 f[MAX_JOINTS];
  Vec3 n[MAX_JOINTS];

  for (u32 j = 0; j < joints; j++) {
    Qt[j] = transpose(temp.rotations[j]);
    // transpose(temp.rotations[j], Qt[j]);
  }

  // note: This is the Newton algorithm in 'Element de robotique' course notes.
  //       Careful because some variables are named differently and uses a slightly different conventions.
  //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
  //       joint turns with respect to the coordinate frame.
  //-- kinematics
  w[0] = vel[0] * manip.joint_axes[0];
  for (u32 j = 1; j < joints; j++)
    w[j] = Qt[j] * w[j - 1] + vel[j] * manip.joint_axes[j];

  Vec3 cdd0 = {0, 0, 9.81f};
  wd[0]     = acc[0] * manip.joint_axes[0];
  cdd[0]    = Qt[0] * cdd0 + cross(wd[0], manip.cog_offsets[0]) + cross(w[0], cross(w[0], manip.cog_offsets[0]));
  for (u32 j = 1; j < joints; j++) {
    wd[j]  = Qt[j] * wd[j - 1] + acc[j] * manip.joint_axes[j] + vel[j] * cross(Qt[j] * w[j - 1], manip.joint_axes[j]);
    cdd[j] = Qt[j] * cdd[j - 1] + cross(wd[j], manip.cog_offsets[j]) + cross(w[j], cross(w[j], manip.cog_offsets[j])) - Qt[j] * cross(wd[j - 1], manip.cog_from_next_joint[j - 1]) - Qt[j] * cross(w[j - 1], cross(w[j - 1], manip.cog_from_next_joint[j - 1]));
  }

  //-- dynamics
  if (manip.has_tool || manip.has_payload) {
    // include tool and payload in final link tensor
    real modified_last_link_mass    = manip.link_masses[manip.n_joints - 1];
    Vec3 modified_last_link_cog     = manip.cog_offsets[manip.n_joints - 1];
    Mat3 modified_last_link_inertia = manip.inertia_tensors[manip.n_joints - 1];

    if (manip.has_tool && manip.has_payload) {
      Vec3 tool_cog_offset_from_last_joint    = manip.tool.position + manip.tool.rotation * manip.tool.cog_offset;
      Vec3 payload_cog_offset_from_last_joint = manip.payload.position + manip.payload.rotation * manip.payload.cog_offset;

      real last_link_mass_with_tool    = manip.link_masses[manip.n_joints - 1];
      Vec3 last_link_cog_with_tool     = manip.cog_offsets[manip.n_joints - 1];
      Mat3 last_link_inertia_with_tool = manip.inertia_tensors[manip.n_joints - 1];

      sum_dynamic_properties(manip.link_masses[manip.n_joints - 1], manip.cog_offsets[manip.n_joints - 1], manip.inertia_tensors[manip.n_joints - 1],
                             manip.tool.mass, tool_cog_offset_from_last_joint, manip.tool.inertia_tensor,
                             last_link_mass_with_tool, last_link_cog_with_tool, last_link_inertia_with_tool);
      sum_dynamic_properties(last_link_mass_with_tool, last_link_cog_with_tool, last_link_inertia_with_tool,
                             manip.payload.mass, payload_cog_offset_from_last_joint, manip.payload.inertia_tensor,
                             modified_last_link_mass, modified_last_link_cog, modified_last_link_inertia);
    } else if (manip.has_tool) {
      Vec3 tool_cog_offset_from_last_joint = manip.tool.position + manip.tool.rotation * manip.tool.cog_offset;

      sum_dynamic_properties(manip.link_masses[manip.n_joints - 1], manip.cog_offsets[manip.n_joints - 1], manip.inertia_tensors[manip.n_joints - 1],
                             manip.tool.mass, tool_cog_offset_from_last_joint, manip.tool.inertia_tensor,
                             modified_last_link_mass, modified_last_link_cog, modified_last_link_inertia);
    } else if (manip.has_payload) { // this if could be a else {}
      Vec3 payload_cog_offset_from_last_joint = manip.payload.position + manip.payload.rotation * manip.payload.cog_offset;

      sum_dynamic_properties(manip.link_masses[manip.n_joints - 1], manip.cog_offsets[manip.n_joints - 1], manip.inertia_tensors[manip.n_joints - 1],
                             manip.payload.mass, payload_cog_offset_from_last_joint, manip.payload.inertia_tensor,
                             modified_last_link_mass, modified_last_link_cog, modified_last_link_inertia);
    }

    Vec3 modified_last_link_cdd = Qt[manip.n_joints - 1] * cdd[manip.n_joints - 1 - 1] + cross(wd[manip.n_joints - 1], modified_last_link_cog) + cross(w[manip.n_joints - 1], cross(w[manip.n_joints - 1], modified_last_link_cog));

    f[joints - 1] = modified_last_link_mass * modified_last_link_cdd;
    n[joints - 1] = modified_last_link_inertia * wd[joints - 1] + cross(w[joints - 1], modified_last_link_inertia * w[joints - 1]) + cross(modified_last_link_cog, f[joints - 1]);
  } else { // no gripper or payload
    f[joints - 1] = manip.link_masses[joints - 1] * cdd[joints - 1];
    n[joints - 1] = manip.inertia_tensors[joints - 1] * wd[joints - 1] + cross(w[joints - 1], manip.inertia_tensors[joints - 1] * w[joints - 1]) + cross(manip.cog_offsets[joints - 1], f[joints - 1]);
  }

  for (int j = (int) joints - 2; j >= 0; j--) {
    f[j] = manip.link_masses[j] * cdd[j] + temp.rotations[j + 1] * f[j + 1];
    n[j] = manip.inertia_tensors[j] * wd[j] + cross(w[j], manip.inertia_tensors[j] * w[j]) + temp.rotations[j + 1] * n[j + 1] + cross(manip.cog_offsets[j], f[j]) - cross(manip.cog_from_next_joint[j], (temp.rotations[j + 1] * f[j + 1]));
  }

  //-- extract torques (last element of each moment vector)
  for (u32 j = 0; j < joints; j++)
    temp.efforts[j] = n[j].z;
}

// todo: move to IK folder?
// inline host_fn double get_error(unsigned int n, const double* x, double* grad, void* data) {
//   Array delta_pose(12);
//   Array current_joint_position;
//   current_joint_position.alias(x, n);
//
//   auto info         = (IK_opt*) data;
//   auto manip        = info->manip;
//   auto desired_pose = info->desired_pose;
//
//   forward_kinematics(manip, current_joint_position);
//   auto current_pose = manip.get_tool_pose();
//
//   auto current_R = rpy2rotation({current_pose[3], current_pose[4], current_pose[5]});
//   auto desired_R = rpy2rotation({desired_pose[3], desired_pose[4], desired_pose[5]});
//
//   for (u32 i = 0; i < 3; i++)
//     delta_pose[i + 9] = current_pose[i] - desired_pose[i];
//
//   for (u32 i = 0; i < 9; i++) {
//     delta_pose[i] = current_R[i] - desired_R[i];
//   }
//
//   return dot(delta_pose, delta_pose);
// }

// inline host_fn Array inverse_kinematics_nlopt(Manipulator manip, Array desired_pose, Array initial_joint_position) {
//   IK_opt info(std::move(manip), std::move(desired_pose));
//
//   auto o      = nlopt_create(nlopt_algorithm::NLOPT_LN_COBYLA, initial_joint_position.size);
//   auto result = nlopt_set_min_objective(o, get_error, &info);
//   Assert(result == NLOPT_SUCCESS);
//   result = nlopt_set_ftol_abs(o, 0.0001);
//   Assert(result == NLOPT_SUCCESS);
//   // result = nlopt_set_xtol_abs(o, 0.0001);
//   Assert(result == NLOPT_SUCCESS);
//   result = nlopt_set_maxtime(o, 50);
//   Assert(result == NLOPT_SUCCESS);
//   result = nlopt_set_maxeval(o, 100000);
//   Assert(result == NLOPT_SUCCESS);
//
//   // launch optimization
//   // note: modifies our copy of initial_joint_position
//   double f;
//   nlopt_optimize(o, initial_joint_position.data, &f);
//
//   nlopt_destroy(o);
//
//   return initial_joint_position;
// }

inline host_fn Manipulator::Manipulator(u32 joint_count, const ManipulatorLimits& limits, const ManipulatorKinematics& kinematics, const ManipulatorDynamics* dynamics, const ManipulatorCapsules* capsules) {
  n_joints = joint_count;
  set_limits(limits);
  set_kinematics(kinematics);
  if (dynamics)
    set_dynamics(*dynamics);
  if (capsules)
    set_capsules(*capsules);
}

inline host_fn void Manipulator::set_limits(const ManipulatorLimits& limits) {
  if (limits.position_max.size) {
    Assert(limits.position_max.size == n_joints);
    std::copy_n(limits.position_max.data, n_joints, position_max.data());
  }
  if (limits.position_min.size) {
    Assert(limits.position_min.size == n_joints);
    std::copy_n(limits.position_min.data, n_joints, position_min.data());
  }
  if (limits.velocity_max.size) {
    Assert(limits.velocity_max.size == n_joints);
    std::copy_n(limits.velocity_max.data, n_joints, velocity_max.data());
  }
  if (limits.acceleration_max.size != 0) {
    Assert(limits.acceleration_max.size == n_joints);
    std::copy_n(limits.acceleration_max.data, n_joints, acceleration_max.data());
  }
  if (limits.jerk_max.size != 0) {
    Assert(limits.jerk_max.size == n_joints);
    std::copy_n(limits.jerk_max.data, n_joints, jerk_max.data());
  }
  if (limits.torque_max.size != 0) {
    Assert(limits.torque_max.size == n_joints);
    std::copy_n(limits.torque_max.data, n_joints, torque_max.data());
  }
  if (limits.tool_speed_max != 0.0) {
    tool_speed_max = limits.tool_speed_max;
  }
}

inline host_fn void Manipulator::set_kinematics(const ManipulatorKinematics& kinematics) {
  Assert(kinematics.joint_offsets.size() >= n_joints);
  Assert(kinematics.joint_axes.size() >= n_joints);
  joint_offsets        = kinematics.joint_offsets;
  joint_axes           = kinematics.joint_axes;
  static_rotations     = kinematics.static_rotations;
  base_position        = kinematics.base_position;
  base_rotation        = kinematics.base_rotation;
  first_joint_position = kinematics.first_joint_position;
}

inline host_fn void Manipulator::set_dynamics(const ManipulatorDynamics& dynamics) {
  Assert(dynamics.link_masses.size() >= n_joints);
  Assert(dynamics.inertia_tensors.size() >= n_joints);
  Assert(dynamics.cog_offsets.size() >= n_joints);
  link_masses     = dynamics.link_masses;
  inertia_tensors = dynamics.inertia_tensors;
  cog_offsets     = dynamics.cog_offsets;
  for (u32 j = 0; j < n_joints; j++) {
    cog_from_next_joint[j] = {-joint_offsets[j] + cog_offsets[j]};
  }
}

inline host_fn void Manipulator::set_capsules(const ManipulatorCapsules& capsules) {
  _n_caps = 0;
  for (auto& cap: capsules.capsule_list) {
    Assert(cap.joint_frame < n_joints + 1);
    _collision_model[_n_caps++] = cap;
  }

  if (capsules.collision_base.size != 0 && capsules.collision_matrix.size != 0) { // todo: fix this... (we need a fail safe)
    _collision_base.resize(_n_caps);
    _collision_base = capsules.collision_base;
    _collision_matrix.resize(_n_caps, _n_caps);
    _collision_matrix      = capsules.collision_matrix;
    _n_internal_collisions = 0;
    for (int i = 0; i < _n_caps; ++i) {
      for (int j = i + 1; j < _n_caps; ++j) {
        _n_internal_collisions += (_collision_matrix(j, i) != 0);
      }
      _n_internal_collisions += (_collision_base[i] != 0);
    }
    _base_sphere = capsules.base_sphere;
    _base_sphere.center += base_position;
  }
}

inline void compute_collision_model(const Manipulator& manip, ManipulatorTempData& manip_data) {
  for (u32 i = 0; i < manip._n_caps; ++i) {
    manip_data.capsule_list[i] = {
            {manip_data.p_j[manip._collision_model[i].joint_frame] + manip_data.rotations_mult[manip._collision_model[i].joint_frame] * manip._collision_model[i].p1},
            {manip_data.p_j[manip._collision_model[i].joint_frame] + manip_data.rotations_mult[manip._collision_model[i].joint_frame] * manip._collision_model[i].p2},
            manip._collision_model[i].radius};
  }

  if (manip.has_tool) {
    Mat3 collision_model_global_rotation = manip_data.rotations_mult[manip.n_joints - 1] * manip_data.tool_rotation;
    Vec3 collision_model_global_position = manip_data.tool_position + collision_model_global_rotation * manip.tool.collision_model.position;
    for (int i = 0; i < manip.tool.collision_model.points.size(); i++) {
      manip_data.tool_collision_model[i] = collision_model_global_position + collision_model_global_rotation * manip.tool.collision_model.points[i];
    }
  }
  if (manip.has_payload) {
    Mat3 collision_model_global_rotation = manip_data.rotations_mult[manip.n_joints - 1] * manip_data.payload_rotation;
    Vec3 collision_model_global_position = manip_data.payload_position + collision_model_global_rotation * manip.payload.collision_model.position;
    for (int i = 0; i < manip.payload.collision_model.points.size(); i++) {
      manip_data.payload_collision_model[i] = collision_model_global_position + collision_model_global_rotation * manip.payload.collision_model.points[i];
    }
  }
}

// todo: reformat for speed. ex: have a list of tuples that contain indices of collidable bodies
inline Array get_internal_collisions(const Manipulator& manip, const ManipulatorTempData& temp) {
  Array distances(manip._n_internal_collisions);
  int   idx = 0;
  for (int i = 0; i < manip._n_caps; ++i) {
    for (int j = i + 1; j < manip._n_caps; ++j) {
      if (manip._collision_matrix(j, i) != 0)
        distances[idx++] = distance(temp.capsule_list[i], temp.capsule_list[j]);
    }
  }
  for (u32 i = 0; i < manip._n_caps; ++i) {
    if (manip._collision_base[i] != 0)
      distances[idx++] = distance(temp.capsule_list[i], manip._base_sphere);
    // if (manip.has_tool && manip._collision_tool[i] != 0)
    //   distances[idx++] = distance(temp.capsule_list[i], temp.tool_collision_model);
    // if (manip.has_payload && manip._collision_payload[i] != 0)
    //   distances[idx++] = distance(temp.capsule_list[i], temp.payload_collision_model);
  }
  return distances;
}

// Attaches a new tool to the manipulator
// Erases the current tool if there is one
inline host_fn void Manipulator::set_tool(const Tool& new_tool) {
  has_tool = true;
  tool     = new_tool;
}

// Attaches a new payload to the manipulator
// Erases the current payload if there is one
inline host_fn void Manipulator::set_payload(const Payload& new_payload) {
  has_payload = true;
  payload     = new_payload;
}

// Remove the current tool from the manipulator
inline host_fn void Manipulator::remove_tool() {
  tool     = {};
  has_tool = false;
}

// Remove the current payload from the manipulator
inline host_fn void Manipulator::remove_payload() {
  has_payload = false;
  payload     = {};
}

inline host_fn real clamped_root(real slope, real h0, real h1) {
  // note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
  real r;
  if (h0 < 0) {
    if (h1 > 0) {
      r = -h0 / slope;
      if (r > 1)
        r = 0.5;
    } else
      r = 1;
  } else
    r = 0;
  return r;
}

inline host_fn void compute_intersection(const real sValue[2], const i32 classify[2], real b, real f00, real f10, i32* edge, real end[][2]) {
  // note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
  constexpr real zero = 0;
  constexpr real half = (real) 0.5;
  constexpr real one  = 1;
  if (classify[0] < 0) {
    edge[0]   = 0;
    end[0][0] = zero;
    end[0][1] = f00 / b;
    if (end[0][1] < zero || end[0][1] > one)
      end[0][1] = half;
    if (classify[1] == 0) {
      edge[1]   = 3;
      end[1][0] = sValue[1];
      end[1][1] = one;
    } else {
      edge[1]   = 1;
      end[1][0] = one;
      end[1][1] = f10 / b;
      if (end[1][1] < zero || end[1][1] > one)
        end[1][1] = half;
    }
  } else if (classify[0] == 0) {
    edge[0]   = 2;
    end[0][0] = sValue[0];
    end[0][1] = zero;
    if (classify[1] < 0) {
      edge[1]   = 0;
      end[1][0] = zero;
      end[1][1] = f00 / b;
      if (end[1][1] < zero || end[1][1] > one)
        end[1][1] = half;
    } else if (classify[1] == 0) {
      edge[1]   = 3;
      end[1][0] = sValue[1];
      end[1][1] = one;
    } else {
      edge[1]   = 1;
      end[1][0] = one;
      end[1][1] = f10 / b;
      if (end[1][1] < zero || end[1][1] > one)
        end[1][1] = half;
    }
  } else {
    edge[0]   = 1;
    end[0][0] = one;
    end[0][1] = f10 / b;
    if (end[0][1] < zero || end[0][1] > one)
      end[0][1] = half;
    if (classify[1] == 0) {
      edge[1]   = 3;
      end[1][0] = sValue[1];
      end[1][1] = one;
    } else {
      edge[1]   = 0;
      end[1][0] = zero;
      end[1][1] = f00 / b;
      if (end[1][1] < zero || end[1][1] > one)
        end[1][1] = half;
    }
  }
}

inline host_fn void compute_minimum_parameters(const i32 edge[2], const real end[][2], real b, real c, real e, real g00, real g10, real g01, real g11, real* parameter) {
  // note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
  constexpr real zero  = 0;
  constexpr real one   = 1;
  real const     delta = end[1][1] - end[0][1];
  if (real h0 = delta * (-b * end[0][0] + c * end[0][1] - e); h0 >= zero) {
    if (edge[0] == 0) {
      parameter[0] = zero;
      parameter[1] = clamped_root(c, g00, g01);
    } else if (edge[0] == 1) {
      parameter[0] = one;
      parameter[1] = clamped_root(c, g10, g11);
    } else {
      parameter[0] = end[0][0];
      parameter[1] = end[0][1];
    }
  } else {
    if (real h1 = delta * (-b * end[1][0] + c * end[1][1] - e); h1 <= zero) {
      if (edge[1] == 0) {
        parameter[0] = zero;
        parameter[1] = clamped_root(c, g00, g01);
      } else if (edge[1] == 1) {
        parameter[0] = one;
        parameter[1] = clamped_root(c, g10, g11);
      } else {
        parameter[0] = end[1][0];
        parameter[1] = end[1][1];
      }
    } else {
      real z       = clamp(h0 / (h0 - h1), 0, 1);
      real omz     = one - z;
      parameter[0] = omz * end[0][0] + z * end[1][0];
      parameter[1] = omz * end[0][1] + z * end[1][1];
    }
  }
}

inline host_fn real two_segment_distance_sqr(Vec3 P0, Vec3 P1, Vec3 Q0, Vec3 Q1) {
  // note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
  auto const P1mP0        = P1 - P0;
  auto const Q1mQ0        = Q1 - Q0;
  auto const P0mQ0        = P0 - Q0;
  real       a            = dot(P1mP0, P1mP0);
  real       b            = dot(P1mP0, Q1mQ0);
  real       c            = dot(Q1mQ0, Q1mQ0);
  real       d            = dot(P1mP0, P0mQ0);
  real       e            = dot(Q1mQ0, P0mQ0);
  real       f00          = d;
  real       f10          = f00 + a;
  real       f01          = f00 - b;
  real       f11          = f10 - b;
  real       g00          = -e;
  real       g10          = g00 - b;
  real       g01          = g00 + c;
  real       g11          = g10 + c;
  real       parameter[2] = {0, 0};
  if (a > 0 && c > 0) {
    real sValue[2] = {
            clamped_root(a, f00, f10),
            clamped_root(a, f01, f11)};
    i32 classify[2] = {0, 0};
    for (size_t i = 0; i < 2; ++i) {
      if (sValue[i] <= 0)
        classify[i] = -1;
      else if (sValue[i] >= 1)
        classify[i] = 1;
      else
        classify[i] = 0;
    }
    if (classify[0] == -1 && classify[1] == -1) {
      parameter[0] = 0;
      parameter[1] = clamped_root(c, g00, g01);
    } else if (classify[0] == +1 && classify[1] == +1) {
      parameter[0] = 1;
      parameter[1] = clamped_root(c, g10, g11);
    } else {
      i32  edge[2] = {0, 0};
      real end[2][2];
      compute_intersection(sValue, classify, b, f00, f10, edge, end);
      compute_minimum_parameters(edge, end, b, c, e, g00, g10, g01, g11, parameter);
    }
  } else {
    if (a > 0) {
      parameter[0] = clamped_root(a, f00, f10);
      parameter[1] = 0;
    } else if (c > 0) {
      parameter[0] = 0;
      parameter[1] = clamped_root(c, g00, g01);
    } else {
      parameter[0] = 0;
      parameter[1] = 0;
    }
  }
  Vec3 closest0 = P0 + parameter[0] * P1mP0;
  Vec3 closest1 = Q0 + parameter[1] * Q1mQ0;
  Vec3 diff     = closest0 - closest1;

  // auto result = sqrt(dot(diff, diff));
  return dot(diff, diff);
}


} // namespace blast
