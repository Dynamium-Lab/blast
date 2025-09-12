#pragma once

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
  r_tool[manip.n_joints - 1] = manip.dv[manip.n_joints - 1];
  for (int i = (int) manip.n_joints - 2; i >= 0; i--) {
    r_tool[i] = manip.dv[i] + temp.rotations[i + 1] * r_tool[i + 1];
  }

  for (int i = 0; i < manip.n_joints; i++) {
    r_tool[i] = temp.rotations_mult[i] * r_tool[i];
  }

  Matrix J_tool(6, manip.n_joints);
  for (int i = 0; i < manip.n_joints; i++) {
    const Vec3 e       = temp.rotations_mult[i] * manip.ev[i]; // replaced e directly in function to skip copy
    const Vec3 cr_tool = cross(e, r_tool[i]);
    // J_tool(0, i) = e.x;
    // J_tool(1, i) = e.y;
    // J_tool(2, i) = e.z;
    // J_tool(3, i) = cr_tool.x;
    // J_tool(4, i) = cr_tool.y;
    // J_tool(5, i) = cr_tool.z;
    J_tool(3, i) = e.x;
    J_tool(4, i) = e.y;
    J_tool(5, i) = e.z;
    J_tool(0, i) = cr_tool.x;
    J_tool(1, i) = cr_tool.y;
    J_tool(2, i) = cr_tool.z;
  }

  return J_tool;
}
//
// inline host_fn void forward_kinematics(Manipulator& manip, const Array& joint_pos) {
//   manip.compute_rotation_matrices(joint_pos);
//
//   manip._rotations_mult[0] = manip.Q_base * manip._rotations[0];
//   for (u32 j = 1; j < manip._rotations_mult.size(); j++) {
//     manip._rotations_mult[j] = manip._rotations_mult[j - 1] * manip._rotations[j];
//   }
//
//   manip._p_j[0] = manip.p_base + manip.Q_base * manip.p_j0;
//   for (u32 j = 1; j < manip._p_j.size(); j++) {
//     manip._p_j[j] = manip._p_j[j - 1] + manip._rotations_mult[j - 1] * manip.dv[j - 1];
//   }
// }

inline host_fn void forward_kinematics(const Manipulator& manip, ManipulatorTempData& temp, const Array& joint_pos) {
  const auto n_joints = joint_pos.size;

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

  for (u32 j = 0; j < n_joints; ++j) {
    auto       s      = sin(joint_pos[j]);
    auto       c      = cos(joint_pos[j]);
    const Mat3 temp_Q = {c, s, 0, -s, c, 0, 0, 0, 1};
    temp.rotations[j] = manip.Q_static[j] * temp_Q;
  }

  temp.rotations_mult[0] = manip.Q_base * temp.rotations[0];
  temp.p_j[0]            = manip.p_base + manip.Q_base * manip.p_j0;
  for (u32 j = 1; j < n_joints; j++) {
    temp.rotations_mult[j] = temp.rotations_mult[j - 1] * temp.rotations[j]; // note: add this to *= temp._rotations_mult ?
    temp.p_j[j]            = temp.p_j[j - 1] + temp.rotations_mult[j - 1] * manip.dv[j - 1];
  }
  temp.p_j[n_joints] = temp.p_j[n_joints - 1] + temp.rotations_mult[n_joints - 1] * manip.dv[n_joints - 1];

  // todo: Handle end-effector
}

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
  w[0] = vel[0] * manip.ev[0];
  for (u32 j = 1; j < joints; j++)
    w[j] = Qt[j] * w[j - 1] + vel[j] * manip.ev[j];

  Vec3 cdd0 = {0, 0, 9.81f};
  wd[0]     = acc[0] * manip.ev[0];
  cdd[0]    = Qt[0] * cdd0 + cross(wd[0], manip.av[0]) + cross(w[0], cross(w[0], manip.av[0]));
  for (u32 j = 1; j < joints; j++) {
    wd[j]  = Qt[j] * wd[j - 1] + acc[j] * manip.ev[j] + vel[j] * cross(Qt[j] * w[j - 1], manip.ev[j]);
    cdd[j] = Qt[j] * cdd[j - 1] + cross(wd[j], manip.av[j]) + cross(w[j], cross(w[j], manip.av[j])) - Qt[j] * cross(wd[j - 1], manip.sv[j - 1]) - Qt[j] * cross(w[j - 1], cross(w[j - 1], manip.sv[j - 1]));
  }

  //-- dynamics
  f[joints - 1] = manip.m[joints - 1] * cdd[joints - 1];
  n[joints - 1] = manip.I[joints - 1] * wd[joints - 1] + cross(w[joints - 1], manip.I[joints - 1] * w[joints - 1]) + cross(manip.av[joints - 1], f[joints - 1]);

  for (int j = (int) joints - 2; j >= 0; j--) {
    f[j] = manip.m[j] * cdd[j] + temp.rotations[j + 1] * f[j + 1];
    n[j] = manip.I[j] * wd[j] + cross(w[j], manip.I[j] * w[j]) + temp.rotations[j + 1] * n[j + 1] + cross(manip.av[j], f[j]) - cross(manip.sv[j], (temp.rotations[j + 1] * f[j + 1]));
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
  if (limits.pmax.size) {
    Assert(limits.pmax.size == n_joints);
    std::copy_n(limits.pmax.data, n_joints, pmax.data());
  }
  if (limits.pmin.size) {
    Assert(limits.pmin.size == n_joints);
    std::copy_n(limits.pmin.data, n_joints, pmin.data());
  }
  if (limits.vmax.size) {
    Assert(limits.vmax.size == n_joints);
    std::copy_n(limits.vmax.data, n_joints, vmax.data());
  }
  if (limits.amax.size != 0) {
    Assert(limits.amax.size == n_joints);
    std::copy_n(limits.amax.data, n_joints, amax.data());
  }
  if (limits.tau_max.size != 0) {
    Assert(limits.tau_max.size == n_joints);
    std::copy_n(limits.tau_max.data, n_joints, tau_max.data());
  }
  if (limits.tcp_max != 0.0) {
    tcp_max = limits.tcp_max;
  }
}

inline host_fn void Manipulator::set_kinematics(const ManipulatorKinematics& kinematics) {
  Assert(kinematics.dv.size() >= n_joints);
  Assert(kinematics.ev.size() >= n_joints);
  dv       = kinematics.dv;
  ev       = kinematics.ev;
  Q_static = kinematics.Q_static;
  p_base   = kinematics.p_base;
  Q_base   = kinematics.Q_base;
  p_j0     = kinematics.p_j0;
}

inline host_fn void Manipulator::set_dynamics(const ManipulatorDynamics& dynamics) {
  Assert(dynamics.m.size() >= n_joints);
  Assert(dynamics.I.size() >= n_joints);
  Assert(dynamics.av.size() >= n_joints);
  m  = dynamics.m;
  I  = dynamics.I;
  av = dynamics.av;
  for (u32 j = 0; j < n_joints; j++) {
    sv[j] = {-dv[j] + av[j]};
  }
}

inline host_fn void Manipulator::set_capsules(const ManipulatorCapsules& capsules) {
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
    _base_sphere.c += p_base;
    // _n_internal_collisions = 0;
  }
}

inline void compute_capsules(const Manipulator& manip, ManipulatorTempData& manip_data) {
  for (u32 i = 0; i < manip._n_caps; ++i) {
    manip_data.capsule_list[i] = {
            {manip_data.p_j[manip._collision_model[i].joint_frame] + manip_data.rotations_mult[manip._collision_model[i].joint_frame] * manip._collision_model[i].p1},
            {manip_data.p_j[manip._collision_model[i].joint_frame] + manip_data.rotations_mult[manip._collision_model[i].joint_frame] * manip._collision_model[i].p2},
            manip._collision_model[i].r};
  }
}


// todo: reformat for speed. ex: have a list of tuples that contain indices of collidable bodies
inline Array get_internal_collisions(const Manipulator& manip, const ManipulatorTempData& temp) {
  Array distances(manip._n_internal_collisions);
  int   idx = 0;
  for (int i = 0; i < manip._collision_matrix.cols; ++i) {
    for (int j = i + 1; j < manip._collision_matrix.rows; ++j) {
      if (manip._collision_matrix(j, i) != 0)
        distances[idx++] = distance(temp.capsule_list[i], temp.capsule_list[j]);
    }
  }
  for (u32 i = 0; i < manip._collision_base.size; ++i) {
    if (manip._collision_base[i] != 0)
      distances[idx++] = distance(temp.capsule_list[i], manip._base_sphere);
  }
  return distances;
}

// inline host_fn Array Manipulator::get_tool_pose() const {
//   Mat3 tmp             = _rotations_mult.back();
//   Vec3 current_rpy     = rotation2rpy(tmp);
//   Vec3 current_cartpos = _p_j.back();
//   return {current_cartpos.x, current_cartpos.y, current_cartpos.z,
//           current_rpy.x, current_rpy.y, current_rpy.z};
// }

inline host_fn void Manipulator::add_end_effector(const EndEffector& ee) {
  dv_ee        = ee.dv_ee;
  Q_ee         = ee.Q_ee;
  m_ee         = ee.m_ee;
  I_ee         = ee.I_ee;
  av_ee        = ee.av_ee;
  end_effector = true;

  dv.back() += Q_ee * ee.dv_ee;
  Q_static.back() *= Q_ee;
  real m_new    = m.back() + ee.m_ee;
  Vec3 av_new   = (m.back() * av.back() + ee.m_ee * ee.av_ee) / m_new;
  Vec3 delta_av = av_new - av.back();

  I.back()(0, 0) += m.back() * delta_av.x * delta_av.x + ee.m_ee * ee.av_ee.x * ee.av_ee.x;
  I.back()(1, 1) += m.back() * delta_av.y * delta_av.y + ee.m_ee * ee.av_ee.y * ee.av_ee.y;
  I.back()(2, 2) += m.back() * delta_av.z * delta_av.z + ee.m_ee * ee.av_ee.z * ee.av_ee.z;
  I.back() += ee.I_ee;

  m.back()  = m_new;
  av.back() = av_new;
  sv.back() = {-dv.back() + av.back()};
}

inline host_fn void Manipulator::set_payload(real m_payload, Vec3 cg_payload, Mat3 I_payload) {
  Vec3 av_payload = cg_payload;
  real m_new      = m.back() + m_payload;
  Vec3 av_new     = (m.back() * av.back() + m_payload * av_payload) / m_new;
  Vec3 delta_av   = av_new - av.back();
  Vec3 av_to_mass = av_payload - av_new;

  I.back()(0, 0) += m.back() * delta_av.x * delta_av.x + m_payload * av_to_mass.x * av_to_mass.x;
  I.back()(1, 1) += m.back() * delta_av.y * delta_av.y + m_payload * av_to_mass.y * av_to_mass.y;
  I.back()(2, 2) += m.back() * delta_av.z * delta_av.z + m_payload * av_to_mass.z * av_to_mass.z;
  I.back() += I_payload;

  m.back()  = m_new;
  av.back() = av_new;
  sv.back() = {-dv.back() + av.back()};
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
