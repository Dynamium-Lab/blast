#pragma once

#include <blast>
#include <optional>

namespace blast {

// Defined Types
struct CollisionModelCapsule;
struct ManipulatorLimits;
struct ManipulatorKinematics;
struct ManipulatorDynamics;
struct ManipulatorCapsules;
struct EndEffector;
struct Manipulator;



/**
 * @struct CollisionModelCapsule
 * @brief Simple capsule primitive for collision checking.
 *
 * @var p1           First endpoint of the capsule (Vec3).
 * @var joint_frame  Index of the joint frame to which this capsule is attached.
 * @var p2           Second endpoint of the capsule (Vec3).
 * @var r            Radius of the capsule.
 */
struct CollisionModelCapsule {
  Vec3 p1;
  u32  joint_frame;
  Vec3 p2;
  real r;
};

/**
 * @struct ManipulatorLimits
 * @brief Actuator limits for a manipulator's joints and TCP.
 *
 * @var pmax    Maximum joint positions.
 * @var pmin    Minimum joint positions.
 * @var vmax    Maximum joint velocities.
 * @var vmin    Minimum joint velocities.
 * @var amax    Maximum joint accelerations.
 * @var amin    Minimum joint accelerations.
 * @var tau_max Maximum joint torques.
 * @var tau_min Minimum joint torques.
 * @var tcp_max Maximum tool‐center‐point speed.
 */
struct ManipulatorLimits {
  Array pmax;    // max joint position
  Array pmin;    // min joint position
  Array vmax;    // max joint velocity
  Array vmin;    // min joint velocity
  Array amax;    // max joint acceleration
  Array amin;    // min joint acceleration
  Array tau_max; // max joint torque
  Array tau_min; // min joint torque

  real tcp_max;  // max tcp speed
};

/**
 * @struct ManipulatorKinematics
 * @brief Geometric configuration of manipulator joints and base.
 *
 * @var dv       Vectors from each joint to the next (in joint frames).
 * @var ev       Unit axis directions for each joint.
 * @var Q_static Static rotation matrices from joint to next (Mat3).
 * @var p_j0     Position of the first joint relative to the base.
 * @var p_base   Position of the manipulator base in the world frame.
 * @var Q_base   Orientation of the manipulator base in the world frame.
 */
struct ManipulatorKinematics {
  std::vector<Vec3> dv;                                   // vector to next joint
  std::vector<Vec3> ev;                                   // direction vectors of joint
  std::vector<Mat3> Q_static;                             // static rotation to next joint
  Vec3              p_j0   = {0.0, 0.0, 0.0};             // position of the first joint relative to base
  Vec3              p_base = {0.0, 0.0, 0.0};             // base position in workspace
  Mat3              Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // base orientation in workspace
};

/**
 * @struct ManipulatorDynamics
 * @brief Dynamic properties of each manipulator link.
 *
 * @var m    Mass of each link.
 * @var I    Inertia tensor of each link (Mat3).
 * @var av   Center of mass offset for each link (Vec3).
 */
struct ManipulatorDynamics {
  std::vector<real> m;  // link mass
  std::vector<Mat3> I;  // inertial tensors
  std::vector<Vec3> av; // center of mass
};

/**
 * @struct ManipulatorCapsules
 * @brief Collision geometry grouping for manipulator.
 *
 * @var base_sphere      Collision sphere around the base.
 * @var capsule_list     List of link capsules for collision.
 * @var collision_matrix Matrix indicating collisions between capsules.
 * @var collision_base   Array indicating collisions with base sphere.
 */
struct ManipulatorCapsules {
  Sphere                             base_sphere;
  std::vector<CollisionModelCapsule> capsule_list;
  ObjMatrix<u8>                      collision_matrix = {}; // which capsules collide
  Array                              collision_base;        // collisions with base sphere
};

/**
 * @struct EndEffector
 * @brief Payload and collision properties of an end effector.
 *
 * @var dv_ee         Offset vector from TCP to payload.
 * @var Q_ee          Rotation from TCP to payload.
 * @var m_ee          Payload mass.
 * @var I_ee          Payload inertia tensor.
 * @var av_ee         Payload center of mass offset.
 * @var capsule_list  Collision capsules for the end effector.
 * @var collision_matrix Collision matrix for payload capsules.
 */
struct EndEffector {
  Vec3 dv_ee = {0.0, 0.0, 0.0};
  Mat3 Q_ee  = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  real m_ee  = 0.0;
  Mat3 I_ee;
  Vec3 av_ee = {0.0, 0.0, 0.0};

  // Collision for end effector
  std::vector<CollisionModelCapsule> capsule_list;
  ObjMatrix<u8>                      collision_matrix = {};
};

/**
 * @class Manipulator
 * @brief High‐level interface combining limits, kinematics, dynamics,
 *        and collision geometry for a manipulator.
 *
 * @note Use the constructor to fully specify at least limits and kinematics.
 */
struct Manipulator {
  u32 joints = 0;

  // Manipulator limits
  Array pmax;
  Array pmin;
  Array vmax;
  Array vmin;
  Array amax;
  Array amin;
  Array tau_max;
  Array tau_min;

  real tcp_max = 0.0; // max TCP speed

  // World‐frame base pose
  Vec3 p_base = {0, 0, 0};
  Mat3 Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  // End effector state
  bool end_effector = false;
  Vec3 dv_ee        = {0, 0, 0};
  Mat3 Q_ee         = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  real m_ee         = 0.0;
  Mat3 I_ee         = eye();
  Vec3 av_ee        = {0, 0, 0};

  // Joint kinematics
  Vec3              p_j0 = {0, 0, 0};
  std::vector<Vec3> dv;       // to next joint
  std::vector<Vec3> ev;       // axis dirs
  std::vector<Mat3> Q_static; // static rotations

  // Link dynamics
  std::vector<real> m;
  std::vector<Mat3> I;
  std::vector<Vec3> av;
  std::vector<Vec3> sv; // COM from next joint

  // Internal caches and collision
  Array                              _efforts;
  std::vector<Mat3>                  _Q;
  std::vector<Mat3>                  _Q_mult;
  std::vector<Vec3>                  _p_j;
  std::vector<CollisionModelCapsule> _collision_model;
  Sphere                             _base_sphere;
  std::vector<Capsule>               _capsule_list;
  ObjMatrix<u8>                      _collision_matrix;
  Array                              _collision_base;
  u32                                _n_caps                = 0;
  u32                                _n_internal_collisions = 0;

  Manipulator() = delete;

  /**
   * @brief Construct a GenericManipulator with required data.
   * @param joint_count  Number of joints (excluding fixed).
   * @param limits       Mandatory ManipulatorLimits.
   * @param kinematics   Mandatory ManipulatorKinematics.
   * @param dynamics     Optional pointer to ManipulatorDynamics.
   * @param capsules     Optional pointer to ManipulatorCapsules.
   */
  host_fn Manipulator(u32 joint_count, const ManipulatorLimits& limits, const ManipulatorKinematics& kinematics, const ManipulatorDynamics* dynamics = nullptr, const ManipulatorCapsules* capsules = nullptr);

  /**
   * @brief Apply provided limits to this manipulator.
   * @param limits  Reference to ManipulatorLimits.
   */
  host_fn void set_limits(const ManipulatorLimits& limits);

  /**
   * @brief Apply provided kinematics to this manipulator.
   * @param kinematics  Reference to ManipulatorKinematics.
   */
  host_fn void set_kinematics(const ManipulatorKinematics& kinematics);

  /**
   * @brief Apply provided dynamics to this manipulator.
   * @param dynamics  Reference to ManipulatorDynamics.
   */
  host_fn void set_dynamics(const ManipulatorDynamics& dynamics);

  /**
   * @brief Apply provided collision capsules to this manipulator.
   * @param capsules  Reference to ManipulatorCapsules.
   */
  host_fn void set_capsules(const ManipulatorCapsules& capsules);

  /**
   * @brief Attach an end effector to this manipulator.
   * @param ee  EndEffector structure with payload and geometry.
   */
  host_fn void add_end_effector(const EndEffector& ee);

  /**
   * @brief Add a payload to the last link, adjusting mass and inertia.
   * @param m_payload    Payload mass.
   * @param cg_payload   Payload center of mass offset.
   * @param I_payload    Payload inertia tensor.
   */
  host_fn void set_payload(real m_payload, Vec3 cg_payload, Mat3 I_payload);

  /**
   * @brief Compute joint rotation matrices for current joint positions (store internally).
   * @param joint_position  Array of joint angles.
   */
  host_fn void compute_rotation_matrices(const Array& joint_position);

  /**
   * @brief Update collision capsules in world frame (store internally).
   */
  host_fn void compute_capsules();

  /**
   * @brief Compute distances for all internal collisions.
   * @return Array of distances between colliding primitives.
   */
  host_fn Array get_internal_collisions() const;

  /**
   * @brief Get the current tool pose as [x,y,z,roll,pitch,yaw].
   * @return Array of six values: position and RPY orientation.
   */
  host_fn Array get_tool_pose() const;
};


inline host_fn Manipulator::Manipulator(u32 joint_count, const ManipulatorLimits& limits, const ManipulatorKinematics& kinematics, const ManipulatorDynamics* dynamics, const ManipulatorCapsules* capsules) {
  joints = joint_count;
  set_limits(limits);
  set_kinematics(kinematics);
  if (dynamics)
    set_dynamics(*dynamics);
  if (capsules)
    set_capsules(*capsules);
}

inline host_fn void Manipulator::set_limits(const ManipulatorLimits& limits) {
  if (limits.pmax.size != 0) {
    Assert(limits.pmax.size == joints);
    pmax = limits.pmax;
  }
  if (limits.pmin.size != 0) {
    Assert(limits.pmin.size == joints);
    pmin = limits.pmin;
  }
  if (limits.vmax.size != 0) {
    Assert(limits.vmax.size == joints);
    vmax = limits.vmax;
  }
  if (limits.vmin.size != 0) {
    Assert(limits.vmin.size == joints);
    vmin = limits.vmin;
  }
  if (limits.amax.size != 0) {
    Assert(limits.amax.size == joints);
    amax = limits.amax;
  }
  if (limits.amin.size != 0) {
    Assert(limits.amin.size == joints);
    amin = limits.amin;
  }
  if (limits.tau_max.size != 0) {
    Assert(limits.tau_max.size == joints);
    tau_max = limits.tau_max;
  }
  if (limits.tau_min.size != 0) {
    Assert(limits.tau_min.size == joints);
    tau_min = limits.tau_min;
  }
  if (limits.tcp_max) {
    tcp_max = limits.tcp_max;
  }
}

inline host_fn void Manipulator::set_kinematics(const ManipulatorKinematics& kinematics) {
  Assert(kinematics.dv.size() == joints);
  Assert(kinematics.ev.size() == joints);
  dv       = kinematics.dv;
  ev       = kinematics.ev;
  Q_static = kinematics.Q_static;
  p_base   = kinematics.p_base;
  Q_base   = kinematics.Q_base;
  p_j0     = kinematics.p_j0;
  _Q.resize(joints);
  _Q_mult.resize(joints);
  _p_j.resize(joints + 1);
}

inline host_fn void Manipulator::set_dynamics(const ManipulatorDynamics& dynamics) {
  Assert(dynamics.m.size() == joints);
  Assert(dynamics.I.size() == joints);
  Assert(dynamics.av.size() == joints);
  m  = dynamics.m;
  I  = dynamics.I;
  av = dynamics.av;
  _efforts.resize(joints);
  sv.resize(joints);
  for (u32 j = 0; j < joints; j++) {
    sv[j] = {-dv[j] + av[j]};
  }
}

inline host_fn void Manipulator::set_capsules(const ManipulatorCapsules& capsules) {
  for (auto& cap: capsules.capsule_list) {
    Assert(cap.joint_frame < _p_j.size());
    _collision_model.push_back(cap);
  }
  _n_caps = _collision_model.size();
  _capsule_list.resize(_n_caps);
  if (capsules.collision_base.size != 0 && capsules.collision_matrix.size != 0) { // todo: fix this... (we need a fail safe)
    _collision_base.resize(_n_caps);
    _collision_base = capsules.collision_base;
    _collision_matrix.resize(_n_caps, _n_caps);
    _collision_matrix = capsules.collision_matrix;
    for (u32 i = 0; i < _n_caps; ++i) {
      for (u32 j = i + 1; j < _n_caps; ++j) {
        _n_internal_collisions += (_collision_matrix(j, i) != 0);
      }
      _n_internal_collisions += (_collision_base[i] != 0);
    }
    _base_sphere = capsules.base_sphere;
    _base_sphere.c += p_base;
    // _n_internal_collisions = 0;
  }
}

inline host_fn void Manipulator::compute_rotation_matrices(const Array& joint_position) {
  Assert(joint_position.size == joints);
  Array s(joints), c(joints);
  sincos(joint_position, s, c);
  for (u32 j = 0; j < joints; ++j) {
    Mat3 temp = {c[j], s[j], 0, -s[j], c[j], 0, 0, 0, 1};
    _Q[j]     = Q_static[j] * temp;
  }
  for (u32 j = joints; j < _Q.size(); ++j) {
    _Q[j] = Q_static[j];
  }
}

inline host_fn void Manipulator::compute_capsules() {
  for (u32 i = 0; i < _collision_model.size(); ++i) {
    _capsule_list[i] = {
            {_p_j[_collision_model[i].joint_frame] + _Q_mult[_collision_model[i].joint_frame] * _collision_model[i].p1},
            {_p_j[_collision_model[i].joint_frame] + _Q_mult[_collision_model[i].joint_frame] * _collision_model[i].p2},
            _collision_model[i].r};
  }
}

inline host_fn Array Manipulator::get_internal_collisions() const {
  Array distances(_n_internal_collisions);
  u32   idx = 0;
  for (u32 i = 0; i < _collision_matrix.cols; ++i) {
    for (u32 j = i + 1; j < _collision_matrix.rows; ++j) {
      if (_collision_matrix(j, i) != 0)
        distances[idx++] = distance(_capsule_list[i], _capsule_list[j]);
    }
  }
  for (u32 i = 0; i < _collision_base.size; ++i) {
    if (_collision_base[i] != 0)
      distances[idx++] = distance(_capsule_list[i], _base_sphere);
  }
  return distances;
}

inline host_fn Array Manipulator::get_tool_pose() const {
  Mat3 tmp             = _Q_mult.back();
  Vec3 current_rpy     = rotation2rpy(tmp);
  Vec3 current_cartpos = _p_j.back();
  return {current_cartpos.x, current_cartpos.y, current_cartpos.z,
          current_rpy.x, current_rpy.y, current_rpy.z};
}

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


} // namespace blast
