#pragma once

#include <blast>
#include "blast_world.hpp"

namespace blast {
// Defined Types
struct Manipulator;
struct ManipulatorLimits;
struct ManipulatorKinematics;
struct ManipulatorDynamics;
struct ManipulatorCapsules;
struct CollisionModelCapsule;
struct EndEffector;
struct ManipulatorTempData;

// todo: which of these need to be exposed to the user?
inline host_fn void forward_kinematics(const Manipulator& manip, ManipulatorTempData& temp, const Array& joint_pos);

inline host_fn Matrix jacobian(const Manipulator& manip, const ManipulatorTempData& temp);
inline host_fn void   dynamics(const Manipulator& manip, ManipulatorTempData& temp, Array& vel, Array& acc);
// double                get_error(unsigned int n, const double* x, double* grad, void* data);
// inline Array          inverse_kinematics_nlopt(Manipulator manip, Array desired_pose, Array initial_joint_position);

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
  std::array<Vec3, MAX_JOINTS>     dv;                                   // vector to next joint
  std::array<Vec3, MAX_JOINTS>     ev;                                   // direction vectors of joint
  std::array<Mat3, MAX_JOINTS + 1> Q_static;                             // static rotation to next joint
  Vec3                             p_base = {0.0, 0.0, 0.0};             // base position in workspace
  Vec3                             p_j0   = {0.0, 0.0, 0.0};             // position of the first joint relative to base
  Mat3                             Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // base orientation in workspace
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
  std::array<real, MAX_JOINTS> m{}; // link mass
  std::array<Mat3, MAX_JOINTS> I{}; // inertial tensors
  std::array<Vec3, MAX_JOINTS> av{};
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


struct ManipulatorTempData {
  std::array<real, MAX_JOINTS>      efforts{};
  std::array<Mat3, MAX_JOINTS>      rotations{};
  std::array<Mat3, MAX_JOINTS>      rotations_mult{};
  std::array<Vec3, MAX_JOINTS + 1>  p_j{};
  std::array<Capsule, CAPSULE_SIZE> capsule_list{};
};

/**
 * @class Manipulator
 * @brief High‐level interface combining limits, kinematics, dynamics,
 *        and collision geometry for a manipulator.
 *
 * @note Use the constructor to fully specify at least limits and kinematics.
 */
struct Manipulator {
  u32 n_joints = 0;

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
  Vec3                             p_j0 = {0, 0, 0};
  std::array<Vec3, MAX_JOINTS>     dv{};       // to next joint
  std::array<Vec3, MAX_JOINTS>     ev{};       // axis dirs
  std::array<Mat3, MAX_JOINTS + 1> Q_static{}; // static rotations

  // Link dynamics
  std::array<real, MAX_JOINTS> m{};
  std::array<Mat3, MAX_JOINTS> I{};
  std::array<Vec3, MAX_JOINTS> av{};
  std::array<Vec3, MAX_JOINTS> sv{}; // COM from next joint

  // Internal caches and collision
  int                                             _n_caps                = 0;
  int                                             _n_internal_collisions = 0;
  std::array<CollisionModelCapsule, CAPSULE_SIZE> _collision_model{};
  ObjMatrix<u8>                                   _collision_matrix{};
  Array                                           _collision_base{};
  Sphere                                          _base_sphere{};

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
   * @brief Update collision capsules in world frame (store internally).
   */
  // host_fn void compute_capsules();

  /**
   * @brief Compute distances for all internal collisions.
   * @return Array of distances between colliding primitives.
   */
  // host_fn Array get_internal_collisions() const;

  /**
   * @brief Get the current tool pose as [x,y,z,roll,pitch,yaw].
   * @return Array of six values: position and RPY orientation.
   */
  host_fn Array get_tool_pose() const;
};

} // namespace blast

#include "manipulator/manipulator.hpp"
// #include "manipulator/Gen3.hpp"
// #include "manipulator/Link6.hpp"
// #include "manipulator/R2.hpp"
