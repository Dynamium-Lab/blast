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
struct Tool;
struct ManipulatorTempData;

// todo: which of these need to be exposed to the user?
inline host_fn void forward_kinematics(const Manipulator& manip, ManipulatorTempData& temp, const Array& joint_pos);

inline host_fn Matrix jacobian(const Manipulator& manip, const ManipulatorTempData& temp);
inline host_fn void   dynamics(const Manipulator& manip, ManipulatorTempData& temp, const Array& vel, const Array& acc);
// double                get_error(unsigned int n, const double* x, double* grad, void* data);
// inline Array          inverse_kinematics_nlopt(Manipulator manip, Array desired_pose, Array initial_joint_position);

/**
 * @struct ManipulatorLimits
 * @brief Actuator limits for a manipulator's joints and Tool.
 *
 * @var position_max    Maximum joint positions.
 * @var position_min    Minimum joint positions.
 * @var velocity_max    Maximum joint velocities.
 * @var acceleration_max Maximum joint accelerations.
 * @var torque_max      Maximum joint torques.
 * @var tool_speed_max   Maximum tool‐center‐point speed.
 */
struct ManipulatorLimits {
  Array position_max;     // max joint position
  Array position_min;     // min joint position
  Array velocity_max;     // max joint velocity
  Array acceleration_max; // max joint acceleration
  Array jerk_max;         // max joint jerk
  Array torque_max;       // max joint torque

  real tool_speed_max;    // max tool speed
};

/**
 * @struct ManipulatorKinematics
 * @brief Geometric configuration of manipulator joints and base.
 *
 * @var joint_offsets     Vectors from each joint to the next (in joint frames).
 * @var joint_axes        Unit axis directions for each joint.
 * @var static_rotations  Static rotation matrices from joint to next (Mat3).
 * @var first_joint_position Position of the first joint relative to the base.
 * @var base_position     Position of the manipulator base in the world frame.
 * @var base_rotation     Orientation of the manipulator base in the world frame.
 */
struct ManipulatorKinematics {
  std::array<Vec3, MAX_JOINTS>     joint_offsets;                                      // vector to next joint
  std::array<Vec3, MAX_JOINTS>     joint_axes;                                         // direction vectors of joint
  std::array<Mat3, MAX_JOINTS + 1> static_rotations;                                   // static rotation to next joint
  Vec3                             base_position        = {0.0, 0.0, 0.0};             // base position in workspace
  Vec3                             first_joint_position = {0.0, 0.0, 0.0};             // position of the first joint relative to base
  Mat3                             base_rotation        = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // base orientation in workspace
};

/**
 * @struct ManipulatorDynamics
 * @brief Dynamic properties of each manipulator link.
 *
 * @var link_masses       Mass of each link.
 * @var inertia_tensors   Inertia tensor of each link (Mat3).
 * @var cog_offsets       Center of mass offset for each link (Vec3).
 */
struct ManipulatorDynamics {
  std::array<real, MAX_JOINTS> link_masses{};     // link mass
  std::array<Mat3, MAX_JOINTS> inertia_tensors{}; // inertial tensors
  std::array<Vec3, MAX_JOINTS> cog_offsets{};     // center of gravity offset
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
 * @struct Payload
 * @brief Kinematic, Dynamic and Collision properties of a payload.
 *
 * @var position        Position vector from Tool to payload.
 * @var rotation        Rotation from Tool to payload.
 * @var mass            Payload mass.
 * @var inertia_tensor  Payload inertia tensor.
 * @var cog_position    Payload center of mass position in local frame.
 * @var collision_model Collision model for the payload.
 */
struct Payload {
  Vec3 position = {0.0, 0.0, 0.0};
  Mat3 rotation = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  real mass     = 0.0;
  Mat3 inertia_tensor;
  Vec3 cog_position = {0.0, 0.0, 0.0}; // center of gravity position

  CollisionModel collision_model;
};

/**
 * @struct Tool
 * @brief Kinematic, Dynamic and Collision properties of a tool.
 *
 * @var position        Position vector from last joint to tool.
 * @var rotation        Rotation from last joint to tool.
 * @var mass            Tool mass.
 * @var inertia_tensor  Tool inertia tensor.
 * @var cog_position    Tool center of mass position in local frame.
 * @var collision_model Collision model for the tool.
 * @var payload Payload being manipulated.
 */
struct Tool {
  Vec3 position = {0.0, 0.0, 0.0};
  Mat3 rotation = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  real mass     = 0.0;
  Mat3 inertia_tensor;
  Vec3 cog_offset = {0.0, 0.0, 0.0}; // center of gravity offset

  CollisionModel collision_model;

  bool    has_payload = false;
  Payload payload;

  inline host_fn void set_payload(Payload payload) {
    has_payload = true;
    payload     = payload;
  }
};

struct ManipulatorTempData {
  std::array<real, MAX_JOINTS>      efforts{};
  std::array<Mat3, MAX_JOINTS>      rotations{};
  std::array<Mat3, MAX_JOINTS>      rotations_mult{};
  std::array<Vec3, MAX_JOINTS + 1>  p_j{};
  std::array<Capsule, MAX_CAPSULES> capsule_list{};
};

/**
 * @class Manipulator
 * @brief High‐level interface combining limits, kinematics, dynamics,
 *        and collision geometry for a manipulator.
 *
 * @note Use the constructor to fully specify at least limits and kinematics.
 */
// todo: Add templating with number of joints instead of a global MAX_JOINTS
struct Manipulator {
  u32 n_joints = 0;

  // Manipulator limits
  std::array<real, MAX_JOINTS> position_max{};
  std::array<real, MAX_JOINTS> position_min{};
  std::array<real, MAX_JOINTS> velocity_max{};
  std::array<real, MAX_JOINTS> acceleration_max{};
  std::array<real, MAX_JOINTS> jerk_max{};
  std::array<real, MAX_JOINTS> torque_max{};

  real tool_speed_max = 0.0; // max Tool speed

  // World‐frame base pose
  Vec3 base_position = {0, 0, 0};
  Mat3 base_rotation = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  // Tool state
  bool has_tool = false;
  Tool tool;

  // Joint kinematics
  Vec3                             first_joint_position = {0, 0, 0};
  std::array<Vec3, MAX_JOINTS>     joint_offsets{};    // to next joint
  std::array<Vec3, MAX_JOINTS>     joint_axes{};       // axis dirs
  std::array<Mat3, MAX_JOINTS + 1> static_rotations{}; // static rotations

  // Link dynamics
  std::array<real, MAX_JOINTS> link_masses{};
  std::array<Mat3, MAX_JOINTS> inertia_tensors{};
  std::array<Vec3, MAX_JOINTS> cog_offsets{};         // center of gravity offset
  std::array<Vec3, MAX_JOINTS> cog_from_next_joint{}; // center of gravity offset

  // Internal caches and collision
  int                                             _n_caps                = 0;
  int                                             _n_internal_collisions = 0;
  std::array<CollisionModelCapsule, MAX_CAPSULES> _collision_model{};
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
   * @brief Set a tool to this manipulator, and adds dynamic and collision contributions to last link.
   * @param tool  Tool structure and geometry.
   */
  host_fn void set_tool(const Tool& tool);
  host_fn void remove_tool();

  /**
   * @brief Set a payload to this manipulator's tool, and adds dynamic and collision contributions to last link.
   * @param payload Payload structure and geometry
   */
  host_fn void set_payload(const Payload& payload);
  host_fn void remove_payload();

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

#include "manipulator/Kinova_Gen3.hpp"
#include "manipulator/Kinova_Link6.hpp"
#include "manipulator/UR5e.hpp"
#include "manipulator/manipulator.hpp"
