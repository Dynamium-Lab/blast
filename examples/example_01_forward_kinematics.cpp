// example_01_forward_kinematics.cpp
//
// This example shows how to:
//   1. Construct a UR5e Manipulator
//   2. Run forward kinematics for a given joint configuration
//   3. Read the resulting joint frame positions (including the TCP)
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_01_forward_kinematics
// Run:
//   ./build/examples/example_01_forward_kinematics

#include <blast>
#include <iostream>
#include "ur5e.hpp"

int main() {
  using namespace blast;

  // -----------------------------------------------------------------------
  // Step 1 — Build the robot model.
  // make_ur5e() returns a fully configured UR5e with limits, kinematics,
  // dynamics, and collision geometry. No external files are needed.
  // -----------------------------------------------------------------------
  Manipulator ur5e = make_ur5e();

  // Optionally reposition the robot in the world (metres, rotation matrix).
  ur5e.base_position = {0.0, 0.0, 0.0};  // world-frame origin
  ur5e.base_rotation = {1, 0, 0,          // identity rotation (robot upright)
                        0, 1, 0,
                        0, 0, 1};

  // -----------------------------------------------------------------------
  // Step 2 — Provide joint positions (radians).
  // ManipulatorTempData is a scratch buffer filled in-place by forward
  // kinematics; it holds frame positions, rotation matrices, etc.
  // -----------------------------------------------------------------------
  ManipulatorTempData temp;

  Array config_home = {0.0, -1.5708, 1.5708, -1.5708, -1.5708, 0.0}; // ~upright

  // -----------------------------------------------------------------------
  // Step 3 — Run forward kinematics.
  // After this call, temp.p_j[i] holds the world-frame position of frame i.
  // Frame 0 is the base; frame n_joints is the TCP (tool centre point).
  // -----------------------------------------------------------------------
  forward_kinematics(ur5e, temp, config_home);

  // -----------------------------------------------------------------------
  // Step 4 — Print results.
  // -----------------------------------------------------------------------
  int tcp_frame = ur5e.n_joints; // index of the TCP frame
  std::cout << "--- Forward kinematics: home configuration ---\n";
  std::cout << "Joint positions (rad): ";
  for (int i = 0; i < ur5e.n_joints; ++i)
    std::cout << config_home[i] << (i + 1 < ur5e.n_joints ? ", " : "\n");

  std::cout << "TCP position (m):\n";
  std::cout << "  x = " << temp.p_j[tcp_frame].x << "\n";
  std::cout << "  y = " << temp.p_j[tcp_frame].y << "\n";
  std::cout << "  z = " << temp.p_j[tcp_frame].z << "\n";

  // Print all intermediate frame positions for inspection.
  std::cout << "\nAll frame positions:\n";
  for (int i = 0; i <= tcp_frame; ++i) {
    std::cout << "  frame[" << i << "]: ("
              << temp.p_j[i].x << ", "
              << temp.p_j[i].y << ", "
              << temp.p_j[i].z << ")\n";
  }

  // Run a second configuration to demonstrate that temp is reusable.
  Array config_stretched = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  forward_kinematics(ur5e, temp, config_stretched);

  std::cout << "\n--- Forward kinematics: all-zeros configuration ---\n";
  std::cout << "TCP position (m):\n";
  std::cout << "  x = " << temp.p_j[tcp_frame].x << "\n";
  std::cout << "  y = " << temp.p_j[tcp_frame].y << "\n";
  std::cout << "  z = " << temp.p_j[tcp_frame].z << "\n";

  return 0;
}
