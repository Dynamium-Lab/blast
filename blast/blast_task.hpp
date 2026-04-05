#pragma once
#include <blast>

namespace blast {


struct TaskBoundary {
  Array position;
  Array velocity;
  Array acceleration;
};

struct Task {
  TaskBoundary start;
  TaskBoundary goal;

  Task() = default;

  explicit Task(u32 n_joints) :
      start{Array(n_joints), Array(n_joints), Array(n_joints)},
      goal{Array(n_joints), Array(n_joints), Array(n_joints)} {
    for (u32 i = 0; i < n_joints; i++) {
      start.position[i] = start.velocity[i] = start.acceleration[i] = 0.0f;
      goal.position[i]  = goal.velocity[i]  = goal.acceleration[i]  = 0.0f;
    }
  }

  static Task stop_to_stop(Array start_pos, Array goal_pos) {
    Task t(start_pos.size);
    t.start.position = std::move(start_pos);
    t.goal.position  = std::move(goal_pos);
    return t;
  }

  Matrix to_matrix() const {
    u32 n = start.position.size;
    Matrix m(n, 6);
    for (u32 i = 0; i < n; i++) {
      m(i, 0) = start.position[i];
      m(i, 1) = start.velocity[i];
      m(i, 2) = start.acceleration[i];
      m(i, 3) = goal.position[i];
      m(i, 4) = goal.velocity[i];
      m(i, 5) = goal.acceleration[i];
    }
    return m;
  }
};


}