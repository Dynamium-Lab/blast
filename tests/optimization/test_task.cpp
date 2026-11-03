#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <blast>
#include <catch2/catch.hpp>

using namespace blast;

TEST_CASE("Task stop_to_stop sets start and goal positions") {
  Array q_start = {(real) 0, -PI / 2, (real) 0, -PI / 2, (real) 0, (real) 0};
  Array q_goal  = {-(real) 1.54, -(real) 1.83, -(real) 2.28, -(real) 0.59, (real) 1.60, (real) 0.023};
  Task  task    = Task::stop_to_stop(q_start, q_goal);

  CHECK(is_close(task.start.position, q_start));
  CHECK(is_close(task.goal.position, q_goal));
  CHECK(is_small(task.start.velocity));
  CHECK(is_small(task.goal.velocity));
  CHECK(is_small(task.start.acceleration));
  CHECK(is_small(task.goal.acceleration));
}

TEST_CASE("Task n_joints constructor initializes to zero") {
  Task task(6);
  for (u32 i = 0; i < 6; i++) {
    CHECK(task.start.position[i] == (real) 0);
    CHECK(task.start.velocity[i] == (real) 0);
    CHECK(task.start.acceleration[i] == (real) 0);
    CHECK(task.goal.position[i] == (real) 0);
    CHECK(task.goal.velocity[i] == (real) 0);
    CHECK(task.goal.acceleration[i] == (real) 0);
  }
}

TEST_CASE("Task round-trips through to_matrix") {
  Array  q_start = {(real) 0, -PI / 2, (real) 0, -PI / 2, (real) 0, (real) 0};
  Array  q_goal  = {-(real) 1.54, -(real) 1.83, -(real) 2.28, -(real) 0.59, (real) 1.60, (real) 0.023};
  Task   task1   = Task::stop_to_stop(q_start, q_goal);
  Matrix m       = task1.to_matrix();
  Task   task2(m);

  CHECK(is_close(task2.start.position, task1.start.position));
  CHECK(is_close(task2.goal.position, task1.goal.position));
  CHECK(is_close(task2.start.velocity, task1.start.velocity));
  CHECK(is_close(task2.goal.velocity, task1.goal.velocity));
  CHECK(is_close(task2.start.acceleration, task1.start.acceleration));
  CHECK(is_close(task2.goal.acceleration, task1.goal.acceleration));
}
