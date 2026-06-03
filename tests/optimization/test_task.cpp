#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <blast>
#include <catch2/catch.hpp>

using namespace blast;

TEST_CASE("Task stop_to_stop sets start and goal positions") {
  Array q_start = {0.0f, -PI / 2, 0.0f, -PI / 2, 0.0f, 0.0f};
  Array q_goal  = {-1.54f, -1.83f, -2.28f, -0.59f, 1.60f, 0.023f};
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
    CHECK(task.start.position[i] == 0.0f);
    CHECK(task.start.velocity[i] == 0.0f);
    CHECK(task.start.acceleration[i] == 0.0f);
    CHECK(task.goal.position[i] == 0.0f);
    CHECK(task.goal.velocity[i] == 0.0f);
    CHECK(task.goal.acceleration[i] == 0.0f);
  }
}

TEST_CASE("Task round-trips through to_matrix") {
  Array  q_start = {0.0f, -PI / 2, 0.0f, -PI / 2, 0.0f, 0.0f};
  Array  q_goal  = {-1.54f, -1.83f, -2.28f, -0.59f, 1.60f, 0.023f};
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
