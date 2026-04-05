#include <blast>

const int total = 10000;
const int units = 10;

void fun() {
  constexpr int time   = total / units; // us
  auto          start  = blast::get_tick_us();
  auto          target = start + time;
  while (blast::get_tick_us() < target) {
  }
}

int main() {
  tf::Executor executor;
  tf::Taskflow taskflow;

  // auto [A, B, C, D] =
  for (int i = 0; i < units; i++) {
    taskflow.emplace(fun);
  }

  // taskflow.dump(std::cout);
  // A.precede(B, C);  // A runs before B and C
  // D.succeed(B, C);  // D runs after  B and C
  {
    auto start = blast::get_tick_us();
    executor.run(taskflow).wait();
    auto stop = blast::get_tick_us();
    std::cout << stop - start << std::endl;
  }
  {
    auto start = blast::get_tick_us();
    for (int i = 0; i < units; i++) {
      fun();
    }
    auto stop = blast::get_tick_us();
    std::cout << stop - start << std::endl;
  }
  return 0;
}
