
#define BLAST_TRACE_LEVEL 2
#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"

int main() {
  using namespace blast;

  // --- Test characteristics ---
  bool position            = true;
  bool velocity            = true;
  bool acceleration        = true;
  bool torque              = true;
  bool tcp_speed           = true;
  bool self_collisions     = true;
  bool external_collisions = true;


  auto link6 = get_generic_Link6();
  auto task  = get_link6_task();

  Optimization opt(link6, task);

  opt.world.add_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1});  // vertical plate (no coll)
  opt.world.add_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, {1, 0, 0, 0, 1, 0, 0, 0, 1}); // table

  opt.guess.type = Guess::custom;
  opt.guess.x0 = Array(opt.bspline.x_len(task), 1.0);

  opt.constraints.position            = position;
  opt.constraints.velocity            = velocity;
  opt.constraints.acceleration        = acceleration;
  opt.constraints.torque              = torque;
  opt.constraints.tcp_speed           = tcp_speed;
  opt.constraints.self_collisions     = self_collisions;
  opt.constraints.external_collisions = external_collisions;

  opt.max_restart = 10;

  auto result = optimize(&opt);

  cout << "Compute time: " << result.compute_time << endl;
  cout << "Trajectory time: " << result.x[result.x.size - 1] << endl;
  cout << "Trajectory success: " << result.success << endl;
  cout << "Number of restarts: " << result.num_restart << endl;

  return 0;
}
