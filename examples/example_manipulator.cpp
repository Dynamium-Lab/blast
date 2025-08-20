
#define BLAST_TRACE_LEVEL 0
// #define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"

int main() {
  using namespace blast;

  Optimization opt(get_generic_gen3(), get_gen3_task());

  opt.world = get_lab_world();

  opt.guess.type   = Guess::random;
  opt.guess.n_shot = 100;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = true;
  opt.constraints.self_collisions     = true;
  opt.constraints.external_collisions = true;

  opt.max_tries = 1;
  opt.success_tolerance = 0.01;

  auto result = optimize(&opt);

  cout << "Compute time:        " << result.compute_time << endl;
  cout << "Trajectory time:     " << result.x[result.x.size - 1] << endl;
  cout << "Trajectory success:  " << result.success << endl;
  cout << "Number of tries:     " << result.num_tries << endl;

  return 0;
}
