
#define BLAST_TRACE_LEVEL 3
#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"

int main() {
  using namespace blast;

  // Optimization opt(get_generic_gen3(), get_gen3_task());
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 100, 5, 6);
  opt.bspline = bspline;

  opt.guess.type = Guess::custom;
  opt.guess.x0   = Array(opt.bspline.x_len(opt.task), 2.0);
  // opt.guess.n_shot = 100;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = false;
  opt.constraints.self_collisions     = false;
  opt.constraints.external_collisions = false;

  opt.max_tries         = 1;
  opt.success_tolerance = 0.01;
  // Sleep(100);

  auto result = optimize(&opt);

  cout << "Compute time:        " << result.compute_time << endl;
  cout << "Function evaluations:" << result.num_eval << endl;
  cout << "Trajectory time:     " << result.x[result.x.size - 1] << endl;
  cout << "Success:             " << result.success << endl;
  cout << "False success:       " << result.success_false << endl;
  cout << "Number of tries:     " << result.num_tries << endl;
  cout << "Stopping criteria:   " << result.nlopt_exit_criteria << endl;

  return 0;
}
