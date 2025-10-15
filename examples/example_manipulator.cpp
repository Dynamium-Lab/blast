
#define BLAST_TRACE_LEVEL 0
#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"

int main() {
  using namespace blast;

  // Optimization opt(get_generic_gen3(), get_gen3_task());
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  // Sleep(200);

  Bspline bspline(16, 110, 5, 6);
  opt.bspline = bspline;

  // opt.guess.type = Guess::custom;
  // opt.guess.x0   = Array(opt.bspline.x_len(opt.task), 1.5);
  // opt.guess.type = Guess::random;
  // opt.guess.n_shot = 100;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = true;
  opt.constraints.self_collisions     = true;
  opt.constraints.external_collisions = true;

  opt.constraints.n_collision_constraints = 1;
  opt.max_tries                           = 1;
  opt.success_tolerance                   = 0.01;
  // Sleep(100);

  // int  max_x       = 2;
  // real increment_x = 0.1;
  int max_iter = 10;

  opt.guess.type = Guess::custom;

  auto   t1 = get_tick_us();
  Result result(&opt);
  // for (real i = increment_x; i < max_x; i += increment_x) {
  for (real i = 0; i < max_iter; i++) {

    // opt.guess.x0        = random_array(opt.bspline.x_len(opt.task), 1);
    // opt.guess.x0.back() = 0.5;
    // result = optimize(&opt);
    result = optimize_with_segments(&opt);

    cout << "Compute time:        " << result.compute_time << endl;
    cout << "Function evaluations:" << result.num_eval << endl;
    cout << "Trajectory time:     " << result.x[result.x.size - 1] << endl;
    cout << "Success:             " << result.success << endl;
    cout << "False success:       " << result.success_false << endl;
    cout << "Number of tries:     " << result.num_tries << endl;
    cout << "Stopping criteria:   " << result.nlopt_exit_criteria << endl;
    cout << endl;
  }
  auto t2 = get_tick_us();
  cout << "Compute time:        " << (double) (t2 - t1) / 1000. << endl;
  cout << "Average compute time:        " << (double) (t2 - t1) / 1000. / (real) max_iter << endl;
  // cout << "Initial guess time:        " << time_for_one_guess << endl;


  return 0;
}
