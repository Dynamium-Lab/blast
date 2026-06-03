//
// Created by nikos on 2025-09-09.
//
#include <blast>

using namespace blast;

int main() {
  int p      = 5;
  int points = 110;
  int joints = 6;
  int n_ctrl = 16;

  Bspline spline(n_ctrl, points, p, joints);
  
  Task  task(joints);
  Array x(spline.x_len(task.to_matrix()));                                            // note: only works this way if task has no NaN values later on

  task.start.position = {2.04549, 0.675615, 0.137326, -0.814844, 0.449157, -3.08683}; // 4
  task.goal.position  = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};          // 0
  x                   = {2.012367432705342,
                         1.9183101026572604,
                         1.7497289098956235,
                         1.7056521658986308,
                         1.7223101992224683,
                         1.823821209145915,
                         2.029475755747117,
                         2.2359775337758414,
                         2.4616270467628607,
                         2.5472728000018465,
                         0.7179444714026678,
                         0.6747835328459124,
                         0.5905813724701455,
                         0.4732084367682011,
                         0.3877419627152058,
                         0.2679800405036057,
                         0.20628662643489376,
                         0.26175257983114314,
                         0.10244647176615736,
                         0.13085517853642165,
                         0.14688471958586236,
                         0.07543553241423792,
                         -0.060942179699366515,
                         -0.22408329710890498,
                         -0.4691395740014129,
                         -0.6189306230897347,
                         -0.6778822187674548,
                         -0.6447207604978721,
                         -0.5187297319175941,
                         -0.41538574895597546,
                         -0.7804945384768247,
                         -0.6930394133807527,
                         -0.4766589490896912,
                         -0.208616172275154,
                         0.03676937270637994,
                         0.285355337026231,
                         0.4269786278073888,
                         0.48795758655073995,
                         0.4505959960586352,
                         0.3575488334455236,
                         0.43912758699716037,
                         0.38082781706673297,
                         0.18780296837181842,
                         0.10228164417466619,
                         0.10092979189901205,
                         0.1966984370981407,
                         0.37117795848532414,
                         0.6580317779330959,
                         0.8703471092778267,
                         0.9605537399961706,
                         -3.0677587708249896,
                         -3.0310307373237357,
                         -2.8258541900851637,
                         -2.820323386173435,
                         -2.710678680133922,
                         -2.7566085960780518,
                         -2.803503611759136,
                         -2.9844253539900607,
                         -3.1311235621554494,
                         -3.163875191630662,
                         0.893};


  print_to_csv(transpose(spline.basis_p), "../../../examples/bsplines_task32/basis_function_p.csv");
  print_to_csv(transpose(spline.basis_v), "../../../examples/bsplines_task32/basis_function_v.csv");
  print_to_csv(transpose(spline.basis_a), "../../../examples/bsplines_task32/basis_function_a.csv");

  spline.compute_trajectory(x, task.to_matrix());
  print_to_csv(transpose(spline.traj.pos), "../../../examples/bsplines_task32/bspline_p.csv");
  print_to_csv(transpose(spline.traj.vel), "../../../examples/bsplines_task32/bspline_v.csv");
  print_to_csv(transpose(spline.traj.acc), "../../../examples/bsplines_task32/bspline_a.csv");
}
