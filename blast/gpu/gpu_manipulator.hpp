#pragma once
#include <math_constants.h>
#include "blast_math.hpp"
#include "gpu/trajectory.cuh"

namespace blast {
__global__ void compute_constraints_kernel();

struct cuGen3MultiTraj {
  Vec3 p_base;
  real m[7];
  Mat3 I[7];
  Vec3 av[7];
  Vec3 dv[7];
  Vec3 sv[7];
  Vec3 ev[7];

  real pmax[7];
  real vmax[7];
  real tau_max[7];

  u32  ntrajectories = 0;
  u32  npoints       = 0;
  bool is_init       = false;
  bool is_alias      = false;

  real* dev_constraints;
  real* host_constraints;

  // HOST
  //----------

  host_fn cuGen3MultiTraj(u32 points, u32 trajs = 1, real mass = 0.0);

  // note: copy constructor creates an alias
  host_fn cuGen3MultiTraj(const cuGen3MultiTraj& other) {
    *this    = other;
    is_alias = true;
  }

  host_fn void compute_constraints(const Array& x, const Matrix& task, cuMultiBsplines& gpu_bsplines);
  host_fn void fetch_constraints();
  host_fn void _init(real mass);
  host_fn ~cuGen3MultiTraj();

  // GENERAL
  //----------

  // Get the start id of the constraint for a given point on a trajectory
  blast_fn u32 con_id(u32 point, u32 traj);
  // get the total number of constraints
  blast_fn u32 nconstraints();

  // DEVICE
  //----------

  // compute the constraints based on the current trajectory point
  dev_fn void compute_constraints_point(const real pos[7], const real vel[7], const real acc[7], real* con);
};

// note: must be a global variable because it is __constant__
// note: much faster if it's __constant__ because all threads access the same location (heavy use of broadcast operations)
const int         MANIP_BROADCAST_ARENA_SIZE = 256;
__constant__ real manip_broadcast_arena[MANIP_BROADCAST_ARENA_SIZE];

//------ HOST FUNCTIONS ------------------------------------------------------------------------------------

host_fn cuGen3MultiTraj::cuGen3MultiTraj(u32 _points, u32 _trajs, real mass) :
    npoints(_points),
    ntrajectories(_trajs) {

  assert_buffer_size<cuGen3MultiTraj, sizeof(manip_broadcast_arena)>(); // compile time check
  cuda_check(cudaMalloc(&dev_constraints, nconstraints() * sizeof(real)));
  cuda_check(cudaMallocHost(&host_constraints, nconstraints() * sizeof(real)));

  this->_init(mass);

  // note: this must be done last
  cuda_check(cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0));
}

host_fn void cuGen3MultiTraj::_init(real mass) {

  // position of the first joint with respect to the table in the center of the base
  p_base = {0, 0, (real) 0.1564};

  // mass
  m[0] = (real) 1.377;
  m[1] = (real) 1.1636;
  m[2] = (real) 1.1636;
  m[3] = (real) 0.93;
  m[4] = (real) 0.678;
  m[5] = (real) 0.678;
  m[6] = (real) 0.364 + (real) 0.921;

  // inertial tensor
  I[0] = {(real) 0.004570, (real) 0.000001, (real) 0.000002, (real) 0.000001, (real) 0.004831, (real) 0.000448, (real) 0.000002, (real) 0.000448, (real) 0.001409};
  I[1] = {(real) 0.011088, (real) 0.000005, (real) 0, (real) 0.000005, (real) 0.001072, (real) 0.000691, (real) 0, (real) 0.000691, (real) 0.011255};
  I[2] = {(real) 0.010932, (real) 0, (real) 0.000007, (real) 0, (real) 0.011127, (real) 0.000606, (real) 0.000007, (real) 0.000606, (real) 0.001043};
  I[3] = {(real) 0.008147, (real) 0.000001, (real) 0, (real) 0.000001, (real) 0.000631, (real) 0.000500, (real) 0, (real) 0.000500, (real) 0.008316};
  I[4] = {(real) 0.001596, (real) 0, (real) 0, (real) 0, (real) 0.001607, (real) 0.000256, (real) 0, (real) 0.000256, (real) 0.000399};
  I[5] = {(real) 0.001641, (real) 0, (real) 0, (real) 0, (real) 0.000410, (real) 0.000278, (real) 0, (real) 0.000278, (real) 0.001641};
  I[6] = {(real) 0.000214, (real) 0, (real) 0.000001, (real) 0, (real) 0.000223, (real) 0.000002, (real) 0.000001, (real) 0.000002, (real) 0.000240};

  // center of mass
  av[0] = {-(real) 0.000023, -(real) 0.010364, -(real) 0.073360};
  av[1] = {-(real) 0.000044, -(real) 0.099580, -(real) 0.013278};
  av[2] = {-(real) 0.000044, -(real) 0.006641, -(real) 0.117892};
  av[3] = {-(real) 0.000018, -(real) 0.075478, -(real) 0.015006};
  av[4] = {(real) 0.000001, -(real) 0.009432, -(real) 0.063883};
  av[5] = {(real) 0.000001, -(real) 0.045483, -(real) 0.009650};
  av[6] = {-(real) 0.000093, (real) 0.000132, -(real) 0.022905};
  Vec3 av_tool(0, 0, -(real) 0.06 - (real) 0.0615);
  av[6] = ((real) 0.364 * av[6] + (real) 0.921 * av_tool) * (1 / ((real) 0.364 + (real) 0.921));

  // vector to next joint
  dv[0] = {(real) 0, (real) 0.0054, -(real) 0.1284};
  dv[1] = {(real) 0, -(real) 0.2104, -(real) 0.0064};
  dv[2] = {(real) 0, -(real) 0.0064, -(real) 0.2104};
  dv[3] = {(real) 0, -(real) 0.2084, -(real) 0.0064};
  dv[4] = {(real) 0, (real) 0, -(real) 0.1059};
  dv[5] = {(real) 0, -(real) 0.1059, (real) 0};
  dv[6] = {(real) 0, (real) 0, -(real) 0.0615 - (real) 0.164};

  // center of mass (from next joint)
  sv[0] = dv[0] - av[0];
  sv[1] = dv[1] - av[1];
  sv[2] = dv[2] - av[2];
  sv[3] = dv[3] - av[3];
  sv[4] = dv[4] - av[4];
  sv[5] = dv[5] - av[5];
  sv[6] = dv[6] - av[6];

  // unit joint direction
  ev[0] = {0, 0, 1};
  ev[1] = {0, 0, 1};
  ev[2] = {0, 0, 1};
  ev[3] = {0, 0, 1};
  ev[4] = {0, 0, 1};
  ev[5] = {0, 0, 1};
  ev[6] = {0, 0, 1};

  // kinematic and dynamic constraints
  // position (rad)
  pmax[0] = INF_REAL;
  pmax[1] = INF_REAL;
  pmax[2] = INF_REAL;
  pmax[3] = (real) 2.58;
  pmax[4] = INF_REAL;
  pmax[5] = (real) 2.1;
  pmax[6] = INF_REAL;

  // velocity (rad/s)
  vmax[0] = (real) 1.745;
  vmax[1] = (real) 1.745;
  vmax[2] = (real) 1.745;
  vmax[3] = (real) 1.745;
  vmax[4] = (real) 2.443;
  vmax[5] = (real) 2.443;
  vmax[6] = (real) 2.443;

  // torque (Nm)
  tau_max[0] = 52;
  tau_max[1] = 52;
  tau_max[2] = 52;
  tau_max[3] = 52;
  tau_max[4] = 17;
  tau_max[5] = 17;
  tau_max[6] = 17;

  is_init = true;
}

host_fn void cuGen3MultiTraj::compute_constraints(const Array& x, const Matrix& task, cuMultiBsplines& gpu_bsplines) {
  gpu_bsplines.compute_control(x, task);
  compute_constraints_kernel<<<ntrajectories, npoints>>>();
  cuda_check_kernel;
  cuda_check(cudaDeviceSynchronize());
}

host_fn void cuGen3MultiTraj::fetch_constraints() {
  cuda_check(cudaMemcpy(host_constraints, dev_constraints, 21 * npoints * sizeof(real), cudaMemcpyDeviceToHost));
}

host_fn cuGen3MultiTraj::~cuGen3MultiTraj() {
  if (!is_alias) {
    cuda_check(cudaFree(dev_constraints));
    dev_constraints = nullptr;
    cuda_check(cudaFreeHost(host_constraints));
    host_constraints = nullptr;
    is_init          = false;
    cuda_check(cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0));
  }
}

//------ GENERAL FUNCTIONS ------------------------------------------------------------------------------------

blast_fn u32 cuGen3MultiTraj::con_id(u32 point, u32 traj) {
  return 21 * (traj * npoints + point);
}

blast_fn u32 cuGen3MultiTraj::nconstraints() {
  return 21 * ntrajectories * npoints;
}

//------ DEVICE FUNCTIONS ------------------------------------------------------------------------------------

__global__ void compute_constraints_kernel() {
  cuMultiBsplines* bspline = (cuMultiBsplines*) bspline_broadcast_arena;
  Assert(!bspline->is_alias);
  cuGen3MultiTraj* manip = (cuGen3MultiTraj*) manip_broadcast_arena;
  Assert(manip->is_init);

  const auto point = threadIdx.x;
  const auto traj  = blockIdx.x;

  // todo: this function accesses global memory, check performance difference with local variables
  bspline->compute_trajectory_point(point, traj);
  manip->compute_constraints_point(bspline->dev_pos, bspline->dev_vel, bspline->dev_acc, &manip->dev_constraints[manip->con_id(point, traj)]);
}

dev_fn void cuGen3MultiTraj::compute_constraints_point(const real pos[7], const real vel[7], const real acc[7], real* con) {
#if BLAST_USE_DOUBLES
  real s[7];
  real c[7];
  ::sincos(pos[0], &s[0], &c[0]);
  ::sincos(pos[1], &s[1], &c[1]);
  ::sincos(pos[2], &s[2], &c[2]);
  ::sincos(pos[3], &s[3], &c[3]);
  ::sincos(pos[4], &s[4], &c[4]);
  ::sincos(pos[5], &s[5], &c[5]);
  ::sincos(pos[6], &s[6], &c[6]);
#else
  real s[7];
  real c[7];
  sincosf(pos[0], &s[0], &c[0]);
  sincosf(pos[1], &s[1], &c[1]);
  sincosf(pos[2], &s[2], &c[2]);
  sincosf(pos[3], &s[3], &c[3]);
  sincosf(pos[4], &s[4], &c[4]);
  sincosf(pos[5], &s[5], &c[5]);
  sincosf(pos[6], &s[6], &c[6]);
#endif

  const Mat3 Q1(c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1);
  const Mat3 Q2(c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0);
  const Mat3 Q3(c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0);
  const Mat3 Q4(c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0);
  const Mat3 Q5(c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0);
  const Mat3 Q6(c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0);
  const Mat3 Q7(c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0);
  const Mat3 Q1t(transpose(Q1));
  const Mat3 Q2t(transpose(Q2));
  const Mat3 Q3t(transpose(Q3));
  const Mat3 Q4t(transpose(Q4));
  const Mat3 Q5t(transpose(Q5));
  const Mat3 Q6t(transpose(Q6));
  const Mat3 Q7t(transpose(Q7));

  //-- Collision constraints
  Vec3 p_orig(0, 0, 0);
  Vec3 p_j2;
  Vec3 p_j3;
  Vec3 p_j4;
  Vec3 p_j6;
  Vec3 p_ee;
  auto p_tmp = p_base;
  auto Q_tmp = Q1;
  p_tmp += Q_tmp * dv[0];
  p_j2 = p_tmp;
  p_tmp += (Q_tmp *= Q2) * dv[1];
  p_j3 = p_tmp;
  p_tmp += (Q_tmp *= Q3) * dv[2];
  p_j4 = p_tmp;
  p_tmp += (Q_tmp *= Q4) * dv[3];
  p_tmp += (Q_tmp *= Q5) * dv[4];
  p_j6 = p_tmp;
  p_tmp += (Q_tmp *= Q6) * dv[5];
  p_tmp += (Q_tmp *= Q7) * dv[6];
  p_ee              = p_tmp;
  const real r1_sqr = (real) 0.09 * (real) 0.09;
  const real r2_sqr = (real) 0.09 * (real) 0.09;
  // Self collisions sqr
  real dist1sqr = two_segment_distance_sqr(p_orig, p_j2, p_j6, p_ee) - r1_sqr;
  real dist2sqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r2_sqr;
  // Collision with table sqr
  const real r4table = (real) 0.05;
  const real r6table = (real) 0.04;
  const Vec3 p_table(0, 0, -(real) 0.0025);
  real       distTJ4 = p_j4.z - p_table.z - r4table;
  real       distTJ6 = p_j6.z - p_table.z - r6table;
  real       distTEE = p_ee.z - p_table.z - r6table;

  con[0] = -dist1sqr;
  con[1] = -dist2sqr;
  con[2] = -distTJ4;
  con[3] = -distTJ6;
  con[4] = -distTEE;

  //-- position constraints
  con[5] = (abs(pos[3]) - pmax[3]) / pmax[3];
  con[6] = (abs(pos[5]) - pmax[5]) / pmax[5];

  //-- velocity constraints
  auto current_result = &con[7];
  for (u32 i = 0; i < 7; i++)
    current_result[i] = (abs(vel[i]) - vmax[i]) / vmax[i];
  current_result += 7;

  //-- Dynamic constraints
  // note: This is the Newton algorithm in 'Element de robotique' course notes.
  //       Careful because some variables are named differently and uses a slightly different conventions.
  //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
  //       joint turns with respect to the coordinate frame.
  Vec3 w1, w2, w3, w4, w5, w6, w7;
  Vec3 wd1, wd2, wd3, wd4, wd5, wd6, wd7;
  Vec3 cdd0 = {0, 0, (real) 9.81};
  Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6, cdd7;
  Vec3 f1, f2, f3, f4, f5, f6, f7;
  Vec3 n1, n2, n3, n4, n5, n6, n7;

  w1 = vel[0] * ev[0];
  w2 = Q2t * w1 + vel[1] * ev[1];
  w3 = Q3t * w2 + vel[2] * ev[2];
  w4 = Q4t * w3 + vel[3] * ev[3];
  w5 = Q5t * w4 + vel[4] * ev[4];
  w6 = Q6t * w5 + vel[5] * ev[5];
  w7 = Q7t * w6 + vel[6] * ev[6];

  wd1  = acc[0] * ev[0];
  cdd1 = Q1t * cdd0 + cross(wd1, av[0]) + cross(w1, cross(w1, av[0]));
  wd2  = Q2t * wd1 + acc[1] * ev[1] + vel[1] * cross(Q2t * w1, ev[1]);
  cdd2 = Q2t * cdd1 + cross(wd2, av[1]) + cross(w2, cross(w2, av[1])) - Q2t * cross(wd1, sv[0]) - Q2t * cross(w1, cross(w1, sv[0]));
  wd3  = Q3t * wd2 + acc[2] * ev[2] + vel[2] * cross(Q3t * w2, ev[2]);
  cdd3 = Q3t * cdd2 + cross(wd3, av[2]) + cross(w3, cross(w3, av[2])) - Q3t * cross(wd2, sv[1]) - Q3t * cross(w2, cross(w2, sv[1]));
  wd4  = Q4t * wd3 + acc[3] * ev[3] + vel[3] * cross(Q4t * w3, ev[3]);
  cdd4 = Q4t * cdd3 + cross(wd4, av[3]) + cross(w4, cross(w4, av[3])) - Q4t * cross(wd3, sv[2]) - Q4t * cross(w3, cross(w3, sv[2]));
  wd5  = Q5t * wd4 + acc[4] * ev[4] + vel[4] * cross(Q5t * w4, ev[4]);
  cdd5 = Q5t * cdd4 + cross(wd5, av[4]) + cross(w5, cross(w5, av[4])) - Q5t * cross(wd4, sv[3]) - Q5t * cross(w4, cross(w4, sv[3]));
  wd6  = Q6t * wd5 + acc[5] * ev[5] + vel[5] * cross(Q6t * w5, ev[5]);
  cdd6 = Q6t * cdd5 + cross(wd6, av[5]) + cross(w6, cross(w6, av[5])) - Q6t * cross(wd5, sv[4]) - Q6t * cross(w5, cross(w5, sv[4]));
  wd7  = Q7t * wd6 + acc[6] * ev[6] + vel[6] * cross(Q7t * w6, ev[6]);
  cdd7 = Q7t * cdd6 + cross(wd7, av[6]) + cross(w7, cross(w7, av[6])) - Q7t * cross(wd6, sv[5]) - Q7t * cross(w6, cross(w6, sv[5]));

  //-- dynamics
  f7 = m[6] * cdd7;
  n7 = I[6] * wd7 + cross(w7, I[6] * w7) + cross(av[6], f7);
  f6 = m[5] * cdd6 + Q7 * f7;
  n6 = I[5] * wd6 + cross(w6, I[5] * w6) + Q7 * n7 + cross(av[5], f6) + cross(sv[5], (Q7 * f7));
  f5 = m[4] * cdd5 + Q6 * f6;
  n5 = I[4] * wd5 + cross(w5, I[4] * w5) + Q6 * n6 + cross(av[4], f5) + cross(sv[4], (Q6 * f6));
  f4 = m[3] * cdd4 + Q5 * f5;
  n4 = I[3] * wd4 + cross(w4, I[3] * w4) + Q5 * n5 + cross(av[3], f4) + cross(sv[3], (Q5 * f5));
  f3 = m[2] * cdd3 + Q4 * f4;
  n3 = I[2] * wd3 + cross(w3, I[2] * w3) + Q4 * n4 + cross(av[2], f3) + cross(sv[2], (Q4 * f4));
  f2 = m[1] * cdd2 + Q3 * f3;
  n2 = I[1] * wd2 + cross(w2, I[1] * w2) + Q3 * n3 + cross(av[1], f2) + cross(sv[1], (Q3 * f3));
  f1 = m[0] * cdd1 + Q2 * f2;
  n1 = I[0] * wd1 + cross(w1, I[0] * w1) + Q2 * n2 + cross(av[0], f1) + cross(sv[0], (Q2 * f2));

  //-- extract torques (last element of each moment vector)
  current_result[0] = (abs(n1.z) - tau_max[0]) / tau_max[0];
  current_result[1] = (abs(n2.z) - tau_max[1]) / tau_max[1];
  current_result[2] = (abs(n3.z) - tau_max[2]) / tau_max[2];
  current_result[3] = (abs(n4.z) - tau_max[3]) / tau_max[3];
  current_result[4] = (abs(n5.z) - tau_max[4]) / tau_max[4];
  current_result[5] = (abs(n6.z) - tau_max[5]) / tau_max[5];
  current_result[6] = (abs(n7.z) - tau_max[6]) / tau_max[6];
}

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
TEST_CASE("GpuCpuManipCorrectness", "[Manipulator]") {
  using namespace blast;

  const u32 points        = 256;
  const u32 joints        = 7;
  const u32 p             = 5;
  const u32 ncontrol      = 24;
  const u32 ntrajectories = 80;

  Bspline         host_bspline(ncontrol, points, p, joints);
  cuMultiBsplines gpu_bsplines(points, joints, p, ncontrol, ntrajectories);

  Gen3_7DOF host_manip;
  host_manip.set_payload(4);
  cuGen3MultiTraj gpu_manip(points, ntrajectories, 4);

  // random task
  Matrix task(joints, 6);
  real   amp = 10;
  for (u32 i = 0; i < task.rows; i++)
    for (u32 j = 0; j < task.cols; j++)
      task(i, j) = amp * random_real();

  // random optimization vector
  Array x(joints * (ncontrol - 6) + 1);
  for (u32 i = 0; i < x.size; i++)
    x[i] = amp * random_real();
  x.back() = abs(x.back());

  // copy the 'x' Array 'ntrajectories' times for the gpu version
  const u32 xlen = host_bspline.xlen(task);
  Array     x_multi(ntrajectories * host_bspline.xlen(task));
  for (int i = 0; i < ntrajectories; i++) {
    for (int j = 0; j < xlen; j++) {
      x_multi[i * xlen + j] = x[j];
    }
  }

  host_bspline.compute_trajectory(x, task);
  auto host_constraints = host_manip.constraints(host_bspline.traj);

  gpu_manip.compute_constraints(x_multi, task, gpu_bsplines);
  gpu_manip.fetch_constraints();

  // Test correctness of the constraints
  for (int traj_id = 0; traj_id < ntrajectories; traj_id++) {
    for (int i = 0; i < (int) host_manip.ncon(points); i++) {
      const auto id = traj_id * host_manip.ncon(points) + i;
      REQUIRE(abs(host_constraints[i] - gpu_manip.host_constraints[id]) < 0.01);
    }
  }
}


#endif // BLAST_ENABLE_TESTS
