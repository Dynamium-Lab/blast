#pragma once
#include <blast>
#include <math_constants.h>
#include "gpu/trajectory.cuh"

namespace blast {
__global__ void compute_constraints_kernel();
__global__ void compute_constraints_kernel_ur5e();

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
  // note: blast_fn (not device_fn) so tests can call it directly on the host
  //       to validate compute_constraints_kernel's thread/block dispatch and
  //       __constant__ broadcast wiring, independent from the physics formula.
  blast_fn void compute_constraints_point(const real pos[7], const real vel[7], const real acc[7], real* con);
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
  CUDA_CHECK(cudaMalloc(&dev_constraints, nconstraints() * sizeof(real)));
  CUDA_CHECK(cudaMallocHost(&host_constraints, nconstraints() * sizeof(real)));

  this->_init(mass);

  // note: this must be done last
  CUDA_CHECK(cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0));
}

host_fn void cuGen3MultiTraj::_init(real mass) {

  // position of the first joint with respect to the table in the center of the base
  p_base = {0, 0, 0.1564f};

  // mass
  m[0] = 1.377f;
  m[1] = 1.1636f;
  m[2] = 1.1636f;
  m[3] = 0.93f;
  m[4] = 0.678f;
  m[5] = 0.678f;
  m[6] = 0.364f + 0.921f;

  // inertial tensor
  I[0] = {0.004570f, 0.000001f, 0.000002f, 0.000001f, 0.004831f, 0.000448f, 0.000002f, 0.000448f, 0.001409f};
  I[1] = {0.011088f, 0.000005f, 0.000000f, 0.000005f, 0.001072f, -0.000691f, 0.000000f, -0.000691f, 0.011255f};
  I[2] = {0.010932f, 0.000000f, -0.000007f, 0.000000f, 0.011127f, 0.000606f, -0.000007f, 0.000606f, 0.001043f};
  I[3] = {0.008147f, -0.000001f, 0.000000f, -0.000001f, 0.000631f, -0.000500f, 0.000000f, -0.000500f, 0.008316f};
  I[4] = {0.001596f, 0.000000f, 0.000000f, 0.000000f, 0.001607f, 0.000256f, 0.000000f, 0.000256f, 0.000399f};
  I[5] = {0.001641f, 0.000000f, 0.000000f, 0.000000f, 0.000410f, -0.000278f, 0.000000f, -0.000278f, 0.001641f};
  I[6] = {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, -0.000002f, 0.000001f, -0.000002f, 0.000240f};

  // center of mass
  av[0] = {-0.000023f, -0.010364f, -0.073360f};
  av[1] = {-0.000044f, -0.099580f, -0.013278f};
  av[2] = {-0.000044f, -0.006641f, -0.117892f};
  av[3] = {-0.000018f, -0.075478f, -0.015006f};
  av[4] = {0.000001f, -0.009432f, -0.063883f};
  av[5] = {0.000001f, -0.045483f, -0.009650f};
  av[6] = {-0.000093f, 0.000132f, -0.022905f};
  Vec3 av_tool(0, 0, -0.06f - 0.0615f);
  av[6] = (0.364f * av[6] + 0.921f * av_tool) * (1 / (0.364f + 0.921f));

  // vector to next joint
  dv[0] = {0.0f, 0.0054f, -0.1284f};
  dv[1] = {0.0f, -0.2104f, -0.0064f};
  dv[2] = {0.0f, -0.0064f, -0.2104f};
  dv[3] = {0.0f, -0.2084f, -0.0064f};
  dv[4] = {0.0f, 0.0f, -0.1059f};
  dv[5] = {0.0f, -0.1059f, 0.0f};
  dv[6] = {0.0f, 0.0f, -0.0615f - 0.164f};

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
  pmax[3] = 2.58f;
  pmax[4] = INF_REAL;
  pmax[5] = 2.1f;
  pmax[6] = INF_REAL;

  // velocity (rad/s)
  vmax[0] = 1.745f;
  vmax[1] = 1.745f;
  vmax[2] = 1.745f;
  vmax[3] = 1.745f;
  vmax[4] = 2.443f;
  vmax[5] = 2.443f;
  vmax[6] = 2.443f;

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
  CUDA_CHECK_KERNEL;
  CUDA_CHECK(cudaDeviceSynchronize());
}

host_fn void cuGen3MultiTraj::fetch_constraints() {
  CUDA_CHECK(cudaMemcpy(host_constraints, dev_constraints, nconstraints() * sizeof(real), cudaMemcpyDeviceToHost));
}

host_fn cuGen3MultiTraj::~cuGen3MultiTraj() {
  if (!is_alias) {
    CUDA_CHECK(cudaFree(dev_constraints));
    dev_constraints = nullptr;
    CUDA_CHECK(cudaFreeHost(host_constraints));
    host_constraints = nullptr;
    is_init          = false;
    CUDA_CHECK(cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0));
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
  const auto pva_id = traj * bspline->points * bspline->joints + point * bspline->joints;
  manip->compute_constraints_point(&bspline->dev_pos[pva_id], &bspline->dev_vel[pva_id], &bspline->dev_acc[pva_id], &manip->dev_constraints[manip->con_id(point, traj)]);
}

// ============================================================================
// UR5e (6-DOF) variant. Structurally identical to cuGen3MultiTraj -- same
// broadcast-arena / kernel-launch pattern -- but checks position, velocity,
// and ACCELERATION (18 = 6+6+6 constraints), not torque. UR5e has a real,
// finite acceleration limit (unlike Gen3, which has none -- every joint's
// acceleration_max is INF_REAL in Kinova_Gen3.hpp, which is why Gen3 uses
// torque as its third constraint instead), and this matches blast's own
// canonical UR5e usage (example_02_trajectory_optimization.cpp uses
// position/velocity/acceleration for UR5e too).
//
// This deliberately has no dynamics: acceleration constraints are a direct
// check against the optimizer's own acc[] input (no Newton-Euler recursion,
// no per-joint rotation matrices, no mass/inertia/kinematic-offset data
// needed at all) -- avoiding the hand-derivation risk a torque-based port
// would carry, since UR5e's per-joint rotation matrices would need to be
// independently derived from its kinematics the same way Gen3's were.
//
// Collision is out of scope for the same reason as Gen3: not part of the
// position/velocity/acceleration PSO-vs-SLSQP comparison.
struct cuUR5eMultiTraj {
  real pmax[6];
  real vmax[6];
  real amax[6];

  u32  ntrajectories = 0;
  u32  npoints       = 0;
  bool is_init       = false;
  bool is_alias      = false;

  real* dev_constraints;
  real* host_constraints;

  // HOST
  //----------

  host_fn cuUR5eMultiTraj(u32 points, u32 trajs = 1, real mass = 0.0);

  host_fn cuUR5eMultiTraj(const cuUR5eMultiTraj& other) {
    *this    = other;
    is_alias = true;
  }

  host_fn void compute_constraints(const Array& x, const Matrix& task, cuMultiBsplines& gpu_bsplines);
  host_fn void fetch_constraints();
  host_fn void _init(real mass);
  host_fn ~cuUR5eMultiTraj();

  // GENERAL
  //----------
  blast_fn u32 con_id(u32 point, u32 traj);
  blast_fn u32 nconstraints();

  // DEVICE
  //----------
  blast_fn void compute_constraints_point(const real pos[6], const real vel[6], const real acc[6], real* con);
};

const int         MANIP_UR5E_BROADCAST_ARENA_SIZE = 256;
__constant__ real manip_ur5e_broadcast_arena[MANIP_UR5E_BROADCAST_ARENA_SIZE];

//------ HOST FUNCTIONS ------------------------------------------------------------------------------------

host_fn cuUR5eMultiTraj::cuUR5eMultiTraj(u32 _points, u32 _trajs, real mass) :
    npoints(_points),
    ntrajectories(_trajs) {

  assert_buffer_size<cuUR5eMultiTraj, sizeof(manip_ur5e_broadcast_arena)>();
  CUDA_CHECK(cudaMalloc(&dev_constraints, nconstraints() * sizeof(real)));
  CUDA_CHECK(cudaMallocHost(&host_constraints, nconstraints() * sizeof(real)));

  this->_init(mass);

  CUDA_CHECK(cudaMemcpyToSymbol(manip_ur5e_broadcast_arena, this, sizeof(*this), 0));
}

host_fn void cuUR5eMultiTraj::_init(real mass) {
  // position (rad) -- all 6 finite for UR5e (blast/manipulator/UR5e.hpp)
  pmax[0] = 6.283200f;
  pmax[1] = 6.283200f;
  pmax[2] = 3.141600f;
  pmax[3] = 6.283200f;
  pmax[4] = 6.283200f;
  pmax[5] = 6.283200f;

  // velocity (rad/s)
  vmax[0] = 3.141600f;
  vmax[1] = 3.141600f;
  vmax[2] = 3.141600f;
  vmax[3] = 3.141600f;
  vmax[4] = 3.141600f;
  vmax[5] = 3.141600f;

  // acceleration (rad/s^2)
  amax[0] = 13.96f;
  amax[1] = 13.96f;
  amax[2] = 13.96f;
  amax[3] = 13.96f;
  amax[4] = 13.96f;
  amax[5] = 13.96f;

  is_init = true;
}

host_fn void cuUR5eMultiTraj::compute_constraints(const Array& x, const Matrix& task, cuMultiBsplines& gpu_bsplines) {
  gpu_bsplines.compute_control(x, task);
  compute_constraints_kernel_ur5e<<<ntrajectories, npoints>>>();
  CUDA_CHECK_KERNEL;
  CUDA_CHECK(cudaDeviceSynchronize());
}

host_fn void cuUR5eMultiTraj::fetch_constraints() {
  CUDA_CHECK(cudaMemcpy(host_constraints, dev_constraints, nconstraints() * sizeof(real), cudaMemcpyDeviceToHost));
}

host_fn cuUR5eMultiTraj::~cuUR5eMultiTraj() {
  if (!is_alias) {
    CUDA_CHECK(cudaFree(dev_constraints));
    dev_constraints = nullptr;
    CUDA_CHECK(cudaFreeHost(host_constraints));
    host_constraints = nullptr;
    is_init          = false;
    CUDA_CHECK(cudaMemcpyToSymbol(manip_ur5e_broadcast_arena, this, sizeof(*this), 0));
  }
}

//------ GENERAL FUNCTIONS ------------------------------------------------------------------------------------

blast_fn u32 cuUR5eMultiTraj::con_id(u32 point, u32 traj) {
  return 18 * (traj * npoints + point);
}

blast_fn u32 cuUR5eMultiTraj::nconstraints() {
  return 18 * ntrajectories * npoints;
}

//------ DEVICE FUNCTIONS ------------------------------------------------------------------------------------

blast_fn void cuUR5eMultiTraj::compute_constraints_point(const real pos[6], const real vel[6], const real acc[6], real* con) {
  //-- position constraints
  for (u32 i = 0; i < 6; i++)
    con[i] = (::abs(pos[i]) - pmax[i]) / pmax[i];

  //-- velocity constraints
  for (u32 i = 0; i < 6; i++)
    con[6 + i] = (::abs(vel[i]) - vmax[i]) / vmax[i];

  //-- acceleration constraints
  for (u32 i = 0; i < 6; i++)
    con[12 + i] = (::abs(acc[i]) - amax[i]) / amax[i];
}

__global__ void compute_constraints_kernel_ur5e() {
  cuMultiBsplines* bspline = (cuMultiBsplines*) bspline_broadcast_arena;
  Assert(!bspline->is_alias);
  cuUR5eMultiTraj* manip = (cuUR5eMultiTraj*) manip_ur5e_broadcast_arena;
  Assert(manip->is_init);

  const auto point = threadIdx.x;
  const auto traj  = blockIdx.x;

  bspline->compute_trajectory_point(point, traj);
  const auto pva_id = traj * bspline->points * bspline->joints + point * bspline->joints;
  manip->compute_constraints_point(&bspline->dev_pos[pva_id], &bspline->dev_vel[pva_id], &bspline->dev_acc[pva_id], &manip->dev_constraints[manip->con_id(point, traj)]);
}

blast_fn void cuGen3MultiTraj::compute_constraints_point(const real pos[7], const real vel[7], const real acc[7], real* con) {
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
  const real r1_sqr = 0.09f * 0.09f;
  const real r2_sqr = 0.09f * 0.09f;
  // Self collisions sqr
  real dist1sqr = two_segment_distance_sqr(p_orig, p_j2, p_j6, p_ee) - r1_sqr;
  real dist2sqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r2_sqr;
  // Collision with table sqr
  const real r4table = 0.05f;
  const real r6table = 0.04f;
  const Vec3 p_table(0, 0, -0.0025f);
  real       distTJ4 = p_j4.z - p_table.z - r4table;
  real       distTJ6 = p_j6.z - p_table.z - r6table;
  real       distTEE = p_ee.z - p_table.z - r6table;

  con[0] = -dist1sqr;
  con[1] = -dist2sqr;
  con[2] = -distTJ4;
  con[3] = -distTJ6;
  con[4] = -distTEE;

  //-- position constraints
  con[5] = (::abs(pos[3]) - pmax[3]) / pmax[3];
  con[6] = (::abs(pos[5]) - pmax[5]) / pmax[5];

  //-- velocity constraints
  auto current_result = &con[7];
  for (u32 i = 0; i < 7; i++)
    current_result[i] = (::abs(vel[i]) - vmax[i]) / vmax[i];
  current_result += 7;

  //-- Dynamic constraints
  // note: This is the Newton algorithm in 'Element de robotique' course notes.
  //       Careful because some variables are named differently and uses a slightly different conventions.
  //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
  //       joint turns with respect to the coordinate frame.
  Vec3 w1, w2, w3, w4, w5, w6, w7;
  Vec3 wd1, wd2, wd3, wd4, wd5, wd6, wd7;
  Vec3 cdd0 = {0, 0, 9.81f};
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
  current_result[0] = (::abs(n1.z) - tau_max[0]) / tau_max[0];
  current_result[1] = (::abs(n2.z) - tau_max[1]) / tau_max[1];
  current_result[2] = (::abs(n3.z) - tau_max[2]) / tau_max[2];
  current_result[3] = (::abs(n4.z) - tau_max[3]) / tau_max[3];
  current_result[4] = (::abs(n5.z) - tau_max[4]) / tau_max[4];
  current_result[5] = (::abs(n6.z) - tau_max[5]) / tau_max[5];
  current_result[6] = (::abs(n7.z) - tau_max[6]) / tau_max[6];
}

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
// note: the previous version of this test used a `Gen3` CPU class
// (.ncon()/.internal_constraints()/.set_payload()) that no longer exists —
// the manipulator API was refactored to the generic Manipulator struct. That
// left this test uncompilable, and it wasn't even wired into any build
// target, so nothing caught it.
//
// This version instead calls cuGen3MultiTraj::compute_constraints_point()
// (now blast_fn, host+device callable — see its declaration above) directly
// on the host, using the exact same trajectory data the GPU kernel consumed.
// Since it's the same function running on both sides, this does not
// independently re-derive the physics formula; what it does validate is the
// GPU-port-specific risk: that compute_constraints_kernel's thread/block
// dispatch, __constant__ broadcast memory (bspline_broadcast_arena /
// manip_broadcast_arena), and con_id() indexing correctly route each
// point/trajectory to the right output slot. That's exactly the class of bug
// found and fixed here (fetch_constraints() previously copied back only
// trajectory 0's data).
TEST_CASE("GpuCpuManipCorrectness", "[Manipulator]") {
  using namespace blast;

  const u32 points        = 256;
  const u32 joints        = 7;
  const u32 p             = 5;
  const u32 ncontrol      = 24;
  const u32 ntrajectories = 8;

  cuMultiBsplines gpu_bsplines(points, joints, p, ncontrol, ntrajectories);
  cuGen3MultiTraj gpu_manip(points, ntrajectories, 4);

  // random task
  Matrix task(joints, 6);
  real   amp = 10;
  for (u32 i = 0; i < task.rows; i++)
    for (u32 j = 0; j < task.cols; j++)
      task(i, j) = amp * random_real();

  // distinct random optimization vector per trajectory (catches
  // trajectory-indexing bugs that identical inputs across trajectories would
  // mask)
  const u32 xlen = gpu_bsplines.host->x_len(task);
  Array     x_multi(ntrajectories * xlen);
  for (u32 i = 0; i < ntrajectories; i++) {
    for (u32 j = 0; j < xlen; j++)
      x_multi[i * xlen + j] = amp * random_real();
    x_multi[i * xlen + xlen - 1] = abs(x_multi[i * xlen + xlen - 1]) + 0.1; // positive duration
  }

  gpu_manip.compute_constraints(x_multi, task, gpu_bsplines);
  gpu_bsplines.fetch_trajectories();
  gpu_manip.fetch_constraints();

  for (u32 traj_id = 0; traj_id < ntrajectories; traj_id++) {
    for (u32 point = 0; point < points; point++) {
      real pos[7], vel[7], acc[7];
      for (u32 j = 0; j < joints; j++) {
        const auto id = traj_id * points * joints + point * joints + j;
        pos[j]        = gpu_bsplines.pos[id];
        vel[j]        = gpu_bsplines.vel[id];
        acc[j]        = gpu_bsplines.acc[id];
      }

      real ref_con[21];
      gpu_manip.compute_constraints_point(pos, vel, acc, ref_con);

      const auto base = gpu_manip.con_id(point, traj_id);
      for (u32 k = 0; k < 21; k++) {
        CAPTURE(traj_id, point, k);
        REQUIRE(abs(ref_con[k] - gpu_manip.host_constraints[base + k]) < 1e-3);
      }
    }
  }
}


#endif // BLAST_ENABLE_TESTS
