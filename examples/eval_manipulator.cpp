//
// Created by nikos on 2025-09-04.
//
#define BLAST_TRACE_LEVEL 0
#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"

#include <fcntl.h> // _O_WRONLY
#include <io.h>    // _open, _dup, _dup2, _close

#define WIN32_LEAN_AND_MEAN
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00 // Windows 10+
#endif
#include <algorithm>
#include <mutex>
#include <thread>
#include <unordered_set>
#include <windows.h>

#include <cassert>
#include <fstream>
#include <stdexcept>


using namespace blast;

constexpr int thread_count = 8;
// ----------------------------- UR5e manip + tasks -----------------------------
inline Manipulator get_generic_ur5e() {
  // Manipulator
  blast::u32 joints = 6;

  // limits
  blast::ManipulatorLimits limits;
  limits.pmax = {6.283200, 6.283200, 3.141600, 6.283200, 6.283200, 6.283200};
  limits.pmin = {-6.283200, -6.283200, -3.141600, -6.283200, -6.283200, -6.283200};
  limits.vmax = {3.141600, 3.141600, 3.141600, 3.141600, 3.141600, 3.141600};
  // limits.vmin = {-3.141600, -3.141600, -3.141600, -3.141600, -3.141600, -3.141600};
  limits.amax = {13.96, 13.96, 13.96, 13.96, 13.96, 13.96};
  // limits.amin = {-13.96, -13.96, -13.96, -13.96, -13.96, -13.96};
  limits.tau_max = {150.000000, 150.000000, 150.000000, 28.000000, 28.000000, 28.000000};
  // limits.tau_min = {-150.000000, -150.000000, -150.000000, -28.000000, -28.000000, -28.000000};
  limits.tcp_max = 1.5; // todo: verify this value

  blast::ManipulatorKinematics kinematics;
  kinematics.dv = {
          Vec3{0, 0, 0},
          {-0.425, 0, 0},
          {-0.3922, 0, 0.1333},
          {0, -0.0997, -0},
          {0, 0.0996, -0},
          {0, 0, 0} // to end effector
  };
  kinematics.ev = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}};
  // kinematics.Q_static.resize(6);
  kinematics.Q_static[0] = {-1.000000, 0.000000, 0.000000, -0.000000, -1.000000, 0.000000, 0.000000, -0.000000, 1.000000}; // modified
  kinematics.Q_static[1] = {1.000000, 0.000000, 0.000000, -0.000000, -0.000000, 1.000000, 0.000000, -1.000000, -0.000000};
  kinematics.Q_static[2] = {1.000000, 0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000000, 1.000000};
  kinematics.Q_static[3] = {1.000000, 0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000000, 1.000000};
  kinematics.Q_static[4] = {1.000000, 0.000000, 0.000000, -0.000000, -0.000000, 1.000000, 0.000000, -1.000000, -0.000000};
  kinematics.Q_static[5] = {1.000000, -0.000000, 0.000000, 0.000000, -0.000000, -1.000000, 0.000000, 1.000000, -0.000000};
  // kinematics.Q_static[6] = {1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000};
  kinematics.p_j0 = {0.000000, 0.000000, 0.162500};

  blast::ManipulatorDynamics dynamics;
  dynamics.m = {3.700000, 8.393000, 2.275000, 1.219000, 1.219000, 0.187900};
  dynamics.I = {
          Mat3{0.0103, 0, 0, 0, 0.0103, 0, 0, 0, 0.0067},
          {0.1339, 0, 0, 0, 0.1339, 0, 0, 0, 0.0151},
          {0.0312, 0, 0, 0, 0.0312, 0, 0, 0, 0.0041},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0001, 0, 0, 0, 0.0001, 0, 0, 0, 0.0001}};
  dynamics.av = {
          Vec3{0, 0, 0},
          {-0.2125, 0, 0.138},
          {-0.1961, 0, 0.007},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, -0.0229}};

  blast::ManipulatorCapsules collisions;

  blast::Sphere base;
  base.c                    = {0, 0, 0.0325};
  base.r                    = 0.09;
  collisions.base_sphere    = base;
  collisions.collision_base = {0, 0, 0, 1, 1, 1, 1};

  collisions.collision_matrix.resize(7, 7);
  collisions.collision_matrix(3, 0) = 1;
  collisions.collision_matrix(4, 0) = 1;
  collisions.collision_matrix(5, 0) = 1;
  collisions.collision_matrix(6, 0) = 1;

  collisions.collision_matrix(4, 1) = 1;
  collisions.collision_matrix(5, 1) = 1;
  collisions.collision_matrix(6, 1) = 1;

  collisions.collision_matrix(6, 3) = 1;

  // mirrored just for consistency
  // collisions.collision_matrix(0, 3) = 1;
  // collisions.collision_matrix(0, 4) = 1;
  // collisions.collision_matrix(0, 5) = 1;
  // collisions.collision_matrix(0, 6) = 1;
  //
  // collisions.collision_matrix(1, 4) = 1;
  // collisions.collision_matrix(1, 5) = 1;
  // collisions.collision_matrix(1, 6) = 1;
  //
  // collisions.collision_matrix(3, 6) = 1;

  blast::CollisionModelCapsule capsule;
  capsule.joint_frame = 0;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, -0.15, 0};
  capsule.r           = 0.09;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 1;
  capsule.p1          = {-0.42, 0, 0.1375};
  capsule.p2          = {0, 0, 0.1375};
  capsule.r           = 0.06;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {0, 0, 0.02};
  capsule.p2          = {0, 0, 0.18};
  capsule.r           = 0.065;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {-0.373156, 0, 0.00850418};
  capsule.p2          = {0.000440611, 0.000121443, 0.00850418};
  capsule.r           = 0.05;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, 0, -0.155};
  capsule.r           = 0.0425;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, 0.03, 0};
  capsule.p2          = {0, -0.1, 0};
  capsule.r           = 0.0425;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 5;
  capsule.p1          = {0, 0, -0.14};
  capsule.p2          = {0, 0, -0.01};
  capsule.r           = 0.038;
  collisions.capsule_list.push_back(capsule);

  blast::Manipulator generic_manip(joints, limits, kinematics, &dynamics, &collisions);
  return generic_manip;
}

inline std::vector<Matrix> get_UR5e_tasks() {
  int                n_joints = 6;
  std::vector<Array> start_pos;
  start_pos.push_back({1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883});
  start_pos.push_back({2.4316, 0.0965, 1.1824, 1.8626, -0.8567, 0.0250});
  start_pos.push_back({-0.704868, 2.17507, 0.834305, 0.132228, 2.27799, 0.00328674});
  start_pos.push_back({-0.819005, 2.50686, 0.344732, -2.85083, -2.39379, -2.93992});
  start_pos.push_back({2.04549, 0.675615, 0.137326, -0.814844, 0.449157, -3.08683});
  start_pos.push_back({2.3499, 1.44305, -1.11122, -0.333013, 0.751039, -3.09482});
  start_pos.push_back({2.36232, 1.06656, -0.132829, -0.933749, 0.788944, -3.13519});
  start_pos.push_back({2.22452, 1.15771, -0.467326, -0.690436, 0.641695, -3.12205});

  std::vector<Array> end_pos;
  end_pos.push_back({2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328});
  end_pos.push_back({0.1706, -2.5340, -1.2649, -2.4857, -1.4118, 2.9011});
  end_pos.push_back({0.3572, 2.9674, 0.4629, 2.8527, -1.2170, 3.1408});
  end_pos.push_back({-0.0436368, -2.65335, 0.593745, 2.05947, -1.61882, -3.12196});
  end_pos.push_back({-2.6558, -0.6703, 0.0599, -2.5311, -2.0486, -0.0127});
  end_pos.push_back({3.12394, -0.553571, -0.509807, 1.06314, 1.52902, -3.10058});
  end_pos.push_back({2.4552, -0.1792, -1.0268, -1.9354, -0.8723, -0.0198});
  end_pos.push_back({2.14686, -0.283198, -0.264971, -2.59327, -0.566506, -0.0169222});

  std::vector<Matrix> task_list(64);
  int                 task_idx = 0;
  for (auto& s: start_pos) {
    for (auto& e: end_pos) {
      task_list[task_idx].resize(n_joints, 6);
      for (int k = 0; k < n_joints; k++) {
        task_list[task_idx](k, 0) = s[k];
        task_list[task_idx](k, 3) = e[k];
      }
      task_idx++;
    }
  }
  return task_list;
}

inline World get_kitchen_open_doors() {
  World world;
  world.add_box(Vec3{0.427, -1.0847, 1.1974}, Vec3{0.36029, 0.06, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, 0.48791, 1.1974}, Vec3{0.36029, 0.06, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 1.5298}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 1.2074}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 0.86498}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.72729, -0.29838, 1.2074}, Vec3{0.06, 0.84629, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{-0.20911, 0.86968, 1.2009}, Vec3{0.032523, 0.42256, 0.36242}, Mat3{-0.70706, -0.70716, 0, -0.70716, -0.70706, 0, 0, 0, 1});
  world.add_box(Vec3{-0.20549, -1.4701, 1.2009}, Vec3{0.037648, 0.42257, 0.36244}, Mat3{-0.70706, 0.70716, 0, 0.70716, -0.70706, 0, 0, 0, 1});
  world.add_box(Vec3{-0.47232, -0.77127, -0.48346}, Vec3{0.023735, 0.41378, 0.41552}, Mat3{4.6327e-05, 0, 1, 0, 1, 0, 1, 0, 4.6327e-05});
  world.add_box(Vec3{0.35698, 0.425, -0.066059}, Vec3{0.41378, 0.8425, 0.43736}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -1.1546, -0.066059}, Vec3{0.41378, 0.03, 0.43736}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, 0.31136}, Vec3{0.41378, 0.38356, 0.06}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, -0.055}, Vec3{0.41378, 0.38356, 0.025}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, -0.47336}, Vec3{0.41378, 0.38356, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  return world;
}

// ----------------------------- Gen3 manip + tasks -----------------------------
inline std::vector<Matrix> get_Gen3_paper_tasks() {

  int                n_joints = 7;
  std::vector<Array> start_pos;
  start_pos.push_back({-1.2215, 1.9006, -2.7479, -1.2570, -1.3301, 1.2168, 1.2276});
  start_pos.push_back({-1.2111, 1.0074, 1.5468, 0.7477, 2.8742, -0.5686, 0.8434});
  start_pos.push_back({-1.29472, 1.21316, -1.87968, -1.33023, -1.50246, 0.207459, 1.93951});
  start_pos.push_back({-0.940587, 0.454046, -2.04557, -1.42284, -2.42939, -0.148302, -2.68966});
  start_pos.push_back({-0.322615, 1.99658, -0.690786, -0.5912, -0.0424388, 0.061149, -0.910936});
  start_pos.push_back({-0.448418, 0.321021, -0.730568, 0.474924, -1.66265, -1.24559, 0.673359});
  start_pos.push_back({-0.839155, 0.567741, -1.33936, 0.149858, 2.3022, 1.30233, -1.94128});
  start_pos.push_back({2.4646, -1.3391, -1.6020, 1.2399, -2.8390, 0.5598, -0.0744});
  start_pos.push_back({-0.9387, 1.9515, 2.1091, 0.5064, -0.4756, 0.5181, 2.9568});
  start_pos.push_back({-1.2804, 1.9185, -2.5408, -1.6593, 1.8834, -1.0995, -2.0906});

  std::vector<Array> end_pos;
  end_pos.push_back({-2.2485, -0.2788, -1.0306, -0.7692, 0.1985, -0.6363, 1.6482});
  end_pos.push_back({0.316206, 0.429902, -2.60506, -1.29654, 1.26543, 0.836061, -0.0189057});
  end_pos.push_back({0.574788, 1.94127, -1.95753, 0.656255, 1.67537, -0.079455, -1.14421});
  end_pos.push_back({-2.29169, -0.592198, 2.25955, -0.10461, -0.0522387, 1.29476, -0.870656});
  end_pos.push_back({-2.0435, -0.5774, -1.2565, -1.3726, 2.1209, -0.1030, -0.0170});
  end_pos.push_back({1.1614, 1.3213, -1.4704, 1.2843, -0.0907, -0.1372, -0.1260});
  end_pos.push_back({0.6055, 2.0039, 0.4019, -0.4372, 1.1506, -0.4391, -3.0771});
  end_pos.push_back({-1.4917, -2.2247, 1.8558, 1.3672, 2.3899, -0.3546, -1.7443});
  end_pos.push_back({-1.9901, -1.8203, -0.2981, -1.6681, -1.9974, -1.4377, -0.9739});
  end_pos.push_back({-2.6543, -1.3441, -1.0625, -0.5575, -2.8525, -0.0867, -0.7317});

  std::vector<Matrix> task_list(100);
  int                 task_idx = 0;
  for (auto& s: start_pos) {
    for (auto& e: end_pos) {
      task_list[task_idx].resize(n_joints, 6);
      for (int k = 0; k < n_joints; k++) {
        task_list[task_idx](k, 0) = s[k];
        task_list[task_idx](k, 3) = e[k];
      }
      task_idx++;
    }
  }
  return task_list;
}

inline World get_bookshelf_thin() {
  World world;
  // blast::add_static_capsule({1.1, 0.25, 0.37}, {1.1, 0.25, 0.51},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({0.8, 0.25, 0.37}, {0.8, 0.25, 0.51},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({1.1, -0.25, 0.62}, {1.1, -0.25, 0.76},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({0.8, -0.25, 0.62}, {0.8, -0.25, 0.76},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({1.1, 0.25, 0.87}, {1.1, 0.25, 1.01},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({0.8, 0.25, 0.87}, {0.8, 0.25, 1.01},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({1.1, -0.25, 1.12}, {1.1, -0.25, 1.26},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({0.8, -0.25, 1.12}, {0.8, -0.25, 1.26},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({1.1, 0.25, 1.37}, {1.1, 0.25, 1.51},  0.03, &bookshelf_thin);
  // blast::add_static_capsule({0.8, 0.25, 1.37}, {0.8, 0.25, 1.51},  0.03, &bookshelf_thin);
  world.add_box({1.0, 0.0, 0.35}, {0.31, 0.5, 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1.0, 0, 0.6}, {0.31, 0.5, 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1.0, 0, 0.85}, {0.31, 0.5, 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1, 0, 1.1}, {0.31, 0.5, 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1, 0, 1.35}, {0.31, 0.5, 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1, 0, 0.9}, {0.31, 0.02, 0.7}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1, 0, 1.6}, {0.31, 0.5, 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({0.7, 0.5, 0.82}, {0.02, 0.02, 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({0.7, -0.5, 0.82}, {0.02, 0.02, 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1.3, 0.5, 0.82}, {0.02, 0.02, 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  world.add_box({1.3, -0.5, 0.82}, {0.02, 0.02, 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  return world;
}

// ------------------------------ Helper functions ------------------------------
struct IOSilencer {
  int             saved_stdout_fd = -1;
  int             null_fd         = -1;
  std::streambuf* old_cout        = nullptr;
  std::streambuf* old_cerr        = nullptr;
  std::ofstream   null_stream;

  IOSilencer() {
    // Flush buffers
    std::cout.flush();
    std::cerr.flush();
    fflush(stdout);

    // Redirect std::cout and std::cerr to /dev/null
    null_stream.open("NUL"); // Windows null device
    old_cout = std::cout.rdbuf();
    old_cerr = std::cerr.rdbuf();
    std::cout.rdbuf(null_stream.rdbuf());
    std::cerr.rdbuf(null_stream.rdbuf());

    // Redirect printf (C stdout)
    saved_stdout_fd = _dup(_fileno(stdout));
    null_fd         = _open("NUL", _O_WRONLY);
    _dup2(null_fd, _fileno(stdout));
  }

  ~IOSilencer() {
    std::cout.rdbuf(old_cout);
    std::cerr.rdbuf(old_cerr);
    fflush(stdout);

    _dup2(saved_stdout_fd, _fileno(stdout));
    _close(null_fd);
    _close(saved_stdout_fd);
  }
};

struct CpuRow {
  ULONG id;    // CpuSet.Id (stable handle for SetThreadSelectedCpuSets)
  WORD  group; // processor group
  ULONG lpi;   // LogicalProcessorIndex (info only)
  ULONG core;  // CoreIndex (same for SMT siblings)
  BYTE  eff;   // EfficiencyClass (0 = most performant if hybrid)
  bool  ok;    // Allocated & !Parked
};

static std::vector<GROUP_AFFINITY> pick_physical_cores(size_t N) {
  std::vector<GROUP_AFFINITY> picks;

  // --- Preferred: CPU set API ---
  DWORD len = 0;
  GetSystemCpuSetInformation(nullptr, 0, &len, nullptr, 0); // note: nullptr process
  if (GetLastError() == ERROR_INSUFFICIENT_BUFFER && len > 0) {
    std::vector<unsigned char> buf(len);
    if (GetSystemCpuSetInformation(reinterpret_cast<SYSTEM_CPU_SET_INFORMATION*>(buf.data()),
                                   (ULONG) buf.size(), &len, nullptr, 0)) {
      std::unordered_set<ULONG> seen_core;
      BYTE                      min_eff = 255, max_eff = 0;

      BYTE* p   = buf.data();
      BYTE* end = p + len;
      struct Row {
        GROUP_AFFINITY ga;
        ULONG          core;
        BYTE           eff;
      };
      std::vector<Row> rows;

      while (p < end) {
        auto* info = reinterpret_cast<SYSTEM_CPU_SET_INFORMATION*>(p);
        if (info->Type == CpuSetInformation) {
          const auto& c = info->CpuSet;
          if (c.Allocated && !c.Parked) {
            Row r{};
            r.ga.Group = c.Group;
            r.ga.Mask  = (KAFFINITY) 1ull << c.LogicalProcessorIndex;
            r.core     = c.CoreIndex;
            r.eff      = c.EfficiencyClass;
            rows.push_back(r);
            if (r.eff < min_eff)
              min_eff = r.eff;
            if (r.eff > max_eff)
              max_eff = r.eff;
          }
        }
        p += info->Size;
      }

      if (!rows.empty()) {
        bool hybrid = (min_eff != max_eff);
        std::sort(rows.begin(), rows.end(),
                  [](const Row& a, const Row& b) { return a.core < b.core; });

        for (const auto& r: rows) {
          if (hybrid && r.eff != min_eff)
            continue; // prefer P-cores
          if (seen_core.insert(r.core).second) {
            picks.push_back(r.ga);
            if (picks.size() == N)
              break;
          }
        }
      }
    }
  }

  if (!picks.empty())
    return picks;

  // --- Fallback: LogicalProcessorInformationEx ---
  DWORD bytes = 0;
  GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &bytes);
  std::vector<unsigned char> buf2(bytes);
  if (GetLogicalProcessorInformationEx(RelationProcessorCore,
                                       reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buf2.data()), &bytes)) {

    BYTE* p2   = buf2.data();
    BYTE* end2 = p2 + bytes;
    while (p2 < end2) {
      auto* ex = reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(p2);
      if (ex->Relationship == RelationProcessorCore) {
        const GROUP_AFFINITY& g = ex->Processor.GroupMask[0];
        if (g.Mask) {
          KAFFINITY      first = g.Mask & (~g.Mask + 1);
          GROUP_AFFINITY ga{};
          ga.Group = g.Group;
          ga.Mask  = first;
          picks.push_back(ga);
          if (picks.size() == N)
            break;
        }
      }
      p2 += ex->Size;
    }
  }

  return picks;
}


static inline void raise_process_priority() {
  SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
}

static inline void configure_current_thread_for_performance() {
  THREAD_POWER_THROTTLING_STATE s{};
  s.Version     = THREAD_POWER_THROTTLING_CURRENT_VERSION;
  s.ControlMask = THREAD_POWER_THROTTLING_EXECUTION_SPEED;
  s.StateMask   = 0; // favor performance (disable EcoQoS)
  SetThreadInformation(GetCurrentThread(), ThreadPowerThrottling, &s, sizeof(s));
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
}

// ===== End helpers =====

struct Config {
  int    n_ctrl     = 0;
  int    n_points   = 0;
  int    n_optim    = 0;
  int    task_idx   = 0;
  int    config_idx = 0;
  Matrix task;
};
constexpr int                   config_size = 64 * 3; // for UR5e 64 tasks * 3 configs // for Gen3 100 tasks * 3 configs
std::array<Config, config_size> config_list;
constexpr int                   _n_optim     = 1;
int                             config_index = 0;
std::mutex                      mut_config;

// Base results
std::vector<Result> result_list;
// Results for each of the four accelerations
std::vector<Result> result_acc1_list;
std::vector<Result> result_acc3_list;
// Results for our final per_segment integration
std::vector<Result> result_segments_list;
std::vector<u32>    task_id;
std::vector<u32>    config_ids;
int                 result_index = 0;
std::mutex          mut_result;

using nlohmann::json;
// helper: convert one Result to the JSON shape you’ve been using
static inline json result_to_json(const Result& r) {
  json j;
  j["success"]       = r.success;
  j["success_false"] = r.success_false;
  j["compute time"]  = r.compute_time;

  auto& jx = j["x"] = json::array();
  jx.get_ref<json::array_t&>().reserve(r.x.size);
  for (u32 k = 0; k < r.x.size; ++k)
    jx.push_back(r.x[k]);

  auto& jx0 = j["x0"] = json::array();
  jx0.get_ref<json::array_t&>().reserve(r.x0.size);
  for (u32 k = 0; k < r.x0.size; ++k)
    jx0.push_back(r.x0[k]);

  return j;
}
struct VariantGroup {
  const char*                label; // "base" | "acc1" | "acc3" | "segments"
  const std::vector<Result>* list;  // pointer to results (same N and task_id alignment)
};

void write_results_json_simple(
        const std::string&                     path,
        const std::vector<u32>&                target_config_ids,
        const std::vector<u32>&                target_task_ids,
        const std::array<Config, config_size>& target_config_list,
        const std::vector<VariantGroup>&       variants) {
  using nlohmann::json;
  json arr = json::array();

  for (const auto& v: variants) {
    const auto& list = *v.list;
    cout << "Variant type : " << v.label << " variant size : " << list.size() << endl;
    for (size_t i = 0; i < list.size(); ++i) {
      const Config& cfg = target_config_list[target_config_ids[i]];
      const Result& r   = list[i];

      json row;
      row["variant"]       = v.label;
      row["task_idx"]      = target_task_ids[i];
      row["n_ctrl"]        = cfg.n_ctrl;
      row["n_points"]      = cfg.n_points;
      row["success"]       = r.success;
      row["success_false"] = r.success_false;
      row["compute_time"]  = r.compute_time;
      if (r.x.size)
        row["trajectory_time"] = r.x.back();

      // dump arrays x and x0
      row["x"] = json::array();
      for (u32 k = 0; k < r.x.size; ++k)
        row["x"].push_back(r.x[k]);

      row["x0"] = json::array();
      for (u32 k = 0; k < r.x0.size; ++k)
        row["x0"].push_back(r.x0[k]);

      arr.push_back(std::move(row));
    }
  }

  std::ofstream f(path);
  f << "{\n\"results\": " << arr.dump(2) << "\n}\n";
}


// We create a config list that will be used for all benchmarks
inline void fill_config_list(std::array<Config, config_size>& target_config_list) {

  // These give rounded number of n_points_per_segment (With higher n_ctrl, we can allow for fewer points per segments)
  std::array<std::tuple<int, int>, 3> bspline_config_list = {std::make_tuple(12, 15), // 15 points per segment
                                                             std::make_tuple(16, 10), // 10 points per segments
                                                             std::make_tuple(20, 10)};

  std::vector<Matrix> tasks = get_Link6_demo1_tasks(); // 13 tasks

  int config_idx = 0;
  int task_idx   = 0;

  for (auto& task: tasks) {
    for (auto& bspline_config: bspline_config_list) {
      auto [n_ctrl, n_points_per_segments] = bspline_config;
      auto n_points                        = n_points_per_segments * (n_ctrl - 5);

      target_config_list[config_idx].n_ctrl     = n_ctrl;
      target_config_list[config_idx].n_points   = n_points;
      target_config_list[config_idx].n_optim    = 100;
      target_config_list[config_idx].task_idx   = task_idx;
      target_config_list[config_idx].config_idx = config_idx;
      target_config_list[config_idx++].task     = task;
    }
    task_idx++;
  }
}

inline void fill_config_list_UR5e(std::array<Config, config_size>& target_config_list) {
  // These give rounded number of n_points_per_segment (With higher n_ctrl, we can allow for fewer points per segments)
  std::array<std::tuple<int, int>, 3> bspline_config_list = {std::make_tuple(12, 15), // 15 points per segment
                                                             std::make_tuple(16, 10), // 10 points per segments
                                                             std::make_tuple(20, 10)};
  // std::array<std::tuple<int, int>, 1> bspline_config_list = {std::make_tuple(16, 15)};
  std::vector<Matrix> tasks = get_UR5e_tasks();

  int config_idx = 0;
  int task_idx   = 0;
  // auto task       = tasks[0];
  for (auto& task: tasks) {
    for (auto& bspline_config: bspline_config_list) {
      auto [n_ctrl, n_points_per_segments] = bspline_config;
      auto n_points                        = n_points_per_segments * (n_ctrl - 5);

      target_config_list[config_idx].n_ctrl     = n_ctrl;
      target_config_list[config_idx].n_points   = n_points;
      target_config_list[config_idx].n_optim    = 100;
      target_config_list[config_idx].task_idx   = task_idx;
      target_config_list[config_idx].config_idx = config_idx;
      target_config_list[config_idx++].task     = task;
    }
    task_idx++;
  }
}

inline void fill_config_list_Gen3(std::array<Config, config_size>& target_config_list) {
  // These give rounded number of n_points_per_segment (With higher n_ctrl, we can allow for fewer points per segments)
  // std::array<std::tuple<int, int>, 3> bspline_config_list = {std::make_tuple(12, 15), // 15 points per segment
  //                                                            std::make_tuple(16, 10), // 10 points per segments
  //                                                            std::make_tuple(20, 10)};
  std::array<std::tuple<int, int>, 1> bspline_config_list = {std::make_tuple(16, 20)};
  std::vector<Matrix>                 tasks               = get_Gen3_paper_tasks();

  int config_idx = 0;
  int task_idx   = 0;
  // auto task       = tasks[0];
  for (auto& task: tasks) {
    for (auto& bspline_config: bspline_config_list) {
      auto [n_ctrl, n_points_per_segments] = bspline_config;
      auto n_points                        = n_points_per_segments * (n_ctrl - 5);

      target_config_list[config_idx].n_ctrl     = n_ctrl;
      target_config_list[config_idx].n_points   = n_points;
      target_config_list[config_idx].n_optim    = 1;
      target_config_list[config_idx].task_idx   = task_idx;
      target_config_list[config_idx].config_idx = config_idx;
      target_config_list[config_idx++].task     = task;
    }
    task_idx++;
  }
}

// ------------------------- Accelerated functions --------------------------------

// acc 1
inline u32 ncon_acc(const Optimization* opt, const int x_idx) {
  const u32 n_points            = opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]; // + 1 since ub is inclusive
  const u32 n_joints            = opt->manip.n_joints;
  const u32 n_constraints_basic = n_points * n_joints;
  u32       n_constraints       = 0;

  if (opt->constraints.position)
    n_constraints += n_constraints_basic;
  if (opt->constraints.velocity)
    n_constraints += n_constraints_basic;
  if (opt->constraints.acceleration)
    n_constraints += n_constraints_basic;
  if (opt->constraints.torque)
    n_constraints += n_constraints_basic;

  if (opt->constraints.tcp_speed)
    n_constraints += n_points;

  if (opt->constraints.external_collisions)
    n_constraints += n_points;
  if (opt->constraints.self_collisions) {
    n_constraints += n_points;
  }

  for (auto& n_con: opt->constraints.n_extra_constraints)
    n_constraints += n_con;

  return n_constraints;
}

inline void compute_constraints_grad1(double* result, const Array& x, const u32 x_idx, Optimization* opt) {
  double* moving_result = result;

  opt->bspline.compute_trajectory(x, opt->task);

  // ObjMatrix<Capsule> capsules(opt->manip._n_caps, (u32) (opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx] + 1) / opt->constraints.n_collision_skip);
  for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    // todo: reset frame???
    // rotations_computed -> bool
    // forward_kinematics_computed -> bool

    ManipulatorTempData manip_data;


    // This is calculated every loop
    auto pos = opt->bspline.traj.pos.col(i);
    Assert(pos.is_alias);

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Capsules");
#endif
      // todo: this cleaner
      if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
        compute_capsules(opt->manip, manip_data);
        // if (opt->constraints.external_collisions) {
        //   if (i % opt->constraints.n_collision_skip == 0 && i / opt->constraints.n_collision_skip < capsules.cols) {
        //     for (int j = 0; j < opt->manip._n_caps; j++) {
        //       capsules(j, i / (int) opt->constraints.n_collision_skip) = manip_data.capsule_list[j];
        //     }
        //   }
        // }
      }
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = bound_constraint(opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        moving_result++;
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        // moving_result[0] = abs_constraint(opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
        moving_result[0] = std::abs(opt->bspline.traj.vel(j, i)) / opt->manip.vmax[j] - 1.0;
        moving_result++;
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = std::abs(opt->bspline.traj.acc(j, i)) / opt->manip.amax[j] - 1.0;
        moving_result++;
      }
    }

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      auto vel = opt->bspline.traj.vel.col(i);
      auto acc = opt->bspline.traj.acc.col(i);
      dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = std::abs(manip_data.efforts[j]) / opt->manip.tau_max[j] - 1.0;
        moving_result++;
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
      moving_result[0] = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
      moving_result++;
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = max(get_internal_collisions(opt->manip, manip_data));
      // for (u32 j = 0; j < tmp_coll.size; j++)
      moving_result[0] = -tmp_coll; //*std::abs(tmp_coll[j]);
      moving_result++;
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisions");
#endif
      auto collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
      // for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
      moving_result[0] = -collisions; //*std::abs(collisions[j]);
      // }
      moving_result++;
    }
  }


  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("CustomConstraints");
#endif
    for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) { // todo: split in extra_constraint per point or once
      auto extra_constraint = opt->constraints.extra_constraints[i];
      extra_constraint(moving_result, opt);
      moving_result += opt->constraints.n_extra_constraints[i];
    }
  }
}

inline void nlopt_constraints_acc1(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, x_len);
  compute_constraints(result, xv, opt);

  // gradients calculation
  if (grad) {
    memset(grad, 0, m * x_len * sizeof(real)); // zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps = 1e-5;
    Array          x_plus(x_len);

    // lb & ub automatically calculated by bspline
    u32 x_per_joint = (x_len - 1) / opt->manip.n_joints; // = n_ctrl - 6 (skip first and last 3)
    u32 mult        = 0;

    for (u32 j = 0; j < x_len - 1; j++) {     // last one is T
      mult      = j / x_per_joint > mult ? mult + 1 : mult;
      u32 x_idx = j - mult * x_per_joint + 3; // todo: fix for NaN in PVA of task

      memcpy(x_plus.data, x, x_len * sizeof(real));
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, x_len);

      // compute constraints
      auto  n_con = ncon_acc(opt, x_idx);
      Array r_plus(n_con);
      compute_constraints_grad1(r_plus.data, x_plus, x_idx, opt);
      auto n_con_lb = ncon_lb_acc(opt, x_idx);

      // todo: fix this loop
      // for (u32 i = 0; i < m; i++)
      for (u32 i = 0; i < n_con; i++)
        grad[(n_con_lb + i) * x_len + j] = (r_plus[i] - result[n_con_lb + i]) / eps;
    }
    {
      // last point T
      u32 j = x_len - 1;
      memcpy(x_plus.data, x, x_len * sizeof(real));
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, x_len);

      // compute constraints
      Array r_plus(m);
      compute_constraints(r_plus.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * x_len + j] = (r_plus[i] - result[i]) / eps;
    }
  }
}

inline Result optimize_acc1(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, 0.01);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = 1000;
  stop.maxtime    = 5;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_acc1;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc1, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 5);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 1000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
      // note: can we initialize grad to 0 here
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      is_valid = max(constraints_points) < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      is_valid_more = max(constraints_more_points) < opt_val_more.success_tolerance;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// acc 2
inline u32 ncon_acc2(const Optimization* opt, const int x_idx) {
  const u32 n_points            = opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]; // + 1 since ub is inclusive
  const u32 n_joints            = opt->manip.n_joints;
  const u32 n_constraints_basic = n_points * n_joints;
  u32       n_constraints       = 0;

  if (opt->constraints.position)
    n_constraints += n_points;
  if (opt->constraints.velocity)
    n_constraints += n_points;
  if (opt->constraints.acceleration)
    n_constraints += n_points;
  if (opt->constraints.torque)
    n_constraints += n_constraints_basic;

  if (opt->constraints.tcp_speed)
    n_constraints += n_points;

  if (opt->constraints.external_collisions)
    n_constraints += n_points;
  if (opt->constraints.self_collisions) {
    n_constraints += n_points;
  }

  for (auto& n_con: opt->constraints.n_extra_constraints)
    n_constraints += n_con;

  return n_constraints;
}

struct ConstraintPerPoint {
  Array  pos_constraint;
  Array  vel_constraint;
  Array  acc_constraint;
  Matrix tor_constraint;
  Array  tcp_constraint;
  Array  collision_constraint;
  Array  self_collision_constraint;

  ConstraintPerPoint(int joints, int points) {
    pos_constraint.resize(points);
    vel_constraint.resize(points);
    acc_constraint.resize(points);
    tor_constraint.resize(joints, points);
    tcp_constraint.resize(points);
    collision_constraint.resize(points);
    self_collision_constraint.resize(points);
  }
  ~ConstraintPerPoint() = default;
};

inline void compute_constraints_grad2(ConstraintPerPoint& constraints, const Array& x, const u32 x_idx, const u32 joint_idx, Optimization* opt) {

  opt->bspline.compute_trajectory(x, opt->task);

  auto n_joints = opt->manip.n_joints;

  // ObjMatrix<Capsule> capsules(opt->manip._n_caps, (u32) (opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx] + 1) / opt->constraints.n_collision_skip);
  for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    // todo: reset frame???
    // rotations_computed -> bool
    // forward_kinematics_computed -> bool

    ManipulatorTempData manip_data;


    // This is calculated every loop
    auto pos = opt->bspline.traj.pos.col(i);
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    Assert(pos.is_alias);

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Capsules");
#endif
      // todo: this cleaner
      // if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
      compute_capsules(opt->manip, manip_data);
      // }
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.pos_constraint[i - opt->bspline.lb[x_idx]] = bound_constraint(pos[joint_idx], opt->manip.pmin[joint_idx], opt->manip.pmax[joint_idx]);

      // }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.vel_constraint[i - opt->bspline.lb[x_idx]] = std::abs(vel[joint_idx]) / opt->manip.vmax[joint_idx] - 1.0;
      // }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.acc_constraint[i - opt->bspline.lb[x_idx]] = std::abs(acc[joint_idx]) / opt->manip.amax[joint_idx] - 1.0;

      // }
    }

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

      for (int j = 0; j < opt->manip.n_joints; j++) {
        constraints.tor_constraint(j, i - opt->bspline.lb[x_idx]) = std::abs(manip_data.efforts[j]) / opt->manip.tau_max[j] - 1.0;
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool                                            = get_J_tool(opt, manip_data);
      real tcp_speed                                         = norm(J_tool * opt->bspline.traj.vel.col(i));
      constraints.tcp_constraint[i - opt->bspline.lb[x_idx]] = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = max(-get_internal_collisions(opt->manip, manip_data));
      // for (u32 j = 0; j < tmp_coll.size; j++)
      constraints.self_collision_constraint[i - opt->bspline.lb[x_idx]] = tmp_coll; //*std::abs(tmp_coll[j]);
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisions");
#endif
      auto collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
      // for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
      constraints.collision_constraint[i - opt->bspline.lb[x_idx]] = -collisions; //*std::abs(collisions[j]);
    }
  }
}

inline void nlopt_constraints_acc2(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);
  compute_constraints(result, xv, opt);

  // gradients calculation
  if (grad) {
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps = 1e-5;

    u32   n_con_lb = 0;
    Array x_plus(xlen);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx       = 0;
    u32 constraint_idx = 0;
    u32 x_idx          = 0;

    auto n_joints = opt->manip.n_joints;

    memcpy(x_plus.data, x, xlen * sizeof(real)); // todo: is this necessary ? done once, and x_plus is modified to original at the end of each iteration

    for (u32 j = 0; j < xlen - 1; j++) {
      // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
      // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
      joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
      x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

      x_plus[j] += eps;                                    // todo: add this is for extra constraints (tcp, collisions)
      opt->bspline.compute_trajectory(x_plus, opt->task);

      auto n_con_per_point = ncon_acc(opt, x_idx) / (opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]);

      auto               n_points = (opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]);
      ConstraintPerPoint constraint(n_joints, n_points);
      compute_constraints_grad2(constraint, x_plus, x_idx, joint, opt);

      n_con_lb = ncon_lb_acc(opt, x_idx); // find the amount of constraints before the current point

      grad_idx = n_con_lb * xlen + j;     // gradients are stored column-wise xlen * npoints
      for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) {
        // lb & ub are inclusive
        grad_idx += joint * xlen;

        // position
        grad[grad_idx] = (constraint.pos_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + joint]) / eps;
        grad_idx += n_joints * xlen;

        // velocity
        grad[grad_idx] = (constraint.vel_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + n_joints + joint]) / eps;
        grad_idx += n_joints * xlen;

        // acceleration
        grad[grad_idx] = (constraint.acc_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + 2 * n_joints + joint]) / eps;
        grad_idx += (n_joints - joint) * xlen;

        // torque
        for (int k = 0; k < n_joints; k++) {
          grad[grad_idx] = (constraint.tor_constraint(k, i - opt->bspline.lb[x_idx]) - result[i * n_con_per_point + 3 * n_joints + k]) / eps;
          grad_idx += xlen;
        }

        // tcp
        grad[grad_idx] = (constraint.tcp_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + 4 * n_joints]) / eps;
        grad_idx += xlen;

        // self col
        grad[grad_idx] = (constraint.self_collision_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + 4 * n_joints + 1]) / eps;
        grad_idx += xlen;

        // col
        grad[grad_idx] = (constraint.collision_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + 4 * n_joints + 2]) / eps;
        grad_idx += xlen;
      }

      x_plus[j] = x[j];
    }

    {
      // last point T
      u32 j = xlen - 1;
      memcpy(x_plus.data, x, xlen * sizeof(real));
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, xlen);

      // compute constraints
      Array r_plus(m);
      compute_constraints(r_plus.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * xlen + j] = (r_plus[i] - result[i]) / eps;
    }
  }
}

inline Result optimize_acc2(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = 1000;
  stop.maxtime    = 30;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_acc2;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc2, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 5);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 1000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
      // note: can we initialize grad to 0 here
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      auto max_con = max(constraints_points);
      // cout << "max_con = " << max_con << endl;
      is_valid = max_con < opt->success_tolerance * 2;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more = max(constraints_more_points);
      // cout << "max_con_more = " << max_con_more << endl;
      is_valid_more = max_con_more < opt->success_tolerance * 2;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// acc 3

template<bool is_grad>
blast_fn void compute_constraints_acc3(double* result, Array& gradient_coeffs, Matrix& dtau_dp, Matrix& dtau_dv, Matrix& dtau_da, /*Array& gradient_coeffs_collisions,*/ const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  ZoneScoped;
#endif

  double* moving_result   = result;
  u32     grad_idx        = 0;
  u32     grad_torque_idx = 0;
  // Array   external_collisions(opt->bspline.n_points);
  // Array   external_collisions_plus(opt->bspline.n_points);

  auto joints = opt->manip.n_joints;

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // Lambda to process common bound constraint operations
  auto process_bound = [&](real value, real bound_max) {
    auto [constraint, gradient_coeff] = abs_constraint_dev<is_grad>(value, bound_max);
    *moving_result++                  = constraint;
    gradient_coeffs[grad_idx++]       = gradient_coeff;
  };

  for (u32 i = 0; i < opt->bspline.n_points; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    ManipulatorTempData manip_data;
    // Array               torque_constraint_plus(joints);

    constexpr real eps = 1e-5;

    auto pos = opt->bspline.traj.pos.col(i);
    Assert(pos.is_alias);

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos);
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      for (int j = 0; j < joints; j++) {
        auto [constraint, gradient_coeff] = bound_constraint_dev<is_grad>(opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        *moving_result++                  = constraint;
        gradient_coeffs[grad_idx++]       = gradient_coeff;
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      for (int j = 0; j < joints; j++) {
        process_bound(opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      for (int j = 0; j < joints; j++) {
        process_bound(opt->bspline.traj.acc(j, i), opt->manip.amax[j]);
      }
    }


    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      auto vel = opt->bspline.traj.vel.col(i);
      auto acc = opt->bspline.traj.acc.col(i);
      dynamics(opt->manip, manip_data, vel, acc);

      Array pos_plus(joints);
      Array vel_plus(joints);
      Array acc_plus(joints);
      Array torque_constraint(joints);

      for (int j = 0; j < joints; j++) {
        torque_constraint[j] = abs_constraint(manip_data.efforts[j], opt->manip.tau_max[j]);
        *moving_result++     = torque_constraint[j];
        grad_idx++; // todo: add analytical gradients
      }

      pos_plus = pos;
      vel_plus = vel;
      acc_plus = acc;

      // grad_coeffs pos
      for (u32 j = 0; j < joints; j++) {
        pos_plus[j] += eps;
        forward_kinematics(opt->manip, manip_data, pos_plus);
        dynamics(opt->manip, manip_data, vel_plus, acc_plus);
        pos_plus[j] -= eps;

        for (u32 k = 0; k < joints; k++) {
          auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
          dtau_dp(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
        }
      }

      forward_kinematics(opt->manip, manip_data, pos_plus);
      // grad_coeffs vel
      for (u32 j = 0; j < joints; j++) {
        vel_plus[j] += eps;
        dynamics(opt->manip, manip_data, vel_plus, acc_plus);
        vel_plus[j] -= eps;

        for (u32 k = 0; k < joints; k++) {
          auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
          dtau_dv(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
        }
      }

      // grad_coeffs acc
      for (u32 j = 0; j < joints; j++) {
        acc_plus[j] += eps;
        dynamics(opt->manip, manip_data, vel_plus, acc_plus);
        acc_plus[j] -= eps;

        for (u32 k = 0; k < joints; k++) {
          auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
          dtau_da(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
        }
      }
    }

    if (opt->constraints.tcp_speed) { // todo: no grad_coeffs yet
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
      *moving_result++ = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
      grad_idx++; // todo: add analytical gradients
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Capsules");
#endif
      if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
        compute_capsules(opt->manip, manip_data);
      }
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
      for (u32 j = 0; j < tmp_coll.size; j++) {
        moving_result[j] = -tmp_coll[j] + 0.01;
      }
      moving_result += tmp_coll.size;
      grad_idx += tmp_coll.size; // todo: add analytical gradients
    }

    //     if (opt->constraints.external_collisions) {
    // #if BLAST_TRACE_LEVEL >= 3
    //       ZoneScopedN("ExternalCollisionsCalculate");
    // #endif
    //       external_collisions[i] = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
    //
    //       for (u32 j = 0; j < joints; j++) {
    //         pos[j] += eps;
    //         forward_kinematics(opt->manip, manip_data, pos);
    //         compute_capsules(opt->manip, manip_data);
    //         external_collisions_plus[i] = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
    //       }
    //     }
  }

  //   if (opt->constraints.external_collisions) {
  // #if BLAST_TRACE_LEVEL >= 3
  //     ZoneScopedN("ExternalCollisionsSort");
  // #endif
  //     Array dist_min(opt->constraints.n_collision_constraints, INF_REAL);
  //     Array dist_min_plus(opt->constraints.n_collision_constraints, INF_REAL);
  //     for (u32 i = 0; i < opt->bspline.n_points; i++) {
  //       for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
  //         if (external_collisions[i] < dist_min[j]) {
  //           for (int k = (int) opt->constraints.n_collision_constraints - 1; k > j; k--) {
  //             dist_min[k]      = dist_min[k - 1];
  //             dist_min_plus[k] = dist_min_plus[k - 1];
  //           }
  //           dist_min[j]      = external_collisions[i];
  //           dist_min_plus[j] = external_collisions_plus[i];
  //           break;
  //         }
  //       }
  //     }
  //     for (u32 i = 0; i < opt->constraints.n_collision_constraints; i++) {
  //       moving_result[i] = -dist_min[i] + 0.01; //*std::abs(collisions[j]);
  //     }
  //     moving_result += opt->constraints.n_collision_constraints;

  //     for (u32 i = 0; i < opt->constraints.n_collision_constraints; i++) {
  //       gradient_coeffs_collisions
  //     }
  //     grad_idx += opt->constraints.n_collision_constraints; // todo: add analytical gradients
  //   }
  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("CustomConstraints");
#endif
    for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) {
      auto extra_constraint = opt->constraints.extra_constraints[i];
      extra_constraint(moving_result, opt);
      moving_result += opt->constraints.n_extra_constraints[i];
      grad_idx += opt->constraints.n_extra_constraints[i]; // todo: add analytical gradients
    }
  }
}

inline blast_fn void compute_constraints_per_point_acc3(double* result, real& external_collisions, const u32 i, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 3
  ZoneScopedN("ConstraintsPerPoint");
#endif
  double* moving_result = result;

  ManipulatorTempData manip_data;

  auto pos = opt->bspline.traj.pos.col(i);
  Assert(pos.is_alias);

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Kinematics");
#endif
    forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
  }

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Capsules");
#endif
    if (opt->constraints.self_collisions || opt->constraints.external_collisions) {
      compute_capsules(opt->manip, manip_data);
    }
  }

  if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Dynamics");
#endif
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

    for (int j = 0; j < opt->manip.n_joints; j++) {
      moving_result[0] = (std::abs(manip_data.efforts[j]) - opt->manip.tau_max[j]) / opt->manip.tau_max[j];
      moving_result++;
    }
  }

  if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("TCPSpeed");
#endif
    auto J_tool      = get_J_tool(opt, manip_data);
    real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
    moving_result[0] = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
    moving_result++;
  }

  if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("SelfCollisions");
#endif
    auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
    for (u32 j = 0; j < tmp_coll.size; j++)
      moving_result[j] = -tmp_coll[j] + 0.01; //*std::abs(tmp_coll[j]);
    moving_result += tmp_coll.size;
  }

  if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ExternalCollisionsCalculate");
#endif
    external_collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
  }
}

inline blast_fn void nlopt_constraints_acc3(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);

  Array  gradient_coeffs(m);                                                        // todo: check performance
  Matrix dtau_dp(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points); // todo: check performance
  Matrix dtau_dv(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points); // todo: check performance
  Matrix dtau_da(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points); // todo: check performance
  if (!grad) {
    compute_constraints_acc3<false>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, xv, opt);
  }
  if (grad) {
    compute_constraints_acc3<true>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, xv, opt);
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps      = 1e-5;
    auto           n_con    = ncon_extras(opt);
    u32            n_con_lb = 0;
    Array          x_plus(xlen);
    Array          r_plus(n_con);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx       = 0;
    u32 constraint_idx = 0;
    u32 x_idx          = 0;

    u32 n_torque         = 0;
    u32 n_tcp_speed      = 0;
    u32 n_self_collision = 0;

    memcpy(x_plus.data, x, xlen * sizeof(real)); // todo: is this necessary ? done once, and x_plus is modified to original at the end of each iteration

    for (u32 j = 0; j < xlen - 1; j++) {         // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
      // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
      joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
      x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

      // x_plus[j] += eps;                                    // todo: add this is for extra constraints (tcp, collisions)
      // opt->bspline.compute_trajectory(x_plus, opt->task);

      n_con_lb = ncon_lb_acc(opt, x_idx); // find the amount of constraints before the current point

      Array external_collisions(opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx]);
      // todo: create alias matrix that points to grad
      // todo: can we change the order in which we store the gradients ?
      grad_idx       = n_con_lb * xlen + j;                                    // gradients are stored column-wise xlen * npoints
      constraint_idx = n_con_lb;
      for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) { // lb & ub are inclusive
        grad_idx += joint * xlen;
        constraint_idx += joint;

        if (opt->constraints.position) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * gradient_coeffs[constraint_idx];
          grad_idx += opt->manip.n_joints * xlen; // increase index by the amount of joints * xlen
          constraint_idx += opt->manip.n_joints;
        }

        if (opt->constraints.velocity) {          // todo: basis_v / t once before this loop
          grad[grad_idx] = opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1] * gradient_coeffs[constraint_idx];
          grad_idx += opt->manip.n_joints * xlen; // increase index by the amount of joints * xlen
          constraint_idx += opt->manip.n_joints;
        }

        if (opt->constraints.acceleration) {                // todo: basis_a / t / t once before this loop
          grad[grad_idx] = opt->bspline.basis_a(x_idx, i) / (opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1] * opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1]) * gradient_coeffs[constraint_idx];
          grad_idx += (opt->manip.n_joints - joint) * xlen; // increase index by the amount of (joints - current joint) * xlen
          constraint_idx += (opt->manip.n_joints - joint);  // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }

        if (opt->constraints.torque) { // todo: add analytical gradients for torque

          for (u32 k = 0; k < opt->manip.n_joints; k++) {
            grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dtau_dp(k, joint + opt->manip.n_joints * i) +
                             opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * dtau_dv(k, joint + opt->manip.n_joints * i) +
                             opt->bspline.basis_a(x_idx, i) / (opt->bspline.traj.t.back() * opt->bspline.traj.t.back()) * dtau_da(k, joint + opt->manip.n_joints * i);
            grad_idx += xlen;
            constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
          }
        }

        if (opt->constraints.tcp_speed || opt->constraints.self_collisions || opt->constraints.external_collisions) {
          compute_constraints_per_point_acc3(r_plus.data, external_collisions[i], i, opt); // todo: not working ...
          // n_torque         = 0;
          n_tcp_speed      = 0;
          n_self_collision = 0;

          if (opt->constraints.tcp_speed) { // todo: not working
            n_tcp_speed = 1;

            grad[grad_idx] = (r_plus.data[0] - result[constraint_idx]) / eps;
            grad_idx += xlen;
            constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
          }

          if (opt->constraints.self_collisions) { // todo: not working
            n_self_collision = opt->manip._n_internal_collisions;
            for (u32 k = n_tcp_speed; k < n_tcp_speed + n_self_collision; k++) {
              grad[grad_idx] = (r_plus.data[k] - result[constraint_idx]) / eps;
              grad_idx += xlen;
              constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
            }
          }
        }
      }

      // x_plus[j] -= eps;
    }

    {
      // last point T
      u32 j = xlen - 1;
      // memcpy(x_plus.data, x, xlen * sizeof(real)); todo: check if necessary
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, xlen);

      // compute constraints
      Array r_plus_T(m);
      compute_constraints(r_plus_T.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * xlen + j] = (r_plus_T[i] - result[i]) / eps;
    }

    if (opt->constraints.show_info) { // when more info is needed per iteration
      Matrix gradients(xlen, m);
      Array  constraints(m);
      for (u32 j = 0; j < xlen; j++) {
        for (u32 i = 0; i < m; i++) {
          gradients(j, i) = grad[i * xlen + j];
          constraints[i]  = result[i];
        }
      }
      opt->constraints.grad_list.push_back(gradients);
      opt->constraints.constr_list.push_back(constraints);
      opt->constraints.x_list.push_back(xv);
    }
  }
}

inline Result optimize_acc3(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.000001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.0001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = 100000;
  stop.maxtime    = 5000;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_acc3;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc3, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      is_valid = max(constraints_points) < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      is_valid_more = max(constraints_more_points) < opt->success_tolerance;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// with_segments (_final)
inline blast_fn void constraints_and_gradients_final(const Array& x, Optimization& opt, Array& constraints, Matrix& grad) {
  // constraints (p,v,a,tor) for each joint, for each segment
  // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]

  // basis: n_ctrl x n_points
  Assert(constraints.is_alias);
  if (grad.size) {
    Assert(grad.is_alias);
    Assert(grad.rows == x.size);
    Assert(grad.cols == constraints.size);
  }
  const int n_segments                = (int) opt.bspline.n_ctrl - (int) opt.bspline.p;
  const int n_points_per_segment      = (int) opt.bspline.n_points / n_segments; // todo: check if fine?
  const int n_joints                  = (int) opt.manip.n_joints;
  const int n_ctrl                    = (int) opt.bspline.n_ctrl;
  const int x_len                     = (int) x.size;
  const int n_capsules                = (int) opt.manip._n_caps;
  const int n_constraints_per_segment = (n_joints * 4); //  todo: remove hard-coded 4 (p.v.a.tau)
  Assert(constraints.size == n_segments * n_constraints_per_segment);

  // limits
  auto& pmax    = opt.manip.pmax;
  auto& pmin    = opt.manip.pmin;
  auto& vmax    = opt.manip.vmax;
  auto& amax    = opt.manip.amax;
  auto& tau_max = opt.manip.tau_max;
  for (int j = 0; j < n_joints; j++) {
    // todo: document the current behaviour in the API
    //        (doesn't currently work if one is inf and the other is not)
    if (pmax[j] == INF_REAL)  // note: replace INF_REAL with huge value
      pmax[j] = 1e300;
    if (pmin[j] == -INF_REAL) // note: replace -INF_REAL with huge negative value
      pmin[j] = -1e300;
  }

  ManipulatorTempData        manip_data;
  std::array<u8, MAX_JOINTS> max_pos_indices{};
  std::array<u8, MAX_JOINTS> max_vel_indices{};
  std::array<u8, MAX_JOINTS> max_acc_indices{};
  std::array<u8, MAX_JOINTS> max_tor_indices{};

  opt.bspline.compute_trajectory(x, opt.task);


  for (int segment = 0; segment < n_segments; segment++) {
#if BLAST_TRACE_LEVEL >= 2
    ZoneScopedN("All Segment Constraints");
#endif
    const int first_affected_control_point = std::max(3, segment);
    const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.p);
    const int n_affected_control_points    = last_affected_control_point - first_affected_control_point + 1; // note: affected_control_points are inclusive, so when we have last = 5, first = 3, we want 3 (5 - 3) + 1
    const int start_point_for_segment      = segment * n_points_per_segment;                                 // note:
    Assert(n_affected_control_points >= 3);
    Assert(n_affected_control_points <= 6);

    Matrix bp(&opt.bspline.basis_p(0, start_point_for_segment), n_ctrl, n_points_per_segment);    // note:
    Matrix bv(&opt.bspline.basis_v(0, start_point_for_segment), n_ctrl, n_points_per_segment);    // note:
    Matrix ba(&opt.bspline.basis_a(0, start_point_for_segment), n_ctrl, n_points_per_segment);    // note:

    Array max_pos_constraints(n_joints, -INF_REAL);                                               // note:
    Array max_vel_constraints(n_joints, -INF_REAL);                                               // note:
    Array max_acc_constraints(n_joints, -INF_REAL);                                               // note:
    Array max_tor_constraints(n_joints, -INF_REAL);                                               // note:

    for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) { // note:
      auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment);              // note:
      auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment);              // note:
      auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment);              // note:

      forward_kinematics(opt.manip, manip_data, p);
      dynamics(opt.manip, manip_data, v, a);

      for (int j = 0; j < n_joints; j++) {
        // position
        if (auto c = bound_constraint(p[j], pmin[j], pmax[j]); c > max_pos_constraints[j]) {
          max_pos_constraints[j] = c;
          max_pos_indices[j]     = point_in_segment; // note:
        }
        // velocity
        if (auto c = std::abs(v[j]) / vmax[j] - 1.0;
            c > max_vel_constraints[j]) {
          max_vel_constraints[j] = c;
          max_vel_indices[j]     = point_in_segment; // note:
        }
        // acceleration
        if (auto c = std::abs(a[j]) / amax[j] - 1.0;
            c > max_acc_constraints[j]) {
          max_acc_constraints[j] = c;
          max_acc_indices[j]     = point_in_segment; // note:
        }
        // torque
        if (auto c = std::abs(manip_data.efforts[j]) / tau_max[j] - 1.0;
            c > max_tor_constraints[j]) {
          max_tor_constraints[j] = c;
          max_tor_indices[j]     = point_in_segment; // note:
        }
      }
    }

    // at this point we have max constraints for pos, vel, acc, tor, for this segment

    // fill in the constraints for the current segment
    // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
    auto fill_idx = segment * n_constraints_per_segment;
    std::copy_n(max_pos_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints; // note (andre): we can use the comma operator because we don't need the result of copy_n()
    std::copy_n(max_vel_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    std::copy_n(max_acc_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    std::copy_n(max_tor_constraints.data, n_joints, &constraints[fill_idx]);


    // The gradient should be a (x_len)x(n_constraints) matrix that looks like this:
    // [dp1/dx1, dp2/dx1, ..., dv1/dx1, dv2/dx1, ..., da1/dx1, da2/dx1, ..., dt1/dx1, dt2/dx1]
    // [dp1/dx2, dp2/dx2, ..., dv1/dx2, dv2/dx2, ..., da1/dx2, da2/dx2, ..., dt1/dx2, dt2/dx2]
    // [dp1/dx3, dp2/dx3, ..., dv1/dx3, dv2/dx3, ..., da1/dx3, da2/dx3, ..., dt1/dx3, dt2/dx3]
    // [dp1/dx4, dp2/dx4, ..., dv1/dx4, dv2/dx4, ..., da1/dx4, da2/dx4, ..., dt1/dx4, dt2/dx4]
    // [dp1/dx5, dp2/dx5, ..., dv1/dx5, dv2/dx5, ..., da1/dx5, da2/dx5, ..., dt1/dx5, dt2/dx5]
    // [dp1/dx6, dp2/dx6, ..., dv1/dx6, dv2/dx6, ..., da1/dx6, da2/dx6, ..., dt1/dx6, dt2/dx6]
    // [.....................]
    // [.....................]
    // [.....................]
    // [.....................]
    // [dp1/dT=0,dp2/dT=0,..., dv1/dT , dv1/dT , ..., da1/dT , da2/dT , ..., dt1/dT , dt2/dT ]
    // *** per segment ***
    // where x is the optimization vector
    if (grad.size) {
      // Matrix (alias) in which we can insert the gradient for the current segment
      Matrix grad_segment(&grad(0, segment * n_constraints_per_segment), x_len, n_constraints_per_segment);
      Assert(grad_segment.is_alias);

      int con = 0;

      // positions
      for (int joint = 0; joint < n_joints; joint++) {

        real pos = opt.bspline.traj.pos(joint, max_pos_indices[joint]);

        // Array of the column where to put the gradient for the current constraint
        Array fill_column = grad_segment.col(con);
        Assert(con == joint); // note:
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // Which values in 'x' affect the current joint's position
        auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
        // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
        x_idx += first_affected_control_point - 3;

        // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
        Array bp_to_use(&bp(first_affected_control_point, max_pos_indices[joint]), n_affected_control_points);
        Assert(bp_to_use.is_alias);

        real coeff = 2.0 * std::abs(opt.bspline.traj.pos(joint, start_point_for_segment + max_pos_indices[joint])) / (pmax[joint] - pmin[joint]); // note:

        for (int i = 0; i < n_affected_control_points; i++) {
          fill_column[x_idx++] = bp_to_use[i] * coeff;
        }

        // note: dp/dT == 0
        con++;
      }

      // velocities
      auto one_over_T = 1 / opt.bspline.traj.t.back();
      for (int joint = 0; joint < n_joints; joint++) {
        // Array of the column where to put the gradient for the current constraint
        Array fill_column = grad_segment.col(con);
        Assert(con == n_joints + joint); // note:
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // Which values in 'x' affect the current joint's position
        auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
        // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
        x_idx += first_affected_control_point - 3;

        Array bv_to_use(&bv(first_affected_control_point, max_vel_indices[joint]), n_affected_control_points);
        Assert(bv_to_use.is_alias);

        real coeff = sign(opt.bspline.traj.vel(joint, start_point_for_segment + max_vel_indices[joint])) / vmax[joint] * one_over_T; // note:

        for (int i = 0; i < n_affected_control_points; i++) {
          fill_column[x_idx++] = bv_to_use[i] * coeff;
        }

        // dvj/dT = - (Cv + 1) / T
        fill_column.back() = -(max_vel_constraints[joint] + 1) * one_over_T;

        con++;
      }

      // accelerations
      auto one_over_T2 = one_over_T * one_over_T;
      for (int joint = 0; joint < n_joints; joint++) {
        // Array of the column where to put the gradient for the current constraint
        Array fill_column = grad_segment.col(con);
        Assert(con == 2 * n_joints + joint); // note:
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // Which values in 'x' affect the current joint's position
        auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
        // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
        x_idx += first_affected_control_point - 3;

        Array ba_to_use(&ba(first_affected_control_point, max_acc_indices[joint]), n_affected_control_points);
        Assert(ba_to_use.is_alias);

        real coeff = sign(opt.bspline.traj.acc(joint, start_point_for_segment + max_acc_indices[joint])) / amax[joint] * one_over_T2; // note:

        for (int i = 0; i < n_affected_control_points; i++) {
          fill_column[x_idx++] = ba_to_use[i] * coeff;
        }

        // daj/dT = -2 * (Ca + 1) / T
        fill_column.back() = -2 * (constraints[con] + 1) * one_over_T;

        con++;
      }

      // torque
      // [dt0/dp0, dt1/dp0, ..., dt4/dp0, dt5/dp0]
      // [dt0/dp1, dt1/dp1, ..., dt4/dp1, dt5/dp1]
      // [dt0/dp2, dt1/dp2, ..., dt4/dp2, dt5/dp2]
      // [dt0/dp3, dt1/dp3, ..., dt4/dp3, dt5/dp3]
      // [.
      // [.
      // [.
      for (int joint = 0; joint < n_joints; joint++) {
        Array fill_column = grad_segment.col(con);
        Assert(con == 3 * n_joints + joint); // note: Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        constexpr real eps            = 1e-5;
        const auto     point          = max_tor_indices[joint];
        auto           p              = opt.bspline.traj.pos.col(start_point_for_segment + point); // note:
        auto           v              = opt.bspline.traj.vel.col(start_point_for_segment + point); // note:
        auto           a              = opt.bspline.traj.acc.col(start_point_for_segment + point); // note:
        auto           old_constraint = max_tor_constraints[joint];                                // note:
        auto           tau_max_now    = tau_max[joint];                                            // note:

        auto p_plus = p;
        auto v_plus = v;
        auto a_plus = a;

        Array bp_to_use(&bp(first_affected_control_point, point), n_affected_control_points);
        Array bv_to_use(&bv(first_affected_control_point, point), n_affected_control_points);
        Array ba_to_use(&ba(first_affected_control_point, point), n_affected_control_points);

        auto x_idx      = first_affected_control_point - 3;
        auto x_idx_skip = n_ctrl - 6;

        for (int j = 0; j < n_joints; j++) {

          // partial derivative of torque constraints w.r.t. position
          // finite difference on position
          p_plus[j] += eps;
          // compute the derivative of constraint(joint) w.r.t. theta(j). (remember, joint != j)
          forward_kinematics(opt.manip, manip_data, p_plus);
          dynamics(opt.manip, manip_data, v_plus, a_plus);
          const real new_constraint_p = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // note: todo: remove 1.01 and use for validation only
          const real dtau_dp          = (new_constraint_p - old_constraint) / eps;
          // reset finite difference
          p_plus[j] = p[j];

          // note: reset forward kinematics because 'v' and 'a' don't change it
          forward_kinematics(opt.manip, manip_data, p); // todo: precompute once

          // partial derivative of torque constraints w.r.t. velocity
          // finite difference on velocity
          v_plus[j] += eps;
          dynamics(opt.manip, manip_data, v_plus, a_plus);
          const real new_constraint_v = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // note: todo: remove 1.01 and use for validation only
          const real dtau_dv          = (new_constraint_v - old_constraint) / eps;
          // reset finite difference
          v_plus[j] = v[j];

          // partial derivative of torque constraints w.r.t. acceleration
          // finite difference on acceleration
          a_plus[j] += eps;
          dynamics(opt.manip, manip_data, v_plus, a_plus);
          const real new_constraint_a = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // note: todo: remove 1.01 and use for validation only
          const real dtau_da          = (new_constraint_a - old_constraint) / eps;
          // reset finite difference
          a_plus[j] = a[j];

          // insert into the gradient
          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx + i] = bp_to_use[i] * dtau_dp +
                                     bv_to_use[i] * dtau_dv * one_over_T +
                                     ba_to_use[i] * dtau_da * one_over_T2;
          }
          x_idx += x_idx_skip;

          // gradient w.r.t. T
          fill_column.back() += dtau_dv * (-v[j] * one_over_T) + dtau_da * (-2 * a[j] * one_over_T); // note:
        }
        con++;                                                                                       // note: moved out of the inner j loop (only changes at the end of all j torques per joint)
      }
    }
  }
}

inline blast_fn void compute_constraints_final(const Array& x, Optimization& opt, Array& constraints) {
  // constraints (p,v,a,tor) for each joint, for each segment
  // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]

  // basis: n_ctrl x n_points
  // Assert(constraints.is_alias); todo: fix for optimization validation

  const int n_segments                = (int) opt.bspline.n_ctrl - (int) opt.bspline.p;
  const int n_points_per_segment      = (int) opt.bspline.n_points / n_segments; // todo: check if fine? -> (Nikos) constructor that sets n_points from n_points_per_segment
  const int n_joints                  = (int) opt.manip.n_joints;
  const int n_ctrl                    = (int) opt.bspline.n_ctrl;
  const int x_len                     = (int) x.size;
  const int n_constraints_per_segment = (n_joints * 4); //  todo: remove hard-coded 4
  Assert(constraints.size == n_segments * n_constraints_per_segment);

  // limits
  auto pmax    = opt.manip.pmax;
  auto pmin    = opt.manip.pmin;
  auto vmax    = opt.manip.vmax;
  auto amax    = opt.manip.amax;
  auto tau_max = opt.manip.tau_max;

  ManipulatorTempData manip_data;
  std::vector<u8>     max_pos_indices(n_joints);
  std::vector<u8>     max_vel_indices(n_joints);
  std::vector<u8>     max_acc_indices(n_joints);
  std::vector<u8>     max_tor_indices(n_joints);

  opt.bspline.compute_trajectory(x, opt.task);


  for (int segment = 0; segment < n_segments; segment++) {
#if BLAST_TRACE_LEVEL >= 2
    ZoneScopedN("All Segment Constraints");
#endif
    const int first_affected_control_point = std::max(3, segment);
    const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.p);
    const int n_affected_control_points    = last_affected_control_point - first_affected_control_point + 1; // note: affected_control_points are inclusive, so when we have last = 5, first = 3, we want 3 (5 - 3) + 1
    const int start_point_for_segment      = segment * n_points_per_segment;                                 // note:
    Assert(n_affected_control_points >= 3);
    Assert(n_affected_control_points <= 6);

    Matrix bp(&opt.bspline.basis_p(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
    Matrix bv(&opt.bspline.basis_v(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
    Matrix ba(&opt.bspline.basis_a(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:

    // Initialize max_constraints with -INF_REAL
    Array max_pos_constraints(n_joints, -INF_REAL);                                               // note:
    Array max_vel_constraints(n_joints, -INF_REAL);                                               // note:
    Array max_acc_constraints(n_joints, -INF_REAL);                                               // note:
    Array max_tor_constraints(n_joints, -INF_REAL);                                               // note:

    for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) { // note:
      auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment);              // note:
      auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment);              // note:
      auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment);              // note:

      forward_kinematics(opt.manip, manip_data, p);
      dynamics(opt.manip, manip_data, v, a);

      for (int j = 0; j < n_joints; j++) {
        {
          // todo (andre): do this when initializing the optimization, not here
          // todo: document the current behaviour in the API
          //        (doesn't currently work if one is inf and the other is not)
          if (pmax[j] == INF_REAL)  // note: replace INF_REAL with huge value
            pmax[j] = 1e300;
          if (pmin[j] == -INF_REAL) // note: replace -INF_REAL with huge negative value
            pmin[j] = -1e300;
        }
        // position
        if (auto c = bound_constraint(p[j], pmin[j], pmax[j]); // note: todo: remove 1.01 and use for validation only
            c > max_pos_constraints[j]) {
          max_pos_constraints[j] = c;
          max_pos_indices[j]     = point_in_segment; // note:
        }
        // velocity
        if (auto c = std::abs(v[j]) / vmax[j] - 1.0; // note: todo: remove 1.01 and use for validation only
            c > max_vel_constraints[j]) {
          max_vel_constraints[j] = c;
          max_vel_indices[j]     = point_in_segment; // note:
        }
        // acceleration
        if (auto c = std::abs(a[j]) / amax[j] - 1.0; // note: todo: remove 1.01 and use for validation only
            c > max_acc_constraints[j]) {
          max_acc_constraints[j] = c;
          max_acc_indices[j]     = point_in_segment; // note:
        }
        // torque
        if (auto c = std::abs(manip_data.efforts[j]) / tau_max[j] - 1.0; // note: todo: remove 1.01 and use for validation only
            c > max_tor_constraints[j]) {
          max_tor_constraints[j] = c;
          max_tor_indices[j]     = point_in_segment; // note:
        }
      }
    }

    // at this point we have max constraints for pos, vel, acc, tor, for this segment

    // fill in the constraints for the current segment
    // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
    auto fill_idx = segment * n_constraints_per_segment;
    std::copy_n(max_pos_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints; // note (andre): we can use the comma operator because we don't need the result of copy_n()
    std::copy_n(max_vel_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    std::copy_n(max_acc_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    std::copy_n(max_tor_constraints.data, n_joints, &constraints[fill_idx]);
  }
}

inline blast_fn void nlopt_constraints_final(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, x_len);

  Array constraints;
  constraints.alias(result, m);
  Matrix gradients;

  if (grad) {
    memset(grad, 0, m * x_len * sizeof(real));
    gradients.alias(grad, x_len, m);
  }

  constraints_and_gradients_final(xv, *opt, constraints, gradients);
}

inline void initialize_optimization_final(Optimization* opt) {
  // Constraints
  if (!opt->constraints.external_collisions) {
    opt->constraints.n_collision_constraints = 0;
  }

  // todo: swap INF for large value

  opt->lb        = Array(opt->x_len(), -HUGE_VAL);
  opt->ub        = Array(opt->x_len(), -HUGE_VAL);
  opt->lb.back() = 0.1;
  opt->ub.back() = 60.0;

  // Task
  auto task = opt->task;
  for (int i = 0; i < opt->manip.n_joints; i++) {
    while (std::abs(task(i, 0) - task(i, 3)) > PI) { // PI since the highest angle between two points on a circle is halfway through, so PI
      // Try updating start value
      real tmp_value = task(i, 0);
      if (task(i, 0) < task(i, 3)) {
        tmp_value += 2 * PI;
      } else {
        tmp_value -= 2 * PI;
      }
      if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
        task(i, 0) = tmp_value;
      } else {
        // Try updating end value
        tmp_value = task(i, 3);
        if (task(i, 3) < task(i, 0)) {
          tmp_value += 2 * PI;
        } else {
          tmp_value -= 2 * PI;
        }
        if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
          task(i, 0) = tmp_value;
        } else {
          break; // Nothing to be done
        }
      }
    }
  }
  opt->task = task;
}

inline void n_con_final(Optimization* opt) {
  opt->constraints.n_constraints = (opt->manip.n_joints * 4) * ((int) opt->bspline.n_ctrl - (int) opt->bspline.p);
}

inline Result optimize_final(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization_final(opt);
  n_con_final(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = 1000;
  stop.maxtime    = 5;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_final;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_finals, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 5);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 1000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
      // note: can we initialize grad to 0 here
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints_final(x, *opt, constraints_points);
      is_valid = max(constraints_points) < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con_final(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints_final(x, opt_val_more, constraints_more_points);
      is_valid_more = max(constraints_more_points) < opt->success_tolerance;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// acc 3 + extras (so acc 4 ?)
template<bool is_grad> // note: n_collision_skip must be 1 for this to work !!!
blast_fn void compute_constraints_acc4(double* result, Array& gradient_coeffs, Matrix& dtau_dp, Matrix& dtau_dv, Matrix& dtau_da, Matrix& dtcp_dp, Matrix& dtcp_dv, Matrix& dselfcol_dp, Matrix& dcol_dp, const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  ZoneScoped;
#endif

  double* moving_result = result;
  u32     grad_idx      = 0;

  constexpr real eps = 1e-5;

  auto joints = opt->manip.n_joints;

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // Lambda to process common bound constraint operations
  auto process_bound = [&](real value, real bound_max) {
    auto [constraint, gradient_coeff] = abs_constraint_dev<is_grad>(value, bound_max);
    *moving_result++                  = constraint;
    gradient_coeffs[grad_idx++]       = gradient_coeff;
  };

  for (u32 i = 0; i < opt->bspline.n_points; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    ManipulatorTempData manip_data;
    // Array               torque_constraint_plus(joints);

    auto pos = opt->bspline.traj.pos.col(i);
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    Assert(pos.is_alias);

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos);
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      for (int j = 0; j < joints; j++) {
        auto [constraint, gradient_coeff] = bound_constraint_dev<is_grad>(opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        *moving_result++                  = constraint;
        gradient_coeffs[grad_idx++]       = gradient_coeff;
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      for (int j = 0; j < joints; j++) {
        process_bound(opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      for (int j = 0; j < joints; j++) {
        process_bound(opt->bspline.traj.acc(j, i), opt->manip.amax[j]);
      }
    }

    Array torque_constraint(joints);
    real  tcp_constraint;
    real  self_collision_constraint;
    real  external_collisions_constraint;

    compute_capsules(opt->manip, manip_data);

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      dynamics(opt->manip, manip_data, vel, acc);
      for (int j = 0; j < joints; j++) {
        torque_constraint[j] = abs_constraint(manip_data.efforts[j], opt->manip.tau_max[j]);
        *moving_result++     = torque_constraint[j];
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * vel);
      tcp_constraint   = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
      *moving_result++ = tcp_constraint;
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      self_collision_constraint = max(-get_internal_collisions(opt->manip, manip_data));
      *moving_result++          = self_collision_constraint;
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisionsCalculate");
#endif
      external_collisions_constraint = -test_collisions_per_point(manip_data.capsule_list, &(opt->world));
      *moving_result++               = external_collisions_constraint;
    }

    if (is_grad) {
      Array pos_plus(joints);
      Array vel_plus(joints);
      Array acc_plus(joints);
      pos_plus = pos;
      vel_plus = vel;
      acc_plus = acc;

      // grad_coeffs pos
      for (u32 j = 0; j < joints; j++) {
        pos_plus[j] += eps;
        forward_kinematics(opt->manip, manip_data, pos_plus);

        if (opt->constraints.torque) {
          dynamics(opt->manip, manip_data, vel_plus, acc_plus);
          for (u32 k = 0; k < joints; k++) {
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
            dtau_dp(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        if (opt->constraints.tcp_speed) {
          auto J_tool_plus         = get_J_tool(opt, manip_data);
          real tcp_speed_plus      = norm(J_tool_plus * vel_plus);
          auto tcp_constraint_plus = bound_constraint(tcp_speed_plus, 0.0, opt->manip.tcp_max);
          dtcp_dp(j, i)            = (tcp_constraint_plus - tcp_constraint) / eps;
        }

        compute_capsules(opt->manip, manip_data);
        if (opt->constraints.self_collisions) {
          auto self_collision_constraint_plus = max(-get_internal_collisions(opt->manip, manip_data));
          dselfcol_dp(j, i)                   = (self_collision_constraint_plus - self_collision_constraint) / eps;
        }

        if (opt->constraints.external_collisions) {
          auto external_collisions_plus = -test_collisions_per_point(manip_data.capsule_list, &(opt->world));
          dcol_dp(j, i)                 = (external_collisions_plus - external_collisions_constraint) / eps;
        }

        pos_plus[j] = pos[j];
      }

      forward_kinematics(opt->manip, manip_data, pos_plus);
      // grad_coeffs vel
      for (u32 j = 0; j < joints; j++) {
        vel_plus[j] += eps;

        if (opt->constraints.torque) {
          dynamics(opt->manip, manip_data, vel_plus, acc_plus);
          for (u32 k = 0; k < joints; k++) {
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
            dtau_dv(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        if (opt->constraints.tcp_speed) {
          auto J_tool_plus         = get_J_tool(opt, manip_data);
          real tcp_speed_plus      = norm(J_tool_plus * vel_plus);
          auto tcp_constraint_plus = bound_constraint(tcp_speed_plus, 0.0, opt->manip.tcp_max);
          dtcp_dv(j, i)            = (tcp_constraint_plus - tcp_constraint) / eps;
        }

        vel_plus[j] = vel[j];
      }

      // grad_coeffs acc
      for (u32 j = 0; j < joints; j++) {
        acc_plus[j] += eps;

        if (opt->constraints.torque) {
          dynamics(opt->manip, manip_data, vel_plus, acc_plus);
          for (u32 k = 0; k < joints; k++) {
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
            dtau_da(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        acc_plus[j] = acc[j];
      }
    }


    //     if (opt->constraints.tcp_speed) {
    // #if BLAST_TRACE_LEVEL >= 3
    //       ZoneScopedN("TCPSpeed");
    // #endif
    //       auto vel = opt->bspline.traj.vel.col(i);
    //
    //       auto J_tool         = get_J_tool(opt, manip_data);
    //       real tcp_speed      = norm(J_tool * vel);
    //       auto tcp_constraint = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
    //       *moving_result++    = tcp_constraint;
    //
    //       if (is_grad) {
    //         Array pos_plus(joints);
    //         Array vel_plus(joints);
    //         pos_plus = pos;
    //         vel_plus = vel;
    //
    //         // grad_coeffs pos
    //         for (u32 j = 0; j < joints; j++) {
    //           pos_plus[j] += eps;
    //           forward_kinematics(opt->manip, manip_data, pos_plus);
    //
    //           auto J_tool_plus         = get_J_tool(opt, manip_data);
    //           real tcp_speed_plus      = norm(J_tool_plus * vel_plus);
    //           auto tcp_constraint_plus = bound_constraint(tcp_speed_plus, 0.0, opt->manip.tcp_max);
    //           dtcp_dp(j, i)            = (tcp_constraint_plus - tcp_constraint) / eps;
    //           pos_plus[j]              = pos[j];
    //         }
    //
    //         forward_kinematics(opt->manip, manip_data, pos_plus);
    //         // grad_coeffs vel
    //         for (u32 j = 0; j < joints; j++) {
    //           vel_plus[j] += eps;
    //
    //           auto J_tool_plus         = get_J_tool(opt, manip_data);
    //           real tcp_speed_plus      = norm(J_tool_plus * vel_plus);
    //           auto tcp_constraint_plus = bound_constraint(tcp_speed_plus, 0.0, opt->manip.tcp_max);
    //           dtcp_dv(j, i)            = (tcp_constraint_plus - tcp_constraint) / eps;
    //
    //           vel_plus[j] = vel[j];
    //         }
    //       }
    //     }
    //
    //     if (opt->constraints.self_collisions) {
    // #if BLAST_TRACE_LEVEL >= 3
    //       ZoneScopedN("SelfCollisions");
    // #endif
    //       compute_capsules(opt->manip, manip_data);
    //       auto self_collision_constraint = -max(get_internal_collisions(opt->manip, manip_data));
    //       *moving_result++               = self_collision_constraint;
    //
    //       if (is_grad) {
    //         Array pos_plus(joints);
    //
    //         for (u32 j = 0; j < joints; j++) {
    //           pos_plus[j] += eps;
    //           forward_kinematics(opt->manip, manip_data, pos_plus);
    //           compute_capsules(opt->manip, manip_data);
    //           auto self_collision_constraint_plus = -max(get_internal_collisions(opt->manip, manip_data));
    //
    //           dselfcol_dp(j, i) = (self_collision_constraint_plus - self_collision_constraint) / eps;
    //
    //           pos_plus[j] = pos[j];
    //         }
    //       }
    //     }
    //
    //     if (opt->constraints.external_collisions) {
    // #if BLAST_TRACE_LEVEL >= 3
    //       ZoneScopedN("ExternalCollisionsCalculate");
    // #endif
    //       compute_capsules(opt->manip, manip_data);
    //       auto external_collisions = -test_collisions_per_point(manip_data.capsule_list, &(opt->world));
    //       *moving_result++         = external_collisions;
    //       if (is_grad) {
    //
    //         Array pos_plus(joints);
    //
    //         for (u32 j = 0; j < joints; j++) {
    //           pos_plus[j] += eps;
    //           forward_kinematics(opt->manip, manip_data, pos_plus);
    //           compute_capsules(opt->manip, manip_data);
    //           auto external_collisions_plus = -test_collisions_per_point(manip_data.capsule_list, &(opt->world));
    //
    //           dcol_dp(j, i) = (external_collisions_plus - external_collisions) / eps;
    //
    //           pos_plus[j] = pos[j];
    //         }
    //       }
    //     }
  }

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("CustomConstraints");
#endif
    for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) {
      auto extra_constraint = opt->constraints.extra_constraints[i];
      extra_constraint(moving_result, opt);
      moving_result += opt->constraints.n_extra_constraints[i];
      grad_idx += opt->constraints.n_extra_constraints[i]; // todo: add analytical gradients
    }
  }
}

inline blast_fn void nlopt_constraints_acc4(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);

  int n_active_constraints = 0;
  if (opt->constraints.position)
    n_active_constraints++;
  if (opt->constraints.velocity)
    n_active_constraints++;
  if (opt->constraints.acceleration)
    n_active_constraints++;

  Array  gradient_coeffs(opt->manip.n_joints * n_active_constraints * opt->bspline.n_points); // todo: check performance
  Matrix dtau_dp(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix dtau_dv(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix dtau_da(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix dtcp_dp(opt->manip.n_joints, opt->bspline.n_points);
  Matrix dtcp_dv(opt->manip.n_joints, opt->bspline.n_points);
  Matrix dselfcol_dp(opt->manip.n_joints, opt->bspline.n_points);
  Matrix dcol_dp(opt->manip.n_joints, opt->bspline.n_points);
  if (!grad) {
    compute_constraints_acc4<false>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtcp_dp, dtcp_dv, dselfcol_dp, dcol_dp, xv, opt);
  }
  if (grad) {
    compute_constraints_acc4<true>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtcp_dp, dtcp_dv, dselfcol_dp, dcol_dp, xv, opt);
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps      = 1e-5;
    u32            n_con_lb = 0;
    Array          x_plus(xlen);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx       = 0;
    u32 constraint_idx = 0;
    u32 x_idx          = 0;

    for (u32 j = 0; j < xlen - 1; j++) { // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
      // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
      joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
      x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

      n_con_lb = ncon_lb_acc(opt, x_idx);                  // find the amount of constraints before the current point

      // todo: create alias matrix that points to grad
      // todo: can we change the order in which we store the gradients ?
      grad_idx       = n_con_lb * xlen + j;                                    // gradients are stored column-wise xlen * npoints
      constraint_idx = opt->manip.n_joints * n_active_constraints * opt->bspline.lb[x_idx];
      for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) { // lb & ub are inclusive
        grad_idx += joint * xlen;
        constraint_idx += joint;

        if (opt->constraints.position) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * gradient_coeffs[constraint_idx];
          grad_idx += opt->manip.n_joints * xlen; // increase index by the amount of joints * xlen
          constraint_idx += opt->manip.n_joints;
        }

        if (opt->constraints.velocity) {          // todo: basis_v / t once before this loop
          grad[grad_idx] = opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1] * gradient_coeffs[constraint_idx];
          grad_idx += opt->manip.n_joints * xlen; // increase index by the amount of joints * xlen
          constraint_idx += opt->manip.n_joints;
        }

        if (opt->constraints.acceleration) {                // todo: basis_a / t / t once before this loop
          grad[grad_idx] = opt->bspline.basis_a(x_idx, i) / (opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1] * opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1]) * gradient_coeffs[constraint_idx];
          grad_idx += (opt->manip.n_joints - joint) * xlen; // increase index by the amount of (joints - current joint) * xlen
          constraint_idx += (opt->manip.n_joints - joint);
        }

        if (opt->constraints.torque) {

          for (u32 k = 0; k < opt->manip.n_joints; k++) {
            grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dtau_dp(k, joint + opt->manip.n_joints * i) +
                             opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * dtau_dv(k, joint + opt->manip.n_joints * i) +
                             opt->bspline.basis_a(x_idx, i) / (opt->bspline.traj.t.back() * opt->bspline.traj.t.back()) * dtau_da(k, joint + opt->manip.n_joints * i);
            grad_idx += xlen;
            // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
          }
        }

        if (opt->constraints.tcp_speed) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dtcp_dp(joint, i) +
                           opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * dtcp_dv(joint, i);
          grad_idx += xlen;
          // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }

        if (opt->constraints.self_collisions) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dselfcol_dp(joint, i);
          grad_idx += xlen;
          // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }

        if (opt->constraints.external_collisions) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dcol_dp(joint, i);
          grad_idx += xlen;
          // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }
      }
    }

    {
      // last point T
      u32 j = xlen - 1;

      auto   n_joints             = opt->manip.n_joints;
      auto   constraint_per_point = m / opt->bspline.n_points;
      auto   one_over_T           = 1 / opt->bspline.traj.t.back();
      Matrix gradients;
      gradients.alias(grad, xlen, m);
      Array constraints;
      constraints.alias(result, m);

      for (int i = 1; i < opt->bspline.n_points; i++) {
        int    constraint_in_point_idx = 0;
        auto   vel                     = opt->bspline.traj.vel.col(i);
        auto   acc                     = opt->bspline.traj.acc.col(i);
        Matrix grad_point(&gradients(0, i * constraint_per_point), xlen, constraint_per_point);
        Array  constraint_point(&constraints[i * constraint_per_point], constraint_per_point);

        // dp/dT == 0
        constraint_in_point_idx += (int) n_joints;
        // dv/dT
        for (int k = 0; k < n_joints; k++) {
          grad_point.data[constraint_in_point_idx * xlen + j] = -(constraint_point[constraint_in_point_idx] + 1) * one_over_T;
          constraint_in_point_idx++;
        }
        // da_dT
        for (int k = 0; k < n_joints; k++) {
          grad_point.data[constraint_in_point_idx * xlen + j] = -2 * (constraint_point[constraint_in_point_idx] + 1) * one_over_T;
          constraint_in_point_idx++;
        }
        // dtau_dT
        for (int k = 0; k < n_joints; k++) {
          for (int l = 0; l < n_joints; l++) {
            grad_point.data[constraint_in_point_idx * xlen + j] += dtau_dv(k, l + n_joints * i) * (-vel[l] * one_over_T) + dtau_da(k, l + n_joints * i) * (-2 * acc[l] * one_over_T);
          }
          constraint_in_point_idx++;
        }
        // dtcp_dT
        {
          for (int k = 0; k < n_joints; k++) {
            grad_point.data[constraint_in_point_idx * xlen + j] += dtcp_dv(k, i) * (-vel[k] * one_over_T);
          }
          constraint_in_point_idx++; // unused, but added for uniformity
        }
      }
    }

    if (opt->constraints.show_info) { // when more info is needed per iteration
      Matrix gradients(xlen, m);
      Array  constraints(m);
      for (u32 j = 0; j < xlen; j++) {
        for (u32 i = 0; i < m; i++) {
          gradients(j, i) = grad[i * xlen + j];
          constraints[i]  = result[i];
        }
      }
      opt->constraints.grad_list.push_back(gradients);
      opt->constraints.constr_list.push_back(constraints);
      opt->constraints.x_list.push_back(xv);
    }
  }
}

inline Result optimize_acc4(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = 100000;
  stop.maxtime    = 30;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_acc4;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc4, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      auto max_con = max(constraints_points);
      // cout << "max_con = " << max_con << endl;
      is_valid = max_con < opt->success_tolerance * 2;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more = max(constraints_more_points);
      // cout << "max_con_more = " << max_con_more << endl;
      is_valid_more = max_con_more < opt->success_tolerance * 2;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// ------------------------------ Eval functions ------------------------------

void eval_function_Link6() {
  for (;;) {
    mut_config.lock();
    if (config_index >= config_list.size()) {
      mut_config.unlock();
      break;
    }
    auto config = config_list[config_index];
    std::cout << "Worker " << std::this_thread::get_id()
              << ", working on unit " << config_index + 1
              << " of " << config_list.size() << std::endl;
    config_index++;
    mut_config.unlock();

    auto  n_ctrl   = config.n_ctrl;
    auto  n_points = config.n_points;
    auto  n_optim  = config.n_optim;
    auto& task     = config.task;

    auto         manip = get_generic_Link6();
    Optimization opt(manip, task);
    Bspline      bspline(n_ctrl, n_points, 5, manip.n_joints);
    opt.bspline = bspline;

    opt.world = get_lab_world();

    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tcp_speed           = true;
    opt.constraints.self_collisions     = true;
    opt.constraints.external_collisions = true;

    opt.max_tries                    = 1;
    opt.success_tolerance            = 0.01;
    opt.constraints.n_collision_skip = 1;

    opt.guess.type = Guess::custom;


    // temp vectors
    std::vector<Result> tmp_result_list, tmp_result_acc1_list, tmp_result_acc3_list,
            tmp_result_segments_list;
    std::vector<u32> tmp_task_id, tmp_config_id;

    tmp_result_list.reserve(n_optim);
    tmp_result_acc1_list.reserve(n_optim);
    tmp_result_acc3_list.reserve(n_optim);
    tmp_result_segments_list.reserve(n_optim);
    tmp_task_id.reserve(n_optim);
    tmp_config_id.reserve(n_optim);

    for (int iter = 0; iter < n_optim; ++iter) {
      opt.guess.x0 = guess_random((opt.bspline), opt.task);
      // opt.guess.x0.back() = 0.5;
      // opt.guess.x0        = random_array(opt.bspline.x_len(opt.task), 1);
      // opt.guess.x0.back() = 0.5;

      auto result = optimize(&opt);
      tmp_result_list.push_back(result);

      auto result_acc1 = optimize_acc2(&opt);
      tmp_result_acc1_list.push_back(result_acc1);

      auto result_acc3 = optimize_acc4(&opt);
      tmp_result_acc3_list.push_back(result_acc3);

      // auto result_segments = optimize_final(&opt);
      auto result_segments = optimize_with_segments(&opt);
      tmp_result_segments_list.push_back(result_segments);

      tmp_task_id.push_back(config.task_idx);
      tmp_config_id.push_back(config.config_idx);
    }

    mut_result.lock();
    result_list.insert(result_list.end(), tmp_result_list.begin(), tmp_result_list.end());
    result_acc1_list.insert(result_acc1_list.end(), tmp_result_acc1_list.begin(), tmp_result_acc1_list.end());
    result_acc3_list.insert(result_acc3_list.end(), tmp_result_acc3_list.begin(), tmp_result_acc3_list.end());
    result_segments_list.insert(result_segments_list.end(), tmp_result_segments_list.begin(), tmp_result_segments_list.end());
    task_id.insert(task_id.end(), tmp_task_id.begin(), tmp_task_id.end());
    config_ids.insert(config_ids.end(), tmp_config_id.begin(), tmp_config_id.end());
    mut_result.unlock();
  }
}

void eval_function_UR5e() {
  for (;;) {
    mut_config.lock();
    if (config_index >= config_list.size()) {
      mut_config.unlock();
      break;
    }
    auto config = config_list[config_index];
    std::cout << "Worker " << std::this_thread::get_id()
              << ", working on unit " << config_index + 1
              << " of " << config_list.size() << std::endl;
    config_index++;
    mut_config.unlock();

    auto  n_ctrl   = config.n_ctrl;
    auto  n_points = config.n_points;
    auto  n_optim  = config.n_optim;
    auto& task     = config.task;

    auto manip   = get_generic_ur5e();
    manip.p_base = {-0.5, -0.3, 0.35};
    manip.Q_base = {-1, 0, 0, 0, -1, 0, 0, 0, 1};

    Optimization opt(manip, task);
    Bspline      bspline(n_ctrl, n_points, 5, manip.n_joints);
    opt.bspline = bspline;

    opt.world = get_kitchen_open_doors();

    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tcp_speed           = true;
    opt.constraints.self_collisions     = true;
    opt.constraints.external_collisions = true;

    opt.max_tries                    = 1;
    opt.success_tolerance            = 0.01;
    opt.constraints.n_collision_skip = 1;


    // temp vectors
    std::vector<Result> tmp_result_list, tmp_result_acc1_list, tmp_result_acc3_list,
            tmp_result_segments_list;
    std::vector<u32> tmp_task_id, tmp_config_id;

    tmp_result_list.reserve(n_optim);
    tmp_result_acc1_list.reserve(n_optim);
    tmp_result_acc3_list.reserve(n_optim);
    tmp_result_segments_list.reserve(n_optim);
    tmp_task_id.reserve(n_optim);
    tmp_config_id.reserve(n_optim);

    for (int iter = 0; iter < n_optim; ++iter) {
      opt.guess.x0 = guess_random((opt.bspline), opt.task);
      // opt.guess.x0.back() = 1.0;
      // opt.guess.x0        = random_array(opt.bspline.x_len(opt.task), 1);
      // opt.guess.x0.back() = 0.5;
      // opt.guess.type = Guess::random;
      initialize_optimization(&opt);
      n_con(&opt);
      auto x         = guess_shot_mean(&opt);
      opt.guess.type = Guess::custom;
      opt.guess.x0   = x;

      // cout << "Base" << endl;
      auto result = optimize(&opt);
      tmp_result_list.push_back(result);

      // cout << "acc1" << endl;
      auto result_acc1 = optimize_acc2(&opt);
      tmp_result_acc1_list.push_back(result_acc1);

      // cout << "acc2" << endl;
      auto result_acc3 = optimize_acc4(&opt);
      tmp_result_acc3_list.push_back(result_acc3);

      // auto result_segments = optimize_final(&opt);
      // cout << "Segments" << endl;
      auto result_segments = optimize_with_segments(&opt);
      tmp_result_segments_list.push_back(result_segments);

      tmp_task_id.push_back(config.task_idx);
      tmp_config_id.push_back(config.config_idx);
    }

    mut_result.lock();
    result_list.insert(result_list.end(), tmp_result_list.begin(), tmp_result_list.end());
    result_acc1_list.insert(result_acc1_list.end(), tmp_result_acc1_list.begin(), tmp_result_acc1_list.end());
    result_acc3_list.insert(result_acc3_list.end(), tmp_result_acc3_list.begin(), tmp_result_acc3_list.end());
    result_segments_list.insert(result_segments_list.end(), tmp_result_segments_list.begin(), tmp_result_segments_list.end());
    task_id.insert(task_id.end(), tmp_task_id.begin(), tmp_task_id.end());
    config_ids.insert(config_ids.end(), tmp_config_id.begin(), tmp_config_id.end());
    mut_result.unlock();
  }
}

void eval_function_Gen3() {
  for (;;) {
    mut_config.lock();
    if (config_index >= config_list.size()) {
      mut_config.unlock();
      break;
    }
    auto config = config_list[config_index];
    std::cout << "Worker " << std::this_thread::get_id()
              << ", working on unit " << config_index + 1
              << " of " << config_list.size() << std::endl;
    config_index++;
    mut_config.unlock();

    auto  n_ctrl   = config.n_ctrl;
    auto  n_points = config.n_points;
    auto  n_optim  = config.n_optim;
    auto& task     = config.task;

    auto manip   = get_generic_gen3();
    manip.p_base = {0.4, 0, 0.6};
    manip.Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    Optimization opt(manip, task);
    Bspline      bspline(n_ctrl, n_points, 5, manip.n_joints);
    opt.bspline = bspline;

    opt.world = get_bookshelf_thin();

    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tcp_speed           = true;
    opt.constraints.self_collisions     = true;
    opt.constraints.external_collisions = true;

    opt.max_tries                    = 1;
    opt.success_tolerance            = 0.01;
    opt.constraints.n_collision_skip = 1;


    // temp vectors
    std::vector<Result> tmp_result_list, tmp_result_acc1_list, tmp_result_acc3_list,
            tmp_result_segments_list;
    std::vector<u32> tmp_task_id, tmp_config_id;

    tmp_result_list.reserve(n_optim);
    tmp_result_acc1_list.reserve(n_optim);
    tmp_result_acc3_list.reserve(n_optim);
    tmp_result_segments_list.reserve(n_optim);
    tmp_task_id.reserve(n_optim);
    tmp_config_id.reserve(n_optim);

    for (int iter = 0; iter < n_optim; ++iter) {
      opt.guess.x0 = guess_random((opt.bspline), opt.task);
      // opt.guess.x0.back() = 1.0;
      // opt.guess.x0        = random_array(opt.bspline.x_len(opt.task), 1);
      // opt.guess.x0.back() = 0.5;
      // opt.guess.type = Guess::random;
      initialize_optimization(&opt);
      n_con(&opt);
      auto x         = guess_shot_mean(&opt);
      opt.guess.type = Guess::custom;
      opt.guess.x0   = x;

      // cout << "Base" << endl;
      auto result = optimize(&opt);
      tmp_result_list.push_back(result);

      // cout << "acc1" << endl;
      auto result_acc1 = optimize_acc2(&opt);
      tmp_result_acc1_list.push_back(result_acc1);

      // cout << "acc2" << endl;
      auto result_acc3 = optimize_acc4(&opt);
      tmp_result_acc3_list.push_back(result_acc3);

      // auto result_segments = optimize_final(&opt);
      // cout << "Segments" << endl;
      auto result_segments = optimize_with_segments(&opt);
      tmp_result_segments_list.push_back(result_segments);

      tmp_task_id.push_back(config.task_idx);
      tmp_config_id.push_back(config.config_idx);
    }

    mut_result.lock();
    result_list.insert(result_list.end(), tmp_result_list.begin(), tmp_result_list.end());
    result_acc1_list.insert(result_acc1_list.end(), tmp_result_acc1_list.begin(), tmp_result_acc1_list.end());
    result_acc3_list.insert(result_acc3_list.end(), tmp_result_acc3_list.begin(), tmp_result_acc3_list.end());
    result_segments_list.insert(result_segments_list.end(), tmp_result_segments_list.begin(), tmp_result_segments_list.end());
    task_id.insert(task_id.end(), tmp_task_id.begin(), tmp_task_id.end());
    config_ids.insert(config_ids.end(), tmp_config_id.begin(), tmp_config_id.end());
    mut_result.unlock();
  }
}

int main() {
  using namespace blast;

  // fill_config_list(config_list); note : uses Link6
  // fill_config_list_UR5e(config_list);
  fill_config_list_Gen3(config_list);
  result_list.reserve(config_list.size() * _n_optim);
  result_acc1_list.reserve(config_list.size() * _n_optim);
  // result_acc2_list.reserve(config_list.size() * _n_optim);
  result_acc3_list.reserve(config_list.size() * _n_optim);
  result_segments_list.reserve(config_list.size() * _n_optim);
  task_id.reserve(config_list.size() * _n_optim);
  config_ids.reserve(config_list.size() * _n_optim);

  raise_process_priority();

  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  auto cores = pick_physical_cores(thread_count);
  for (int t = 0; t < thread_count; ++t) {
    workers.emplace_back([&, t] {
      configure_current_thread_for_performance();

      // eval_function_Link6();
      eval_function_UR5e();
      // eval_function_Gen3();
    });

    // Pin the thread immediately after creation
    if (t < (int) cores.size()) {
      SetThreadGroupAffinity(
              (HANDLE) workers.back().native_handle(),
              &cores[t],
              nullptr);
    }
  }

  for (auto& th: workers)
    th.join();

  cout << "Workers are done" << endl;

  cout << "Crunching the numbers" << endl;

  const std::vector<VariantGroup> variants = {
          {"base", &result_list},
          {"acc1", &result_acc1_list},
          {"acc2", &result_acc3_list},
          {"segments", &result_segments_list},
  };

  write_results_json_simple(
          "results_all_Gen3_new.json",
          config_ids,
          task_id,
          config_list,
          variants);


  return 0;
}
