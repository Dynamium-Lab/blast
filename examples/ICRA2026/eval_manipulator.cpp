//
// Created by nikos on 2025-09-04.
//
#define BLAST_TRACE_LEVEL 0
#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

#include "../../tests/test_helper/test_helper.hpp"

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
  limits.tcp_max = 2.0; // todo: verify this value

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

// ============================================================================
// Cross-platform performance/thread helpers
// ============================================================================
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00
#endif
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <vector>
#include <windows.h>

static std::vector<GROUP_AFFINITY> pick_physical_cores(size_t N) {
  std::vector<GROUP_AFFINITY> picks;

  DWORD len = 0;
  GetSystemCpuSetInformation(nullptr, 0, &len, nullptr, 0);
  if (GetLastError() == ERROR_INSUFFICIENT_BUFFER && len > 0) {
    std::vector<unsigned char> buf(len);
    if (GetSystemCpuSetInformation(reinterpret_cast<SYSTEM_CPU_SET_INFORMATION*>(buf.data()),
                                   (ULONG) buf.size(), &len, nullptr, 0)) {
      std::unordered_set<ULONG> seen_core;
      BYTE                      min_eff = 255, max_eff = 0;
      struct Row {
        GROUP_AFFINITY ga;
        ULONG          core;
        BYTE           eff;
      };
      std::vector<Row> rows;

      BYTE* p   = buf.data();
      BYTE* end = p + len;
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
            min_eff = std::min(min_eff, r.eff);
            max_eff = std::max(max_eff, r.eff);
          }
        }
        p += info->Size;
      }

      if (!rows.empty()) {
        bool hybrid = (min_eff != max_eff);
        std::sort(rows.begin(), rows.end(), [](auto& a, auto& b) { return a.core < b.core; });
        for (auto& r: rows) {
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

  // fallback if empty
  if (picks.empty()) {
    DWORD bytes = 0;
    GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &bytes);
    std::vector<unsigned char> buf2(bytes);
    if (GetLogicalProcessorInformationEx(RelationProcessorCore,
                                         reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buf2.data()), &bytes)) {
      BYTE* p   = buf2.data();
      BYTE* end = p + bytes;
      while (p < end) {
        auto* ex = reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(p);
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
        p += ex->Size;
      }
    }
  }

  return picks;
}

inline void raise_process_priority() {
  SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
}

inline void configure_current_thread_for_performance() {
  THREAD_POWER_THROTTLING_STATE s{};
  s.Version     = THREAD_POWER_THROTTLING_CURRENT_VERSION;
  s.ControlMask = THREAD_POWER_THROTTLING_EXECUTION_SPEED;
  s.StateMask   = 0; // favor performance
  SetThreadInformation(GetCurrentThread(), ThreadPowerThrottling, &s, sizeof(s));
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
}

#else // =============================== LINUX =================================

#include <cstring>
#include <iostream>
#include <pthread.h>
#include <sched.h>
#include <sys/resource.h>
#include <unistd.h>
#include <vector>

struct GROUP_AFFINITY {}; // dummy placeholder for portability

// Simple round-robin CPU pinning on Linux
static std::vector<int> pick_physical_cores(size_t N) {
  std::vector<int> cores;
  int              ncpu = std::thread::hardware_concurrency();
  for (int i = 0; i < (int) N && i < ncpu; ++i)
    cores.push_back(i);
  return cores;
}

inline void raise_process_priority() {
  if (setpriority(PRIO_PROCESS, 0, -10) != 0)
    std::cerr << "Warning: cannot raise process priority (" << strerror(errno) << ")\n";
}

inline void configure_current_thread_for_performance() {
  pthread_t this_thread = pthread_self();
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);

  static thread_local int next_core = 0;
  int                     ncpu      = std::thread::hardware_concurrency();
  int                     core_id   = next_core++ % ncpu;

  CPU_SET(core_id, &cpuset);
  if (pthread_setaffinity_np(this_thread, sizeof(cpu_set_t), &cpuset) != 0)
    std::cerr << "Warning: cannot set CPU affinity (" << strerror(errno) << ")\n";

  // Optional: try to lower nice value further
  setpriority(PRIO_PROCESS, 0, -10);
}

#endif
// ============================================================================


// ===== End helpers =====

struct Config {
  int    n_ctrl     = 0;
  int    n_points   = 0;
  int    n_optim    = 0;
  int    task_idx   = 0;
  int    config_idx = 0;
  Matrix task;
};
constexpr int                   config_size = 64; // for UR5e 64 tasks
std::array<Config, config_size> config_list;
constexpr int                   _n_optim     = 1;
int                             config_index = 0;
std::mutex                      mut_config;

// Base results
std::vector<Result> result_base_list;
// Results for each of the four accelerations
std::vector<Result> result_acc1_list;
std::vector<Result> result_acc2_list;
// Results for our final per_segment integration
std::vector<Result> result_acc3_list;
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

  cout << "Crunching the numbers" << endl;

  for (const auto& v: variants) {
    const auto& list = *v.list;
    cout << "Variant type : " << v.label << " variant size : " << list.size() << endl;
    for (size_t i = 0; i < list.size(); ++i) {
      const Config& cfg = target_config_list[target_config_ids[i]];
      const Result& r   = list[i];

      json row;
      row["variant"]                          = v.label;
      row["task_idx"]                         = target_task_ids[i];
      row["n_ctrl"]                           = cfg.n_ctrl;
      row["n_points"]                         = cfg.n_points;
      row["success"]                          = r.success;
      row["success_false"]                    = r.success_false;
      row["max_constraint_value"]             = r.max_constraint_value;
      row["max_constraint_idx"]               = r.max_constraint_idx;
      row["max_constraint_more_points_value"] = r.max_constraint_more_points_value;
      row["max_constraint_more_points_idx"]   = r.max_constraint_more_points_idx;
      row["nlopt_exit_criteria"]              = r.nlopt_exit_criteria;
      row["num_eval"]                         = r.num_eval;
      row["compute_time"]                     = r.compute_time;
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
inline void fill_config_list_UR5e(std::array<Config, config_size>& target_config_list) {
  // These give rounded number of n_points_per_segment (With higher n_ctrl, we can allow for fewer points per segments)
  // std::array<std::tuple<int, int>, 3> bspline_config_list = {std::make_tuple(12, 15), // 15 points per segment
  //                                                            std::make_tuple(16, 10), // 10 points per segments
  //                                                            std::make_tuple(20, 10)};
  std::array<std::tuple<int, int>, 1> bspline_config_list = {std::make_tuple(16, 10)};
  std::vector<Matrix>                 tasks               = get_UR5e_tasks();

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

  if (opt->constraints.self_collisions) {
    n_constraints += n_points;
  }
  if (opt->constraints.external_collisions)
    n_constraints += n_points * opt->manip._n_caps;

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
  Matrix collision_constraint;
  Array  self_collision_constraint;

  ConstraintPerPoint(int joints, int points, int n_capsules) {
    pos_constraint.resize(points);
    vel_constraint.resize(points);
    acc_constraint.resize(points);
    tor_constraint.resize(joints, points);
    tcp_constraint.resize(points);
    self_collision_constraint.resize(points);
    collision_constraint.resize(n_capsules, points);
  }
  ~ConstraintPerPoint() = default;
};

inline void compute_constraints_grad1(ConstraintPerPoint& constraints, const Array& x, const u32 x_idx, const u32 joint_idx, Optimization* opt) {

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
      // auto collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
      Array max_col_constraints(opt->manip._n_caps, -INF_REAL);
      for (int capsule_id = 0; capsule_id < opt->manip._n_caps; capsule_id++) {
        real       dist_min = INF_REAL;
        const auto capsule  = manip_data.capsule_list[capsule_id];

        // check against boxes
        int count = 0;
        for (const auto& box: opt->world.boxes) {
          if (const auto dist = distance(capsule, box);
              dist < dist_min) {
            dist_min = dist;
          }
          count++;
        }

        // check against capsules
        count = 0;
        for (const auto caps: opt->world.capsules) {
          if (const auto dist = distance(capsule, caps);
              dist < dist_min) {
            dist_min = dist;
          }
          count++;
        }

        // check against spheres
        count = 0;
        for (const auto sphere: opt->world.spheres) {
          if (const auto dist = distance(capsule, sphere);
              dist < dist_min) {
            dist_min = dist;
          }
          count++;
        }

        dist_min = -dist_min; // negative distance is positive constraint

        // update worst position for the current capsule if necessary
        if (dist_min > max_col_constraints[capsule_id]) {
          max_col_constraints[capsule_id] = dist_min;
        }
      }
      for (int k = 0; k < opt->manip._n_caps; k++)
        constraints.collision_constraint(k, i - opt->bspline.lb[x_idx]) = max_col_constraints[k];
    }
  }
}

inline void nlopt_constraints_acc1(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);
  {

#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Constraints");
#endif
    compute_constraints(result, xv, opt);
  }

  // gradients calculation
  if (grad) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Grad");
#endif
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps = 1e-5;

    u32   n_con_lb = 0;
    Array x_plus(xlen);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx = 0;
    u32 x_idx    = 0;

    auto n_joints   = opt->manip.n_joints;
    auto n_capsules = opt->manip._n_caps;

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
      ConstraintPerPoint constraint(n_joints, n_points, n_capsules);
      compute_constraints_grad1(constraint, x_plus, x_idx, joint, opt);

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
        for (int k = 0; k < opt->manip._n_caps; k++) {
          grad[grad_idx] = (constraint.collision_constraint(k, i - opt->bspline.lb[x_idx]) - result[i * n_con_per_point + 4 * n_joints + 2 + k]) / eps;
          grad_idx += xlen;
        }
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

inline Result optimize_acc1(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  // initialize_optimization(opt);
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
  ub.back() = 30.0;
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
  nlopt_res = nlopt_set_maxtime(o, 30);
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
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance * 2;
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
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance * 2;

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
template<bool is_grad> // note: n_collision_skip must be 1 for this to work !!!
blast_fn void compute_constraints_acc2(double* result, Array& gradient_coeffs, Matrix& dtau_dp, Matrix& dtau_dv, Matrix& dtau_da, Matrix& dtcp_dp, Matrix& dtcp_dv, Matrix& dselfcol_dp, std::vector<Matrix>& dcol_dp, const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  ZoneScoped;
#endif

  enum class CollisionObjectType {
    box,
    sphere,
    capsule,
  };
  struct CollisionEntities {
    // object in the world
    CollisionObjectType other_object_type = CollisionObjectType::box;
    union {
      Box     box{};
      Sphere  sphere;
      Capsule capsule;
    };

    // which point in time to look up the capsule
    int point_in_segment = 0;
  };

  double* moving_result = result;
  u32     grad_idx      = 0;

  constexpr real eps = 1e-5;

  auto joints     = opt->manip.n_joints;
  auto n_capsules = opt->manip._n_caps;

  std::array<CollisionEntities, MAX_CAPSULES> max_collision_entities{};

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
#if BLAST_TRACE_LEVEL >= 2
    ZoneScopedN("ConstraintSinglePoint");
#endif

    ManipulatorTempData manip_data;
    // Array               torque_constraint_plus(joints);

    auto pos = opt->bspline.traj.pos.col(i);
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    Assert(pos.is_alias);

    Array torque_constraint(joints);
    real  tcp_constraint;
    real  self_collision_constraint;
    Array max_col_constraints(n_capsules, -INF_REAL);

    {
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("Constraints");
#endif


      {
#if BLAST_TRACE_LEVEL >= 3
        ZoneScopedN("Constraints per point");
#endif


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
          // check every capsule with world
          for (int capsule_id = 0; capsule_id < n_capsules; capsule_id++) {
            real       dist_min = INF_REAL;
            const auto capsule  = manip_data.capsule_list[capsule_id];

            CollisionEntities collision_objects{};

            // check against boxes
            int count = 0;
            for (const auto& box: opt->world.boxes) {
              if (const auto dist = distance(capsule, box);
                  dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::box;
                collision_objects.box               = box;
                collision_objects.point_in_segment  = i;
              }
              count++;
            }

            // check against capsules
            count = 0;
            for (const auto caps: opt->world.capsules) {
              if (const auto dist = distance(capsule, caps);
                  dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::capsule;
                collision_objects.capsule           = capsule;
                collision_objects.point_in_segment  = i;
              }
              count++;
            }

            // check against spheres
            count = 0;
            for (const auto sphere: opt->world.spheres) {
              if (const auto dist = distance(capsule, sphere);
                  dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::sphere;
                collision_objects.sphere            = sphere;
                collision_objects.point_in_segment  = i;
              }
              count++;
            }

            dist_min = -dist_min; // negative distance is positive constraint

            // update worst position for the current capsule if necessary
            if (dist_min > max_col_constraints[capsule_id]) {
              max_col_constraints[capsule_id]    = dist_min;
              max_collision_entities[capsule_id] = collision_objects;
            }
          }
          for (int k = 0; k < n_capsules; k++)
            *moving_result++ = max_col_constraints[k];
        }
      }
    }
    if (is_grad) {
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("Grad");
#endif
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
          for (int capsule_id = 0; capsule_id < n_capsules; capsule_id++) {
            const auto  capsule = manip_data.capsule_list[capsule_id];
            real        distance_plus;
            const auto& objects = max_collision_entities[capsule_id];
            switch (objects.other_object_type) {
              case CollisionObjectType::box: {
                distance_plus = distance(capsule, objects.box);
                break;
              }
              case CollisionObjectType::capsule: {
                distance_plus = distance(capsule, objects.capsule);
                break;
              }
              case CollisionObjectType::sphere: {
                distance_plus = distance(capsule, objects.sphere);
                break;
              }
            }
            distance_plus = -distance_plus;
            // auto external_collisions_plus = -test_collisions_per_point(manip_data.capsule_list, &(opt->world));
            dcol_dp[i].resize(n_capsules, joints);
            dcol_dp[i](capsule_id, j) = (distance_plus - max_col_constraints[capsule_id]) / eps;
          }
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

inline blast_fn void nlopt_constraints_acc2(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
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

  Array               gradient_coeffs(opt->manip.n_joints * n_active_constraints * opt->bspline.n_points); // todo: check performance
  Matrix              dtau_dp(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix              dtau_dv(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix              dtau_da(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix              dtcp_dp(opt->manip.n_joints, opt->bspline.n_points);
  Matrix              dtcp_dv(opt->manip.n_joints, opt->bspline.n_points);
  Matrix              dselfcol_dp(opt->manip.n_joints, opt->bspline.n_points);
  std::vector<Matrix> dcol_dp; // (n_caps, n_joints)
  dcol_dp.resize(opt->bspline.n_points);
  // (opt->manip.n_joints, opt->bspline.n_points);
  if (!grad) {

    compute_constraints_acc2<false>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtcp_dp, dtcp_dv, dselfcol_dp, dcol_dp, xv, opt);
  }
  if (grad) {
    compute_constraints_acc2<true>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtcp_dp, dtcp_dv, dselfcol_dp, dcol_dp, xv, opt);
    {
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("Grad Fill");
#endif
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
            for (int k = 0; k < opt->manip._n_caps; k++) {
              grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dcol_dp[i](k, joint);
              grad_idx += xlen;
            }
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
}

inline Result optimize_acc2(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  // initialize_optimization(opt);
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
  ub.back() = 30.0;
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
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance * 2;
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
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance * 2;

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

    opt.guess.type = Guess::custom;

    // temp vectors
    std::vector<Result> tmp_result_list, tmp_result_acc1_list, tmp_result_acc2_list, tmp_result_acc3_list;
    std::vector<u32>    tmp_task_id, tmp_config_id;

    tmp_result_list.reserve(n_optim);
    tmp_result_acc1_list.reserve(n_optim);
    tmp_result_acc2_list.reserve(n_optim);
    tmp_result_acc3_list.reserve(n_optim);
    tmp_task_id.reserve(n_optim);
    tmp_config_id.reserve(n_optim);

    for (int iter = 0; iter < n_optim; ++iter) {
      opt.guess.x0 = guess_random((opt.bspline), opt.task);
      // opt.guess.x0.back() = 0.5;
      auto result_base = optimize(&opt);
      tmp_result_list.push_back(result_base);

      auto result_acc1 = optimize_acc1(&opt);
      tmp_result_acc1_list.push_back(result_acc1);

      auto result_acc2 = optimize_acc2(&opt);
      tmp_result_acc2_list.push_back(result_acc2);

      auto result_acc3 = optimize_with_segments(&opt);
      tmp_result_acc3_list.push_back(result_acc3);

      tmp_task_id.push_back(config.task_idx);
      tmp_config_id.push_back(config.config_idx);
    }

    mut_result.lock();
    result_base_list.insert(result_base_list.end(), tmp_result_list.begin(), tmp_result_list.end());
    result_acc1_list.insert(result_acc1_list.end(), tmp_result_acc1_list.begin(), tmp_result_acc1_list.end());
    result_acc2_list.insert(result_acc2_list.end(), tmp_result_acc2_list.begin(), tmp_result_acc2_list.end());
    result_acc3_list.insert(result_acc3_list.end(), tmp_result_acc3_list.begin(), tmp_result_acc3_list.end());
    task_id.insert(task_id.end(), tmp_task_id.begin(), tmp_task_id.end());
    config_ids.insert(config_ids.end(), tmp_config_id.begin(), tmp_config_id.end());
    mut_result.unlock();
  }
}

int main() {
  using namespace blast;

  fill_config_list_UR5e(config_list);

  result_base_list.reserve(config_list.size() * _n_optim);
  result_acc1_list.reserve(config_list.size() * _n_optim);
  result_acc2_list.reserve(config_list.size() * _n_optim);
  result_acc3_list.reserve(config_list.size() * _n_optim);
  task_id.reserve(config_list.size() * _n_optim);
  config_ids.reserve(config_list.size() * _n_optim);

  raise_process_priority();

  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  auto cores = pick_physical_cores(thread_count);
  for (int t = 0; t < thread_count; ++t) {
    workers.emplace_back([&, t] {
      configure_current_thread_for_performance();
      eval_function_UR5e();
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

  const std::vector<VariantGroup> variants = {
          {"base", &result_base_list},
          {"acc1", &result_acc1_list},
          {"acc2", &result_acc2_list},
          {"segments", &result_acc3_list},
  };

  write_results_json_simple(
          "../../../examples/ICRA2026/results_16x110_tmp.json",
          config_ids,
          task_id,
          config_list,
          variants);


  return 0;
}
