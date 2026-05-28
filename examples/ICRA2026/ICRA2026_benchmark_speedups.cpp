#define BLAST_TRACE_LEVEL 0
#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>
#include <mutex>
#include <thread>

#include <cassert>
#include <fstream>
#include <stdexcept>

using namespace blast;

constexpr int thread_count = 8;
// ----------------------------- UR5e tasks + world -----------------------------
inline std::vector<Task> get_UR5e_tasks() {
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

  std::vector<Task> task_list;
  task_list.reserve(64);
  for (auto& s: start_pos)
    for (auto& e: end_pos)
      task_list.push_back(Task::stop_to_stop(s, e));
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
  int  n_ctrl     = 0;
  int  n_points   = 0;
  int  n_optim    = 0;
  int  task_idx   = 0;
  int  config_idx = 0;
  Task task;
};

constexpr int                   config_size = 64; // for UR5e 64 tasks
std::array<Config, config_size> config_list;
constexpr int                   _n_optim     = 100;
int                             config_index = 0;
std::mutex                      mut_config;

// Base results
std::vector<Result> result_base_list;
// Results for each of the four speedups
std::vector<Result> result_spdup1_list;
std::vector<Result> result_spdup2_list;
// Results for our final per_segment integration
std::vector<Result> result_spdup3_list;
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
  const char*                label; // "base" | "spdup1" | "spdup2" | "spdup3"
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
  std::vector<Task>                   tasks               = get_UR5e_tasks();

  int config_idx = 0;
  int task_idx   = 0;
  // auto task       = tasks[0];
  for (auto& task: tasks) {
    for (auto& bspline_config: bspline_config_list) {
      auto [n_ctrl, n_points_per_segments] = bspline_config;
      auto n_points                        = n_points_per_segments * (n_ctrl - 5);

      target_config_list[config_idx].n_ctrl     = n_ctrl;
      target_config_list[config_idx].n_points   = n_points;
      target_config_list[config_idx].n_optim    = _n_optim;
      target_config_list[config_idx].task_idx   = task_idx;
      target_config_list[config_idx].config_idx = config_idx;
      target_config_list[config_idx++].task     = task;
    }
    task_idx++;
  }
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

    auto manip          = make_UR5e();
    manip.base_position = {-0.5, -0.3, 0.35};
    manip.base_rotation = {-1, 0, 0, 0, -1, 0, 0, 0, 1};

    Optimization opt(manip, task);
    Bspline      bspline(n_ctrl, n_points, 5, manip.n_joints);
    opt.bspline = bspline;

    opt.world = get_kitchen_open_doors();

    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tool_speed          = true;
    opt.constraints.self_collisions     = true;
    opt.constraints.external_collisions = true;

    opt.max_tries                    = 1;
    opt.success_tolerance            = 0.01;
    opt.constraints.n_collision_skip = 1;

    opt.guess.type = Guess::custom;

    // temp vectors
    std::vector<Result> tmp_result_list, tmp_result_spdup1_list, tmp_result_spdup2_list, tmp_result_spdup3_list;
    std::vector<u32>    tmp_task_id, tmp_config_id;

    tmp_result_list.reserve(n_optim);
    tmp_result_spdup1_list.reserve(n_optim);
    tmp_result_spdup2_list.reserve(n_optim);
    tmp_result_spdup3_list.reserve(n_optim);
    tmp_task_id.reserve(n_optim);
    tmp_config_id.reserve(n_optim);

    for (int iter = 0; iter < n_optim; ++iter) {
      opt.guess.initial_x = guess_random((opt.bspline), opt.task);
      // opt.guess.initial_x.back() = 0.5;
      auto result_base = optimize(&opt, blast::OptimizationMethod::baseline);
      tmp_result_list.push_back(result_base);

      auto result_spdup1 = optimize(&opt, blast::OptimizationMethod::with_analytical_pva);
      tmp_result_spdup1_list.push_back(result_spdup1);

      auto result_spdup2 = optimize(&opt, blast::OptimizationMethod::with_analytical_dynamics);
      tmp_result_spdup2_list.push_back(result_spdup2);

      auto result_spdup3 = optimize(&opt, blast::OptimizationMethod::with_segments);
      tmp_result_spdup3_list.push_back(result_spdup3);

      tmp_task_id.push_back(config.task_idx);
      tmp_config_id.push_back(config.config_idx);
    }

    mut_result.lock();
    result_base_list.insert(result_base_list.end(), tmp_result_list.begin(), tmp_result_list.end());
    result_spdup1_list.insert(result_spdup1_list.end(), tmp_result_spdup1_list.begin(), tmp_result_spdup1_list.end());
    result_spdup2_list.insert(result_spdup2_list.end(), tmp_result_spdup2_list.begin(), tmp_result_spdup2_list.end());
    result_spdup3_list.insert(result_spdup3_list.end(), tmp_result_spdup3_list.begin(), tmp_result_spdup3_list.end());
    task_id.insert(task_id.end(), tmp_task_id.begin(), tmp_task_id.end());
    config_ids.insert(config_ids.end(), tmp_config_id.begin(), tmp_config_id.end());
    mut_result.unlock();
  }
}

int main() {
  using namespace blast;

  fill_config_list_UR5e(config_list);

  result_base_list.reserve(config_list.size() * _n_optim);
  result_spdup1_list.reserve(config_list.size() * _n_optim);
  result_spdup2_list.reserve(config_list.size() * _n_optim);
  result_spdup3_list.reserve(config_list.size() * _n_optim);
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

#ifdef _WIN32
    // Only available on Windows — pin threads to P-cores
    if (t < (int) cores.size()) {
      SetThreadGroupAffinity(
              (HANDLE) workers.back().native_handle(),
              &cores[t],
              nullptr);
    }
#endif
  }


  for (auto& th: workers)
    th.join();

  cout << "Workers are done" << endl;

  const std::vector<VariantGroup> variants = {
          {"base", &result_base_list},
          {"spdup1", &result_spdup1_list},
          {"spdup2", &result_spdup2_list},
          {"spdup3", &result_spdup3_list},
  };

  write_results_json_simple(
          "../../../../examples/ICRA2026/results_tmp.json",
          config_ids,
          task_id,
          config_list,
          variants);


  return 0;
}
