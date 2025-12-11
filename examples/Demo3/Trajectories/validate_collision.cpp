
#define BLAST_TRACE_LEVEL 0
#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"

using namespace blast;

int main() {
  auto link6   = get_generic_Link6_fixed();
  link6.p_base = {0, 0, 0.833}; // link6
  link6.Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  auto world   = get_demo2_world();

  auto inflated_gen3   = get_generic_gen3_fixed();
  inflated_gen3.p_base = {1.28, 0.025, 0.73};
  inflated_gen3.Q_base = {-1, 0, 0, 0, -1, 0, 0, 0, 1};
//   for (auto& caps: inflated_gen3._collision_model) { // todo: fix std::array
//     caps.r += 0.05; // 5 cm buffer
//   }

  auto trajectory_gen3 = read_csv_trajectory_no_header("../../../examples/Demo3/Trajectories/trajectory_full_gen3.csv", ",");
  auto gen3_trajectory_pos = trajectory_gen3.pos;
  
  auto trajectory_link6 = read_csv_trajectory_no_header("../../../examples/Demo3/Trajectories/trajectory_full_link6.csv", ",");
  auto link6_trajectory_pos = trajectory_link6.pos;

  // Add padding
  int max_pt = std::max(link6_trajectory_pos.cols, gen3_trajectory_pos.cols);
  bool gen3_is_max = max_pt == gen3_trajectory_pos.cols;

  if (gen3_is_max) {
    int original_size = link6_trajectory_pos.cols;
    link6_trajectory_pos.resize(link6_trajectory_pos.rows, max_pt);
    for (int i = original_size; i < max_pt; i++) {
        for (int j = 0; j < link6_trajectory_pos.rows; j++) {
            link6_trajectory_pos(j, i) = link6_trajectory_pos(j, original_size);
        }
    }
  } else {
    int original_size = gen3_trajectory_pos.cols;
    gen3_trajectory_pos.resize(gen3_trajectory_pos.rows, max_pt);
    for (int i = original_size; i < max_pt; i++) {
        for (int j = 0; j < gen3_trajectory_pos.rows; j++) {
            gen3_trajectory_pos(j, i) = gen3_trajectory_pos(j, original_size);
        }
    }
  }

  // test collisions
  ManipulatorTempData gen3_tmp_data;
  ManipulatorTempData link6_tmp_data;

  Array dist_min(max_pt, INF_REAL);
  
  for (int i = 0; i < max_pt; i++) {
    forward_kinematics(link6, link6_tmp_data, link6_trajectory_pos.col(i));
    compute_capsules(link6, link6_tmp_data);
    
    if (i == 6386) {
        int stop = 0;
    }
    forward_kinematics(inflated_gen3, gen3_tmp_data, gen3_trajectory_pos.col(i));
    compute_capsules(inflated_gen3, gen3_tmp_data);


    for (int j = 0; j < link6._n_caps; j++) {
        for (int k = 0; k < inflated_gen3._n_caps; k++) {
            real dist = distance(link6_tmp_data.capsule_list[j], gen3_tmp_data.capsule_list[k]);
            if (dist < dist_min[i]) {
                dist_min[i] = dist;
            }
        }
    }
  }

  // output results
  std::cout << "Minimal distance is d = " << min(dist_min) << ", found at position " << argmin(dist_min) << std::endl;
}