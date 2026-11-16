#pragma once
#include <blast>

namespace blast {
// Extracted scenes from MotionBenchMaker (github : https://github.com/KavrakiLab/motion_bench_maker)

inline host_fn World get_bookshelf_small() {
  World bookshelf_small;
  bookshelf_small.add_capsule({(real) 0.9, 0, (real) 1.01}, {(real) 0.9, 0, (real) 1.15}, (real) 0.03);
  bookshelf_small.add_capsule({(real) 0.7, 0, (real) 1.01}, {(real) 0.7, 0, (real) 1.15}, (real) 0.03);
  bookshelf_small.add_capsule({(real) 0.5, 0, (real) 1.01}, {(real) 0.5, 0, (real) 1.15}, (real) 0.03);
  bookshelf_small.add_box({1, 0, 1}, {(real) 0.6, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_small.add_box({1, -(real) 0.5, (real) 1.15}, {(real) 0.6, (real) 0.02, (real) 0.17}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_small.add_box({1, (real) 0.5, (real) 1.15}, {(real) 0.6, (real) 0.02, (real) 0.17}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_small.add_box({1, 0, (real) 1.3}, {(real) 0.6, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  return bookshelf_small;
}

inline host_fn World get_bookshelf_tall() {
  World bookshelf_tall;
  bookshelf_tall.add_capsule({(real) 0.9, 0, (real) 1.31}, {(real) 0.9, 0, (real) 1.45}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.7, 0, (real) 1.31}, {(real) 0.7, 0, (real) 1.45}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.5, 0, (real) 1.31}, {(real) 0.5, 0, (real) 1.45}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.9, 0, (real) 0.71}, {(real) 0.9, 0, (real) 0.85}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.7, 0, (real) 0.71}, {(real) 0.7, 0, (real) 0.85}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.5, 0, (real) 0.71}, {(real) 0.5, 0, (real) 0.85}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.9, 0, (real) 1.01}, {(real) 0.9, 0, (real) 1.15}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.7, 0, (real) 1.01}, {(real) 0.7, 0, (real) 1.15}, (real) 0.03);
  bookshelf_tall.add_capsule({(real) 0.5, 0, (real) 1.01}, {(real) 0.5, 0, (real) 1.15}, (real) 0.03);
  bookshelf_tall.add_box({1, 0, (real) 0.7}, {(real) 0.6, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_tall.add_box({1, 0, 1}, {(real) 0.6, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_tall.add_box({1, 0, (real) 1.3}, {(real) 0.6, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_tall.add_box({1, 0, (real) 1.6}, {(real) 0.6, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_tall.add_box({1, -(real) 0.5, (real) 0.82}, {(real) 0.65, (real) 0.02, (real) 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_tall.add_box({1, (real) 0.5, (real) 0.82}, {(real) 0.65, (real) 0.02, (real) 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  return bookshelf_tall;
}

inline host_fn World get_bookshelf_thin() {
  World bookshelf_thin;
  bookshelf_thin.add_capsule({(real) 1.1, (real) 0.25, (real) 0.37}, {(real) 1.1, (real) 0.25, (real) 0.51}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 0.8, (real) 0.25, (real) 0.37}, {(real) 0.8, (real) 0.25, (real) 0.51}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 1.1, -(real) 0.25, (real) 0.62}, {(real) 1.1, -(real) 0.25, (real) 0.76}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 0.8, -(real) 0.25, (real) 0.62}, {(real) 0.8, -(real) 0.25, (real) 0.76}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 1.1, (real) 0.25, (real) 0.87}, {(real) 1.1, (real) 0.25, (real) 1.01}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 0.8, (real) 0.25, (real) 0.87}, {(real) 0.8, (real) 0.25, (real) 1.01}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 1.1, -(real) 0.25, (real) 1.12}, {(real) 1.1, -(real) 0.25, (real) 1.26}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 0.8, -(real) 0.25, (real) 1.12}, {(real) 0.8, -(real) 0.25, (real) 1.26}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 1.1, (real) 0.25, (real) 1.37}, {(real) 1.1, (real) 0.25, (real) 1.51}, (real) 0.03);
  bookshelf_thin.add_capsule({(real) 0.8, (real) 0.25, (real) 1.37}, {(real) 0.8, (real) 0.25, (real) 1.51}, (real) 0.03);
  bookshelf_thin.add_box({1, 0, (real) 0.35}, {(real) 0.31, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({1, 0, (real) 0.6}, {(real) 0.31, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({1, 0, (real) 0.85}, {(real) 0.31, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({1, 0, (real) 1.1}, {(real) 0.31, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({1, 0, (real) 1.35}, {(real) 0.31, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({1, 0, (real) 0.9}, {(real) 0.31, (real) 0.02, (real) 0.7}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({1, 0, (real) 1.6}, {(real) 0.31, (real) 0.5, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({(real) 0.7, (real) 0.5, (real) 0.82}, {(real) 0.02, (real) 0.02, (real) 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({(real) 0.7, -(real) 0.5, (real) 0.82}, {(real) 0.02, (real) 0.02, (real) 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({(real) 1.3, (real) 0.5, (real) 0.82}, {(real) 0.02, (real) 0.02, (real) 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  bookshelf_thin.add_box({(real) 1.3, -(real) 0.5, (real) 0.82}, {(real) 0.02, (real) 0.02, (real) 0.82}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  return bookshelf_thin;
}

inline host_fn World get_scene_box() {
  World box;
  box.add_capsule({(real) 0.8, 0, (real) 0.48}, {(real) 0.8, 0, (real) 0.62}, (real) 0.03);
  box.add_box({(real) 0.8, 0, (real) 0.44}, {(real) 0.35, (real) 0.35, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  box.add_box({(real) 0.8, -(real) 0.35, (real) 0.8}, {(real) 0.35, (real) 0.02, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  box.add_box({(real) 0.8, (real) 0.35, (real) 0.8}, {(real) 0.35, (real) 0.02, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  box.add_box({(real) 0.45, 0, (real) 0.7}, {(real) 0.02, (real) 0.35, (real) 0.3}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  box.add_box({(real) 0.9, 0, (real) 1.35}, {(real) 0.35, (real) 0.35, (real) 0.02}, {-(real) 0.707107, 0, (real) 0.707107, 0, -1, 0, (real) 0.707107, 0, (real) 0.707107});
  box.add_box({(real) 1.15, 0, (real) 0.8}, {(real) 0.02, (real) 0.35, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  return box;
}

inline host_fn World get_scene_cage() {
  World cage;
  cage.add_box({(real) 0.8, 0, (real) 0.52}, {(real) 0.035, (real) 0.035, (real) 0.035}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  cage.add_box({(real) 0.8, 0, (real) 0.44}, {(real) 0.35, (real) 0.35, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  cage.add_box({(real) 0.8, -(real) 0.35, (real) 0.8}, {(real) 0.35, (real) 0.02, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  cage.add_box({(real) 0.8, (real) 0.35, (real) 0.8}, {(real) 0.35, (real) 0.02, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  cage.add_box({(real) 0.45, 0, (real) 0.9}, {(real) 0.02, (real) 0.35, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  cage.add_box({(real) 0.45, 0, (real) 0.6}, {(real) 0.02, (real) 0.35, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  cage.add_box({(real) 0.85, 0, (real) 1.15}, {(real) 0.325, (real) 0.325, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  cage.add_box({(real) 1.15, 0, (real) 0.8}, {(real) 0.02, (real) 0.35, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  return cage;
}

inline host_fn World get_scene_table() {
  World table;
  table.add_capsule({(real) 0.85, 0, (real) 0.74}, {(real) 0.85, 0, (real) 0.86}, (real) 0.03);
  table.add_box({(real) 0.75, (real) 0.4, (real) 0.85}, {(real) 0.125, (real) 0.125, (real) 0.125}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 1.5, (real) 0.85, (real) 0.35}, {(real) 0.025, (real) 0.025, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 0.6, (real) 0.85, (real) 0.35}, {(real) 0.025, (real) 0.025, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 1.5, -(real) 0.85, (real) 0.35}, {(real) 0.025, (real) 0.025, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 0.6, -(real) 0.85, (real) 0.35}, {(real) 0.025, (real) 0.025, (real) 0.35}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 1.05, 0, (real) 0.7}, {(real) 0.6, 1, (real) 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_capsule({(real) 1.35, 0, (real) 0.675}, {(real) 1.35, 0, (real) 1.025}, (real) 0.05);
  table.add_box({(real) 1.05, -(real) 0.2, (real) 0.9}, {(real) 0.1, (real) 0.01, (real) 0.2}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 0.65, (real) 0.2, (real) 0.9}, {(real) 0.01, (real) 0.1, (real) 0.2}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 0.65, -(real) 0.2, (real) 0.9}, {(real) 0.1, (real) 0.025, (real) 0.175}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  table.add_box({(real) 1.05, (real) 0.2, (real) 0.9}, {(real) 0.1, (real) 0.025, (real) 0.175}, {-1, 0, 0, 0, -1, 0, 0, 0, 1});
  return table;
}

inline host_fn World get_kitchen_no_doors() {
  World world;
  world.add_box({(real) 0.427, -(real) 1.0847, (real) 1.1974}, {(real) 0.36029, (real) 0.06, (real) 0.36242}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.427, (real) 0.48791, (real) 1.1974}, {(real) 0.36029, (real) 0.06, (real) 0.36242}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.427, -(real) 0.29838, (real) 1.5298}, {(real) 0.36029, (real) 0.84629, (real) 0.03}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.427, -(real) 0.29838, (real) 1.2074}, {(real) 0.36029, (real) 0.84629, (real) 0.03}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.427, -(real) 0.29838, (real) 0.86498}, {(real) 0.36029, (real) 0.84629, (real) 0.03}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.72729, -(real) 0.29838, (real) 1.2074}, {(real) 0.06, (real) 0.84629, (real) 0.36242}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.35698, (real) 0.425, -(real) 0.066059}, {(real) 0.41378, (real) 0.8425, (real) 0.43736}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.35698, -(real) 1.1546, -(real) 0.066059}, {(real) 0.41378, (real) 0.03, (real) 0.43736}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.35698, -(real) 0.80106, (real) 0.31136}, {(real) 0.41378, (real) 0.38356, (real) 0.06}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.35698, -(real) 0.80106, -(real) 0.055}, {(real) 0.41378, (real) 0.38356, (real) 0.025}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box({(real) 0.35698, -(real) 0.80106, -(real) 0.47336}, {(real) 0.41378, (real) 0.38356, (real) 0.03}, {1, 0, 0, 0, 1, 0, 0, 0, 1});
  return world;
}

inline host_fn World get_demo2_world() {
  World result;
  result.add_box({(real) 0.7, 0, (real) 0.381}, {1, (real) 0.75, (real) 0.381}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});                       // table 76 cm high
  result.add_box({(real) 0.67, -(real) 0.1475, (real) 0.96}, {(real) 0.35, (real) 0.025, (real) 0.2}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // vertical plate (no coll)
  // result.add_box({0.7, 0, -0.55}, {1.0, 0.75, 0.4}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});                 // table 76 cm high
  // result.add_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // vertical plate (no coll)
  return result;
}

inline host_fn World get_lab_world() {
  World world;
  // Lab environment
  // add_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); //table
  world.add_box({(real) 0.7, 0, -(real) 0.55}, {1, (real) 0.75, (real) 0.4}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});                                 // table
  world.add_box({-(real) 0.5654, -(real) 0.8145, (real) 0.3248}, {(real) 0.381, (real) 0.635, (real) 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // computer and screens next to link6
  world.add_box({(real) 0.4506, -(real) 1.3479, (real) 0.3248}, {(real) 0.635, (real) 0.381, (real) 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // computer and screens next to gen3lite
  world.add_box({0, 0, -(real) 0.4091}, {(real) 0.1459, (real) 0.2018, (real) 0.4091}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));                       // link6 base
  world.add_box({(real) 0.6665, (real) 1.1286, 0}, {(real) 0.508, (real) 0.3175, (real) 1.0457}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));             // UR5
  world.add_box({0, -(real) 0.75, 1}, {(real) 0.10, (real) 0.15, 2}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));                                         // cables

  // DEMO 1 bin 1
  // world.add_box({0.4114, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // bottom
  // world.add_box({0.2114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // front
  // world.add_box({0.6114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // back
  // world.add_box({0.4114, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // left
  // world.add_box({0.4114, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // right

  // // DEMO 1 bin 2
  // world.add_box({0.8716, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // bottom
  // world.add_box({0.6716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // front
  // world.add_box({1.0716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // back
  // world.add_box({0.8716, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // left
  // world.add_box({0.8716, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // right

  // DEMO 1 obstacle
  world.add_box({(real) 0.67, -(real) 0.1475, -(real) 0.0562}, {(real) 0.35, (real) 0.025, (real) 0.4}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // vertical plate (no coll)

  return world;
}
} // namespace blast
