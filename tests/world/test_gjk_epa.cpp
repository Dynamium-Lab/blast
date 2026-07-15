#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include <blast>
#include "catch2/catch.hpp"

#include "test_helper/test_dynamic_objects.hpp"
#include "test_helper/test_point_clouds.hpp"

using namespace blast;

TEST_CASE("Test - GJK for box vs capsule") {
  World     world;
  Capsule   capsule;
  gjkresult result;
  capsule.radius = 1;
  add_dynamic_obstacles(world);
  int i = 0;

  for (int j = 0; j < world.dynamic_boxes.size(); j++) {
    auto dynamic_box = world.dynamic_boxes[j];
    for (int i = 0; i < dynamic_box.trajectory.size(); i++) {
      Box  box   = dynamic_box.trajectory[i];
      real a     = 2;
      capsule.p1 = {a * random_real(), a * random_real(), a * random_real()};
      capsule.p2 = {a * random_real(), a * random_real(), a * random_real()};
      CHECK(is_close(distance_GJK(capsule, box), distance(capsule, box)));
      if (!is_close(distance_GJK(capsule, box), distance(capsule, box))) {
        std::cout << "FAIL clean : " << distance_GJK(capsule, box) << "!=" << distance(capsule, box) << std::endl;
        printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
      }
      if (!is_close(distance_GJK_simple(capsule, box), distance(capsule, box))) {
        std::cout << "FAIL messy: " << distance_GJK_simple(capsule, box) << "!=" << distance(capsule, box) << std::endl;
        printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
        std::cout << "Box ID: " << i << "\nDynamicBox ID: " << j << std::endl;
      }
      CHECK(is_close(distance_GJK_simple(capsule, box), distance(capsule, box)));
      // if (i % 41 == 0) {
      //   BENCHMARK("distance()") {
      //     real dist = distance(capsule, box);
      //   };
      //   BENCHMARK("distance_GJK_OG") {
      //     result = GJK_solve_gjk_simple(capsule, box);
      //   };
      //   BENCHMARK("distance_GJK()") {
      //     real dist = distance_GJK(capsule, box);
      //   };
      //   std::cout << std::endl;
      // }
      // i++;
    }
  }
}

std::vector<std::pair<Capsule, std::pair<int, int>>> debugging_GJK = {
        // should all fail in GJK simple
        {{{-1.7853, -1.4976, -1.3274}, {-0.4435, 1.1769, -1.8909}, 1.00}, {16, 1}},
        {{{-0.5725, -1.1564, -1.8425}, {-0.6842, 0.7947, 1.4433}, 1.00}, {23, 1}},
        {{{-0.3146, 1.5202, -1.7772}, {1.9790, 1.8273, 1.8606}, 1.00}, {44, 0}},
        // {{{0.7349, -0.1461, 0.2158}, {0.5890, 1.6033, 1.8372}, 1.00}, {58, 1}}, // should enter EPA -> no mistake when using Ericson EPA
        {{{-1.6627, -0.5169, 0.9055}, {1.7553, 1.4242, 0.1141}, 1.00}, {68, 1}},
        {{{-0.4764, -1.6489, 1.6199}, {1.4919, 1.0736, 0.4002}, 1.00}, {77, 1}},
        // {{{1.1418, 1.8783, 1.4663}, {-0.8472, -0.1556, -0.9498}, 1.00}, {95, 1}}, // should enter EPA
        {{{-0.1298, -1.3521, 0.3215}, {1.9320, -1.2444, 1.9411}, 1.00}, {0, 0}}, // clean failed on this one!! FAIL clean : 0.639496!=-0.659484
};

// TEST_CASE("Debug - distance_GJK() vs distance() - box vs capsule") {
//   World     world;
//   gjkresult result;
//   add_dynamic_obstacles(world);
//   int i = 0;

//   for (int i = 0; i < debugging_GJK.size(); i++) {
//     DynamicBox dynamic_box = world.dynamic_boxes[debugging_GJK[i].second.second];
//     Capsule    capsule     = debugging_GJK[i].first;
//     Box        box         = dynamic_box.trajectory[debugging_GJK[i].second.first];

//     result = GJK_solve_gjk_simple(capsule, box);
//     CHECK(is_close(distance_GJK(capsule, box), distance(capsule, box)));
//     CHECK(is_close(result.minimal_distance, distance(capsule, box)));

//     // if (!is_close(distance_GJK(capsule, box), distance(capsule, box))) {
//     //   std::cout << "FAIL clean : " << result.minimal_distance << "!=" << distance(capsule, box) << std::endl;
//     //   printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//     // }
//     // result = GJK_solve_gjk_simple(capsule, box);
//     if (!is_close(result.minimal_distance, distance(capsule, box))) {
//       std::cout << "FAIL messy: " << result.minimal_distance << "!=" << distance(capsule, box) << std::endl;
//       printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//       std::cout << "Box ID: " << debugging_GJK[i].second.first << /* "\nDynamicBox ID: " << j <<*/ std::endl;
//     }
//   }
// }

std::vector<std::pair<World, std::string>> worlds = {
        {get_bookshelf_small(), "get_bookshelf_small()"},
        {get_bookshelf_tall(), "get_bookshelf_tall()"},
        {get_bookshelf_thin(), "get_bookshelf_thin()"},
        {get_scene_box(), "get_scene_box()"},
        {get_scene_cage(), "get_scene_cage()"},
        {get_scene_table(), "get_scene_table()"},
        {get_kitchen_no_doors(), "get_kitchen_no_doors()"},
        {get_lab_world(), "get_lab_world()"},
};

TEST_CASE("Test - GJK for point clouds") {
  Capsule           capsule;
  std::vector<Vec3> caps_pts;
  int               k = 1;
  for (int i = 0; i < worlds.size(); i++) {
    std::cout << "\nTest - " << worlds[i].second << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    auto tests = get_boxes_and_point_clouds(worlds[i].first);
    // for (int j = 0; j < 1; j++) {
    for (int j = 0; j < tests.size(); j++) {
      Box        box   = tests[j].first;
      PointCloud cloud = tests[j].second;

      real a         = 2;
      capsule.radius = 1;
      capsule.p1     = {a * random_real(), a * random_real(), a * random_real()};
      capsule.p2     = {a * random_real(), a * random_real(), a * random_real()};
      // capsule = {{0.1459, -1.2412, 0.0733}, {1.9317, -0.9945, 1.7238}, 1.00}; // doesn't enter EPA
      // capsule = {{-0.9016, -1.8075, 1.6098}, {-1.6268, 1.0098, 0.2212}, 1.00}; // fails somewhere...

      caps_pts = {capsule.p1, capsule.p2};

      // auto result = solve_general_GJK(caps_pts, cloud.set);
      CHECK(is_close(distance(capsule, box), general_GJK(caps_pts, cloud.set) - capsule.radius));
      CHECK(is_close(distance(capsule, box), solve_general_GJK(caps_pts, cloud.set).minimal_distance - capsule.radius));

      if (!is_close(distance(capsule, box), general_GJK(caps_pts, cloud.set) - capsule.radius)) {
        real dist = general_GJK(caps_pts, cloud.set) - capsule.radius;
        std::cout << "FAIL clean : " << general_GJK(caps_pts, cloud.set) - capsule.radius << "!=" << distance(capsule, box) << std::endl;
        printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
        std::cout << "Box ID: " << j << std::endl;
        std::cout << "Failure no. " << k << std::endl;
        k++;
      }
      if (!is_close(distance(capsule, box), solve_general_GJK(caps_pts, cloud.set).minimal_distance - capsule.radius)) {
        real dist = solve_general_GJK(caps_pts, cloud.set).minimal_distance - capsule.radius;
        std::cout << "FAIL messy: " << solve_general_GJK(caps_pts, cloud.set).minimal_distance - capsule.radius << "!=" << distance(capsule, box) << std::endl;
        printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
        std::cout << "Box ID: " << j << std::endl;
        std::cout << "Failure no. " << k << std::endl;
        k++;
      }
    }
  }
}

// TEST_CASE("Test - GJK for point clouds") {
//   Capsule           capsule;
//   std::vector<Vec3> caps_pts;
//   for (auto pair: worlds) {
//     std::cout << "\nTest - " << pair.second << std::endl;
//     std::cout << "------------------------------------------------------------" << std::endl;
//     auto tests = get_boxes_and_point_clouds(pair.first);

//     std::cout << "[";
//     for (auto pt: tests[0].second.set) {
//       std::cout << pt.x << ", " << pt.y << ", " << pt.z << ";\n";
//     }
//     std::cout << "];\n";
//   }
// }
