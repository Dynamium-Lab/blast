#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include <blast>
#include "catch2/catch.hpp"

#include "test_helper/test_dynamic_objects.hpp"
#include "test_helper/test_point_clouds.hpp"

using namespace blast;

// TEST_CASE("Test - GJK for box vs capsule") {
//   World     world;
//   Capsule   capsule;
//   gjkresult result;
//   capsule.radius = 1;
//   add_dynamic_obstacles(world);
//   int i = 0;

//   for (int j = 0; j < world.dynamic_boxes.size(); j++) {
//     auto dynamic_box = world.dynamic_boxes[j];
//     for (int i = 0; i < dynamic_box.trajectory.size(); i++) {
//       Box  box   = dynamic_box.trajectory[i];
//       real a     = 2;
//       capsule.p1 = {a * random_real(), a * random_real(), a * random_real()};
//       capsule.p2 = {a * random_real(), a * random_real(), a * random_real()};
//       CHECK(is_close(distance_GJK(capsule, box), distance(capsule, box)));
//       if (!is_close(distance_GJK(capsule, box), distance(capsule, box))) {
//         std::cout << "FAIL clean : " << distance_GJK(capsule, box) << "!=" << distance(capsule, box) << std::endl;
//         printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//       }
//       if (!is_close(distance_GJK_simple(capsule, box), distance(capsule, box))) {
//         std::cout << "FAIL messy: " << distance_GJK_simple(capsule, box) << "!=" << distance(capsule, box) << std::endl;
//         printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//         std::cout << "Box ID: " << i << "\nDynamicBox ID: " << j << std::endl;
//       }
//       CHECK(is_close(distance_GJK_simple(capsule, box), distance(capsule, box)));
//       // if (i % 41 == 0) {
//       //   BENCHMARK("distance()") {
//       //     real dist = distance(capsule, box);
//       //   };
//       //   BENCHMARK("distance_GJK_OG") {
//       //     result = GJK_solve_gjk_simple(capsule, box);
//       //   };
//       //   BENCHMARK("distance_GJK()") {
//       //     real dist = distance_GJK(capsule, box);
//       //   };
//       //   std::cout << std::endl;
//       // }
//       // i++;
//     }
//   }
// }

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
//   World world;
//   // gjkresult result;
//   add_dynamic_obstacles(world);
//   int i = 0;

//   for (int i = 0; i < debugging_GJK.size(); i++) {
//     DynamicBox dynamic_box = world.dynamic_boxes[debugging_GJK[i].second.second];
//     Capsule    capsule     = debugging_GJK[i].first;
//     Box        box         = dynamic_box.trajectory[debugging_GJK[i].second.first];


//     CHECK(is_close(distance_GJK(capsule, box), distance(capsule, box)));
//     CHECK(is_close(distance_GJK_simple(capsule, box), distance(capsule, box)));

//     if (!is_close(distance_GJK(capsule, box), distance(capsule, box))) {
//       std::cout << "FAIL clean : " << distance_GJK_simple(capsule, box) << "!=" << distance(capsule, box) << std::endl;
//       printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//     }
//     // result = GJK_solve_gjk_simple(capsule, box);
//     if (!is_close(distance_GJK_simple(capsule, box), distance(capsule, box))) {
//       std::cout << "FAIL messy: " << distance_GJK_simple(capsule, box) << "!=" << distance(capsule, box) << std::endl;
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

// TEST_CASE("Test - GJK for point clouds") {
//   Capsule           capsule;
//   std::vector<Vec3> caps_pts;
//   int               k = 1;
//   // for (int i = 0; i < 1; i++) {
//   for (int i = 0; i < worlds.size(); i++) {
//     std::cout << "\nTest - " << worlds[i].second << std::endl;
//     std::cout << "------------------------------------------------------------" << std::endl;
//     auto tests = get_boxes_and_point_clouds(worlds[i].first);
//     // for (int j = 0; j < 1; j++) {
//     for (int j = 0; j < tests.size(); j++) {
//       Box        box   = tests[j].first;
//       PointCloud cloud = tests[j].second;

//       // box   = {{0.45, 0, 0.6}, {0.02, 0.35, 0.02}, {-1, 0, 0, 0, -1, 0, 0, 0, 1}};
//       // cloud = point_cloud_from_box(box);

//       real a         = 2;
//       capsule.radius = 1;
//       capsule.p1     = {a * random_real(), a * random_real(), a * random_real()};
//       capsule.p2     = {a * random_real(), a * random_real(), a * random_real()};
//       // capsule = {{0.1459, -1.2412, 0.0733}, {1.9317, -0.9945, 1.7238}, 1.00}; // doesn't enter EPA
//       // capsule = {{-0.9016, -1.8075, 1.6098}, {-1.6268, 1.0098, 0.2212}, 1.00}; // fails somewhere...
//       // capsule = {{-0.1298, -1.3521, 0.3215}, {1.9320, -1.2444, 1.9411}, 1.00};
//       // capsule = {{-1.4225, -0.8594, 1.4846}, {1.4531, 0.0481, 1.5521}, 1.00};

//       caps_pts = {capsule.p1, capsule.p2};

//       // auto result = solve_general_GJK(caps_pts, cloud.points);
//       // CHECK(is_close(distance(capsule, box), general_GJK(caps_pts, cloud.points) - capsule.radius));
//       CHECK(is_close(distance(capsule, box), distance_GJK(capsule, box)));
//       // CHECK(is_close(distance(capsule, box), solve_general_GJK(caps_pts, cloud.points).minimal_distance - capsule.radius));
//       CHECK(is_close(distance(capsule, box), distance_GJK_simple(capsule, box)));

//       if (!is_close(distance(capsule, box), general_GJK(caps_pts, cloud.points) - capsule.radius)) {
//         real dist = general_GJK(caps_pts, cloud.points) - capsule.radius;
//         std::cout << "FAIL clean : " << general_GJK(caps_pts, cloud.points) - capsule.radius << "!=" << distance(capsule, box) << std::endl;
//         printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//         std::cout << "Box ID: " << j << std::endl;
//         std::cout << "Failure no. " << k << std::endl;
//         k++;
//       }
//       if (!is_close(distance(capsule, box), solve_general_GJK(caps_pts, cloud.points).minimal_distance - capsule.radius)) {
//         real dist = solve_general_GJK(caps_pts, cloud.points).minimal_distance - capsule.radius;
//         std::cout << "FAIL messy: " << solve_general_GJK(caps_pts, cloud.points).minimal_distance - capsule.radius << "!=" << distance(capsule, box) << std::endl;
//         printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//         std::cout << "Box ID: " << j << std::endl;
//         std::cout << "Failure no. " << k << std::endl;
//         k++;
//       }
//     }
//   }
// }

// TEST_CASE("Test - GJK for point clouds") {
//   Capsule           capsule;
//   std::vector<Vec3> caps_pts;
//   for (auto pair: worlds) {
//     std::cout << "\nTest - " << pair.second << std::endl;
//     std::cout << "------------------------------------------------------------" << std::endl;
//     auto tests = get_boxes_and_point_clouds(pair.first);

//     std::cout << "[";
//     for (auto pt: tests[0].second.points) {
//       std::cout << pt.x << ", " << pt.y << ", " << pt.z << ";\n";
//     }
//     std::cout << "];\n";
//   }
// }

// Capsule capsule = {{1.4442, -0.2837, 1.7364}, {-1.3301, -0.9690, -0.3969}, 1.00};
// Box box = {{0.9, 0, 1.35}, {0.35, 0.35, 0.02}, {-0.708, 0.0, 0.708, 0.0, -1.001, 0.0, 0.708, 0.0, 0.707}};

// TEST_CASE("Test - GJK for point clouds") {
//   std::vector<Vec3> caps_pts;
//   int               k = 1;

//   Box        box   = {{0.9, 0, 1.35}, {0.35, 0.35, 0.02}, {-0.70686, 0.0, 0.70736, 0.0, -1.0, 0.0, 0.70736, 0.0, 0.70686}};
//   PointCloud cloud = point_cloud_from_box(box);

//   Capsule capsule = {{1.4442, -0.2837, 1.7364}, {-1.3301, -0.9690, -0.3969}, 1.00};

//   // real a         = 2;
//   // capsule.radius = 1;
//   // capsule.p1     = {a * random_real(), a * random_real(), a * random_real()};
//   // capsule.p2     = {a * random_real(), a * random_real(), a * random_real()};

//   caps_pts = {capsule.p1, capsule.p2};

//   // No point in checking.. we know they fail for now
//   real real_dist = distance(capsule, box);
//   CHECK(is_close(distance(capsule, box), general_GJK(caps_pts, cloud.points) - capsule.radius));
//   CHECK(is_close(distance(capsule, box), distance_GJK(capsule, box)));
//   CHECK(is_close(distance(capsule, box), solve_general_GJK(caps_pts, cloud.points).minimal_distance - capsule.radius));
//   CHECK(is_close(distance(capsule, box), distance_GJK_simple(capsule, box)));

//   {
//     real dist = distance_GJK(capsule, box);
//     if (!is_close(dist, real_dist)) {
//       std::cout << "FAIL ericson : " << dist << "!=" << real_dist << std::endl;
//       // printf("capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//       printf("capsule = {{%.6f,%.6f,%.6f},{%.6f,%.6f,%.6f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//     }
//   }
//   {
//     real dist = distance_GJK_simple(capsule, box);
//     if (!is_close(dist, real_dist)) {
//       std::cout << "FAIL simple : " << dist << "!=" << real_dist << std::endl;
//       // printf("Capsule = {{%.4f,%.4f,%.4f},{%.4f,%.4f,%.4f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//       printf("capsule = {{%.6f,%.6f,%.6f},{%.6f,%.6f,%.6f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
//     }
//   }
// }

real determinant(Mat3 R) {
  return R(0, 0) * (R(1, 1) * R(2, 2) - R(1, 2) * R(2, 1)) - R(0, 1) * (R(1, 0) * R(2, 2) - R(1, 2) * R(2, 0)) + R(0, 2) * (R(1, 0) * R(2, 1) - R(1, 1) * R(2, 0));
}

Mat3 random_rotation_matrix() {
  real alpha = 3.14159;
  real beta  = 3.14159;
  real gamma = 3.14159;

  alpha *= random_real();
  beta *= random_real();
  gamma *= random_real();

  Mat3 Rx = {1, 0, 0, 0, cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha)};
  Mat3 Ry = {cos(beta), 0, -sin(beta), 0, 1, 0, sin(beta), 0, cos(beta)};
  Mat3 Rz = {cos(gamma), sin(gamma), 0, -sin(gamma), cos(gamma), 0, 0, 0, 1};

  if (!is_close(determinant(Rz * Ry * Rx), 1.0))
    std::cout << "Not orthogonal" << std::endl;

  return (Rz * Ry * Rx);
}

TEST_CASE("Test - GJK Box vs Capsule (random generation)") {
  Capsule           capsule;
  Box               box;
  PointCloud        cloud;
  std::vector<Vec3> line(2);
  Vec3              center, extents, rpy;
  Mat3              rotation;


  int                num_tests = 1;
  std::vector<float> epa(num_tests);
  for (int i = 0; i < num_tests; i++) {
    // Random capsule
    real a         = 3.0;
    capsule.radius = 0.5;
    capsule.p1     = {a * random_real(), a * random_real(), a * random_real()};
    capsule.p2     = {a * random_real(), a * random_real(), a * random_real()};

    // Random box
    real b  = 3.0;
    real c  = 5.0;
    center  = {b * random_real(), b * random_real(), b * random_real()};
    extents = {c * random_real(), c * random_real(), c * random_real()};
    // center -= Vec3{10, 10, 10};
    for (int i = 0; i < 3; i++) {
      extents[i] = std::abs(extents[i]) + 0.1; // enforce positive non zero extents
    }
    rpy      = {180 * random_real(), 180 * random_real(), 180 * random_real()};
    rotation = rpy2rotation(rpy);
    // rotation = random_rotation_matrix();
    box = {center, extents, rotation};

    capsule = {{-1.961081, -2.771850, -1.583547}, {2.475456, -2.403686, 1.102632}, 0.50};
    box     = {{-2.507173, 2.988569, 1.571170}, {3.736067, 3.489249, 2.199382}, {-0.822818, -0.175611, -0.540492, -0.369301, 0.888106, 0.273651, 0.431958, 0.424769, -0.795603}};

    // Get point sets
    cloud = point_cloud_from_box(box, true);
    line  = {capsule.p1, capsule.p2};

    // Checks
    real real_dist = distance(capsule, box);
    {
      real dist = distance_GJK(capsule, box);
      CHECK(is_close(real_dist, dist));
      // CHECK(is_close(real_dist, solve_general_GJK(line, cloud.points).minimal_distance - capsule.radius));
      if (!is_close(dist, real_dist)) {
        std::cout << "FAIL ericson : " << dist << "!=" << real_dist << std::endl;
        printf("capsule = {{%.6f,%.6f,%.6f},{%.6f,%.6f,%.6f},%.2f}\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
        printf("box = {{%.6f,%.6f,%.6f},{%.6f,%.6f,%.6f},{%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f}}\n",
               box.center.x, box.center.y, box.center.z,
               box.extents.x, box.extents.y, box.extents.z,
               box.rotation[0], box.rotation[1], box.rotation[2],
               box.rotation[3], box.rotation[4], box.rotation[5],
               box.rotation[6], box.rotation[7], box.rotation[8]);
      }
    }
    {
      real dist = distance_GJK_simple(capsule, box);
      CHECK(is_close(real_dist, dist));
      // CHECK(is_close(real_dist, general_GJK(line, cloud.points) - capsule.radius));
      if (!is_close(dist, real_dist)) {
        std::cout << "FAIL simple : " << dist << "!=" << real_dist << std::endl;
      }
    }

    if (real_dist + capsule.radius < 0) {
      epa[i] = 1.0f;
    } else {
      epa[i] = 0.0f;
    }

    // {
    //   ZoneScopedN("sleep_for");
    //   std::this_thread::sleep_for(std::chrono::microseconds(100));
    // }

    // Benchmarking
    // if (real_dist < capsule.radius)
    //   std::cout << "EPA reached" << std::endl;
    // else
    //   std::cout << "GJK only" << std::endl;
    // BENCHMARK("distance") {
    //   real dist = distance(capsule, box);
    // };
    // BENCHMARK("distance_GJK (ericson)") {
    //   real dist = distance_GJK(capsule, box);
    // };
    // BENCHMARK("general_GJK (ericson)") {
    //   real dist = general_GJK(line, cloud.points) - capsule.radius;
    // };
    // BENCHMARK("distance_GJK_simple") {
    //   real dist = distance_GJK_simple(capsule, box);
    // };
    // BENCHMARK("solve_general_GJK") {
    //   real dist = solve_general_GJK(line, cloud.points).minimal_distance - capsule.radius;
    // };
    // std::cout << "\n--------------------------------------------" << std::endl;
  }

  real sum = 0.0;
  for (int i = 0; i < epa.size(); i++) {
    sum += epa[i];
  }

  std::cout << "Entering EPA algorithm " << sum / epa.size() * 100.0 << " % of the time" << std::endl;
}
