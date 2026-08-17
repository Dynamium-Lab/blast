#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include <blast>
#include "catch2/catch.hpp"

#include "test_helper/test_dynamic_objects.hpp"
#include "test_helper/test_point_clouds.hpp"

using namespace blast;

TEST_CASE("Test - GJK Box vs Capsule (random generation)") {
  Capsule           capsule;
  Box               box;
  PointCloud        cloud;
  std::vector<Vec3> line(2);
  Vec3              center, extents, rpy;
  Mat3              rotation;

  // Run GJKs to warm up (for benchmarking purposes)
  for (int i = 0; i < 20; i++) {
    real a         = 3.0;
    capsule.radius = 0.5;
    capsule.p1     = {a * random_real(), a * random_real(), a * random_real()};
    capsule.p2     = {a * random_real(), a * random_real(), a * random_real()};

    // Random box
    real b  = 3.0;
    real c  = 5.0;
    center  = {b * random_real(), b * random_real(), b * random_real()};
    extents = {c * random_real(), c * random_real(), c * random_real()};
    // center -= Vec3{20, 20, 20}; // enforce never reaching EPA
    for (int j = 0; j < 3; j++) {
      extents[j] = std::abs(extents[j]) + 0.1; // enforce positive non zero extents
    }
    rpy      = {180 * random_real(), 180 * random_real(), 180 * random_real()};
    rotation = rpy2rotation(rpy);
    box      = {center, extents, rotation};

    distance(capsule, box);
    distance_GJK(capsule, box);
    distance_GJK_simple(capsule, box);
  }

  // Tests
  int                num_tests = 1e3;
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
    // center -= Vec3{20, 20, 20}; // enforce never reaching EPA
    for (int j = 0; j < 3; j++) {
      extents[j] = std::abs(extents[j]) + 0.1; // enforce positive non zero extents
    }
    rpy      = {180 * random_real(), 180 * random_real(), 180 * random_real()};
    rotation = rpy2rotation(rpy);
    box      = {center, extents, rotation};

    // Get point sets
    cloud = point_cloud_from_box(box, 0);
    line  = {capsule.p1, capsule.p2};

    // Checks (in random order)
    real real_dist = distance(capsule, box);

    std::vector<std::pair<
            std::function<real(Capsule, Box)>,
            std::string>>
                       functions = {{distance_GJK_simple, "Messy GJK"}, {distance_GJK, "Ericson GJK"}};
    std::random_device rd;
    std::mt19937       g(rd());
    std::shuffle(functions.begin(), functions.end(), g);

    for (const auto& function: functions) {
      real dist = function.first(capsule, box);
      CHECK(is_close(real_dist, dist));
      if (!is_close(real_dist, dist, 1e-4)) {
        std::cout << "Fail - " << function.second << ": " << std::fixed << std::setprecision(10) << dist << "!=" << real_dist << std::endl;
        printf("capsule = {{%.10f,%.10f,%.10f},{%.10f,%.10f,%.10f},%.2f};\n", capsule.p1.x, capsule.p1.y, capsule.p1.z, capsule.p2.x, capsule.p2.y, capsule.p2.z, capsule.radius);
        printf("box = {{%.10f,%.10f,%.10f},{%.10f,%.10f,%.10f},{%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f}};\n",
               box.center.x, box.center.y, box.center.z,
               box.extents.x, box.extents.y, box.extents.z,
               box.rotation[0], box.rotation[1], box.rotation[2],
               box.rotation[3], box.rotation[4], box.rotation[5],
               box.rotation[6], box.rotation[7], box.rotation[8]);
      }
    }

    if (real_dist + capsule.radius < 0) {
      epa[i] = 1.0f;
    } else {
      epa[i] = 0.0f;
    }

#if TRACY_ACTIVE
    {
      ZoneScopedN("sleep_for");
      std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
#endif

    // Benchmarking
    // Only run benchmarking with much lower test count
    // if (real_dist + capsule.radius < 0)
    //   std::cout << "EPA reached" << std::endl;
    // else
    //   std::cout << "GJK only" << std::endl;
    // BENCHMARK("distance") {
    //   real dist = distance(capsule, box);
    //   return;
    // };
    // BENCHMARK("distance_GJK (ericson)") {
    //   real dist = distance_GJK(capsule, box);
    //   return;
    // };
    // BENCHMARK("general_GJK (ericson)") {
    //   real dist = general_GJK(&line[0], line.size(), &cloud.points[0], cloud.points.size()) - capsule.radius;
    //   return;
    // };
    // BENCHMARK("distance_GJK_simple") {
    //   real dist = distance_GJK_simple(capsule, box);
    //   return;
    // };
    // BENCHMARK("solve_general_GJK") {
    //   real dist = solve_general_GJK(&line[0], line.size(), &cloud.points[0], cloud.points.size()) - capsule.radius;
    //   return;
    // };
    // std::cout << "\n--------------------------------------------" << std::endl;
  }

  real sum = 0.0;
  for (int i = 0; i < epa.size(); i++) {
    sum += epa[i];
  }

  std::cout << "Entering EPA algorithm " << sum / epa.size() * 100.0 << " % of the time" << std::endl;
}

TEST_CASE("Benchmark - GJK computation time vs point cloud size") {
  Capsule           capsule;
  Box               box;
  PointCloud        cloud;
  std::vector<Vec3> line(2);
  Vec3              center, extents, rpy;
  Mat3              rotation;


  int                num_tests = 2;
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
    // center -= Vec3{20, 20, 20}; // enforce never reaching EPA
    for (int j = 0; j < 3; j++) {
      extents[j] = std::abs(extents[j]) + 0.1; // enforce positive non zero extents
    }
    rpy      = {180 * random_real(), 180 * random_real(), 180 * random_real()};
    rotation = rpy2rotation(rpy);
    box      = {center, extents, rotation};

    real real_dist = distance(capsule, box);
    if (real_dist + capsule.radius < 0) {
      epa[i] = 1.0f;
    } else {
      epa[i] = 0.0f;
    }

    // Benchmarking
    if (real_dist < capsule.radius)
      std::cout << "EPA reached" << std::endl;
    else
      std::cout << "GJK only" << std::endl;
    for (int i = 0; i < 8; i++) {
      cloud = point_cloud_from_box(box, i * 10);
      std::cout << "\nPoints: " << 8 + i * 10 << std::endl;
      BENCHMARK("general_GJK()") {
        real dist = general_GJK(&line[0], line.size(), &cloud.points[0], cloud.points.size()) - capsule.radius;
        return dist;
      };
    }
    for (int i = 0; i < 8; i++) {
      cloud = point_cloud_from_box(box, i * 10);
      std::cout << "\nPoints: " << 8 + i * 10 << std::endl;
      BENCHMARK("solve_general_GJK()") {
        real dist = solve_general_GJK(&line[0], line.size(), &cloud.points[0], cloud.points.size()) - capsule.radius;
        return dist;
      };
    }
    std::cout << "\n--------------------------------------------" << std::endl;
  }

  real sum = 0.0;
  for (int i = 0; i < epa.size(); i++) {
    sum += epa[i];
  }

  std::cout << "Entering EPA algorithm " << sum / epa.size() * 100.0 << " % of the time" << std::endl;
}
