#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include <blast>
#include "catch2/catch.hpp"

using namespace blast;

// Capsule - Sphere collision tests
struct collision_test_sph {
  Sphere  sphere;
  Capsule capsule;
  real    expected_dist;
};
collision_test_sph test[] = {
        {{{4.27, 1.73, 3.74}, 1}, {{4, 4.2, 4.3}, {0.5, 0.4, 0.9}, 0.5}, 0.42},
        {{{2.75, 1.48, 3.44}, 1}, {{4, 4.2, 4.3}, {0.5, 0.4, 0.9}, 0.5}, -0.25},
        {{{5.72, 6.47, 3.78}, 1}, {{4, 4.2, 4.3}, {0.5, 0.4, 0.9}, 0.5}, 1.4},
        {{{4.74, 5.26, 4.52}, 1}, {{4, 4.2, 4.3}, {0.5, 0.4, 0.9}, 0.5}, -0.18},
};

// Capsule - Capsule collision tests
struct collision_test_caps {
  Capsule capsule1;
  Capsule capsule2;
  real    expected_dist;
};
collision_test_caps test_caps[] = {
        {{{0.5, 0.95, 6.73}, {4, 4.9, 3.51}, 0.5}, {{4, 4.2, 4.3}, {0.5, 0.4, 0.9}, 0.5}, -0.54},
        {{{5.16, 4.9, 4.37}, {3.55, 0.95, 8.84}, 0.5}, {{4, 4.2, 4.3}, {0.5, 0.4, 0.9}, 0.5}, 0.16},
        {{{5.16, 4.9, 4.37}, {3.55, 0.95, 8.84}, 0.5}, {{5.45, 4.78, 4.57}, {6.57, 1, -0.2}, 0.5}, -0.66},
        {{{2.59, 1.88, 4.22}, {18.85, 10.07, 6.1}, 1}, {{0.55, 21.32, 3.08}, {5.07, 3.62, 4.19}, 1}, -1.44},
        {{{-9.26, 12.81, 7.98}, {6.63, 4.62, 4.04}, 1}, {{0.55, 21.32, 3.08}, {5.07, 3.62, 4.19}, 1}, -1.52},
};

// Capsule - Box collision tests
struct collision_test_box {
  Box     box;
  Capsule capsule;
  real    expected_dist = 0;
};
collision_test_box test_obb[] = {
        /*Test1*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{-1.21, -5.18, 18.05}, {-3.89, 8.59, 6.3}, 1}, 1.63659624},
        /*Test2*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{10.46, 0.97, 3.7}, {7.79, 14.74, -8.04}, 1}, 1.50169942},
        /*Test3*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{10.05, 1.11, 12.87}, {7.37, 14.89, 1.13}, 1}, 0.33397901},
        /*Test4*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{4.82, 6.9, 12.84}, {4.7, -7.14, 24.59}, 1}, 4.00838054},
        /*Test5*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{-3.06, 5.73, 1.67}, {-0.39, -8.05, 13.41}, 1}, 0.89513819},
        /*Test6*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{4.97, 2.79, 12.23}, {2.29, 16.57, 0.48}, 1}, 1.69981481},
        /*Test7*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{11.49, -0.06, 8.32}, {14.17, -13.84, 20.07}, 1}, 0.49000000},
        /*Test8*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{6.76, -1.55, 8}, {9.44, -15.33, 19.74}, 1}, 0.45000000},
        /*Test9*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{-1.15, 13.92, -7.48}, {1.53, 0.14, 4.27}, 1}, 0.73046237},
        /*Test10*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{4.54, 0.1, 10.59}, {1.86, 13.87, -1.15}, 1}, -1.00000000},
        /*Test11*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{10.76, 0.28, 8.08}, {13.44, -13.5, 19.83}, 1}, -0.21961454},
        /*Test12*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{4.96, 0.43, 8.3}, {7.64, -13.35, 20.05}, 1}, -1.53000000},
        /*Test13*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{4.48, -4.07, 0.76}, {8.64, -4.41, 18.58}, 1}, 3.06923695},
        /*Test14*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{17.48, 2.95, 13.77}, {-0.82, 3.11, 13.77}, 1}, 2.41054264},
        /*Test15*/ {{{5, 0, 9}, {0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{11.18, 8.56, 4.82}, {11.04, -6.44, 15.29}, 1}, 0.09912546},
        // Box.R and capsule.r changed
        /*Test16*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{2.94, -6.06, 5.74}, {4.01, -9.86, 0.98}, 0.5}, 1.00474115},
        /*Test17*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{1.64, 9.52, 11.7}, {2.71, 5.72, 6.94}, 0.5}, 0.22669533},
        /*Test18*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{3, 5.69, 6.2}, {4.07, 1.89, 1.44}, 0.5}, 0.42102531},
        /*Test19*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{3.1, 4.23, 6.22}, {4.17, 0.43, 1.46}, 0.5}, 0.10695205},
        /*Test20*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{6.35, 8.57, 12.85}, {7.42, 4.77, 8.09}, 0.5}, 0.85902522},
        /*Test21*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{1.88, -2.47, 7.07}, {2.95, -6.27, 2.31}, 0.5}, 0.37892454},
        /*Test22*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{1.8, 3.83, 14.57}, {2.87, 0.03, 9.81}, 0.5}, 1.47889394},
        /*Test23*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{4.28, 6.32, 8.33}, {3.21, 10.12, 13.09}, 0.5}, 0.82000000},
        /*Test24*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{2.36, 2.3, 6.31}, {3.43, -1.5, 1.55}, 0.5}, 0.26887914},
        /*Test25*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{2.93, 7.6, 12.29}, {4, 3.8, 7.53}, 0.5}, -0.93234019},
        /*Test26*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{7.55, 0.92, 10.24}, {6.48, 4.72, 15}, 0.5}, -0.32852683},
        /*Test27*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{6.79, 5, 13.29}, {7.86, 1.2, 8.52}, 0.5}, -0.40212621},
        /*Test28*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{5.66, -0.92, 7.44}, {8.04, 4.16, 9.96}, 0.5}, 0.87078210},
        /*Test29*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{8.94, 3.63, 11.01}, {8.94, -2.55, 11.01}, 0.5}, 1.24844065},
        /*Test30*/ {{{5, 0, 9}, {0.1, 5, 3}, {0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707}}, {{5.91, 5.9, 7.39}, {4.14, 5.9, 13.32}, 0.5}, 0.40000000},

        /*Test31*/ {{{0, 0, 0}, {1.421459, 0.796365, 0.574752}, {1, 0, 0, 0, 1, 0, 0, 0, 1}}, {{2.447825, 6.30643376, -5.224298}, {-5.4486154, -5.149872, 5.1825604}, 0}, -0.05877295},
        /*Test32*/ {{{0, 0, 0}, {1.7236252662212059, 1.6611691154387149, 0.20442662553572166}, {1, 0, 0, 0, 1, 0, 0, 0, 1}}, {{-4.6785661670641119, 1.5869561202501099, -4.7958108710825522}, {6.9379851583657839, -2.4864207406657761, 4.5569750611707658}, 0}, -0.438675},
        /*Test33*/ {{{0, 0, 0}, {1.7865788377999734, 1.4615853859974308, 1.7957946652540584}, {0.16899262598095025, 0.69208966716438747, -0.70175022976009782, -0.73366055668912455, 0.56377691511991157, 0.37933860538637498, 0.65816690886329632, 0.45074103716231956, 0.60303303184416857}}, {{-0.050034305603475548, 0.40002630593895194, 2.3444724393109695}, {1.3614165219982919, 0.62636518062708535, -1.2242779633088015}, 0}, -1.3557069221709168},
};

TEST_CASE("Primitive distance functions test", "[World]") {
  using namespace blast;

  real TEST_EPSILON = 1e-2;

  std::vector<Vec3> v1 = {{-4.6607382574035014, -3.4686124971115451, 0.89349474531514073},
                          {6.4517093458638977, 6.9310459212756506, -0.24683715666348305}};
  std::vector<Vec3> v2 = {{2.4373332611708229, 2.5910471347625212, -0.04489407117738407},
                          {2.5547766616001675, 4.0525983620611852, 0.19409068785959360},
                          {1.8335195262912714, 2.7351427068670993, -0.62940465909289878},
                          {1.9509629267206159, 4.1966939341657632, -0.39041990005592114},
                          {1.6790895426140411, 2.5264996514086859, 0.72247776841327949},
                          {1.7965329430433856, 3.9880508787073494, 0.96146252745025707},
                          {1.0752758077344895, 2.6705952235132639, 0.13796718049776469},
                          {1.1927192081638340, 4.1321464508119279, 0.37695193953474238}};
  // real dist = general_gjk(v1, v2);

  // Capsule - Sphere test
  for (auto t: test) {
    real dist = distance(t.capsule, t.sphere);
    CHECK(is_close(dist, t.expected_dist, TEST_EPSILON));
  }

  // Capsule - Capsule test
  for (auto t: test_caps) {
    real dist = distance(t.capsule1, t.capsule2);
    CHECK(is_close(dist, t.expected_dist, TEST_EPSILON));
  }

  // Capsule - OBB test (vectors accelerated)
  for (const auto& [box, capsule, expected_dist]: test_obb) {
    real dist = distance(capsule, box);
    CHECK(is_close(dist, expected_dist, TEST_EPSILON));
    auto new_caps = capsule;
    new_caps.p1   = capsule.p2;
    new_caps.p2   = capsule.p1;
    dist          = distance(new_caps, box);
    CHECK(is_close(dist, expected_dist, TEST_EPSILON));
  }
}

TEST_CASE("Collision method benchmarks", "[World]") {
  using namespace blast;

  // todo: remove?
  // BENCHMARK("box - capsule Box test with vectors") {
  //     real dist;
  //     for (auto t : test_obb) {
  //         dist = distance_old(t.box, t.capsule);
  //         // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
  //     }
  //     return dist;
  // };

  BENCHMARK("box - capsule Box test with vectors accelerated") {
    real dist = 0;
    for (const auto& t: test_obb) {
      dist = distance(t.capsule, t.box);
      // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
    }
    return dist;
  };

  // todo: remove?
  // BENCHMARK("box - capsule Box test with GJK") {
  //     real dist;
  //     for (auto t : test_obb) {
  //         dist = gjk_box_capsule(t.capsule, t.box);
  //         // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
  //     }
  //     return dist;
  // };
}


// todo: remove?
// TEST_CASE("Collision method comparison exhaustive (Box-capsules)", "[World]") {
//     using namespace blast;

//     real TEST_EPSILON = 1e-2;

//     vector<Box> obb_list;
//     vector<capsule> caps_list;
//     int n = 1000;
//     // bool error_info = false;

//     Array s(3);
//     Array c(3);
//     Array angles(3);

//     Array obb_e(3);
//     Array obb_c(3);

//     Array seg_p1(3);
//     Array seg_p2(3);

//     for(u32 i = 0; i < n; i++) {
//         fill_random(angles, PI);
//         blast::sincos(angles, s, c);

//         Mat3 Rz = {c[0], -s[0], 0, s[0], c[0], 0, 0, 0, 1};
//         Mat3 Ry = {c[1], 0, s[1], 0, 1, 0, -s[1], 0, c[1]};
//         Mat3 Rx = {1, 0, 0, 0, c[2], -s[2], 0, s[2], c[2]};
//         auto R = Rz*Ry*Rx;
//         CHECK(abs(1.0 - sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2])) < TEST_EPSILON);
//         CHECK(abs(1.0 - sqrt(R[3]*R[3] + R[4]*R[4] + R[5]*R[5])) < TEST_EPSILON);
//         CHECK(abs(1.0 - sqrt(R[6]*R[6] + R[7]*R[7] + R[8]*R[8])) < TEST_EPSILON);

//         // Generate n random Box
//         fill_random(obb_e, 2);
//         obb_e = array_abs(obb_e);
//         fill_random(obb_c, 5);
//         obb_c = array_abs(obb_c);
//         obb_list.push_back({{obb_c[0], obb_c[1], obb_c[2]}, {obb_e[0], obb_e[1], obb_e[2]}, R});

//         // Generate n random capsule
//         fill_random(seg_p1, 7);
//         fill_random(seg_p2, 7);
//         auto r = 0.0;
//         caps_list.push_back({{seg_p1[0], seg_p1[1], seg_p1[2]}, {seg_p2[0], seg_p2[1], seg_p2[2]}, r});
//     }

//     vector<capsule> caps_failed;
//     vector<Box> obb_failed;
//     for (int cap = 0; cap < caps_list.size(); cap++) {
//         // Box collisions
//         for (int i = 0; i < obb_list.size(); i++) {
//             auto dist_min_vec = distance_old(obb_list[i], caps_list[cap]);
//             auto dist_min_vector_acc = distance(obb_list[i], caps_list[cap]);
//             // auto dist_min_new = distmin_hull(obb_list[i], caps_list[cap]);
//             // auto dist_min_gjk = gjk_box_capsule(caps_list[cap], obb_list[i]);

//             CHECK(abs(dist_min_vec - dist_min_vector_acc) < TEST_EPSILON);
//             if (abs(dist_min_vec - dist_min_vector_acc) > TEST_EPSILON) {
//                 // save and test capsule and obb in future
//                 caps_failed.push_back(caps_list[cap]);
//                 obb_failed.push_back(obb_list[i]);
//             }
//         }
//     }
//     // real total_error = 0;
//     // real max_error = 0;
//     // bool error_when_neg = false;
//     // bool error_when_pos = false;
//     // std::vector<Vec3> error_distmin_distminnew(caps_failed.size());
//     // for (u32 i = 0; i < caps_failed.size(); i++) {
//     //     auto dist_min_vec = distance_old(obb_failed[i], caps_failed[i]);
//     //     auto dist_min_new = distmin_hull(obb_failed[i], caps_failed[i]);
//     //     auto dist_min_gjk = gjk_box_capsule(caps_failed[i], obb_failed[i]);
//     //     auto dist_min_pierce = distmin_pierce(obb_failed[i], caps_failed[i]);
//     //     real error = abs(dist_min_new - dist_min_vec) * 100 / dist_min;
//     //     total_error += error;
//     //     int dist = 0;
//     //     max_error = error > max_error ? error : max_error;
//     //     if (error > 100)
//     //         real found_error = 0;
//     //     if (dist_min_vec >=0)
//     //         error_when_pos = true;
//     //     if (dist_min_vec < 0) {
//     //         error_when_neg = true;
//     //         real check = dist_min_new + dist_min_vec;
//     //     }
//     //     if (error_info) {
//     //         Vec3 error_point = { error, dist_min_vec, dist_min_new };
//     //         error_distmin_distminnew.push_back(error_point);
//     //     }
//     // }
//     // real avg_error = total_error / caps_failed.size();
//     // real avg_error_over_all = total_error / (n*n);
//     // real percent_error = caps_failed.size() / (n*n);
//     // std::cout << caps_failed.size() << " failed out of " << n*n << " tests ( " << percent_error << " % ) \n";
//     // std::cout << "Average error of " << avg_error << "% and a max error of " << max_error << "%. \n";
//     // std::cout << "Counting the tests which passed, we find an average error of: " << avg_error_over_all << " %.  \n";
//     // std::cout << "Error when positive dist (true or false): " << error_when_pos << " \n";
//     // std::cout << "Error when negative dist (true or false): " << error_when_neg << " \n";
//     // if (error_info) {
//     //     #include <fstream>
//     //     // Open a CSV file for writing
//     //     std::ofstream csvFile("C:\\Users\\thoma\\Desktop\\example.csv");
//     //     // Check if the file is open
//     //     if (!csvFile.is_open()) {
//     //         std::cerr << "Error opening file!" << std::endl;
//     //     }
//     //     // Write headers to the CSV file
//     //     csvFile << "error (%), dist_min_vec, dist_min_new" << std::endl;
//     //     // Write data
//     //     for (int i = 0; i < error_distmin_distminnew.size(); i++) {
//     //         Vec3 error_pt = error_distmin_distminnew[i];
//     //         csvFile << (error_pt).x << ";" << (error_pt).y << ";" << (error_pt).z << std::endl;
//     //     }
//     //     // Close the file
//     //     csvFile.close();
//     //     std::cout << "CSV file created successfully." << std::endl;
//     // }
// }
