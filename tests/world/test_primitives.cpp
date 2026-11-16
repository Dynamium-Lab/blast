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
        {{{(real) 4.27, (real) 1.73, (real) 3.74}, 1}, {{4, (real) 4.2, (real) 4.3}, {(real) 0.5, (real) 0.4, (real) 0.9}, (real) 0.5}, (real) 0.42},
        {{{(real) 2.75, (real) 1.48, (real) 3.44}, 1}, {{4, (real) 4.2, (real) 4.3}, {(real) 0.5, (real) 0.4, (real) 0.9}, (real) 0.5}, -(real) 0.25},
        {{{(real) 5.72, (real) 6.47, (real) 3.78}, 1}, {{4, (real) 4.2, (real) 4.3}, {(real) 0.5, (real) 0.4, (real) 0.9}, (real) 0.5}, (real) 1.4},
        {{{(real) 4.74, (real) 5.26, (real) 4.52}, 1}, {{4, (real) 4.2, (real) 4.3}, {(real) 0.5, (real) 0.4, (real) 0.9}, (real) 0.5}, -(real) 0.18},
};

// Capsule - Capsule collision tests
struct collision_test_caps {
  Capsule capsule1;
  Capsule capsule2;
  real    expected_dist;
};
collision_test_caps test_caps[] = {
        {{{(real) 0.5, (real) 0.95, (real) 6.73}, {4, (real) 4.9, (real) 3.51}, (real) 0.5}, {{4, (real) 4.2, (real) 4.3}, {(real) 0.5, (real) 0.4, (real) 0.9}, (real) 0.5}, -(real) 0.54},
        {{{(real) 5.16, (real) 4.9, (real) 4.37}, {(real) 3.55, (real) 0.95, (real) 8.84}, (real) 0.5}, {{4, (real) 4.2, (real) 4.3}, {(real) 0.5, (real) 0.4, (real) 0.9}, (real) 0.5}, (real) 0.16},
        {{{(real) 5.16, (real) 4.9, (real) 4.37}, {(real) 3.55, (real) 0.95, (real) 8.84}, (real) 0.5}, {{(real) 5.45, (real) 4.78, (real) 4.57}, {(real) 6.57, 1, -(real) 0.2}, (real) 0.5}, -(real) 0.66},
        {{{(real) 2.59, (real) 1.88, (real) 4.22}, {(real) 18.85, (real) 10.07, (real) 6.1}, 1}, {{(real) 0.55, (real) 21.32, (real) 3.08}, {(real) 5.07, (real) 3.62, (real) 4.19}, 1}, -(real) 1.44},
        {{{-(real) 9.26, (real) 12.81, (real) 7.98}, {(real) 6.63, (real) 4.62, (real) 4.04}, 1}, {{(real) 0.55, (real) 21.32, (real) 3.08}, {(real) 5.07, (real) 3.62, (real) 4.19}, 1}, -(real) 1.52},
};

// Capsule - Box collision tests
struct collision_test_box {
  Box     box;
  Capsule capsule;
  real    expected_dist = 0;
};
collision_test_box test_obb[] = {
        /*Test1*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{-(real) 1.21, -(real) 5.18, (real) 18.05}, {-(real) 3.89, (real) 8.59, (real) 6.3}, 1}, (real) 1.63659624},
        /*Test2*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 10.46, (real) 0.97, (real) 3.7}, {(real) 7.79, (real) 14.74, -(real) 8.04}, 1}, (real) 1.50169942},
        /*Test3*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 10.05, (real) 1.11, (real) 12.87}, {(real) 7.37, (real) 14.89, (real) 1.13}, 1}, (real) 0.33397901},
        /*Test4*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 4.82, (real) 6.9, (real) 12.84}, {(real) 4.7, -(real) 7.14, (real) 24.59}, 1}, (real) 4.00838054},
        /*Test5*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{-(real) 3.06, (real) 5.73, (real) 1.67}, {-(real) 0.39, -(real) 8.05, (real) 13.41}, 1}, (real) 0.89513819},
        /*Test6*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 4.97, (real) 2.79, (real) 12.23}, {(real) 2.29, (real) 16.57, (real) 0.48}, 1}, (real) 1.69981481},
        /*Test7*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 11.49, -(real) 0.06, (real) 8.32}, {(real) 14.17, -(real) 13.84, (real) 20.07}, 1}, (real) 0.49000000},
        /*Test8*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 6.76, -(real) 1.55, 8}, {(real) 9.44, -(real) 15.33, (real) 19.74}, 1}, (real) 0.45000000},
        /*Test9*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{-(real) 1.15, (real) 13.92, -(real) 7.48}, {(real) 1.53, (real) 0.14, (real) 4.27}, 1}, (real) 0.73046237},
        /*Test10*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 4.54, (real) 0.1, (real) 10.59}, {(real) 1.86, (real) 13.87, -(real) 1.15}, 1}, -1},
        /*Test11*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 10.76, (real) 0.28, (real) 8.08}, {(real) 13.44, -(real) 13.5, (real) 19.83}, 1}, -(real) 0.21961454},
        /*Test12*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 4.96, (real) 0.43, (real) 8.3}, {(real) 7.64, -(real) 13.35, (real) 20.05}, 1}, -(real) 1.53000000},
        /*Test13*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 4.48, -(real) 4.07, (real) 0.76}, {(real) 8.64, -(real) 4.41, (real) 18.58}, 1}, (real) 3.06923695},
        /*Test14*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 17.48, (real) 2.95, (real) 13.77}, {-(real) 0.82, (real) 3.11, (real) 13.77}, 1}, (real) 2.41054264},
        /*Test15*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {0, -1, 0, 1, 0, 0, 0, 0, 1}}, {{(real) 11.18, (real) 8.56, (real) 4.82}, {(real) 11.04, -(real) 6.44, (real) 15.29}, 1}, (real) 0.09912546},
        // Box.R and capsule.r changed
        /*Test16*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 2.94, -(real) 6.06, (real) 5.74}, {(real) 4.01, -(real) 9.86, (real) 0.98}, (real) 0.5}, (real) 1.00474115},
        /*Test17*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 1.64, (real) 9.52, (real) 11.7}, {(real) 2.71, (real) 5.72, (real) 6.94}, (real) 0.5}, (real) 0.22669533},
        /*Test18*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{3, (real) 5.69, (real) 6.2}, {(real) 4.07, (real) 1.89, (real) 1.44}, (real) 0.5}, (real) 0.42102531},
        /*Test19*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 3.1, (real) 4.23, (real) 6.22}, {(real) 4.17, (real) 0.43, (real) 1.46}, (real) 0.5}, (real) 0.10695205},
        /*Test20*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 6.35, (real) 8.57, (real) 12.85}, {(real) 7.42, (real) 4.77, (real) 8.09}, (real) 0.5}, (real) 0.85902522},
        /*Test21*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 1.88, -(real) 2.47, (real) 7.07}, {(real) 2.95, -(real) 6.27, (real) 2.31}, (real) 0.5}, (real) 0.37892454},
        /*Test22*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 1.8, (real) 3.83, (real) 14.57}, {(real) 2.87, (real) 0.03, (real) 9.81}, (real) 0.5}, (real) 1.47889394},
        /*Test23*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 4.28, (real) 6.32, (real) 8.33}, {(real) 3.21, (real) 10.12, (real) 13.09}, (real) 0.5}, (real) 0.82000000},
        /*Test24*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 2.36, (real) 2.3, (real) 6.31}, {(real) 3.43, -(real) 1.5, (real) 1.55}, (real) 0.5}, (real) 0.26887914},
        /*Test25*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 2.93, (real) 7.6, (real) 12.29}, {4, (real) 3.8, (real) 7.53}, (real) 0.5}, -(real) 0.93234019},
        /*Test26*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 7.55, (real) 0.92, (real) 10.24}, {(real) 6.48, (real) 4.72, 15}, (real) 0.5}, -(real) 0.32852683},
        /*Test27*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 6.79, 5, (real) 13.29}, {(real) 7.86, (real) 1.2, (real) 8.52}, (real) 0.5}, -(real) 0.40212621},
        /*Test28*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 5.66, -(real) 0.92, (real) 7.44}, {(real) 8.04, (real) 4.16, (real) 9.96}, (real) 0.5}, (real) 0.87078210},
        /*Test29*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 8.94, (real) 3.63, (real) 11.01}, {(real) 8.94, -(real) 2.55, (real) 11.01}, (real) 0.5}, (real) 1.24844065},
        /*Test30*/ {{{5, 0, 9}, {(real) 0.1, 5, 3}, {(real) 0.707, 0, -(real) 0.707, 0, 1, 0, (real) 0.707, 0, (real) 0.707}}, {{(real) 5.91, (real) 5.9, (real) 7.39}, {(real) 4.14, (real) 5.9, (real) 13.32}, (real) 0.5}, (real) 0.40000000},

        /*Test31*/ {{{0, 0, 0}, {(real) 1.421459, (real) 0.796365, (real) 0.574752}, {1, 0, 0, 0, 1, 0, 0, 0, 1}}, {{(real) 2.447825, (real) 6.30643376, -(real) 5.224298}, {-(real) 5.4486154, -(real) 5.149872, (real) 5.1825604}, 0}, -(real) 0.05877295},
        /*Test32*/ {{{0, 0, 0}, {(real) 1.7236252662212059, (real) 1.6611691154387149, (real) 0.20442662553572166}, {1, 0, 0, 0, 1, 0, 0, 0, 1}}, {{-(real) 4.6785661670641119, (real) 1.5869561202501099, -(real) 4.7958108710825522}, {(real) 6.9379851583657839, -(real) 2.4864207406657761, (real) 4.5569750611707658}, 0}, -(real) 0.438675},
        /*Test33*/ {{{0, 0, 0}, {(real) 1.7865788377999734, (real) 1.4615853859974308, (real) 1.7957946652540584}, {(real) 0.16899262598095025, (real) 0.69208966716438747, -(real) 0.70175022976009782, -(real) 0.73366055668912455, (real) 0.56377691511991157, (real) 0.37933860538637498, (real) 0.65816690886329632, (real) 0.45074103716231956, (real) 0.60303303184416857}}, {{-(real) 0.050034305603475548, (real) 0.40002630593895194, (real) 2.3444724393109695}, {(real) 1.3614165219982919, (real) 0.62636518062708535, -(real) 1.2242779633088015}, 0}, -(real) 1.3557069221709168},
};

TEST_CASE("Primitive distance functions test", "[World]") {
  using namespace blast;

  real TEST_EPSILON = (real) 1e-2;

  std::vector<Vec3> v1 = {{-(real) 4.6607382574035014, -(real) 3.4686124971115451, (real) 0.89349474531514073},
                          {(real) 6.4517093458638977, (real) 6.9310459212756506, -(real) 0.24683715666348305}};
  std::vector<Vec3> v2 = {{(real) 2.4373332611708229, (real) 2.5910471347625212, -(real) 0.04489407117738407},
                          {(real) 2.5547766616001675, (real) 4.0525983620611852, (real) 0.19409068785959360},
                          {(real) 1.8335195262912714, (real) 2.7351427068670993, -(real) 0.62940465909289878},
                          {(real) 1.9509629267206159, (real) 4.1966939341657632, -(real) 0.39041990005592114},
                          {(real) 1.6790895426140411, (real) 2.5264996514086859, (real) 0.72247776841327949},
                          {(real) 1.7965329430433856, (real) 3.9880508787073494, (real) 0.96146252745025707},
                          {(real) 1.0752758077344895, (real) 2.6705952235132639, (real) 0.13796718049776469},
                          {(real) 1.1927192081638340, (real) 4.1321464508119279, (real) 0.37695193953474238}};
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
