#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <blast>
#include "test_helper/test_functions.hpp"

TEST_CASE("ObjMatrix Copy Constructor", "[ObjMatrix]") {
  blast::ObjMatrix<int> mat1(2, 3);
  mat1(0, 0) = 1;
  mat1(0, 1) = 2;
  mat1(0, 2) = 3;
  mat1(1, 0) = 4;
  mat1(1, 1) = 5;
  mat1(1, 2) = 6;

  blast::ObjMatrix<int> mat2 = mat1;

  REQUIRE(mat2(0, 0) == 1);
  REQUIRE(mat2(0, 1) == 2);
  REQUIRE(mat2(0, 2) == 3);
  REQUIRE(mat2(1, 0) == 4);
  REQUIRE(mat2(1, 1) == 5);
  REQUIRE(mat2(1, 2) == 6);
}

TEST_CASE("ObjMatrix Move Constructor", "[ObjMatrix]") {
  blast::ObjMatrix<int> mat1(2, 3);
  mat1(0, 0) = 1;
  mat1(0, 1) = 2;
  mat1(0, 2) = 3;
  mat1(1, 0) = 4;
  mat1(1, 1) = 5;
  mat1(1, 2) = 6;

  blast::ObjMatrix<int> mat2(std::move(mat1));

  REQUIRE(mat2(0, 0) == 1);
  REQUIRE(mat2(0, 1) == 2);
  REQUIRE(mat2(0, 2) == 3);
  REQUIRE(mat2(1, 0) == 4);
  REQUIRE(mat2(1, 1) == 5);
  REQUIRE(mat2(1, 2) == 6);

  REQUIRE(mat1.cols == 0);
  REQUIRE(mat1.rows == 0);
  REQUIRE(mat1.size == 0);
  REQUIRE(mat1.data.size() == 0);
}

TEST_CASE("ObjMatrix Copy Assignment", "[ObjMatrix]") {
  blast::ObjMatrix<int> mat1(2, 3);
  mat1(0, 0) = 1;
  mat1(0, 1) = 2;
  mat1(0, 2) = 3;
  mat1(1, 0) = 4;
  mat1(1, 1) = 5;
  mat1(1, 2) = 6;

  blast::ObjMatrix<int> mat2(2, 3);
  mat2 = mat1;

  REQUIRE(mat2(0, 0) == 1);
  REQUIRE(mat2(0, 1) == 2);
  REQUIRE(mat2(0, 2) == 3);
  REQUIRE(mat2(1, 0) == 4);
  REQUIRE(mat2(1, 1) == 5);
  REQUIRE(mat2(1, 2) == 6);
}

TEST_CASE("ObjMatrix Move Assignment", "[ObjMatrix]") {
  blast::ObjMatrix<int> mat1(2, 3);
  mat1(0, 0) = 1;
  mat1(0, 1) = 2;
  mat1(0, 2) = 3;
  mat1(1, 0) = 4;
  mat1(1, 1) = 5;
  mat1(1, 2) = 6;

  blast::ObjMatrix<int> mat2(2, 3);
  mat2 = std::move(mat1);

  REQUIRE(mat2(0, 0) == 1);
  REQUIRE(mat2(0, 1) == 2);
  REQUIRE(mat2(0, 2) == 3);
  REQUIRE(mat2(1, 0) == 4);
  REQUIRE(mat2(1, 1) == 5);
  REQUIRE(mat2(1, 2) == 6);

  REQUIRE(mat1.cols == 0);
  REQUIRE(mat1.rows == 0);
  REQUIRE(mat1.size == 0);
  REQUIRE(mat1.data.size() == 0);
}

TEST_CASE("ObjMatrix Resize", "[ObjMatrix]") {
  blast::ObjMatrix<int> mat(2, 3);
  mat.resize(3, 2);

  REQUIRE(mat.rows == 3);
  REQUIRE(mat.cols == 2);
  REQUIRE(mat.size == 6);
}

TEST_CASE("ObjMatrix Capsules", "[ObjMatrix]") {

  blast::ObjMatrix<blast::Capsule> caps_matrix(7, 1);

  blast::Link6 manip;
  manip._efforts.resize(manip.joints);
  manip._Q.resize(manip.joints);
  manip._Q_mult.resize(manip.joints);
  manip._p_j.resize(manip.joints + 1);
  manip._capsule_list.resize(manip._n_caps);

  blast::Array pos = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  blast::forward_kinematics(manip, pos);

  manip.compute_capsules();

  blast::World world;
  blast::add_box({1.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world);

  for (blast::u32 i = 0; i < manip.joints; i++)
    caps_matrix(i, 0) = manip._capsule_list[i];

  auto distance = blast::test_collisions(caps_matrix, &world, 1, 0, 0);

  auto tmp_dist = 0.5 - manip._p_j[1].x - 0.065;
  REQUIRE(blast::is_close(distance[0], (blast::real) tmp_dist));
}
