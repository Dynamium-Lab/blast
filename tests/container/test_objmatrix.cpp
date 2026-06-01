#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <blast>
#include "blast_utilities.hpp"
#include "manipulator/UR5e.hpp"
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

  blast::Manipulator manip = blast::make_UR5e();

  blast::Array               pos = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  blast::ManipulatorTempData data;

  blast::forward_kinematics(manip, data, pos);
  blast::compute_capsules(manip, data);

  blast::World world;
  world.add_box({1.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {1, 0, 0, 0, 1, 0, 0, 0, 1});

  for (blast::u32 i = 0; i < manip.n_joints; i++)
    caps_matrix(i, 0) = data.capsule_list[i];

  auto distance = blast::test_collisions(caps_matrix, &world, 7, 0, 0);

  auto tmp_dist = 0.5 - data.p_j[1].x - data.capsule_list[0].radius;
  CHECK(blast::is_close(distance[0], (blast::real) tmp_dist));
}
