#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>
#include "test_helper.hpp"
#include "manipulator/UR5e.hpp"

TEST_CASE("Conversions Q & rx, ry, rz", "[Math][Misc]") {
  using namespace blast;
  int                 n_tests = 1;
  auto                manip   = make_UR5e();
  ManipulatorTempData data;

  for (u32 i = 0; i < n_tests; i++) {

    Array joint_pos(6);

    fill_random(joint_pos, PI);

    forward_kinematics(manip, data, joint_pos);
    auto rxyz_rot = rotation2rpy(data.rotations_mult[5]);
    auto rot_rxyz = rpy2rotation(rxyz_rot);

    CHECK(is_close(data.rotations_mult[5], rot_rxyz, 0.001));
  }
}

TEST_CASE("sign() - positive, negative, zero", "[Math][Misc]") {
  using namespace blast;
  CHECK(sign(1.0) == 1.0);
  CHECK(sign(-1.0) == -1.0);
  CHECK(sign(0.0) == 0.0);
  CHECK(sign(3.14) == 1.0);
  CHECK(sign(-0.001) == -1.0);
}

TEST_CASE("sign_no_zero() - zero returns 1", "[Math][Misc]") {
  using namespace blast;
  CHECK(sign_no_zero(1.0) == 1.0);
  CHECK(sign_no_zero(-1.0) == -1.0);
  CHECK(sign_no_zero(0.0) == 1.0);
}

TEST_CASE("wrap2pi() - scalar wrapping", "[Math][Misc]") {
  using namespace blast;
  CHECK(std::abs(wrap2pi(1.0) - 1.0) < blast::test::abs_tol);
  CHECK(std::abs(wrap2pi(-1.0) - (-1.0)) < blast::test::abs_tol);
  CHECK(std::abs(wrap2pi(0.0) - 0.0) < blast::test::abs_tol);
  CHECK(std::abs(wrap2pi(4.0) - (4.0 - 2 * PI)) < blast::test::abs_tol);
  CHECK(std::abs(wrap2pi(-4.0) - (-4.0 + 2 * PI)) < blast::test::abs_tol);
}

TEST_CASE("wrap2pi() - Array overload", "[Math][Misc]") {
  using namespace blast;
  Array a = {0.0, 4.0, -4.0};
  auto  b = wrap2pi(a);
  CHECK(std::abs(b[0] - 0.0) < blast::test::abs_tol);
  CHECK(std::abs(b[1] - (4.0 - 2 * PI)) < blast::test::abs_tol);
  CHECK(std::abs(b[2] - (-4.0 + 2 * PI)) < blast::test::abs_tol);
}

TEST_CASE("wrap_to_180() - real and float overloads", "[Math][Misc]") {
  using namespace blast;
  CHECK(std::abs(wrap_to_180(90.0) - 90.0) < blast::test::abs_tol);
  CHECK(std::abs(wrap_to_180(270.0) - (-90.0)) < blast::test::abs_tol);
  CHECK(std::abs(wrap_to_180(-270.0) - 90.0) < blast::test::abs_tol);
  CHECK(std::abs(wrap_to_180(90.0f) - 90.0f) < 1e-4f);
  CHECK(std::abs(wrap_to_180(270.0f) - (-90.0f)) < 1e-4f);
}

TEST_CASE("deg2rad() / rad2deg() - scalar", "[Math][Misc]") {
  using namespace blast;
  CHECK(std::abs(deg2rad(0.0) - 0.0) < blast::test::abs_tol);
  CHECK(std::abs(deg2rad(180.0) - PI) < blast::test::abs_tol);
  CHECK(std::abs(deg2rad(90.0) - PI / 2) < blast::test::abs_tol);
  CHECK(std::abs(rad2deg(0.0) - 0.0) < blast::test::abs_tol);
  CHECK(std::abs(rad2deg(PI) - 180.0) < blast::test::abs_tol);
  CHECK(std::abs(rad2deg(PI / 2) - 90.0) < blast::test::abs_tol);
  CHECK(std::abs(rad2deg(deg2rad(45.0)) - 45.0) < blast::test::abs_tol);
  CHECK(std::abs(deg2rad(rad2deg(1.0)) - 1.0) < blast::test::abs_tol);
}

TEST_CASE("deg2rad() / rad2deg() - Array overloads", "[Math][Misc]") {
  using namespace blast;
  Array angles_deg = {0.0, 90.0, 180.0, 360.0};
  Array angles_rad = deg2rad(angles_deg);
  CHECK(std::abs(angles_rad[0] - 0.0) < blast::test::abs_tol);
  CHECK(std::abs(angles_rad[1] - PI / 2) < blast::test::abs_tol);
  CHECK(std::abs(angles_rad[2] - PI) < blast::test::abs_tol);
  CHECK(std::abs(angles_rad[3] - 2 * PI) < blast::test::abs_tol);
  Array back = rad2deg(angles_rad);
  CHECK(std::abs(back[0] - 0.0) < blast::test::abs_tol);
  CHECK(std::abs(back[1] - 90.0) < blast::test::abs_tol);
  CHECK(std::abs(back[2] - 180.0) < blast::test::abs_tol);
  CHECK(std::abs(back[3] - 360.0) < blast::test::abs_tol);
}

TEST_CASE("clamp() - below min, in range, above max", "[Math][Misc]") {
  using namespace blast;
  CHECK(clamp(5.0, 0.0, 10.0) == 5.0);
  CHECK(clamp(-1.0, 0.0, 10.0) == 0.0);
  CHECK(clamp(11.0, 0.0, 10.0) == 10.0);
  CHECK(clamp(0.0, 0.0, 10.0) == 0.0);
  CHECK(clamp(10.0, 0.0, 10.0) == 10.0);
}

TEST_CASE("clamp_inplace() - modifies value and returns same reference", "[Math][Misc]") {
  using namespace blast;
  real  v1 = 5.0;
  real& r1 = clamp_inplace(v1, 0.0, 10.0);
  CHECK(v1 == 5.0);
  CHECK(&r1 == &v1);
  real v2 = -1.0;
  clamp_inplace(v2, 0.0, 10.0);
  CHECK(v2 == 0.0);
  real v3 = 11.0;
  clamp_inplace(v3, 0.0, 10.0);
  CHECK(v3 == 10.0);
}

TEST_CASE("rpy2rotation() - identity angles produce identity matrix", "[Math][Misc]") {
  using namespace blast;
  Mat3 R        = rpy2rotation({0.0, 0.0, 0.0});
  Mat3 expected = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  CHECK(is_close(R, expected, blast::test::abs_tol));
}

TEST_CASE("rpy2rotation() / rotation2rpy() - explicit roundtrip", "[Math][Misc]") {
  using namespace blast;
  Vec3 rpy_in  = {0.1, 0.2, 0.3};
  Mat3 R       = rpy2rotation(rpy_in);
  Vec3 rpy_out = rotation2rpy(R);
  CHECK(is_close(rpy_in, rpy_out, blast::test::abs_tol));
}

TEST_CASE("rpy2rotation() / rotation2rpy() - gimbal lock preserves rotation matrix", "[Math][Misc]") {
  using namespace blast;
  // At pitch = pi/2 (gimbal lock), roll/yaw decomposition is degenerate.
  // The round-trip rotation matrix should still be identical.
  Vec3 rpy_in  = {0.0, PI / 2, 0.5};
  Mat3 R       = rpy2rotation(rpy_in);
  Vec3 rpy_out = rotation2rpy(R);
  Mat3 R2      = rpy2rotation(rpy_out);
  CHECK(is_close(R, R2, 1e-9));
}
