#include <blast>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include "conversions.hpp"

namespace nb = nanobind;
using namespace blast;

// Helper: copy std::array<Vec3, N> into a Python list of Vec3
template<size_t N>
static nb::list vec3_array_to_list(const std::array<Vec3, N>& a, size_t n) {
  nb::list lst;
  for (size_t i = 0; i < n; i++)
    lst.append(a[i]);
  return lst;
}

// Helper: copy std::array<Mat3, N> into a Python list of Mat3
template<size_t N>
static nb::list mat3_array_to_list(const std::array<Mat3, N>& a, size_t n) {
  nb::list lst;
  for (size_t i = 0; i < n; i++)
    lst.append(a[i]);
  return lst;
}

void bind_manipulator(nb::module_& m) {
  // ------------------------------------------------------------------
  // CollisionModelCapsule
  // ------------------------------------------------------------------
  nb::class_<CollisionModelCapsule>(m, "CollisionModelCapsule")
          .def(nb::init<>())
          .def_rw("p1", &CollisionModelCapsule::p1,
                  "First endpoint of the capsule (Vec3).")
          .def_rw("p2", &CollisionModelCapsule::p2,
                  "Second endpoint of the capsule (Vec3).")
          .def_rw("joint_frame", &CollisionModelCapsule::joint_frame,
                  "Index of the joint frame this capsule is attached to.")
          .def_rw("radius", &CollisionModelCapsule::radius,
                  "Radius of the capsule.");

  // ------------------------------------------------------------------
  // ManipulatorLimits
  // ------------------------------------------------------------------
  nb::class_<ManipulatorLimits>(m, "ManipulatorLimits")
          .def(nb::init<>())
          // Getters return copies (capsule-owned) → rv_policy::move so nanobind
          // does not also try to apply reference_internal (which would conflict).
          .def_prop_rw("position_max", [](const ManipulatorLimits& l) { return array_to_copy(l.position_max); }, [](ManipulatorLimits& l, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { l.position_max = array_from_numpy(a); }, nb::rv_policy::move, "Maximum joint positions (rad).")
          .def_prop_rw("position_min", [](const ManipulatorLimits& l) { return array_to_copy(l.position_min); }, [](ManipulatorLimits& l, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { l.position_min = array_from_numpy(a); }, nb::rv_policy::move, "Minimum joint positions (rad).")
          .def_prop_rw("velocity_max", [](const ManipulatorLimits& l) { return array_to_copy(l.velocity_max); }, [](ManipulatorLimits& l, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { l.velocity_max = array_from_numpy(a); }, nb::rv_policy::move, "Maximum joint velocities (rad/s).")
          .def_prop_rw("acceleration_max", [](const ManipulatorLimits& l) { return array_to_copy(l.acceleration_max); }, [](ManipulatorLimits& l, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { l.acceleration_max = array_from_numpy(a); }, nb::rv_policy::move, "Maximum joint accelerations (rad/s²).")
          .def_prop_rw("torque_max", [](const ManipulatorLimits& l) { return array_to_copy(l.torque_max); }, [](ManipulatorLimits& l, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { l.torque_max = array_from_numpy(a); }, nb::rv_policy::move, "Maximum joint torques (N·m).")
          .def_rw("tool_speed_max", &ManipulatorLimits::tool_speed_max, "Maximum tool-center-point speed (m/s).");

  // ------------------------------------------------------------------
  // ManipulatorKinematics
  // ------------------------------------------------------------------
  nb::class_<ManipulatorKinematics>(m, "ManipulatorKinematics")
          .def(nb::init<>())
          .def_prop_rw("joint_offsets", [](const ManipulatorKinematics& k) { return vec3_array_to_list(k.joint_offsets, MAX_JOINTS); }, [](ManipulatorKinematics& k, nb::list lst) {
                size_t i = 0;
                for (auto item : lst) {
                    if (i >= (size_t)MAX_JOINTS)
                        throw nb::value_error("Too many joint offsets (max MAX_JOINTS)");
                    k.joint_offsets[i++] = nb::cast<Vec3>(item);
                } }, "Vectors from each joint to the next, in the joint frame (list of Vec3).")
          .def_prop_rw("joint_axes", [](const ManipulatorKinematics& k) { return vec3_array_to_list(k.joint_axes, MAX_JOINTS); }, [](ManipulatorKinematics& k, nb::list lst) {
                size_t i = 0;
                for (auto item : lst) {
                    if (i >= (size_t)MAX_JOINTS)
                        throw nb::value_error("Too many joint axes (max MAX_JOINTS)");
                    k.joint_axes[i++] = nb::cast<Vec3>(item);
                } }, "Unit rotation axes for each joint (list of Vec3).")
          .def_prop_rw("static_rotations", [](const ManipulatorKinematics& k) { return mat3_array_to_list(k.static_rotations, MAX_JOINTS + 1); }, [](ManipulatorKinematics& k, nb::list lst) {
                size_t i = 0;
                for (auto item : lst) {
                    if (i >= (size_t)(MAX_JOINTS + 1))
                        throw nb::value_error("Too many static rotations (max MAX_JOINTS+1)");
                    k.static_rotations[i++] = nb::cast<Mat3>(item);
                } }, "Static rotation matrices between consecutive joint frames (list of Mat3, length MAX_JOINTS+1).")
          .def_rw("base_position", &ManipulatorKinematics::base_position, "Base position in the world frame (Vec3).")
          .def_rw("first_joint_position", &ManipulatorKinematics::first_joint_position, "Position of the first joint relative to the base (Vec3).")
          .def_rw("base_rotation", &ManipulatorKinematics::base_rotation, "Base orientation in the world frame (Mat3).");

  // ------------------------------------------------------------------
  // ManipulatorDynamics
  // ------------------------------------------------------------------
  nb::class_<ManipulatorDynamics>(m, "ManipulatorDynamics")
          .def(nb::init<>())
          .def_prop_rw("link_masses", [](const ManipulatorDynamics& d) { return real_array_to_numpy(d.link_masses, MAX_JOINTS); }, [](ManipulatorDynamics& d, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { real_array_from_numpy(d.link_masses, a); }, nb::rv_policy::move, "Mass of each link (kg). Length must be <= MAX_JOINTS.")
          .def_prop_rw("inertia_tensors", [](const ManipulatorDynamics& d) { return mat3_array_to_list(d.inertia_tensors, MAX_JOINTS); }, [](ManipulatorDynamics& d, nb::list lst) {
                size_t i = 0;
                for (auto item : lst) {
                    if (i >= (size_t)MAX_JOINTS)
                        throw nb::value_error("Too many inertia tensors");
                    d.inertia_tensors[i++] = nb::cast<Mat3>(item);
                } }, "Inertia tensor for each link (list of Mat3).")
          .def_prop_rw("cog_offsets", [](const ManipulatorDynamics& d) { return vec3_array_to_list(d.cog_offsets, MAX_JOINTS); }, [](ManipulatorDynamics& d, nb::list lst) {
                size_t i = 0;
                for (auto item : lst) {
                    if (i >= (size_t)MAX_JOINTS)
                        throw nb::value_error("Too many CoG offsets");
                    d.cog_offsets[i++] = nb::cast<Vec3>(item);
                } }, "Center-of-gravity offset for each link (list of Vec3).");

  // ------------------------------------------------------------------
  // ManipulatorCapsules
  // ------------------------------------------------------------------
  nb::class_<ManipulatorCapsules>(m, "ManipulatorCapsules")
          .def(nb::init<>())
          .def_rw("base_sphere", &ManipulatorCapsules::base_sphere,
                  "Collision sphere around the robot base.")
          .def_rw("capsule_list", &ManipulatorCapsules::capsule_list,
                  "List of CollisionModelCapsule objects.")
          .def_prop_rw("collision_base", [](const ManipulatorCapsules& c) { return array_to_copy(c.collision_base); }, [](ManipulatorCapsules& c, nb::ndarray<nb::numpy, real, nb::ndim<1>, nb::c_contig> a) { c.collision_base = array_from_numpy(a); }, nb::rv_policy::move, "Boolean-like array: 1 if capsule[i] can collide with the base sphere.")
          .def_prop_rw("collision_matrix", [](const ManipulatorCapsules& c) -> nb::ndarray<nb::numpy, uint8_t> {
                auto& mat = c.collision_matrix;
                int n = mat.rows * mat.cols;
                uint8_t* buf = new uint8_t[n];
                std::memcpy(buf, mat.data.data(), (size_t)n);
                nb::capsule owner(buf, [](void* p) noexcept {
                    delete[] static_cast<uint8_t*>(p);
                });
                size_t  shape[2]   = {(size_t)mat.rows, (size_t)mat.cols};
                // column-major strides (uint8 = 1 byte)
                int64_t strides[2] = {1, (int64_t)mat.rows};
                return nb::ndarray<nb::numpy, uint8_t>(buf, 2, shape, owner, strides); }, [](ManipulatorCapsules& c, nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>> a) {
                int rows = (int)a.shape(0), cols = (int)a.shape(1);
                c.collision_matrix.resize(rows, cols);
                for (int r = 0; r < rows; r++)
                    for (int col = 0; col < cols; col++)
                        c.collision_matrix(r, col) = a(r, col); }, nb::rv_policy::move, "Symmetric adjacency matrix (uint8): collision_matrix[i,j]=1 means capsule i and j can collide.");

  // ------------------------------------------------------------------
  // Tool
  // ------------------------------------------------------------------
  nb::class_<Tool>(m, "Tool")
          .def(nb::init<>())
          .def_rw("tool_offset", &Tool::tool_offset, "Offset vector to the tool (Vec3).")
          .def_rw("tool_rotation", &Tool::tool_rotation, "Rotation to the tool frame (Mat3).")
          .def_rw("mass", &Tool::mass, "Tool mass (kg).")
          .def_rw("inertia_tensor", &Tool::inertia_tensor, "Tool inertia tensor (Mat3).")
          .def_rw("cog_offset", &Tool::cog_offset, "Tool center-of-gravity offset (Vec3).")
          .def_rw("capsule_list", &Tool::capsule_list, "Collision capsules for the tool.")
          .def_prop_rw("collision_matrix", [](const Tool& t) -> nb::ndarray<nb::numpy, uint8_t> {
                auto& mat = t.collision_matrix;
                int n = mat.rows * mat.cols;
                uint8_t* buf = new uint8_t[n];
                std::memcpy(buf, mat.data.data(), (size_t)n);
                nb::capsule owner(buf, [](void* p) noexcept {
                    delete[] static_cast<uint8_t*>(p);
                });
                size_t  shape[2]   = {(size_t)mat.rows, (size_t)mat.cols};
                int64_t strides[2] = {1, (int64_t)mat.rows};
                return nb::ndarray<nb::numpy, uint8_t>(buf, 2, shape, owner, strides); }, [](Tool& t, nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>> a) {
                int rows = (int)a.shape(0), cols = (int)a.shape(1);
                t.collision_matrix.resize(rows, cols);
                for (int r = 0; r < rows; r++)
                    for (int col = 0; col < cols; col++)
                        t.collision_matrix(r, col) = a(r, col); }, nb::rv_policy::move);

  // ------------------------------------------------------------------
  // Manipulator
  // ------------------------------------------------------------------
  nb::class_<Manipulator>(m, "Manipulator")
          .def(nb::init<u32, const ManipulatorLimits&, const ManipulatorKinematics&>(),
               nb::arg("n_joints"), nb::arg("limits"), nb::arg("kinematics"),
               "Construct a manipulator with kinematics only (no dynamics or collision model).")
          .def("__init__", [](Manipulator* self, u32 n, const ManipulatorLimits& l, const ManipulatorKinematics& k, const ManipulatorDynamics& d) { new (self) Manipulator(n, l, k, &d, nullptr); }, nb::arg("n_joints"), nb::arg("limits"), nb::arg("kinematics"), nb::arg("dynamics"), "Construct with kinematics and dynamics.")
          .def("__init__", [](Manipulator* self, u32 n, const ManipulatorLimits& l, const ManipulatorKinematics& k, const ManipulatorDynamics& d, const ManipulatorCapsules& c) { new (self) Manipulator(n, l, k, &d, &c); }, nb::arg("n_joints"), nb::arg("limits"), nb::arg("kinematics"), nb::arg("dynamics"), nb::arg("capsules"), "Construct with kinematics, dynamics, and collision capsules.")
          .def_ro("n_joints", &Manipulator::n_joints, "Number of joints.")
          .def("set_limits", &Manipulator::set_limits, nb::arg("limits"))
          .def("set_kinematics", &Manipulator::set_kinematics, nb::arg("kinematics"))
          .def("set_dynamics", &Manipulator::set_dynamics, nb::arg("dynamics"))
          .def("set_capsules", &Manipulator::set_capsules, nb::arg("capsules"))
          .def("add_tool", &Manipulator::add_tool, nb::arg("tool"))
          .def("set_payload", &Manipulator::set_payload, nb::arg("mass"), nb::arg("cog"), nb::arg("inertia"))
          // Expose a few internal fields for inspection
          .def_ro("base_position", &Manipulator::base_position)
          .def_ro("base_rotation", &Manipulator::base_rotation)
          .def_ro("has_tool", &Manipulator::has_tool)
          .def("__repr__", [](const Manipulator& manip) { return "<Manipulator n_joints=" + std::to_string(manip.n_joints) + ">"; });
}
