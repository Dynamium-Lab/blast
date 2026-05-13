#include <blast>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include "conversions.hpp"

namespace nb = nanobind;
using namespace blast;

void bind_world(nb::module_& m) {
  // ------------------------------------------------------------------
  // Box
  // ------------------------------------------------------------------
  nb::class_<Box>(m, "Box")
          .def(nb::init<>())
          .def_rw("center", &Box::center,
                  "Center point of the box (Vec3).")
          .def_rw("extents", &Box::extents,
                  "Positive half-width extents along each axis (Vec3).")
          .def_rw("rotation", &Box::rotation,
                  "Rotation matrix: local x/y/z-axes (Mat3).")
          .def("__repr__", [](const Box& b) {
            return "<Box center=(" + std::to_string(b.center.x) + ", " +
                   std::to_string(b.center.y) + ", " + std::to_string(b.center.z) + ")>";
          });

  // ------------------------------------------------------------------
  // Sphere
  // ------------------------------------------------------------------
  nb::class_<Sphere>(m, "Sphere")
          .def(nb::init<>())
          .def_rw("center", &Sphere::center,
                  "Center of the sphere (Vec3).")
          .def_rw("radius", &Sphere::radius,
                  "Radius of the sphere.")
          .def("__repr__", [](const Sphere& s) {
            return "<Sphere center=(" + std::to_string(s.center.x) + ", " +
                   std::to_string(s.center.y) + ", " + std::to_string(s.center.z) +
                   ") radius=" + std::to_string(s.radius) + ">";
          });

  // ------------------------------------------------------------------
  // Capsule
  // ------------------------------------------------------------------
  nb::class_<Capsule>(m, "Capsule")
          .def(nb::init<>())
          .def_rw("p1", &Capsule::p1, "First endpoint (Vec3).")
          .def_rw("p2", &Capsule::p2, "Second endpoint (Vec3).")
          .def_rw("radius", &Capsule::radius, "Radius of the capsule.")
          .def("__repr__", [](const Capsule& c) {
            return "<Capsule radius=" + std::to_string(c.radius) + ">";
          });

  // ------------------------------------------------------------------
  // DynamicBox
  // ------------------------------------------------------------------
  nb::class_<DynamicBox>(m, "DynamicBox")
          .def(nb::init<>())
          .def_rw("n_points", &DynamicBox::n_points)
          .def_rw("start_time", &DynamicBox::start_time)
          .def_rw("end_time", &DynamicBox::end_time)
          .def_rw("trajectory", &DynamicBox::trajectory,
                  "List of Box poses, one per interpolation point.")
          .def("lookup", &DynamicBox::lookup, nb::arg("time"),
               "Interpolate the box pose at the given time.");

  // ------------------------------------------------------------------
  // DynamicSphere
  // ------------------------------------------------------------------
  nb::class_<DynamicSphere>(m, "DynamicSphere")
          .def(nb::init<>())
          .def_rw("n_points", &DynamicSphere::n_points)
          .def_rw("start_time", &DynamicSphere::start_time)
          .def_rw("end_time", &DynamicSphere::end_time)
          .def_rw("trajectory", &DynamicSphere::trajectory)
          .def("lookup", &DynamicSphere::lookup, nb::arg("time"));

  // ------------------------------------------------------------------
  // DynamicCapsule
  // ------------------------------------------------------------------
  nb::class_<DynamicCapsule>(m, "DynamicCapsule")
          .def(nb::init<>())
          .def_rw("n_points", &DynamicCapsule::n_points)
          .def_rw("start_time", &DynamicCapsule::start_time)
          .def_rw("end_time", &DynamicCapsule::end_time)
          .def_rw("trajectory", &DynamicCapsule::trajectory)
          .def("lookup", &DynamicCapsule::lookup, nb::arg("time"));

  // ------------------------------------------------------------------
  // World
  // ------------------------------------------------------------------
  nb::class_<World>(m, "World")
          .def(nb::init<>())
          .def_ro("size", &World::size, "Total number of static obstacle primitives.")
          .def_ro("boxes", &World::boxes, "List of static Box obstacles.")
          .def_ro("spheres", &World::spheres, "List of static Sphere obstacles.")
          .def_ro("capsules", &World::capsules, "List of static Capsule obstacles.")
          // Static obstacles
          .def("add_box",
               nb::overload_cast<const Box&>(&World::add_box),
               nb::arg("box"), "Add a Box obstacle.")
          .def("add_box",
               nb::overload_cast<Vec3, Vec3, Mat3>(&World::add_box),
               nb::arg("center"), nb::arg("half_extents"), nb::arg("rotation"),
               "Add a box specified by center, half-extents and rotation matrix.")
          .def("add_sphere",
               nb::overload_cast<const Sphere&>(&World::add_sphere),
               nb::arg("sphere"), "Add a Sphere obstacle.")
          .def("add_sphere",
               nb::overload_cast<Vec3, real>(&World::add_sphere),
               nb::arg("center"), nb::arg("radius"),
               "Add a sphere specified by center and radius.")
          .def("add_capsule",
               nb::overload_cast<const Capsule&>(&World::add_capsule),
               nb::arg("capsule"), "Add a Capsule obstacle.")
          .def("add_capsule",
               nb::overload_cast<Vec3, Vec3, real>(&World::add_capsule),
               nb::arg("p1"), nb::arg("p2"), nb::arg("radius"),
               "Add a capsule specified by two endpoints and radius.")
          // Dynamic obstacles
          .def("add_dynamic_box",
               nb::overload_cast<const DynamicBox&>(&World::add_dynamic_box),
               nb::arg("box"))
          .def("add_dynamic_sphere",
               nb::overload_cast<const DynamicSphere&>(&World::add_dynamic_sphere),
               nb::arg("sphere"))
          .def("add_dynamic_capsule",
               nb::overload_cast<const DynamicCapsule&>(&World::add_dynamic_capsule),
               nb::arg("capsule"))
          .def("__repr__", [](const World& w) {
            return "<World boxes=" + std::to_string(w.boxes.size()) +
                   " spheres=" + std::to_string(w.spheres.size()) +
                   " capsules=" + std::to_string(w.capsules.size()) + ">";
          });
}
