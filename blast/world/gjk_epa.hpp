#pragma once

#include <blast>
#include <queue>
#include "tracy/Tracy.hpp"

namespace blast {

#define COLLISION_EPSILON 1e-9

struct EPA_dist_norm {
  Vec3 normal;
  real distance;
};

struct closeface {
  real distance;
  Vec3 normal;
  int  index;
};

// Ericson structs
struct Simplex {
  Vec3 P;
  Vec3 a;
  Vec3 b;
  Vec3 c;
  Vec3 d;
  int  size;
};

struct PolytopeFace {
  Vec3 p1;
  Vec3 p2;
  Vec3 p3;
  Vec3 n;
  real dist = 0.0;
};

host_fn Vec3 ptint(Segment seg, Vec3 point) {
  Vec3 ab = seg.p2 - seg.p1;
  real t  = dot(point - seg.p1, ab) / dot(ab, ab);

  t = clamp(t, 0, 1);

  Vec3 d = seg.p1 + t * ab;
  return d;
}

host_fn EPA_dist_norm distmin_origin(Triangle tri) {
  EPA_dist_norm min_dist;
  Vec3          o  = tri.p1;
  Vec3          v1 = tri.p2 - tri.p1;
  Vec3          v2 = tri.p3 - tri.p1;
  Vec3          n  = cross(v1, v2);
  if (dot(n, o) < 0)
    n = -n;

  real distmin = INF_REAL;
  if (point_in_triangle(v1, v2, o, {0, 0, 0})) {
    min_dist.distance = dot(o, 1 / norm(n) * n);
    min_dist.normal   = n;
  }

  Segment seg[3];

  seg[0].p1 = tri.p1;
  seg[0].p2 = tri.p2;
  seg[1].p1 = tri.p2;
  seg[1].p2 = tri.p3;
  seg[2].p1 = tri.p3;
  seg[2].p2 = tri.p1;

  for (int i = 0; i < 3; i++) {
    Vec3 d = closest_point_origin(seg[i]);
    if (dot(d, d) < distmin * distmin) {
      min_dist.distance = norm(d);
      min_dist.normal   = d;
    }
  }
  return min_dist;
}

host_fn std::vector<real> test_collision(std::vector<Capsule>& caps_list, World* world, int n_var) {
  std::vector<real> dist_min(n_var, INF_REAL);
  real              dist;

  for (int c = 0; c < caps_list.size(); c++) {
    // box
    for (int i = 0; i < world->boxes.size(); i++) {
      dist = distance(caps_list[c], world->boxes[i]);
      for (int j = 0; j < n_var; j++) {
        if (dist < dist_min[j]) {
          for (int k = n_var; k > j; k--) {
            dist_min[k] = dist_min[k - 1];
          }
          dist_min[j] = dist;
          break;
        }
      }
    }

    // capsule
    for (int i = 0; i < world->capsules.size(); i++) {
      dist = distance(caps_list[c], world->capsules[i]);
      for (int j = 0; j < n_var; j++) {
        if (dist < dist_min[j]) {
          for (int k = n_var; k > j; k--) {
            dist_min[k] = dist_min[k - 1];
          }
          dist_min[j] = dist;
          break;
        }
      }
    }

    // sphere
    for (int i = 0; i < world->spheres.size(); i++) {
      dist = distance(caps_list[c], world->spheres[i]);
      for (int j = 0; j < n_var; j++) {
        if (dist < dist_min[j]) {
          for (int k = n_var; k > j; k--) {
            dist_min[k] = dist_min[k - 1];
          }
          dist_min[j] = dist;
          break;
        }
      }
    }
  }
  return dist_min;
}

//    note: this is a simpler version of the collision test which only allows for one minimum distance to be returned.
//    It could easily be adapted to include all distances as well.
//    real dist_min;
//
//     for (int i = 0; i < size(world.OBBlist); i++) {
//         real dist = distmin(caps, world.OBBlist[i]);
//         if (dist < dist_min)
//             dist_min = dist;
//     }
//
//     for (int i = 0; i < size(world.capslist); i++) {
//         real dist = distmin(caps, world.capslist[i]);
//         if (dist < dist_min)
//             dist_min = dist;
//     }
//
//     for (int i = 0; i < size(world.sphlist); i++) {
//         real dist = distmin(caps, world.sphlist[i]);
//         if (dist < dist_min)
//             dist_min = dist;
//     }
//
//     for (int i = 0; i < size(world.cyllist); i++) {
//         real dist = distmin(caps, world.cyllist[i]);
//         if (dist < dist_min)
//             dist_min = dist;
//     }
// real n_pts = size(robot.pts);
// real n_link = n_pts - 1;
// capsule link[n_link];
//
// for (int i = 0; i < n_pts - 1; i++) {
//     link[i].p1 = robot.pts[i];
//     link[i].p2 = robot.pts[i+1];
//     link[i].r = robot.r[i];
// }
//
// for (int i = 0; i < n_link; i++) {
//     for (int j = 2 + i; j < n_link; j++) {
//         dist = distmin(link[i], link[i+j]);
//         if (dist < dist_min)
//             dist_min = dist;
//     }
// }
//
//     return dist_min;
// }

// ======================================
//            GJK algorithm
// ======================================
// This version of the GJK algorithm is adapted from : https://www.youtube.com/watch?v=DGVZYdlw_uo

struct ComplexSimplex {
  Vec3 a;
  Vec3 a1;
  Vec3 a2;

  Vec3 b;
  Vec3 b1;
  Vec3 b2;

  Vec3 c;
  Vec3 c1;
  Vec3 c2;

  Vec3 d;
  Vec3 d1;
  Vec3 d2;

  int count;
};

struct KodaVertices {
  int               count;
  std::vector<Vec3> fixed; // points
};

struct GJKResult {
  bool           intersection;
  ComplexSimplex final_simplex;
  Vec3           A_closept;
  Vec3           B_closept;
  real           minimal_distance;
};

host_fn bool GJK_same_direction(Vec3 a, Vec3 b) {
  return dot(a, b) > 0;
}

host_fn real GJK_triangle_area_2d(real x1, real y1, real x2, real y2, real x3, real y3) {
  return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

host_fn Vec3 GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 c, Vec3 p) {
  Vec3 ba  = b - a;
  Vec3 ca  = c - a;
  Vec3 abc = cross(ba, ca);
  real nu;
  real nv;
  real ood;

  real x = std::abs(abc.x);
  real y = std::abs(abc.y);
  real z = std::abs(abc.z);

  if ((x >= y) && (x >= z)) {
    // x is the largest so project onto the yz plane
    nu  = GJK_triangle_area_2d(p.y, p.z, b.y, b.z, c.y, c.z);
    nv  = GJK_triangle_area_2d(p.y, p.z, c.y, c.z, a.y, a.z);
    ood = 1 / abc.x;
  } else if ((y >= x) && (y >= z)) {
    // y is the largest so project onto the xz plane
    nu  = GJK_triangle_area_2d(p.x, p.z, b.x, b.z, c.x, c.z);
    nv  = GJK_triangle_area_2d(p.x, p.z, c.x, c.z, a.x, a.z);
    ood = 1 / -abc.y;
  } else {
    // z is the largest so project onto the xy plane
    nu  = GJK_triangle_area_2d(p.x, p.y, b.x, b.y, c.x, c.y);
    nv  = GJK_triangle_area_2d(p.x, p.y, c.x, c.y, a.x, a.y);
    ood = 1 / abc.z;
  }

  Assert(isinf(ood) == false);

  Vec3 result;
  result.x = nu * ood;
  result.y = nv * ood;
  result.z = 1 - result.x - result.y;
  return result;
}

host_fn real GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 p) {
  Vec3 ab = b - a;
  Vec3 ap = p - a;

  real t = dot(ap, ab) / dot(ab, ab);
  return t;
}

host_fn Vec3 GJK_convert_cartesian(Vec3 a, Vec3 b, real barycentric) {
  Vec3 result;
  result.x = a.x + barycentric * (b.x - a.x);
  result.y = a.y + barycentric * (b.y - a.y);
  result.z = a.z + barycentric * (b.z - a.z);
  return result;
}

host_fn Vec3 GJK_convert_cartesian(Vec3 a, Vec3 b, Vec3 c, Vec3 barycentric) {
  Vec3 result;
  result.x = a.x * barycentric.x + b.x * barycentric.y + c.x * barycentric.z;
  result.y = a.y * barycentric.x + b.y * barycentric.y + c.y * barycentric.z;
  result.z = a.z * barycentric.x + b.z * barycentric.y + c.z * barycentric.z;
  return result;
}

host_fn TwoPts GJK_get_local_points(ComplexSimplex simplex, Vec3 p) {
  TwoPts loc;
  real   barycentric;
  Vec3   barycentric3;

  switch (simplex.count) {
    case 1:
      loc.p1 = simplex.a1;
      loc.p2 = simplex.a2;
      break;

    case 2:
      barycentric = GJK_convert_barycentric(simplex.a, simplex.b, p);
      loc.p1      = GJK_convert_cartesian(simplex.a1, simplex.b1, barycentric);
      loc.p2      = GJK_convert_cartesian(simplex.a2, simplex.b2, barycentric);
      break;

    case 3:
      barycentric3 = GJK_convert_barycentric(simplex.a, simplex.b, simplex.c, p);
      loc.p1       = GJK_convert_cartesian(simplex.a1, simplex.b1, simplex.c1, barycentric3);
      loc.p2       = GJK_convert_cartesian(simplex.a2, simplex.b2, simplex.c2, barycentric3);
      Assert(isnan(loc.p1.x) == false);
      break;

    default:
      Assert(false /*"invalid simplex dimension for nearest feature resolution."*/);
  }
  Assert(isnan(loc.p1.x) == false);
  return loc;
}

host_fn Vec3 GJK_get_support(const std::vector<Vec3>& vertices, Vec3 direction) {
  real largest_dot    = dot(vertices[0], direction);
  Vec3 largest_vertex = vertices[0];
  for (auto& vertex: vertices) {
    const real current_dot = dot(vertex, direction);
    if (current_dot > largest_dot) {
      largest_vertex = vertex;
      largest_dot    = current_dot;
    }
  }
  Assert(isinf(largest_dot) == false);
  return largest_vertex;
}

host_fn TwoPts GJK_solve_simplex2(ComplexSimplex& simplex) { // finds next search direction (ao) and closest point on line (ab) for 2pt simplex
  Vec3 ab = simplex.b - simplex.a;
  Vec3 ao = -simplex.a;

  const auto d = dot(ab, ao);
  if (d > 0) {
    simplex.count = 2;
    ao            = cross(cross(ab, ao), ab);
    ab            = simplex.a + d / dot(ab, ab) * ab;
  } else {
    simplex.count = 1;
    ab            = simplex.a;
  }

  return {ab, ao};
}

host_fn TwoPts GJK_solve_simplex3(ComplexSimplex& simplex) {
  Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
  Vec3 ac  = simplex.c - simplex.a;
  Vec3 ao  = -simplex.a;

  if (GJK_same_direction(cross(abc, ac), ao)) {
    // the origin is in the direction of the triangle normal
    if (GJK_same_direction(ac, ao)) {
      // the origin is nearest to the line ac
      // simplex c remains as c
      // simplex.b_all = simplex.c_all;
      simplex.b     = simplex.c;
      simplex.b1    = simplex.c1;
      simplex.b2    = simplex.c2;
      simplex.count = 2;
      real t        = dot(ao, ac) / dot(ac, ac);
      return {simplex.a + t * ac, cross(cross(ac, ao), ac)};
    } else {
      Vec3 ab = simplex.b - simplex.a;
      if (GJK_same_direction(ab, ao)) {
        // origin is closest to the line ab
        simplex.count = 2;
        real t        = dot(ao, ab) / dot(ab, ab);
        return {simplex.a + t * ab, cross(cross(ab, ao), ab)};
      } else {
        // the origin is nearest to the point a
        simplex.count = 1;
        return {simplex.a, ao};
      }
    }
  } else {
    Vec3 ab = simplex.b - simplex.a;
    if (GJK_same_direction(cross(ab, abc), ao)) {
      if (GJK_same_direction(ab, ao)) {
        // the origin is closest to line ab
        simplex.count = 2;
        real t        = dot(ao, ab) / dot(ab, ab);
        return {simplex.a + t * ab, cross(cross(ab, ao), ab)};
      } else {
        // the origin is nearest to the point a
        simplex.count = 1;
        return {simplex.a, ao};
      }
    } else {
      if (GJK_same_direction(abc, ao)) {
        // the origin is closest to the triangle abc
        Vec3 bo    = -simplex.b;
        Vec3 co    = -simplex.c;
        real d1    = dot(ab, ao);
        real d2    = dot(ac, ao);
        real d3    = dot(ab, bo);
        real d4    = dot(ac, bo);
        real d5    = dot(ab, co);
        real d6    = dot(ac, co);
        real va    = d3 * d6 - d5 * d4;
        real vb    = d5 * d2 - d1 * d6;
        real vc    = d1 * d4 - d3 * d2;
        real denom = 1 / (va + vb + vc);
        real v     = vb * denom;
        real w     = vc * denom;
        return {simplex.a + ab * v + ac * w, abc};
      } else {
        // the origin is closest to the triangle acb (other side of the triangle)
        // simplex.b_all, simplex.c_all = swap(simplex.b_all, simplex.c_all);
        Vec3 sb  = simplex.c;
        Vec3 sb1 = simplex.c1;
        Vec3 sb2 = simplex.c2;

        simplex.c  = simplex.b;
        simplex.c1 = simplex.b1;
        simplex.c2 = simplex.b2;

        simplex.b  = sb;
        simplex.b1 = sb1;
        simplex.b2 = sb2;

        Vec3 bo    = -simplex.b;
        Vec3 co    = -simplex.c;
        real d1    = dot(ac, ao);
        real d2    = dot(ab, ao);
        real d3    = dot(ac, co);
        real d4    = dot(ab, co);
        real d5    = dot(ac, bo);
        real d6    = dot(ab, bo);
        real va    = d3 * d6 - d5 * d4;
        real vb    = d5 * d2 - d1 * d6;
        real vc    = d1 * d4 - d3 * d2;
        real denom = 1 / (va + vb + vc);
        real v     = vb * denom;
        real w     = vc * denom;
        return {simplex.a + ab * v + ac * w, -abc};
      }
    }
  }
}

host_fn TwoPts GJK_solve_simplex4(ComplexSimplex& simplex) {
  Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
  Vec3 acd = cross(simplex.c - simplex.a, simplex.d - simplex.a);
  Vec3 adb = cross(simplex.d - simplex.a, simplex.b - simplex.a);

  Vec3 ao = -simplex.a;

  real abc_dir = GJK_same_direction(abc, ao);
  real acd_dir = GJK_same_direction(acd, ao);
  real adb_dir = GJK_same_direction(adb, ao);

  if (!abc_dir && !acd_dir && !adb_dir) {
    // the origin is inside the simplex
    return {{0, 0, 0}, {0, 0, 0}};
  }
  if (abc_dir && !acd_dir && !adb_dir) {
    // the origin is near abc
    simplex.count = 3;
    return GJK_solve_simplex3(simplex);
  }
  if (!abc_dir && acd_dir && !adb_dir) {
    // the origin is near acd
    // simplex.b_all = simplex.c_all;
    simplex.b  = simplex.c;
    simplex.b1 = simplex.c1;
    simplex.b2 = simplex.c2;
    // simplex.c_all = simplex.d_all;
    simplex.c  = simplex.d;
    simplex.c1 = simplex.d1;
    simplex.c2 = simplex.d2;

    simplex.count = 3;
    return GJK_solve_simplex3(simplex);
  }
  if (!abc_dir && !acd_dir && adb_dir) {
    // the origin is near adb
    // simplex.c_all = simplex.b_all;
    simplex.c  = simplex.b;
    simplex.c1 = simplex.b1;
    simplex.c2 = simplex.b2;

    // simplex.b_all = simplex.d_all;
    simplex.b  = simplex.d;
    simplex.b1 = simplex.d1;
    simplex.b2 = simplex.d2;

    simplex.count = 3;
    return GJK_solve_simplex3(simplex);
  }

  // the origin potentially falls on multiple triangles
  ComplexSimplex simplex_abc = simplex;
  simplex_abc.count          = 3;

  // There may be a better way to do this.

  ComplexSimplex simplex_acd = simplex;
  // simplex_acd.b_all = simplex_acd.c_all;
  simplex_acd.b  = simplex.c;
  simplex_acd.b1 = simplex.c1;
  simplex_acd.b2 = simplex.c2;
  // simplex_acd.c_all = simplex_acd.d_all;
  simplex_acd.c  = simplex.d;
  simplex_acd.c1 = simplex.d1;
  simplex_acd.c2 = simplex.d2;

  simplex_acd.count = 3;

  ComplexSimplex simplex_adb = simplex;
  // simplex_adb.c_all = simplex_adb.b_all;
  simplex_adb.c  = simplex.b;
  simplex_adb.c1 = simplex.b1;
  simplex_adb.c2 = simplex.b2;
  // simplex_adb.c_all = simplex_adb.d_all;
  simplex_adb.b  = simplex.d;
  simplex_adb.b1 = simplex.d1;
  simplex_adb.b2 = simplex.d2;

  simplex_adb.count = 3;

  TwoPts solved;
  solved       = GJK_solve_simplex3(simplex_abc);
  Vec3 p_abc   = solved.p1;
  Vec3 dir_abc = solved.p2;
  solved       = GJK_solve_simplex3(simplex_adb);
  Vec3 p_adb   = solved.p1;
  Vec3 dir_adb = solved.p2;
  solved       = GJK_solve_simplex3(simplex_acd);
  Vec3 p_acd   = solved.p1;
  Vec3 dir_acd = solved.p2;

  real abc_d2 = dot(p_abc, p_abc);
  real acd_d2 = dot(p_acd, p_acd);
  real adb_d2 = dot(p_adb, p_adb);

  if ((abc_d2 <= acd_d2) && (abc_d2 <= adb_d2)) {
    simplex = simplex_abc;
    return {p_abc, dir_abc};
  } else if ((acd_d2 <= abc_d2) && (acd_d2 <= adb_d2)) {
    simplex = simplex_acd;
    return {p_acd, dir_acd};
  } else if ((adb_d2 <= abc_d2) && (adb_d2 <= acd_d2)) {
    simplex = simplex_adb;
    return {p_adb, dir_adb};
  }

  Assert(false); // This souldn't be reachable

  // the origin isn't outside of any plane, so it's inside the tetrahedron
  return {{0, 0, 0}, {0, 0, 0}};
}

// ======================================
//            EPA algorithm
// ======================================
// This version of the EPA algorithm is adapted from : https://dyn4j.org/2010/05/epa-expanding-polytope-algorithm/

host_fn bool check_same_point(Vec3 p1, Vec3 p2) {
  return ((p1.x - p2.x) * (p1.x - p2.x) < COLLISION_EPSILON && (p1.y - p2.y) * (p1.y - p2.y) < COLLISION_EPSILON && (p1.z - p2.z) * (p1.z - p2.z) < COLLISION_EPSILON);
}

host_fn bool check_same_triangle(Triangle tri1, Triangle tri2) {
  if (check_same_point(tri1.p1, tri2.p1) == 0 && check_same_point(tri1.p1, tri2.p2) == 0 && check_same_point(tri1.p1, tri2.p3) == 0)
    return 0;
  if (check_same_point(tri1.p2, tri2.p1) == 0 && check_same_point(tri1.p2, tri2.p2) == 0 && check_same_point(tri1.p2, tri2.p3) == 0)
    return 0;
  if (check_same_point(tri1.p3, tri2.p1) == 0 && check_same_point(tri1.p3, tri2.p2) == 0 && check_same_point(tri1.p3, tri2.p3) == 0)
    return 0;
  return 1;
}

host_fn closeface find_closest_face(std::vector<Triangle> faces) {
  closeface closest;
  closest.distance = INF_REAL;
  for (int i = 0; i < faces.size(); i++) {
    Triangle current_triangle;
    current_triangle.p1 = faces[i].p1;
    current_triangle.p2 = faces[i].p2;
    current_triangle.p3 = faces[i].p3;

    Vec3 vec1 = current_triangle.p2 - current_triangle.p1;
    Vec3 vec2 = current_triangle.p3 - current_triangle.p1;

    Vec3 origin = {0, 0, 0};
    Vec3 n      = cross(vec1, vec2);
    if (norm(n) < COLLISION_EPSILON)
      continue; // use this for now.. consider switching back to break;
    if (dot(n, current_triangle.p1) < 0)
      n = -n;

    EPA_dist_norm distance = distmin_origin(current_triangle);
    if (distance.distance < closest.distance) {
      closest.distance = distance.distance;
      closest.normal   = distance.normal;
      closest.index    = i;
    }
  }
  return closest;
}

// there is a bug in this EPA !!
// host_fn real solve_EPA(ComplexSimplex simplex, std::vector<Vec3> v1, std::vector<Vec3> v2) {
//   std::vector<Triangle> s{4}; // convex hull (no longer simplex)
//   s[0] = {simplex.a, simplex.b, simplex.c};
//   s[1] = {simplex.a, simplex.b, simplex.d};
//   s[2] = {simplex.a, simplex.c, simplex.d};
//   s[3] = {simplex.b, simplex.c, simplex.d};
//
//   while (true) {
//     // obtain the triangle closest to origin on polytope s
//     closeface e = find_closest_face(s);
//
//     // find support point further in direction of triangle normal
//     Vec3 support1 = GJK_get_support(v1, -e.normal);
//     Vec3 support2 = GJK_get_support(v2, e.normal);
//     Vec3 p        = support2 - support1;
//
//     // check support point (whether it's further than current closest triangle)
//     real d = dot(p, e.normal);
//     if (std::abs(d) - std::abs(e.distance) < COLLISION_EPSILON) { // COLLISION_EPSILON might be too small (was previously 1e-2)
//       // polytope cannot be expanded any further
//       return -std::abs(d);
//     } else {
//       // polytope can be expanded so add support point
//       // delete the triangle which is currently the closest and create three triangles that are inserted in its place
//
//       Vec3 p1 = s[e.index].p1;
//       Vec3 p2 = s[e.index].p2;
//       Vec3 p3 = s[e.index].p3;
//
//       Triangle new_tri1 = {p1, p2, p};
//       Triangle new_tri2 = {p1, p3, p};
//       Triangle new_tri3 = {p2, p3, p};
//
//       // delete old closest triangle
//       s.erase(s.begin() + e.index);
//
//       // in the case where two iterations return the same support point, it can be the case that two
//       // faces are created on top of one another. In this case, the face is always inside the polytope and it should therefore be deleted.
//       for (int i = size(s) - 1; i >= 0; i--) {
//         if (check_same_triangle(s[i], new_tri1)) {
//           s.erase(s.begin() + i); // if the triangle was already in the list, then we must delete it as this face is now inside the simplex
//           break;
//         } else if (i == 0) {
//           s.push_back(new_tri1); // if it is a new triangle, add it to the list
//           break;
//         }
//       }
//
//       for (int i = size(s) - 1; i >= 0; i--) {
//         if (check_same_triangle(s[i], new_tri2)) {
//           s.erase(s.begin() + i);
//           break;
//         } else if (i == 0) {
//           s.push_back(new_tri2); // if it is a new triangle, add it to the list
//           break;
//         }
//       }
//
//       for (int i = size(s) - 1; i >= 0; i--) {
//         if (check_same_triangle(s[i], new_tri3)) {
//           s.erase(s.begin() + i);
//           break;
//         } else if (i == 0) {
//           s.push_back(new_tri3); // if it is a new triangle, add it to the list
//           break;
//         }
//       }
//     }
//   }
// }

host_fn real distmin_origin(PolytopeFace face) {
  Vec3 a = face.p1;
  Vec3 b = face.p2;
  Vec3 c = face.p3;

  Vec3 ab = b - a;
  Vec3 ac = c - a;
  Vec3 bc = c - b;

  // Compute parametric position s for projection P’ of P on AB,
  // P’ = A + s*AB, s = snom/(snom+sdenom)
  float snom = dot(-a, ab), sdenom = dot(-b, a - b);

  // Compute parametric position t for projection P’ of P on AC,
  // P’ = A + t*AC, s = tnom/(tnom+tdenom)
  float tnom = dot(-a, ac), tdenom = dot(-c, a - c);

  if (snom <= 0.0f && tnom <= 0.0f)
    return norm(a); // Vertex region early out

  // Compute parametric position u for projection P’ of P on BC,
  // P’ = B + u*BC, u = unom/(unom+udenom)
  float unom = dot(-b, bc), udenom = dot(-c, b - c);
  if (sdenom <= 0.0f && unom <= 0.0f)
    return norm(b); // Vertex region early out

  if (tdenom <= 0.0f && udenom <= 0.0f)
    return norm(c); // Vertex region early out

  // P is outside (or on) AB if the triple scalar product [N PA PB] <= 0
  // Vec3 n = face.n;
  Vec3  n  = cross(b - a, c - a);
  float vc = dot(n, cross(a, b));

  // If P outside AB and within feature region of AB,
  // return projection of P onto AB
  if (vc <= 0.0f && snom >= 0.0f && sdenom >= 0.0f)
    return norm(a + snom / (snom + sdenom) * ab);

  // P is outside (or on) BC if the triple scalar product [N PB PC] <= 0
  float va = dot(n, cross(b, c));

  // If P outside BC and within feature region of BC,
  // return projection of P onto BC
  if (va <= 0.0f && unom >= 0.0f && udenom >= 0.0f)
    return norm(b + unom / (unom + udenom) * bc);
  // P is outside (or on) CA if the triple scalar product [N PC PA] <= 0
  float vb = dot(n, cross(c, a));

  // If P outside CA and within feature region of CA,
  // return projection of P onto CA
  if (vb <= 0.0f && tnom >= 0.0f && tdenom >= 0.0f)
    return norm(a + tnom / (tnom + tdenom) * ac);
  // P must project inside face region. Compute Q using barycentric coordinates
  float u = va / (va + vb + vc);
  float v = vb / (va + vb + vc);
  float w = 1.0f - u - v; // = vc / (va + vb + vc)
  return norm(u * a + v * b + w * c);
}

real solve_EPA(ComplexSimplex simplex, const std::vector<Vec3>& v1, const std::vector<Vec3>& v2) {
  ZoneScopedN("Messy EPA");
  // create polytope faces
  std::vector<PolytopeFace> faces;
  std::vector<int>          active_faces;
  std::vector<int>          deleted_faces;
  std::vector<TwoPts>       loose_edges;

  int n_faces = v1.size() * v2.size();
  faces.reserve(2 * n_faces);
  active_faces.reserve(n_faces);
  deleted_faces.reserve(8); // probably won't exceed this
  loose_edges.reserve(12);
  int id = 0;

  Vec3 ab = simplex.b - simplex.a;
  Vec3 ac = simplex.c - simplex.a;
  Vec3 ad = simplex.d - simplex.a;
  Vec3 bc = simplex.c - simplex.b;
  Vec3 bd = simplex.d - simplex.b;

  Vec3 n       = cross(ab, ac);
  real dot_a_n = dot(simplex.a, n);
  n            = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n            = (1 / norm(n)) * n;
  real dist    = distmin_origin({simplex.a, simplex.b, simplex.c, n});
  faces.push_back({simplex.a, simplex.b, simplex.c, n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(ab, ad);
  dot_a_n = dot(simplex.a, n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.a, simplex.b, simplex.d, n});
  faces.push_back({simplex.a, simplex.b, simplex.d, n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(ac, ad);
  dot_a_n = dot(simplex.a, n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.a, simplex.c, simplex.d, n});
  faces.push_back({simplex.a, simplex.c, simplex.d, n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(bc, bd);
  dot_a_n = dot(simplex.d, n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.b, simplex.c, simplex.d, n});
  faces.push_back({simplex.b, simplex.c, simplex.d, n, dist});
  active_faces.push_back(id);
  id++;

  real min_dist;
  int  idx;

  while (true) {
    // Find closest face
    min_dist = INF_REAL;
    for (auto active_id: active_faces) {
      if (faces[active_id].dist < min_dist) {
        min_dist = faces[active_id].dist;
        idx      = active_id;
      }
    }

    // obtain a new support point in the direction of the edge normal
    Vec3 support1 = GJK_get_support(v1, -faces[idx].n);
    Vec3 support2 = GJK_get_support(v2, faces[idx].n);
    Vec3 p        = support2 - support1;

    // If the vertex does not expand the polytope in the direction of the normal, the minimum distance
    // is with the closest face (unchanged). Compute and return.
    dist = dot(p - faces[idx].p1, faces[idx].n);
    if (dist < COLLISION_EPSILON) { // previously 1e-2
      break;
    }

    // Get all the faces which are "seen" by the new point, add them to deleted_faces and delete them.
    deleted_faces.clear();
    for (int i = 0; i < active_faces.size();) {
      if (dot(faces[active_faces[i]].n, p - faces[active_faces[i]].p1) > 0) {
        deleted_faces.push_back(active_faces[i]);
        active_faces[i] = active_faces.back();
        active_faces.pop_back();
      } else {
        i++;
      }
    }

    // Create new convex hull by determining which line segments are contained in the deleted faces.
    // These segments must be deleted and the ones which are only contained in one face will be added
    // as parts of the new faces (with the new point forming the remaining two line segments).
    TwoPts current_edge[3];
    loose_edges.clear();
    for (auto idx: deleted_faces) {
      current_edge[0] = {faces[idx].p1, faces[idx].p2};
      current_edge[1] = {faces[idx].p1, faces[idx].p3};
      current_edge[2] = {faces[idx].p2, faces[idx].p3};

      int  loose_edge_idx = 0;
      bool found_edge     = false;
      for (const auto& edge: current_edge) {
        found_edge = false;
        for (int k = 0; k < loose_edges.size(); k++) { // Is the current edge already in loose_edges ?
          if ((edge.p1 == loose_edges[k].p1 && edge.p2 == loose_edges[k].p2) ||
              (edge.p2 == loose_edges[k].p1 && edge.p1 == loose_edges[k].p2)) {
            // edge is already in list
            found_edge     = true;
            loose_edge_idx = k;
            break;
          }
        }

        if (found_edge == false)
          loose_edges.push_back(edge);
        else {
          loose_edges[loose_edge_idx] = std::move(loose_edges.back());
          loose_edges.pop_back();
        }
      }
    }

    // rebuild polytope with new faces
    // PolytopeFace new_face;
    Vec3 n;
    real dot_p1_n;
    for (int i = 0; i < loose_edges.size(); i++) {
      n        = cross(p - loose_edges[i].p1, p - loose_edges[i].p2);
      dot_p1_n = dot(loose_edges[i].p1, n);
      n        = dot_p1_n > 0 ? n : -n;
      n        = (1 / norm(n)) * n;

      dist = distmin_origin({loose_edges[i].p1, loose_edges[i].p2, p, n});
      faces.push_back({loose_edges[i].p1, loose_edges[i].p2, p, n, dist});
      active_faces.push_back(id);
      id++;
    }
  }
  return -std::abs(min_dist);
}

host_fn GJKResult solve_general_GJK(const std::vector<Vec3>& v1, const std::vector<Vec3>& v2) {
  // ZoneScopedN("Messy GJK");
  ComplexSimplex simplex;
  GJKResult      results;

  Vec3 direction = v1[0] - v2[0];
  simplex.a1     = GJK_get_support(v1, -direction);
  simplex.a2     = GJK_get_support(v2, direction);
  simplex.a      = simplex.a2 - simplex.a1; // todo: a1 and a2 necessary? Same for b1, b2, ...
  simplex.count  = 1;

  TwoPts solved;
  Vec3   p;
  int    old_simplex_count;

  while (true) {
    old_simplex_count = simplex.count;

    if (simplex.count == 1) {
      p         = simplex.a;
      direction = -simplex.a;
    } else {
      solved    = (simplex.count == 2) ? GJK_solve_simplex2(simplex) : ((simplex.count == 3) ? GJK_solve_simplex3(simplex) : GJK_solve_simplex4(simplex));
      p         = solved.p1; // current closest point on polytope
      direction = solved.p2;
    }

    if (dot(p, p) < COLLISION_EPSILON && old_simplex_count == 4) {
      results.intersection     = true;
      results.final_simplex    = simplex;
      results.A_closept        = {};
      results.B_closept        = {};
      results.minimal_distance = solve_EPA(simplex, v1, v2);
      Assert(isnan(results.minimal_distance) == false);
      break;
    }

    Vec3 support1 = GJK_get_support(v1, -direction);
    Vec3 support2 = GJK_get_support(v2, direction);

    Vec3 support         = support2 - support1;
    real goal_post       = dot(p, direction);       // current support point in search direction
    real current_highest = dot(support, direction); // new support point in search direction

    // if new support point isn't better, terminate
    if (current_highest - goal_post <= COLLISION_EPSILON) {
      if (simplex.count == 3 && norm(cross(simplex.b - simplex.a, simplex.c - simplex.a)) < 1e-9) { // if 3 pt simplex is line when it should be triangle
        Segment seg_test;
        // Calculate all points in minkowski difference
        // Make a list of all the points that have the same dot product as our point p (closest point)
        std::vector<Vec3> pts_list;
        for (auto& vertex1: v1) {
          for (auto& vertex2: v2) {
            Vec3 new_pt = vertex2 - vertex1;                                      // point on minkowski difference
            if (std::abs(goal_post - dot(new_pt, direction)) < COLLISION_EPSILON) // if new point is close to current support point in search direction
              pts_list.push_back(new_pt);
          }
        }

        // Within this list, find the two points which are the most extreme
        // Create a segment with these two points and test segment-point with the origin
        Vec3 new_dir = simplex.b - simplex.a;
        real max_dot = -INF_REAL;
        real min_dot = INF_REAL;
        for (auto& pts: pts_list) {
          real current_dot = dot(pts, new_dir);
          if (current_dot > max_dot) {
            max_dot     = current_dot;
            seg_test.p1 = pts;
          }
          if (current_dot < min_dot) { // used to be else if
            min_dot     = current_dot;
            seg_test.p2 = pts;
          }
        }
        simplex.a  = seg_test.p1;
        simplex.a1 = GJK_get_support(v1, -simplex.a);
        simplex.a2 = GJK_get_support(v2, simplex.a);

        simplex.b  = seg_test.p2;
        simplex.b1 = GJK_get_support(v1, -simplex.b);
        simplex.b2 = GJK_get_support(v2, simplex.b);

        simplex.count = 2;
      }

      TwoPts local = GJK_get_local_points(simplex, p);

      results.intersection     = false;
      results.A_closept        = local.p1;
      results.B_closept        = local.p2;
      results.final_simplex    = simplex;
      results.minimal_distance = norm(results.A_closept - results.B_closept);
      Assert(isnan(results.minimal_distance) == false);
      break;
    }

    Assert(simplex.count < 4); /*"You cannot have a simplex that's a tetrahedron at this point."*/

    simplex.d  = simplex.c;
    simplex.d1 = simplex.c1;
    simplex.d2 = simplex.c2;

    simplex.c  = simplex.b;
    simplex.c1 = simplex.b1;
    simplex.c2 = simplex.b2;

    simplex.b  = simplex.a;
    simplex.b1 = simplex.a1;
    simplex.b2 = simplex.a2;

    simplex.count += 1;

    simplex.a1 = support1;
    simplex.a2 = support2;
    simplex.a  = support;
  }

  return results;
}

host_fn real distance_GJK_simple(const Capsule& caps, const Box& box) {
  ZoneScopedN("Messy GJK");
  Vec3 size_x_org = {box.extents.x, 0, 0};
  Vec3 size_y_org = {0, box.extents.y, 0};
  Vec3 size_z_org = {0, 0, box.extents.z};
  Vec3 size_x     = box.rotation * size_x_org;
  Vec3 size_y     = box.rotation * size_y_org;
  Vec3 size_z     = box.rotation * size_z_org;

  std::vector<Vec3> v1(2);
  v1[0] = caps.p1;
  v1[1] = caps.p2;

  std::vector<Vec3> v2(8);
  v2[0] = box.center + size_x + size_y + size_z;
  v2[1] = box.center + size_x + size_y - size_z;
  v2[2] = box.center + size_x - size_y + size_z;
  v2[3] = box.center + size_x - size_y - size_z;
  v2[4] = box.center - size_x + size_y + size_z;
  v2[5] = box.center - size_x + size_y - size_z;
  v2[6] = box.center - size_x - size_y + size_z;
  v2[7] = box.center - size_x - size_y - size_z;

  auto result = solve_general_GJK(v1, v2);

  return (result.minimal_distance - caps.radius);
}

// ======================================
//            GJK algorithm
// ======================================
// Adapted from Ericson

host_fn void GJK_solve_simplex2_Ericson(Simplex* simplex) {
  Vec3 ab = (*simplex).b - (*simplex).a;
  real t  = dot(-(*simplex).a, ab) / dot(ab, ab);

  t = clamp(t, 0, 1);

  (*simplex).P = (*simplex).a + t * ab;

  if (t == 1) {
    Vec3 a_temp     = (*simplex).a;
    (*simplex).a    = (*simplex).b;
    (*simplex).b    = a_temp;
    (*simplex).size = 1;
  } else if (t == 0) {
    (*simplex).size = 1;
  }

  return;
}

host_fn void GJK_solve_simplex3_Ericson_old(Simplex* simplex) {
  // Check if P in vertex region outside A
  Vec3 a = (*simplex).a;
  Vec3 b = (*simplex).b;
  Vec3 c = (*simplex).c;

  Vec3 ac = c - a;
  // Vec3 ap = - a; // since p = 0 (originally = p - a)
  Vec3 ab = b - a;
  real d1 = dot(ab, -a);
  real d2 = dot(ac, -a);
  // if (d1 <= 0.0f && d2 <= 0.0f) {                  -> Impossible, or c would have been wrongly chosen
  //     (*simplex).size = 1;
  //     (*simplex).P = a; // barycentric coordinates (1,0,0)
  //     return;
  // }

  // Check if P in vertex region outside B
  // Vec3 bp = - b;  // since p = 0 (originally = p - b)
  real d3 = dot(ab, -b);
  real d4 = dot(ac, -b);
  // if (d3 >= 0.0f && d4 <= d3) {                    -> Impossible, or c would have been wrongly chosen
  //     (*simplex).size = 1;
  //     (*simplex).a = b;
  //     (*simplex).P = b;
  //     return;
  // }

  // // Check if P in edge region of AB, if so return projection of P onto AB
  real vc = d1 * d4 - d3 * d2;
  // if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {    -> Impossible, or c would have been wrongly chosen
  //     real v = d1 / (d1 - d3);
  //     (*simplex).size = 2;
  //     (*simplex).P = a + v * ab ; // barycentric coordinates (1-v,v,0)
  //     return;
  // }
  // Check if P in vertex region outside C
  // Vec3 cp = - c; // since p = 0 (originally = p - c)
  real d5 = dot(ab, -c);
  real d6 = dot(ac, -c);
  if (d6 >= 0.0f && d5 <= d6) {
    (*simplex).a    = c; // barycentric coordinates (0,0,1)
    (*simplex).c    = a;
    (*simplex).size = 1;
    (*simplex).P    = c;
    return;
  }

  // Check if P in edge region of BC, if so return projection of P onto BC
  real va = d3 * d6 - d5 * d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    real w          = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    (*simplex).size = 2;
    (*simplex).a    = b;
    (*simplex).b    = c;
    (*simplex).c    = a;
    (*simplex).P    = b + w * (c - b); // barycentric coordinates (0,1-w,w)
    return;
  }

  // Check if P in edge region of AC, if so return projection of P onto AC
  real vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    real w          = d2 / (d2 - d6);
    (*simplex).size = 2;
    (*simplex).b    = c;
    (*simplex).c    = a;
    (*simplex).P    = a + w * ac; // barycentric coordinates (1-w,0,w)
    return;
  }

  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  real denom   = 1.0f / (va + vb + vc);
  real v       = vb * denom;
  real w       = vc * denom;
  (*simplex).P = a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
  return;
}

host_fn void GJK_solve_simplex3_Ericson(Simplex* simplex) {
  Vec3 a = simplex->a;
  Vec3 b = simplex->b;
  Vec3 c = simplex->c;

  Vec3 ab = b - a;
  Vec3 ac = c - a;

  real d1 = dot(ab, -a);
  real d2 = dot(ac, -a);

  // Check if P in vertex region outside A
  if (d1 <= 0.0f && d2 <= 0.0f) {
    simplex->size = 1;
    // simplex->a    = a;
    simplex->P = a;
    return;
  }

  // Check if P in vertex region outside B
  real d3 = dot(ab, -b);
  real d4 = dot(ac, -b);
  if (d3 >= 0.0f && d4 <= d3) {
    simplex->size = 1;
    simplex->a    = b;
    simplex->P    = b;
    return;
  }

  // Check if P in edge region of AB (return projection of P onto AB)
  real vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
    real v        = d1 / (d1 - d3);
    simplex->size = 2;
    // simplex->a    = a;
    // simplex->b    = b;
    simplex->P = a + v * ab;
    return;
  }

  // Check if P in vertex region outside C
  real d5 = dot(ab, -c);
  real d6 = dot(ac, -c);
  if (d6 >= 0.0f && d5 <= d6) {
    simplex->size = 1;
    simplex->a    = c;
    simplex->P    = c;
    return;
  }

  // Check if P in edge region of BC (return projection of P onto BC)
  real va = d3 * d6 - d5 * d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    real w        = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    simplex->size = 2;
    simplex->a    = b;
    simplex->b    = c;
    simplex->P    = b + w * (c - b);
    return;
  }

  // Check if P in edge region of AC (return projection of P onto AC)
  real vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    real w        = d2 / (d2 - d6);
    simplex->size = 2;
    // simplex->a    = a;
    simplex->b = c;
    simplex->P = a + w * ac;
    return;
  }

  // P is inside face region
  real denom = 1.0f / (va + vb + vc);
  real v     = vb * denom;
  real w     = vc * denom;

  simplex->P = a + ab * v + ac * w;
  return;
}

host_fn int point_outside_plane(Vec3 p, Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
  // The last input (d) is the one that will be tested
  real signp = dot(p - a, cross(b - a, c - a)); // [AP AB AC]
  real signd = dot(d - a, cross(b - a, c - a)); // [AD AB AC]
  // Points on opposite sides if expression signs are opposite
  return signp * signd < 0.0f;
}

host_fn void GJK_solve_simplex4_Ericson(Simplex* simplex) {
  Vec3 a       = (*simplex).a;
  Vec3 b       = (*simplex).b;
  Vec3 c       = (*simplex).c;
  Vec3 d       = (*simplex).d;
  Vec3 p       = {0, 0, 0};
  (*simplex).P = p;

  Simplex simplex_temp = *simplex;

  // Start out assuming point inside all halfspaces, so closest to itself
  Vec3 closestPt  = p;
  real bestSqDist = INF_REAL;

  Vec3 a_temp;
  Vec3 b_temp;
  Vec3 c_temp;
  Vec3 d_temp;
  // If point outside face bdc then compute closest point on bcd
  if (point_outside_plane(p, b, c, d, a)) {
    simplex_temp.a    = b;
    simplex_temp.b    = c;
    simplex_temp.c    = d;
    simplex_temp.size = 3;
    GJK_solve_simplex3_Ericson(&simplex_temp);
    real sqDist = dot(simplex_temp.P, simplex_temp.P);
    if (sqDist < bestSqDist) {
      a_temp = b;
      b_temp = c;
      c_temp = d;
      d_temp = a;

      bestSqDist = sqDist;
      closestPt  = simplex_temp.P;
    }
    simplex_temp = *simplex; // resets to original simplex
  }
  // Repeat test for face acd
  if (point_outside_plane(p, a, c, d, b)) {
    // simplex_temp.a    = a;
    simplex_temp.b    = c;
    simplex_temp.c    = d;
    simplex_temp.size = 3;
    GJK_solve_simplex3_Ericson(&simplex_temp);
    real sqDist = dot(simplex_temp.P, simplex_temp.P);
    if (sqDist < bestSqDist) {
      a_temp = a;
      b_temp = c;
      c_temp = d;
      d_temp = b;

      bestSqDist = sqDist;
      closestPt  = simplex_temp.P;
    }
    // Simplex simplex_temp = *simplex;
    simplex_temp = *simplex;
  }
  // Repeat test for face adb
  if (point_outside_plane(p, a, b, d, c)) {
    // simplex_temp.a    = a;
    // simplex_temp.b    = b;
    simplex_temp.c    = d;
    simplex_temp.size = 3;
    GJK_solve_simplex3_Ericson(&simplex_temp);
    real sqDist = dot(simplex_temp.P, simplex_temp.P);
    if (sqDist < bestSqDist) {
      a_temp = a;
      b_temp = b;
      c_temp = d;
      d_temp = c;

      bestSqDist = sqDist;
      closestPt  = simplex_temp.P;
    }
    // Simplex simplex_temp = *simplex;
    simplex_temp = *simplex;
  }
  // Repeat test for face abc
  if (point_outside_plane(p, a, b, c, d)) {
    // simplex_temp.a    = a;
    // simplex_temp.b    = b;
    // simplex_temp.c    = c;
    simplex_temp.size = 3;
    GJK_solve_simplex3_Ericson(&simplex_temp);
    real sqDist = dot(simplex_temp.P, simplex_temp.P);
    // Update best closest point if (squared) distance is less than current best
    if (sqDist < bestSqDist) {
      a_temp = a;
      b_temp = b;
      c_temp = c;
      d_temp = d;

      bestSqDist = sqDist;
      closestPt  = simplex_temp.P;
    }
    // Simplex simplex_temp = *simplex;
    simplex_temp = *simplex;
  }

  if (dot(closestPt, closestPt) > COLLISION_EPSILON) {
    (*simplex).size = 3;
    (*simplex).a    = a_temp;
    (*simplex).b    = b_temp;
    (*simplex).c    = c_temp;
    (*simplex).d    = d_temp;
    (*simplex).P    = closestPt;

    // GJK_solve_simplex3_Ericson(simplex);
  }
  return;
}

// host_fn real distmin_origin(PolytopeFace face) {
//   Vec3 a = face.p1;
//   Vec3 b = face.p2;
//   Vec3 c = face.p3;
//
//   Vec3 ab = b - a;
//   Vec3 ac = c - a;
//   Vec3 bc = c - b;
//
//   // Compute parametric position s for projection P’ of P on AB,
//   // P’ = A + s*AB, s = snom/(snom+sdenom)
//   float snom = dot(-a, ab), sdenom = dot(-b, a - b);
//
//   // Compute parametric position t for projection P’ of P on AC,
//   // P’ = A + t*AC, s = tnom/(tnom+tdenom)
//   float tnom = dot(-a, ac), tdenom = dot(-c, a - c);
//
//   if (snom <= 0.0f && tnom <= 0.0f)
//     return norm(a); // Vertex region early out
//
//   // Compute parametric position u for projection P’ of P on BC,
//   // P’ = B + u*BC, u = unom/(unom+udenom)
//   float unom = dot(-b, bc), udenom = dot(-c, b - c);
//   if (sdenom <= 0.0f && unom <= 0.0f)
//     return norm(b); // Vertex region early out
//
//   if (tdenom <= 0.0f && udenom <= 0.0f)
//     return norm(c); // Vertex region early out
//
//   // P is outside (or on) AB if the triple scalar product [N PA PB] <= 0
//   // Vec3 n = face.n;
//   Vec3  n  = cross(b - a, c - a);
//   float vc = dot(n, cross(a, b));
//
//   // If P outside AB and within feature region of AB,
//   // return projection of P onto AB
//   if (vc <= 0.0f && snom >= 0.0f && sdenom >= 0.0f)
//     return norm(a + snom / (snom + sdenom) * ab);
//
//   // P is outside (or on) BC if the triple scalar product [N PB PC] <= 0
//   float va = dot(n, cross(b, c));
//
//   // If P outside BC and within feature region of BC,
//   // return projection of P onto BC
//   if (va <= 0.0f && unom >= 0.0f && udenom >= 0.0f)
//     return norm(b + unom / (unom + udenom) * bc);
//   // P is outside (or on) CA if the triple scalar product [N PC PA] <= 0
//   float vb = dot(n, cross(c, a));
//
//   // If P outside CA and within feature region of CA,
//   // return projection of P onto CA
//   if (vb <= 0.0f && tnom >= 0.0f && tdenom >= 0.0f)
//     return norm(a + tnom / (tnom + tdenom) * ac);
//   // P must project inside face region. Compute Q using barycentric coordinates
//   float u = va / (va + vb + vc);
//   float v = vb / (va + vb + vc);
//   float w = 1.0f - u - v; // = vc / (va + vb + vc)
//   return norm(u * a + v * b + w * c);
// }

real solve_EPA(Simplex simplex, const std::vector<Vec3>& v1, const std::vector<Vec3>& v2) {
  ZoneScopedN("Ericson EPA");
  // create a face vector that has three points and a normal
  std::vector<PolytopeFace> faces;
  std::vector<int>          active_faces;
  std::vector<int>          deleted_faces;
  std::vector<TwoPts>       loose_edges;

  int n_faces = v1.size() * v2.size(); // may have to be bigger for active faces
  faces.reserve(2 * n_faces);
  active_faces.reserve(n_faces);
  deleted_faces.reserve(8); // probably won't exceed this
  loose_edges.reserve(16);
  int id = 0;

  Vec3 ab = simplex.b - simplex.a;
  Vec3 ac = simplex.c - simplex.a;
  Vec3 ad = simplex.d - simplex.a;
  Vec3 bc = simplex.c - simplex.b;
  Vec3 bd = simplex.d - simplex.b;

  Vec3 n       = cross(ab, ac);
  real dot_a_n = dot(simplex.a, n);
  n            = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n            = (1 / norm(n)) * n;
  real dist    = distmin_origin({simplex.a, simplex.b, simplex.c, n});
  faces.push_back({simplex.a, simplex.b, simplex.c, n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(ab, ad);
  dot_a_n = dot(simplex.a, n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.a, simplex.b, simplex.d, n});
  faces.push_back({simplex.a, simplex.b, simplex.d, n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(ac, ad);
  dot_a_n = dot(simplex.a, n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.a, simplex.c, simplex.d, n});
  faces.push_back({simplex.a, simplex.c, simplex.d, n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(bc, bd);
  dot_a_n = dot(simplex.d, n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.b, simplex.c, simplex.d, n});
  faces.push_back({simplex.b, simplex.c, simplex.d, n, dist});
  active_faces.push_back(id);
  id++;

  real min_dist;
  int  idx;

  while (true) {
    // Find closest face
    min_dist = INF_REAL;
    for (auto active_id: active_faces) {
      if (faces[active_id].dist < min_dist) {
        min_dist = faces[active_id].dist;
        idx      = active_id;
      }
    }

    // obtain a new support point in the direction of the edge normal
    Vec3 support1 = GJK_get_support(v1, -faces[idx].n);
    Vec3 support2 = GJK_get_support(v2, faces[idx].n);
    Vec3 p        = support2 - support1;

    // If the vertex does not expand the polytope in the direction of the normal, the minimum distance
    // is with the closest face (unchanged). Compute and return.
    dist = dot(p - faces[idx].p1, faces[idx].n);
    if (dist < COLLISION_EPSILON) { // previously 1e-2
      break;
    }

    // Get all the faces which are "seen" by the new point, add them to deleted_faces and delete them.
    deleted_faces.clear();
    for (int i = 0; i < active_faces.size();) {
      if (dot(faces[active_faces[i]].n, p - faces[active_faces[i]].p1) > 0) {
        deleted_faces.push_back(active_faces[i]);
        active_faces[i] = active_faces.back();
        active_faces.pop_back();
      } else {
        i++;
      }
    }

    // Create new convex hull by determining which line segments are contained in the deleted faces.
    // These segments must be deleted and the ones which are only contained in one face will be added
    // as parts of the new faces (with the new point forming the remaining two line segments).
    TwoPts current_edge[3];
    loose_edges.clear();
    for (auto idx: deleted_faces) {
      current_edge[0] = {faces[idx].p1, faces[idx].p2};
      current_edge[1] = {faces[idx].p1, faces[idx].p3};
      current_edge[2] = {faces[idx].p2, faces[idx].p3};

      int  loose_edge_idx = 0;
      bool found_edge     = false;
      for (const auto& edge: current_edge) {
        found_edge = false;
        for (int k = 0; k < loose_edges.size(); k++) { // Is the current edge already in loose_edges ?
          if ((edge.p1 == loose_edges[k].p1 && edge.p2 == loose_edges[k].p2) ||
              (edge.p2 == loose_edges[k].p1 && edge.p1 == loose_edges[k].p2)) {
            // edge is already in list
            found_edge     = true;
            loose_edge_idx = k;
            break;
          }
        }

        if (found_edge == false)
          loose_edges.push_back(edge);
        else {
          loose_edges[loose_edge_idx] = std::move(loose_edges.back());
          loose_edges.pop_back();
        }
      }
    }

    // rebuild simplex with new faces
    // PolytopeFace new_face;
    Vec3 n;
    real dot_p1_n;
    for (int i = 0; i < loose_edges.size(); i++) {
      n        = cross(p - loose_edges[i].p1, p - loose_edges[i].p2);
      dot_p1_n = dot(loose_edges[i].p1, n);
      n        = dot_p1_n > 0 ? n : -n;
      n        = (1 / norm(n)) * n;

      dist = distmin_origin({loose_edges[i].p1, loose_edges[i].p2, p, n});
      faces.push_back({loose_edges[i].p1, loose_edges[i].p2, p, n, dist});
      active_faces.push_back(id);
      id++;
    }
  }
  return -std::abs(min_dist);
}

// host_fn real solve_EPA(Simplex simplex, std::vector<Vec3> v1, std::vector<Vec3> v2) {
//   ZoneScopedN("Ericson EPA");
//   // create a face vector that has three points and a normal
//   std::vector<PolytopeFace> faces;
//   faces.reserve(v1.size() * v2.size());
//
//   Vec3 ab = simplex.b - simplex.a;
//   Vec3 ac = simplex.c - simplex.a;
//   Vec3 ad = simplex.d - simplex.a;
//   Vec3 bc = simplex.c - simplex.b;
//   Vec3 bd = simplex.d - simplex.b;
//
//   Vec3 n1 = cross(ab, ac);
//   Vec3 n2 = cross(ab, ad);
//   Vec3 n3 = cross(ac, ad);
//   Vec3 n4 = cross(bc, bd);
//
//   real dot_a_n1 = dot(simplex.a, n1);
//   real dot_a_n2 = dot(simplex.a, n2);
//   real dot_a_n3 = dot(simplex.a, n3);
//   real dot_d_n4 = dot(simplex.d, n4);
//
//   n1 = dot_a_n1 > 0 ? n1 : (dot_a_n1 < 0 ? -n1 : (dot(n1, simplex.d) >= 0 ? -n1 : n1));
//   n2 = dot_a_n2 > 0 ? n2 : (dot_a_n2 < 0 ? -n2 : (dot(n2, simplex.d) >= 0 ? -n2 : n2));
//   n3 = dot_a_n3 > 0 ? n3 : (dot_a_n3 < 0 ? -n3 : (dot(n3, simplex.d) >= 0 ? -n3 : n3));
//   n4 = dot_d_n4 > 0 ? n4 : (dot_d_n4 < 0 ? -n4 : (dot(n4, simplex.a) >= 0 ? -n4 : n4));
//
//   n1 = (1 / norm(n1)) * n1;
//   n2 = (1 / norm(n2)) * n2;
//   n3 = (1 / norm(n3)) * n3;
//   n4 = (1 / norm(n4)) * n4;
//
//   faces.push_back({simplex.a, simplex.b, simplex.c, n1});
//   faces.push_back({simplex.a, simplex.b, simplex.d, n2});
//   faces.push_back({simplex.a, simplex.c, simplex.d, n3});
//   faces.push_back({simplex.b, simplex.c, simplex.d, n4});
//
//   real min_dist;
//   int  idx;
//   real dist;
//
//   while (true) {
//     // Find closest face
//     min_dist = INF_REAL;
//     for (int i = 0; i < size(faces); i++) {
//       real current_dist = distmin_origin(faces[i]);
//       if (current_dist < min_dist) {
//         min_dist = current_dist;
//         idx      = i;
//       }
//     }
//
//     // obtain a new support point in the direction of the edge normal
//     Vec3 support1 = GJK_get_support(v1, -faces[idx].n);
//     Vec3 support2 = GJK_get_support(v2, faces[idx].n);
//     Vec3 p        = support2 - support1;
//
//     // If the vertex does not expand the polytope in the direction of the normal, the minimum distance
//     // is with the closest face (unchanged). Compute and return.
//     dist = dot(p - faces[idx].p1, faces[idx].n);
//     if (dist < COLLISION_EPSILON) { // previously 1e-2
//       break;
//     }
//
//     // Get all the faces which are "seen" by the new point, add them to deleted_faces and delete them.
//     std::vector<PolytopeFace> deleted_faces;
//     deleted_faces.reserve(faces.size());
//     for (int i = 0; i < size(faces); i++) {
//       if (dot(faces[i].n, p - faces[i].p1) > 0) {
//         deleted_faces.push_back(faces[i]);
//         faces.erase(faces.begin() + i);
//         i--; // bug fix maybe
//       }
//     }
//
//     // Create new convex hull by determining which line segments are contained in the deleted faces.
//     // These segments must be deleted and the ones which are only contained in one face will be added
//     // as parts of the new faces (with the new point forming the remaining two line segments).
//     TwoPts              current_edge[3];
//     std::vector<TwoPts> loose_edges;
//     for (int i = 0; i < size(deleted_faces); i++) { // for all deleted faces found
//       current_edge[0] = {deleted_faces[i].p1, deleted_faces[i].p2};
//       current_edge[1] = {deleted_faces[i].p1, deleted_faces[i].p3};
//       current_edge[2] = {deleted_faces[i].p2, deleted_faces[i].p3};
//
//       int  loose_edge_idx = 0;
//       bool found_edge     = false;
//       for (int j = 0; j < 3; j++) {                    // for all three edges of each face
//         found_edge = false;
//         for (int k = 0; k < loose_edges.size(); k++) { // Is the current edge already in loose_edges ?
//           if ((current_edge[j].p1 == loose_edges[k].p1 && current_edge[j].p2 == loose_edges[k].p2) ||
//               (current_edge[j].p2 == loose_edges[k].p1 && current_edge[j].p1 == loose_edges[k].p2)) {
//             // edge is already in list
//             found_edge     = true;
//             loose_edge_idx = k;
//             break;
//           }
//         }
//
//         if (found_edge == false)
//           loose_edges.push_back(current_edge[j]);
//         else {
//           loose_edges.erase(loose_edges.begin() + loose_edge_idx);
//         }
//       }
//     }
//
//     // rebuild simplex with new faces
//     PolytopeFace new_face;
//     Vec3     n;
//     real     dot_p1_n;
//     for (int i = 0; i < loose_edges.size(); i++) {
//       n        = cross(p - loose_edges[i].p1, p - loose_edges[i].p2);
//       dot_p1_n = dot(loose_edges[i].p1, n);
//       n        = dot_p1_n > 0 ? n : -n;
//       n        = (1 / norm(n)) * n;
//
//       new_face = {loose_edges[i].p1, loose_edges[i].p2, p, n};
//       faces.push_back(new_face);
//     }
//   }
//   return -std::abs(min_dist);
// }

host_fn real general_GJK(const std::vector<Vec3>& set1, const std::vector<Vec3>& set2) {
  // ZoneScopedN("Ericson GJK");
  // Implemented from the basic algorithm described in Collision Detection manual by Ericson.

  // 1. Initializing simplex to a point from a random direction
  Vec3 direction = set1[0] - set2[0];

  Vec3 a1 = GJK_get_support(set1, -direction);
  Vec3 a2 = GJK_get_support(set2, direction);
  Vec3 V  = a2 - a1;

  Simplex simplex;
  simplex.a    = V;
  simplex.size = 1;
  real dot_P;

  while (true) {
    // 2. Computing the point P of minimum norm in CH(simplex)
    switch (simplex.size) {
      case 1:
        simplex.P = simplex.a;
        break;
      case 2:
        GJK_solve_simplex2_Ericson(&simplex);
        break;
      case 3:
        GJK_solve_simplex3_Ericson(&simplex);
        break;
      case 4:
        GJK_solve_simplex4_Ericson(&simplex);
        break;
      default:
        Assert(false);
    }

    dot_P = dot(simplex.P, simplex.P);

    // 3. If P is the origin itself, stop and return A and B as intersecting.
    if (dot_P < COLLISION_EPSILON && simplex.size == 4) {
      real dist = solve_EPA(simplex, set1, set2);
      return dist;
    }

    // 4. Remove any points from simplex not determining the subsimplex of simplex in which P lies.
    // (Done automatically in GJK_solve_simplex functions)

    // 5. Find the next supporting point in direction -P
    a1 = GJK_get_support(set1, simplex.P);
    a2 = GJK_get_support(set2, -simplex.P);
    V  = a2 - a1;

    // 6. If V is no more exremal in direction -P than P itself, stop and return A and B as not intersecting.
    // The length of the vector from the origin to P is the separation distance of A and B.
    real ans1 = dot(V, -simplex.P) / dot_P;

    if (ans1 + 1 <= COLLISION_EPSILON) // No more progress is being made
      break;

    // 7. Add V to simplex and go to 2.
    // todo : optimize this part of code using std::vector instead of Vec3 for a, b, c, d
    Assert(simplex.size <= 3);

    if (simplex.size == 1)
      simplex.b = V;
    if (simplex.size == 2)
      simplex.c = V;
    if (simplex.size == 3)
      simplex.d = V;
    simplex.size += 1;
  }

  return std::sqrt(dot_P);
}

host_fn real distance_GJK(const Capsule& caps, const Box& box) {
  ZoneScopedN("Ericson GJK");
  Assert(is_close(norm(box.rotation.col_copy(0)), 1.0, 1e-6));
  Assert(is_close(norm(box.rotation.col_copy(1)), 1.0, 1e-6));
  Assert(is_close(norm(box.rotation.col_copy(2)), 1.0, 1e-6));
  Assert(is_close(norm(Vec3(box.rotation(0, 0), box.rotation(0, 1), box.rotation(0, 2))), 1.0, 1e-6));
  Assert(is_close(norm(Vec3(box.rotation(1, 0), box.rotation(1, 1), box.rotation(1, 2))), 1.0, 1e-6));
  Assert(is_close(norm(Vec3(box.rotation(2, 0), box.rotation(2, 1), box.rotation(2, 2))), 1.0, 1e-6));

  // Initialization of the eight OBB points
  Vec3 size_x_org = {box.extents.x, 0, 0};
  Vec3 size_y_org = {0, box.extents.y, 0};
  Vec3 size_z_org = {0, 0, box.extents.z};
  Vec3 size_x     = box.rotation * size_x_org;
  Vec3 size_y     = box.rotation * size_y_org;
  Vec3 size_z     = box.rotation * size_z_org;

  std::vector<Vec3> v1(2);
  v1[0] = caps.p1;
  v1[1] = caps.p2;

  std::vector<Vec3> v2(8);
  v2[0] = box.center + size_x + size_y + size_z;
  v2[1] = box.center + size_x + size_y - size_z;
  v2[2] = box.center + size_x - size_y + size_z;
  v2[3] = box.center + size_x - size_y - size_z;
  v2[4] = box.center - size_x + size_y + size_z;
  v2[5] = box.center - size_x + size_y - size_z;
  v2[6] = box.center - size_x - size_y + size_z;
  v2[7] = box.center - size_x - size_y - size_z;

  real dist = general_GJK(v1, v2) - caps.radius;
  return dist;
}

} // namespace blast
