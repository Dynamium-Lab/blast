#pragma once

#include <blast>
#include "tracy/Tracy.hpp"

namespace blast {

#define COLLISION_EPSILON 1e-9
#define TRACY_ACTIVE 2 // 0 nothing, 1 big functions, 2 medium functions, 3 all functions
#define GJK_ACTIVE 1   // 0 nothing, 1 ericson, 2 messy, 3 both
int MAX_ITERATIONS = 64;

// structs
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

struct Simplex {
  Vec3 P;
  Vec3 pts[4]; // (a, b, c, d)
  int  size;
};

struct PolytopeFace {
  Vec3 p1;
  Vec3 p2;
  Vec3 p3;
  Vec3 n;
  real dist = 0.0;
};

// ======================================
//            GJK algorithm
// ======================================
// This version of the GJK algorithm is adapted from : https://www.youtube.com/watch?v=DGVZYdlw_uo

inline host_fn real GJK_triangle_area_2d(real x1, real y1, real x2, real y2, real x3, real y3) {
#if TRACY_ACTIVE >= 3 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

host_fn Vec3 GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 c, Vec3 p) {
#if TRACY_ACTIVE >= 3 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
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

inline host_fn real GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 p) {
#if TRACY_ACTIVE >= 3 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 ab = b - a;
  Vec3 ap = p - a;

  real t = dot(ap, ab) / dot(ab, ab);
  return t;
}

inline host_fn Vec3 GJK_convert_cartesian(Vec3 a, Vec3 b, real barycentric) {
#if TRACY_ACTIVE >= 3 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 result;
  result.x = a.x + barycentric * (b.x - a.x);
  result.y = a.y + barycentric * (b.y - a.y);
  result.z = a.z + barycentric * (b.z - a.z);
  return result;
}

inline host_fn Vec3 GJK_convert_cartesian(Vec3 a, Vec3 b, Vec3 c, Vec3 barycentric) {
#if TRACY_ACTIVE >= 3 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif;
  Vec3 result;
  result.x = a.x * barycentric.x + b.x * barycentric.y + c.x * barycentric.z;
  result.y = a.y * barycentric.x + b.y * barycentric.y + c.y * barycentric.z;
  result.z = a.z * barycentric.x + b.z * barycentric.y + c.z * barycentric.z;
  return result;
}

host_fn TwoPts GJK_get_local_points(const ComplexSimplex& simplex, Vec3 p) {
#if TRACY_ACTIVE >= 3 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
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
      Assert(false);
  }
  Assert(isnan(loc.p1.x) == false);
  return loc;
}

inline host_fn Vec3 GJK_get_support(Vec3* vertices, size_t count, Vec3 direction) {
#if TRACY_ACTIVE >= 2
  ZoneScoped;
#endif
  real   largest_dot = dot(vertices[0], direction);
  size_t largest_idx = 0;
  for (size_t i = 1; i < count; i++) {
    const real current_dot = dot(vertices[i], direction);
    if (current_dot > largest_dot) {
      largest_idx = i;
      largest_dot = current_dot;
    }
  }
  Assert(isinf(largest_dot) == false);
  return vertices[largest_idx];
}

host_fn TwoPts GJK_solve_simplex2(ComplexSimplex& simplex) { // finds next search direction (ao) and closest point on line (ab) for 2pt simplex
#if TRACY_ACTIVE >= 2 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
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
#if TRACY_ACTIVE >= 2 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
  Vec3 ac  = simplex.c - simplex.a;
  Vec3 ao  = -simplex.a;

  if (dot(cross(abc, ac), ao) > 0) {
    // the origin is in the direction of the triangle normal
    if (dot(ac, ao) > 0) {
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
      if (dot(ab, ao) > 0) {
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
    if (dot(cross(ab, abc), ao) > 0) {
      if (dot(ab, ao) > 0) {
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
      if (dot(abc, ao) > 0) {
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
#if TRACY_ACTIVE >= 2 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
  Vec3 acd = cross(simplex.c - simplex.a, simplex.d - simplex.a);
  Vec3 adb = cross(simplex.d - simplex.a, simplex.b - simplex.a);

  Vec3 ao = -simplex.a;

  bool abc_dir = dot(abc, ao) > 0;
  bool acd_dir = dot(acd, ao) > 0;
  bool adb_dir = dot(adb, ao) > 0;

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

host_fn real distmin_origin(const PolytopeFace& face) {
#if TRACY_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 a = face.p1;
  Vec3 b = face.p2;
  Vec3 c = face.p3;

  Vec3 ab = b - a;
  Vec3 ac = c - a;
  Vec3 bc = c - b;

  // Check origin outisde vertex A
  float snom = dot(-a, ab);
  float tnom = dot(-a, ac);
  if (snom <= 0.0f && tnom <= 0.0f)
    return norm(a);

  // Check origin outisde vertex B
  float sdenom = dot(b, ab);
  float unom   = dot(-b, bc);
  if (sdenom <= 0.0f && unom <= 0.0f)
    return norm(b);

  // Check origin outisde vertex C
  float tdenom = dot(c, ac);
  float udenom = dot(c, bc);
  if (tdenom <= 0.0f && udenom <= 0.0f)
    return norm(c);

  Vec3 n = cross(b - a, c - a);

  // Check origin outside edge AB
  float vc = dot(n, cross(a, b));
  if (vc <= 0.0f && snom >= 0.0f && sdenom >= 0.0f)
    return norm(a + snom / (snom + sdenom) * ab);

  // Check origin outside edge BC
  float va = dot(n, cross(b, c));
  if (va <= 0.0f && unom >= 0.0f && udenom >= 0.0f)
    return norm(b + unom / (unom + udenom) * bc);

  // Check origin outside edge AC
  float vb = dot(n, cross(c, a));
  if (vb <= 0.0f && tnom >= 0.0f && tdenom >= 0.0f)
    return norm(a + tnom / (tnom + tdenom) * ac);

  // Origin must be on ABC
  float denom = va + vb + vc;
  float u     = va / denom;
  float v     = vb / denom;
  float w     = 1.0f - u - v; // = vc / (va + vb + vc)
  return norm(u * a + v * b + w * c);
}

real solve_EPA(ComplexSimplex simplex, Vec3* v1, size_t count1, Vec3* v2, size_t count2) {
#if TRACY_ACTIVE >= 1 && GJK_ACTIVE >= 2
  ZoneScopedN("Messy EPA");
#endif
  // create polytope faces
  std::vector<PolytopeFace> faces;
  std::vector<int>          active_faces;
  std::vector<int>          deleted_faces;
  std::vector<TwoPts>       loose_edges;

  int n_faces = count1 * count2;
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
  // n            = (1 / norm(n)) * n;
  real dist = distmin_origin({simplex.a, simplex.b, simplex.c, n});
  faces.push_back({simplex.a, simplex.b, simplex.c, n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(ab, ad);
  dot_a_n = dot(simplex.a, n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
  // n       = (1 / norm(n)) * n;
  dist = distmin_origin({simplex.a, simplex.b, simplex.d, n});
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

  // while (true) {
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    if (i == MAX_ITERATIONS - 1) {
      std::cout << "Max iterations reached - messy EPA" << std::endl;
    }
    // Find closest face
    min_dist = INF_REAL;
    for (auto active_id: active_faces) {
      if (faces[active_id].dist < min_dist) {
        min_dist = faces[active_id].dist;
        idx      = active_id;
      }
    }

    // obtain a new support point in the direction of the edge normal
    Vec3 support1 = GJK_get_support(v1, count1, -faces[idx].n);
    Vec3 support2 = GJK_get_support(v2, count2, faces[idx].n);
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

host_fn real solve_general_GJK(Vec3* v1, size_t count1, Vec3* v2, size_t count2) {
  // #if TRACY_ACTIVE >= 1
  //   ZoneScoped;
  //   // ZoneScopedN("Messy GJK");
  // #endif
  ComplexSimplex simplex;
  real           dist_min;

  Vec3 direction = v1[0] - v2[0];
  simplex.a1     = GJK_get_support(v1, count1, -direction);
  simplex.a2     = GJK_get_support(v2, count2, direction);
  simplex.a      = simplex.a2 - simplex.a1;
  simplex.count  = 1;

  TwoPts solved;
  Vec3   p;
  int    old_simplex_count;

  // while (true) {
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    if (i == MAX_ITERATIONS - 1) {
      std::cout << "Max iterations reached - messy GJK" << std::endl;
    }
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
      dist_min = solve_EPA(simplex, v1, count1, v2, count2);
      Assert(isnan(dist_min) == false);
      break;
    }

    Vec3 support1 = GJK_get_support(v1, count1, -direction);
    Vec3 support2 = GJK_get_support(v2, count2, direction);

    Vec3 support     = support2 - support1;
    real old_support = dot(p, direction); // (in direction)
    real new_support = dot(support, direction);

    // if no improvement, terminate
    if (new_support - old_support <= COLLISION_EPSILON) {
      if (simplex.count == 3 && norm(cross(simplex.b - simplex.a, simplex.c - simplex.a)) < 1e-9) { // if 3 pt simplex is line when it should be triangle
        Segment seg_test;
        // Calculate all points in minkowski difference and find the two points which are the most extreme
        // Create a segment with these two points and test segment-point with the origin
        Vec3 new_dir = simplex.b - simplex.a;
        real max_dot = -INF_REAL;
        real min_dot = INF_REAL;
        int  min_id1 = -1;
        int  min_id2 = -1;
        int  max_id1 = -1;
        int  max_id2 = -1;
        for (int i = 0; i < count1; i++) {
          for (int j = 0; j < count2; j++) {
            Vec3 new_pt = v2[j] - v1[i];
            if (std::abs(old_support - dot(new_pt, direction)) < COLLISION_EPSILON) {
              real current_dot = dot(new_pt, new_dir);
              if (current_dot > max_dot) {
                max_dot     = current_dot;
                seg_test.p1 = new_pt;
                max_id1     = i;
                max_id2     = j;
              }
              if (current_dot < min_dot) {
                min_dot     = current_dot;
                seg_test.p2 = new_pt;
                min_id1     = i;
                max_id2     = j;
              }
            }
          }
        }

        simplex.a  = seg_test.p1;
        simplex.a1 = v1[max_id1];
        simplex.a2 = v2[max_id2];

        simplex.b  = seg_test.p2;
        simplex.b1 = v1[min_id1];
        simplex.b2 = v2[min_id2];

        simplex.count = 2;
      }

      TwoPts local = GJK_get_local_points(simplex, p);

      dist_min = norm(local.p1 - local.p2);
      Assert(isnan(dist_min) == false);
      break;
    }

    Assert(simplex.count < 4); // You cannot have a simplex that's a tetrahedron at this point.

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

  return dist_min;
}

host_fn real distance_GJK_simple(const Capsule& caps, const Box& box) {
#if TRACY_ACTIVE >= 1 && GJK_ACTIVE >= 2
  ZoneScopedN("Messy GJK");
#endif
  Vec3 size_x_org = {box.extents.x, 0, 0};
  Vec3 size_y_org = {0, box.extents.y, 0};
  Vec3 size_z_org = {0, 0, box.extents.z};
  Vec3 size_x     = box.rotation * size_x_org;
  Vec3 size_y     = box.rotation * size_y_org;
  Vec3 size_z     = box.rotation * size_z_org;

  std::array<Vec3, 2> v1;
  v1[0] = caps.p1;
  v1[1] = caps.p2;

  std::array<Vec3, 8> v2;
  v2[0] = box.center + size_x + size_y + size_z;
  v2[1] = box.center + size_x + size_y - size_z;
  v2[2] = box.center + size_x - size_y + size_z;
  v2[3] = box.center + size_x - size_y - size_z;
  v2[4] = box.center - size_x + size_y + size_z;
  v2[5] = box.center - size_x + size_y - size_z;
  v2[6] = box.center - size_x - size_y + size_z;
  v2[7] = box.center - size_x - size_y - size_z;

  auto result = solve_general_GJK(&v1[0], 2, &v2[0], 8);

  return (result - caps.radius);
}

// ======================================
//            GJK algorithm
// ======================================
// Adapted from Ericson

host_fn void GJK_solve_simplex2_Ericson(Simplex* simplex) {
#if TRACY_ACTIVE >= 2 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScoped;
#endif
  Vec3 ab = (*simplex).pts[1] - (*simplex).pts[0];
  real t  = dot(-(*simplex).pts[0], ab) / dot(ab, ab);

  t = clamp(t, 0, 1);

  (*simplex).P = (*simplex).pts[0] + t * ab;

  if (t == 1) {
    Vec3 a_temp       = (*simplex).pts[0];
    (*simplex).pts[0] = (*simplex).pts[1];
    (*simplex).pts[1] = a_temp;
    (*simplex).size   = 1;
  } else if (t == 0) {
    (*simplex).size = 1;
  }

  return;
}

host_fn void GJK_solve_simplex3_Ericson(Simplex* simplex) {
#if TRACY_ACTIVE >= 2 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScoped;
#endif
  Vec3 a = simplex->pts[0];
  Vec3 b = simplex->pts[1];
  Vec3 c = simplex->pts[2];

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
    simplex->size   = 1;
    simplex->pts[0] = b;
    simplex->P      = b;
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
    simplex->size   = 1;
    simplex->pts[0] = c;
    simplex->P      = c;
    return;
  }

  // Check if P in edge region of BC (return projection of P onto BC)
  real va = d3 * d6 - d5 * d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    real w          = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    simplex->size   = 2;
    simplex->pts[0] = b;
    simplex->pts[1] = c;
    simplex->P      = b + w * (c - b);
    return;
  }

  // Check if P in edge region of AC (return projection of P onto AC)
  real vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    real w        = d2 / (d2 - d6);
    simplex->size = 2;
    // simplex->pts[0]    = a;
    simplex->pts[1] = c;
    simplex->P      = a + w * ac;
    return;
  }

  // P is inside face region
  real denom = 1.0f / (va + vb + vc);
  real v     = vb * denom;
  real w     = vc * denom;

  simplex->P = a + ab * v + ac * w;
  return;
}

inline host_fn int point_outside_plane(Vec3 p, Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
#if TRACY_ACTIVE >= 3 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScoped;
#endif
  // The last input (d) is the one that will be tested
  real signp = dot(p - a, cross(b - a, c - a)); // [AP AB AC]
  real signd = dot(d - a, cross(b - a, c - a)); // [AD AB AC]
  // Points on opposite sides if expression signs are opposite
  return signp * signd < 0.0f;
}

// host_fn void GJK_solve_simplex4_Ericson(Simplex* simplex) {
// #if TRACY_ACTIVE >= 2 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
//   ZoneScoped;
// #endif

//   // d is the most recently added point in standard GJK convention,
//   // but here we work relative to the origin O = (0, 0, 0).
//   // Vector AO, BO, CO, DO are simply -a, -b, -c, -d.
//   const Vec3 a = simplex->pts[0];
//   const Vec3 b = simplex->pts[1];
//   const Vec3 c = simplex->pts[2];
//   const Vec3 d = simplex->pts[3]; // Last added point (most reliable search direction)

//   const Vec3 ao  = {-a.x, -a.y, -a.z};
//   const Vec3 bo  = {-b.x, -b.y, -b.z};
//   const Vec3 co  = {-c.x, -c.y, -c.z};
//   const Vec3 do_ = {-d.x, -d.y, -d.z};

//   // Edge vectors pointing away from d
//   const Vec3 da = {a.x - d.x, a.y - d.y, a.z - d.z};
//   const Vec3 db = {b.x - d.x, b.y - d.y, b.z - d.z};
//   const Vec3 dc = {c.x - d.x, c.y - d.y, c.z - d.z};

//   // Triangle face normals (unnormalized) pointing outwards
//   const Vec3 abc_norm = cross({b.x - a.x, b.y - a.y, b.z - a.z}, {c.x - a.x, c.y - a.y, c.z - a.x});
//   const Vec3 dbc_norm = cross(db, dc);
//   const Vec3 dca_norm = cross(dc, da);
//   const Vec3 dab_norm = cross(da, db);

//   // 1. Check Face DBC
//   if (dot(dbc_norm, do_) > 0.0f) {
//     if (dot(cross(db, dbc_norm), do_) > 0.0f) {
//       // Region DB
//       simplex->size   = 2;
//       simplex->pts[0] = d;
//       simplex->pts[1] = b;
//       // Project O onto segment DB
//       real t     = dot(do_, db) / dot(db, db);
//       simplex->P = {d.x + t * db.x, d.y + t * db.y, d.z + t * db.z};
//       return;
//     }
//     if (dot(cross(dbc_norm, dc), do_) > 0.0f) {
//       // Region DC
//       simplex->size   = 2;
//       simplex->pts[0] = d;
//       simplex->pts[1] = c;
//       real t          = dot(do_, dc) / dot(dc, dc);
//       simplex->P      = {d.x + t * dc.x, d.y + t * dc.y, d.z + t * dc.z};
//       return;
//     }
//     // Region DBC (Triangle Face)
//     simplex->size   = 3;
//     simplex->pts[0] = d;
//     simplex->pts[1] = b;
//     simplex->pts[2] = c;
//     real inv_denom  = 1.0f / dot(dbc_norm, dbc_norm);
//     real dist       = dot(dbc_norm, do_) * inv_denom;
//     simplex->P      = {dist * dbc_norm.x, dist * dbc_norm.y, dist * dbc_norm.z};
//     return;
//   }

//   // 2. Check Face DCA
//   if (dot(dca_norm, do_) > 0.0f) {
//     if (dot(cross(dca_norm, da), do_) > 0.0f) {
//       // Region DA
//       simplex->size   = 2;
//       simplex->pts[0] = d;
//       simplex->pts[1] = a;
//       real t          = dot(do_, da) / dot(da, da);
//       simplex->P      = {d.x + t * da.x, d.y + t * da.y, d.z + t * da.z};
//       return;
//     }
//     // Region DCA (Triangle Face)
//     simplex->size   = 3;
//     simplex->pts[0] = d;
//     simplex->pts[1] = c;
//     simplex->pts[2] = a;
//     real inv_denom  = 1.0f / dot(dca_norm, dca_norm);
//     real dist       = dot(dca_norm, do_) * inv_denom;
//     simplex->P      = {dist * dca_norm.x, dist * dca_norm.y, dist * dca_norm.z};
//     return;
//   }

//   // 3. Check Face DAB
//   if (dot(dab_norm, do_) > 0.0f) {
//     // Region DAB (Triangle Face)
//     simplex->size   = 3;
//     simplex->pts[0] = d;
//     simplex->pts[1] = a;
//     simplex->pts[2] = b;
//     real inv_denom  = 1.0f / dot(dab_norm, dab_norm);
//     real dist       = dot(dab_norm, do_) * inv_denom;
//     simplex->P      = {dist * dab_norm.x, dist * dab_norm.y, dist * dab_norm.z};
//     return;
//   }

//   // 4. Origin is inside the tetrahedron (Collision detected)
//   simplex->size = 4;
//   simplex->P    = {0.0f, 0.0f, 0.0f};
// }

host_fn void GJK_solve_simplex4_Ericson(Simplex* simplex) {
#if TRACY_ACTIVE >= 2 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScoped;
#endif
  Vec3 a       = (*simplex).pts[0];
  Vec3 b       = (*simplex).pts[1];
  Vec3 c       = (*simplex).pts[2];
  Vec3 d       = (*simplex).pts[3];
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
    simplex_temp.pts[0] = b;
    simplex_temp.pts[1] = c;
    simplex_temp.pts[2] = d;
    simplex_temp.size   = 3;
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
    // simplex_temp.pts[0]    = a;
    simplex_temp.pts[1] = c;
    simplex_temp.pts[2] = d;
    simplex_temp.size   = 3;
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
    simplex_temp = *simplex;
  }
  // Repeat test for face adb
  if (point_outside_plane(p, a, b, d, c)) {
    // simplex_temp.pts[0]    = a;
    // simplex_temp.pts[1]    = b;
    simplex_temp.pts[2] = d;
    simplex_temp.size   = 3;
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
    // simplex_temp = *simplex;
  }
  // // Repeat test for face abc // should never be this
  // if (point_outside_plane(p, a, b, c, d)) {
  //   // simplex_temp.pts[0]    = a;
  //   // simplex_temp.pts[1]    = b;
  //   // simplex_temp.pts[2]    = c;
  //   simplex_temp.size = 3;
  //   GJK_solve_simplex3_Ericson(&simplex_temp);
  //   real sqDist = dot(simplex_temp.P, simplex_temp.P);
  //   // Update best closest point if (squared) distance is less than current best
  //   if (sqDist < bestSqDist) {
  //     a_temp = a;
  //     b_temp = b;
  //     c_temp = c;
  //     d_temp = d;

  //     bestSqDist = sqDist;
  //     closestPt  = simplex_temp.P;
  //   }
  //   // simplex_temp = *simplex;
  // }

  if (dot(closestPt, closestPt) > COLLISION_EPSILON) {
    (*simplex).size   = 3;
    (*simplex).pts[0] = a_temp;
    (*simplex).pts[1] = b_temp;
    (*simplex).pts[2] = c_temp;
    (*simplex).pts[3] = d_temp;
    (*simplex).P      = closestPt;
  }
  return;
}

real solve_EPA(Simplex simplex, Vec3* v1, size_t count1, Vec3* v2, size_t count2) {
#if TRACY_ACTIVE >= 1 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScopedN("Ericson EPA");
#endif
  // create a face vector that has three points and a normal
  std::vector<PolytopeFace> faces;
  std::vector<int>          active_faces;
  std::vector<int>          deleted_faces;
  std::vector<TwoPts>       loose_edges;

  int n_faces = count1 * count2; // may have to be bigger for active faces
  faces.reserve(2 * n_faces);
  active_faces.reserve(n_faces);
  deleted_faces.reserve(8); // probably won't exceed this
  loose_edges.reserve(16);
  int id = 0;

  Vec3 ab = simplex.pts[1] - simplex.pts[0];
  Vec3 ac = simplex.pts[2] - simplex.pts[0];
  Vec3 ad = simplex.pts[3] - simplex.pts[0];
  Vec3 bc = simplex.pts[2] - simplex.pts[1];
  Vec3 bd = simplex.pts[3] - simplex.pts[1];

  Vec3 n       = cross(ab, ac);
  real dot_a_n = dot(simplex.pts[0], n);
  n            = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[3]) >= 0 ? -n : n));
  n            = (1 / norm(n)) * n;
  real dist    = distmin_origin({simplex.pts[0], simplex.pts[1], simplex.pts[2], n});
  faces.push_back({simplex.pts[0], simplex.pts[1], simplex.pts[2], n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(ab, ad);
  dot_a_n = dot(simplex.pts[0], n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[3]) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.pts[0], simplex.pts[1], simplex.pts[3], n});
  faces.push_back({simplex.pts[0], simplex.pts[1], simplex.pts[3], n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(ac, ad);
  dot_a_n = dot(simplex.pts[0], n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[3]) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.pts[0], simplex.pts[2], simplex.pts[3], n});
  faces.push_back({simplex.pts[0], simplex.pts[2], simplex.pts[3], n, dist});
  active_faces.push_back(id);
  id++;

  n       = cross(bc, bd);
  dot_a_n = dot(simplex.pts[3], n);
  n       = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[3]) >= 0 ? -n : n));
  n       = (1 / norm(n)) * n;
  dist    = distmin_origin({simplex.pts[1], simplex.pts[2], simplex.pts[3], n});
  faces.push_back({simplex.pts[1], simplex.pts[2], simplex.pts[3], n, dist});
  active_faces.push_back(id);
  id++;

  real min_dist;
  int  idx;

  // while (true) {
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    if (i == MAX_ITERATIONS - 1) {
      std::cout << "Max iterations reached - ericson EPA" << std::endl;
    }
    // Find closest face
    min_dist = INF_REAL;
    for (auto active_id: active_faces) {
      if (faces[active_id].dist < min_dist) {
        min_dist = faces[active_id].dist;
        idx      = active_id;
      }
    }

    // obtain a new support point in the direction of the edge normal
    Vec3 support1 = GJK_get_support(v1, count1, -faces[idx].n);
    Vec3 support2 = GJK_get_support(v2, count2, faces[idx].n);
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
//   Vec3 ab = simplex.pts[1] - simplex.pts[0];
//   Vec3 ac = simplex.pts[2] - simplex.pts[0];
//   Vec3 ad = simplex.pts[3] - simplex.pts[0];
//   Vec3 bc = simplex.pts[2] - simplex.pts[1];
//   Vec3 bd = simplex.pts[3] - simplex.pts[1];
//
//   Vec3 n1 = cross(ab, ac);
//   Vec3 n2 = cross(ab, ad);
//   Vec3 n3 = cross(ac, ad);
//   Vec3 n4 = cross(bc, bd);
//
//   real dot_a_n1 = dot(simplex.pts[0], n1);
//   real dot_a_n2 = dot(simplex.pts[0], n2);
//   real dot_a_n3 = dot(simplex.pts[0], n3);
//   real dot_d_n4 = dot(simplex.pts[3], n4);
//
//   n1 = dot_a_n1 > 0 ? n1 : (dot_a_n1 < 0 ? -n1 : (dot(n1, simplex.pts[3]) >= 0 ? -n1 : n1));
//   n2 = dot_a_n2 > 0 ? n2 : (dot_a_n2 < 0 ? -n2 : (dot(n2, simplex.pts[3]) >= 0 ? -n2 : n2));
//   n3 = dot_a_n3 > 0 ? n3 : (dot_a_n3 < 0 ? -n3 : (dot(n3, simplex.pts[3]) >= 0 ? -n3 : n3));
//   n4 = dot_d_n4 > 0 ? n4 : (dot_d_n4 < 0 ? -n4 : (dot(n4, simplex.pts[0]) >= 0 ? -n4 : n4));
//
//   n1 = (1 / norm(n1)) * n1;
//   n2 = (1 / norm(n2)) * n2;
//   n3 = (1 / norm(n3)) * n3;
//   n4 = (1 / norm(n4)) * n4;
//
//   faces.push_back({simplex.pts[0], simplex.pts[1], simplex.pts[2], n1});
//   faces.push_back({simplex.pts[0], simplex.pts[1], simplex.pts[3], n2});
//   faces.push_back({simplex.pts[0], simplex.pts[2], simplex.pts[3], n3});
//   faces.push_back({simplex.pts[1], simplex.pts[2], simplex.pts[3], n4});
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

host_fn real general_GJK(Vec3* set1, size_t count1, Vec3* set2, size_t count2) {
  // #if TRACY_ACtracy
  // ZoneScopedN("Ericson GJK");
  // Implemented from the basic algorithm described in Collision Detection manual by Ericson.

  // 1. Initializing simplex to a point from a random direction
  Vec3 direction = set1[0] - set2[0];

  Vec3 a1 = GJK_get_support(set1, count1, -direction);
  Vec3 a2 = GJK_get_support(set2, count2, direction);
  Vec3 V  = a2 - a1;

  Simplex simplex;
  simplex.pts[0] = V;
  simplex.size   = 1;
  real dot_P;

  // while (true) {
  for (int i = 0; i < MAX_ITERATIONS; i++) { // max iteration count
    if (i == MAX_ITERATIONS - 1) {
      std::cout << "max iterations reached - ericson GJK" << std::endl;
      break;
    }
    // 2. Computing the point P of minimum norm in CH(simplex)
    switch (simplex.size) {
      case 1:
        simplex.P = simplex.pts[0];
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
      real dist = solve_EPA(simplex, set1, count1, set2, count2);
      return dist;
    }

    // 4. Remove any points from simplex not determining the subsimplex of simplex in which P lies.
    // (Done automatically in GJK_solve_simplex functions)

    // 5. Find the next supporting point in direction -P
    a1 = GJK_get_support(set1, count1, simplex.P);
    a2 = GJK_get_support(set2, count2, -simplex.P);
    V  = a2 - a1;

    // 6. If V is no more exremal in direction -P than P itself, stop and return A and B as not intersecting.
    // The length of the vector from the origin to P is the separation distance of A and B.
    real ans1 = dot(V, -simplex.P) / dot_P;
    if (ans1 + 1 <= COLLISION_EPSILON) // No more progress is being made
      break;

    // 7. Add V to simplex and go to 2.
    Assert(simplex.size <= 3);

    simplex.pts[simplex.size] = V;
    simplex.size++;
  }

  return std::sqrt(dot_P);
}

host_fn real distance_GJK(const Capsule& caps, const Box& box) {
#if TRACY_ACTIVE >= 1 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScopedN("Ericson GJK");
#endif
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

  std::array<Vec3, 2> v1;
  v1[0] = caps.p1;
  v1[1] = caps.p2;

  std::array<Vec3, 8> v2;
  v2[0] = box.center + size_x + size_y + size_z;
  v2[1] = box.center + size_x + size_y - size_z;
  v2[2] = box.center + size_x - size_y + size_z;
  v2[3] = box.center + size_x - size_y - size_z;
  v2[4] = box.center - size_x + size_y + size_z;
  v2[5] = box.center - size_x + size_y - size_z;
  v2[6] = box.center - size_x - size_y + size_z;
  v2[7] = box.center - size_x - size_y - size_z;

  real dist = general_GJK(&v1[0], 2, &v2[0], 8) - caps.radius;
  return dist;
}

} // namespace blast
