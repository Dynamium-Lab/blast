#pragma once

#include <blast>
// #include "tracy/Tracy.hpp"

namespace blast {

#define TRACY_ACTIVE 0 // 0 nothing, 1 big functions, 2 medium functions, 3 all functions
#define GJK_ACTIVE 0   // 0 nothing, 1 ericson, 2 messy, 3 both
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
  Vec3 pts[4];
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

// these functions currently not used ////////////////////////////////////////////////////////////
inline host_fn real GJK_triangle_area_2d(real x1, real y1, real x2, real y2, real x3, real y3) {
#if TRACY_ACTIVE >= 3 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

inline host_fn Vec3 GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 c, Vec3 p) {
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

inline host_fn TwoPts GJK_get_local_points(const ComplexSimplex& simplex, Vec3 p) {
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
//////////////////////////////////////////////////////////////////////////////////////////////////

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

inline host_fn TwoPts GJK_solve_simplex2(ComplexSimplex& simplex) { // finds next search direction (ao) and closest point on line (ab) for 2pt simplex
#if TRACY_ACTIVE >= 2 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 ab = simplex.b - simplex.a;
  Vec3 ao = -simplex.a;

  const auto d = dot(ab, ao);
  if (d > 0) {
    simplex.count = 2;
    ao            = cross(cross(ab, ao), ab);
    ab            = simplex.a + (1.0 / dot(ab, ab)) * d * ab;
  } else {
    simplex.count = 1;
    ab            = simplex.a;
  }

  return {ab, ao};
}

inline host_fn TwoPts GJK_solve_simplex3(ComplexSimplex& simplex) {
#if TRACY_ACTIVE >= 2 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 ab  = simplex.b - simplex.a;
  Vec3 ac  = simplex.c - simplex.a;
  Vec3 ao  = -simplex.a;
  Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);

  if (dot(cross(abc, ac), ao) > 0) {
    // the origin is in the direction of the triangle normal
    if (dot(ac, ao) > 0) {
      // the origin is nearest to the line ac
      simplex.b     = simplex.c;
      simplex.b1    = simplex.c1;
      simplex.b2    = simplex.c2;
      simplex.count = 2;
      real t        = dot(ao, ac) / dot(ac, ac);
      return {simplex.a + t * ac, cross(cross(ac, ao), ac)};
    } else {
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

inline host_fn TwoPts GJK_solve_simplex4(ComplexSimplex& simplex) {
#if TRACY_ACTIVE >= 2 && GJK_ACTIVE >= 2
  ZoneScoped;
#endif
  Vec3 ab = simplex.b - simplex.a;
  Vec3 ac = simplex.c - simplex.a;
  Vec3 ad = simplex.d - simplex.a;
  Vec3 ao = -simplex.a;

  // face normals
  Vec3 abc = cross(ab, ac);
  Vec3 acd = cross(ac, ad);
  Vec3 adb = cross(ad, ab);

  int abc_outside = dot(abc, ao) > 0;
  int acd_outside = dot(acd, ao) > 0;
  int adb_outside = dot(adb, ao) > 0;

  // check face abc
  if (abc_outside && !acd_outside && !adb_outside) {
    // drop d
    simplex.count = 3;
    return GJK_solve_simplex3(simplex);
  }

  // 2. check face acd
  if (acd_outside && !abc_outside && !adb_outside) {
    // b = c, c = d
    simplex.b  = simplex.c;
    simplex.b1 = simplex.c1;
    simplex.b2 = simplex.c2;

    simplex.c  = simplex.d;
    simplex.c1 = simplex.d1;
    simplex.c2 = simplex.d2;

    simplex.count = 3;
    return GJK_solve_simplex3(simplex);
  }

  // check face adb
  if (adb_outside && !abc_outside && !acd_outside) {
    // b = d, c = b
    simplex.c  = simplex.b;
    simplex.c1 = simplex.b1;
    simplex.c2 = simplex.b2;

    simplex.b  = simplex.d;
    simplex.b1 = simplex.d1;
    simplex.b2 = simplex.d2;

    simplex.count = 3;
    return GJK_solve_simplex3(simplex);
  }

  // origin inside simplex
  if (!abc_outside && !acd_outside && !adb_outside) {
    return {{0, 0, 0}, {0, 0, 0}};
  }

  // Must check all the faces unfortunately // really really unfortunate
  simplex.count              = 3;
  ComplexSimplex simplex_abc = simplex;
  // simplex_abc.count          = 3;

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

  // simplex_acd.count = 3;

  ComplexSimplex simplex_adb = simplex;
  // simplex_adb.c_all = simplex_adb.b_all;
  simplex_adb.c  = simplex.b;
  simplex_adb.c1 = simplex.b1;
  simplex_adb.c2 = simplex.b2;
  // simplex_adb.c_all = simplex_adb.d_all;
  simplex_adb.b  = simplex.d;
  simplex_adb.b1 = simplex.d1;
  simplex_adb.b2 = simplex.d2;

  // simplex_adb.count = 3;

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

  if (abc_d2 <= acd_d2 && abc_d2 <= adb_d2) {
    simplex = simplex_abc;
    return {p_abc, dir_abc};
  } else if (acd_d2 <= adb_d2) {
    simplex = simplex_acd;
    return {p_acd, dir_acd};
  } else {
    simplex = simplex_adb;
    return {p_adb, dir_adb};
  }
}

// ======================================
//            EPA algorithm
// ======================================

inline host_fn real distmin_origin(const PolytopeFace& face) {
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
  float w     = 1.0f - u - v; // = vc / denom
  return norm(u * a + v * b + w * c);
}

inline host_fn real solve_EPA(ComplexSimplex simplex, Vec3* v1, size_t count1, Vec3* v2, size_t count2) {
#if TRACY_ACTIVE >= 1 && GJK_ACTIVE >= 2
  ZoneScopedN("Messy EPA");
#endif
  // create a face vector that has three points and a normal
  PolytopeFace faces[64];
  int          n_faces = 0;
  TwoPts       loose_edges[32];
  int          n_edges = 0;


  Vec3 ab = simplex.b - simplex.a;
  Vec3 ac = simplex.c - simplex.a;
  Vec3 ad = simplex.d - simplex.a;
  Vec3 bc = simplex.c - simplex.b;
  Vec3 bd = simplex.d - simplex.b;

  Vec3 n;
  real dot_a_n, dist;
  if (simplex.count == 3) {
    n                = cross(ab, ac);
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.a, simplex.b, simplex.c, n});
    faces[n_faces++] = {simplex.a, simplex.b, simplex.c, n, dist};
    faces[n_faces++] = {simplex.a, simplex.c, simplex.b, -n, dist};
  } else {
    n                = cross(ab, ac);
    dot_a_n          = dot(simplex.a, n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.a, simplex.b, simplex.c, n});
    faces[n_faces++] = {simplex.a, simplex.b, simplex.c, n, dist};

    n                = cross(ab, ad);
    dot_a_n          = dot(simplex.a, n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.a, simplex.b, simplex.d, n});
    faces[n_faces++] = {simplex.a, simplex.b, simplex.d, n, dist};

    n                = cross(ac, ad);
    dot_a_n          = dot(simplex.a, n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.d) >= 0 ? -n : n));
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.a, simplex.c, simplex.d, n});
    faces[n_faces++] = {simplex.a, simplex.c, simplex.d, n, dist};

    n                = cross(bc, bd);
    dot_a_n          = dot(simplex.d, n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.c) >= 0 ? -n : n)); // simplex.d -> [2]?
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.b, simplex.c, simplex.d, n});
    faces[n_faces++] = {simplex.b, simplex.c, simplex.d, n, dist};
  }
  real min_dist;
  int  idx;

  for (int j = 0; j < MAX_ITERATIONS; j++) {
    // Find closest face
    min_dist = INF_REAL;
    for (int k = 0; k < n_faces; k++) {
      if (faces[k].dist < min_dist) {
        min_dist = faces[k].dist;
        idx      = k;
      }
    }

    // obtain a new support point in the direction of the edge normal
    Vec3 support1 = GJK_get_support(v1, count1, -faces[idx].n);
    Vec3 support2 = GJK_get_support(v2, count2, faces[idx].n);
    Vec3 p        = support2 - support1;

    // If the vertex does not expand the polytope in the direction of the normal, the minimum distance
    // is with the closest face (unchanged). Compute and return.
    dist = dot(p - faces[idx].p1, faces[idx].n);
    if (dist < BLAST_GJK_EPSILON) {
      break;
    }

    // Get all the faces which are "seen" by the new point

    // Create new convex hull by determining which line segments are contained in the faces "seen" by new points.
    // These segments must be deleted and the ones which are only contained in one face will be added
    // as parts of the new faces (with the new point forming the remaining two line segments).
    TwoPts current_edge[3];
    n_edges = 0;

    for (int i = 0; i < n_faces;) {
      if (dot(faces[i].n, p - faces[i].p1) > 0) {
        int indx = i;

        current_edge[0] = {faces[indx].p1, faces[indx].p2};
        current_edge[1] = {faces[indx].p1, faces[indx].p3};
        current_edge[2] = {faces[indx].p2, faces[indx].p3};

        faces[indx] = std::move(faces[--n_faces]);

        int  loose_edge_idx = 0;
        bool found_edge     = false;
        for (const auto& edge: current_edge) {
          found_edge = false;
          for (int k = 0; k < n_edges; k++) { // Is the current edge already in loose_edges ?
            if ((edge.p1 == loose_edges[k].p1 && edge.p2 == loose_edges[k].p2) ||
                (edge.p2 == loose_edges[k].p1 && edge.p1 == loose_edges[k].p2)) {
              // edge is already in list
              found_edge     = true;
              loose_edge_idx = k;
              break;
            }
          }

          if (found_edge == false)
            loose_edges[n_edges++] = edge;
          else {
            loose_edges[loose_edge_idx] = std::move(loose_edges[--n_edges]);
          }
        }
      } else {
        i++;
      }
    }

    // rebuild simplex with new faces
    Vec3 n;
    real dot_p1_n;
    for (int i = 0; i < n_edges; i++) {
      n        = cross(p - loose_edges[i].p1, p - loose_edges[i].p2);
      dot_p1_n = dot(loose_edges[i].p1, n);
      n        = dot_p1_n > 0 ? n : -n;
      n        = (1 / norm(n)) * n;

      dist             = distmin_origin({loose_edges[i].p1, loose_edges[i].p2, p, n});
      faces[n_faces++] = {loose_edges[i].p1, loose_edges[i].p2, p, n, dist};
      if (n_faces == 64) { // remove?
        break;
      }
    }
    Assert(n_faces > 0);
    if (j == MAX_ITERATIONS - 1) {
      std::cout << "Max iterations reached - messy EPA" << std::endl;
    }
  }
  return -std::abs(min_dist);
}

inline host_fn real solve_general_GJK(Vec3* v1, size_t count1, Vec3* v2, size_t count2) {
  ComplexSimplex simplex;
  real           dist_min;

  Vec3 direction = v2[0] - v1[0];
  simplex.a1     = GJK_get_support(v1, count1, -direction);
  simplex.a2     = GJK_get_support(v2, count2, direction);
  simplex.a      = simplex.a2 - simplex.a1;
  simplex.count  = 1;

  TwoPts solved;
  Vec3   p;
  int    old_simplex_count;

  for (int i = 0; i < MAX_ITERATIONS; i++) {
    old_simplex_count = simplex.count;

    if (simplex.count == 1) {
      p         = simplex.a;
      direction = -simplex.a;
    } else {
      solved    = (simplex.count == 2) ? GJK_solve_simplex2(simplex) : ((simplex.count == 3) ? GJK_solve_simplex3(simplex) : GJK_solve_simplex4(simplex));
      p         = solved.p1; // current closest point on polytope
      direction = solved.p2; // original
      // direction = -solved.p1; // Try this
    }

    real dot_P = dot(p, p); // just for debugging remove later

    if (dot(p, p) < BLAST_GJK_EPSILON_SQ && old_simplex_count == 4) {
      dist_min = solve_EPA(simplex, v1, count1, v2, count2);
      Assert(isnan(dist_min) == false);
      break;
    }

    Vec3 support1 = GJK_get_support(v1, count1, -direction);
    Vec3 support2 = GJK_get_support(v2, count2, direction);

    Vec3 support     = support2 - support1;
    direction        = (1 / norm(direction)) * direction; // necessary in CI-float build
    real old_support = dot(p, direction);
    real new_support = dot(support, direction);

    // if no improvement, terminate
    if (new_support - old_support <= BLAST_GJK_EPSILON) {
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
            if (std::abs(old_support - dot(new_pt, direction)) < BLAST_GJK_EPSILON_SQ) {
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

      dist_min = norm(p);
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
    if (i == MAX_ITERATIONS - 1) {
      std::cout << "Max iterations reached - messy GJK" << std::endl;
      dist_min = norm(p);
    }
  }
  return dist_min;
}

inline host_fn real distance_GJK_simple(const Capsule& caps, const Box& box) {
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

inline host_fn void GJK_solve_simplex2_Ericson(Simplex* simplex) {
#if TRACY_ACTIVE >= 2 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScoped;
#endif
  Vec3 ab = (*simplex).pts[1] - (*simplex).pts[0];
  real t  = dot(-(*simplex).pts[0], ab) / dot(ab, ab);

  t = clamp(t, 0, 1);

  (*simplex).P = (*simplex).pts[0] + t * ab;

  if (t == 1) {
    (*simplex).pts[0] = (*simplex).pts[1];
    (*simplex).size   = 1;
  }

  return;
}

inline host_fn void GJK_solve_simplex3_Ericson(Simplex* simplex) {
#if TRACY_ACTIVE >= 2 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScoped;
#endif
  Vec3 a = simplex->pts[0];
  Vec3 b = simplex->pts[1];
  Vec3 c = simplex->pts[2];

  Vec3 ab = b - a;
  Vec3 ac = c - a;

  // P in vertex region outside A not possible unless coming from simplex4 (hence GJK_solve_simplex3_from4_Ericson())
  real d1 = dot(ab, -a);
  real d2 = dot(ac, -a);

  // P in vertex region outside B not possible unless coming from simplex4
  real d3 = dot(ab, -b);
  real d4 = dot(ac, -b);

  // P in edge region of AB not possible unless coming from simplex4
  real vc = d1 * d4 - d3 * d2;

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
    real w          = d2 / (d2 - d6);
    simplex->size   = 2;
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

inline host_fn void GJK_simplex3_from4_Ericson(Simplex* simplex) {
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
    simplex->P    = a;
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
    simplex->P    = a + v * ab;
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
    real w          = d2 / (d2 - d6);
    simplex->size   = 2;
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

inline host_fn void GJK_solve_simplex4_Ericson(Simplex* simplex) {
#if TRACY_ACTIVE >= 2 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScoped;
#endif
  Vec3 a = simplex->pts[0];
  Vec3 b = simplex->pts[1];
  Vec3 c = simplex->pts[2];
  Vec3 d = simplex->pts[3]; // Most recently added support point in GJK

  Vec3 da  = a - d;
  Vec3 db  = b - d;
  Vec3 dc  = c - d;
  Vec3 do_ = -d; // Vector from d to Origin (p = {0,0,0})

  // 1. Check for duplicate points
  Vec3 dbc    = cross(db, dc);
  real volume = dot(da, dbc);

  if (std::abs(volume) < BLAST_GJK_EPSILON) {
    simplex->size = 3;
    if (dot(da, da) < BLAST_GJK_EPSILON) {
      simplex->pts[0] = b;
      simplex->pts[1] = c;
      simplex->pts[2] = d;
    } else if (dot(db, db) < BLAST_GJK_EPSILON) {
      simplex->pts[1] = c;
      simplex->pts[2] = d;
    } else if (dot(dc, dc) < BLAST_GJK_EPSILON) {
      simplex->pts[2] = d;
    } else {
      // simplex->pts[0] = b; // this is not right
      // simplex->pts[1] = c;
      // simplex->pts[2] = d;

      Simplex simplex_abd = {{0, 0, 0}, {simplex->pts[0], simplex->pts[1], simplex->pts[3]}, 3};
      Simplex simplex_acd = {{0, 0, 0}, {simplex->pts[0], simplex->pts[2], simplex->pts[3]}, 3};
      Simplex simplex_bcd = {{0, 0, 0}, {simplex->pts[1], simplex->pts[2], simplex->pts[3]}, 3};

      GJK_simplex3_from4_Ericson(&simplex_abd);
      GJK_simplex3_from4_Ericson(&simplex_acd);
      GJK_simplex3_from4_Ericson(&simplex_bcd);

      real dotP     = dot(simplex->P, simplex->P);
      real dotP_abd = dot(simplex_abd.P, simplex_abd.P);
      real dotP_acd = dot(simplex_acd.P, simplex_acd.P);
      real dotP_bcd = dot(simplex_bcd.P, simplex_bcd.P);

      if (dotP <= dotP_abd && dotP <= dotP_acd && dotP <= dotP_bcd) {
        simplex->size = 3;
        return;
      } else if (dotP_abd <= dotP_acd && dotP_abd <= dotP_bcd) {
        *simplex = simplex_abd;
        return;
      } else if (dotP_acd <= dotP_bcd) {
        *simplex = simplex_acd;
        return;
      } else {
        *simplex = simplex_bcd;
        return;
      }
    }
    GJK_simplex3_from4_Ericson(simplex);
    return;
  }

  // 2. Compute face normals for faces sharing point d
  Vec3 dca = cross(dc, da);
  Vec3 dab = cross(da, db);

  // Determine if Origin is on the outer side of each face plane relative to the opposite vertex
  int bcd_outside = (dot(dbc, do_) * dot(dbc, a - d)) <= 0; // face bcd
  int acd_outside = (dot(dca, do_) * dot(dca, b - d)) <= 0; // face acd
  int adb_outside = (dot(dab, do_) * dot(dab, c - d)) <= 0; // face abd

  // 3. Outward Evaluation
  // d is most recent point -> the origin cannot lie on the far side of face abc. todo: is this true?

  if (bcd_outside && !acd_outside && !adb_outside) {
    // Origin is outside face bcd
    simplex->pts[0] = b;
    simplex->pts[1] = c;
    simplex->pts[2] = d;
    simplex->size   = 3;
    GJK_simplex3_from4_Ericson(simplex);
    return;
  }

  if (acd_outside && !bcd_outside && !adb_outside) {
    // Origin is outside face acd
    // simplex->pts[0] = a;
    simplex->pts[1] = c;
    simplex->pts[2] = d;
    simplex->size   = 3;
    GJK_simplex3_from4_Ericson(simplex);
    return;
  }

  if (adb_outside && !bcd_outside && !acd_outside) {
    // Origin is outside face adb
    // simplex->pts[0] = a;
    // simplex->pts[1] = b;
    simplex->pts[2] = d;
    simplex->size   = 3;
    GJK_simplex3_from4_Ericson(simplex);
    return;
  }

  if (!bcd_outside && !acd_outside && !adb_outside) {
    // 4. Origin is inside all half-spaces -> intersection
    simplex->size = 4;
    simplex->P    = Vec3{0, 0, 0};
    return;
  }

  Simplex simplex_abd = {{0, 0, 0}, {simplex->pts[0], simplex->pts[1], simplex->pts[3]}, 3};
  Simplex simplex_acd = {{0, 0, 0}, {simplex->pts[0], simplex->pts[2], simplex->pts[3]}, 3};
  Simplex simplex_bcd = {{0, 0, 0}, {simplex->pts[1], simplex->pts[2], simplex->pts[3]}, 3};

  GJK_simplex3_from4_Ericson(&simplex_abd);
  GJK_simplex3_from4_Ericson(&simplex_acd);
  GJK_simplex3_from4_Ericson(&simplex_bcd);

  real dotP     = dot(simplex->P, simplex->P);
  real dotP_abd = dot(simplex_abd.P, simplex_abd.P);
  real dotP_acd = dot(simplex_acd.P, simplex_acd.P);
  real dotP_bcd = dot(simplex_bcd.P, simplex_bcd.P);

  if (dotP <= dotP_abd && dotP <= dotP_acd && dotP <= dotP_bcd) {
    simplex->size = 3;
    return;
  } else if (dotP_abd <= dotP_acd && dotP_abd <= dotP_bcd) {
    *simplex = simplex_abd;
    return;
  } else if (dotP_acd <= dotP_bcd) {
    *simplex = simplex_acd;
    return;
  } else {
    *simplex = simplex_bcd;
    return;
  }
}

inline host_fn real solve_EPA(Simplex simplex, Vec3* v1, size_t count1, Vec3* v2, size_t count2) {
#if TRACY_ACTIVE >= 1 && (GJK_ACTIVE == 1 || GJK_ACTIVE == 3)
  ZoneScopedN("Ericson EPA");
#endif
  // create a face vector that has three points and a normal
  PolytopeFace faces[64]; // change size probably
  int          n_faces = 0;
  TwoPts       loose_edges[32];
  int          n_edges = 0;


  Vec3 ab = simplex.pts[1] - simplex.pts[0];
  Vec3 ac = simplex.pts[2] - simplex.pts[0];
  Vec3 ad = simplex.pts[3] - simplex.pts[0];
  Vec3 bc = simplex.pts[2] - simplex.pts[1];
  Vec3 bd = simplex.pts[3] - simplex.pts[1];

  Vec3 n;
  real dot_a_n, dist;
  if (simplex.size == 3) {
    n                = cross(ab, ac);
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.pts[0], simplex.pts[1], simplex.pts[2], n});
    faces[n_faces++] = {simplex.pts[0], simplex.pts[1], simplex.pts[2], n, dist};
    faces[n_faces++] = {simplex.pts[0], simplex.pts[2], simplex.pts[1], -n, dist};
  } else {

    Vec3 n           = cross(ab, ac);
    real dot_a_n     = dot(simplex.pts[0], n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[3]) >= 0 ? -n : n));
    n                = (1 / norm(n)) * n;
    real dist        = distmin_origin({simplex.pts[0], simplex.pts[1], simplex.pts[2], n});
    faces[n_faces++] = {simplex.pts[0], simplex.pts[1], simplex.pts[2], n, dist};

    n                = cross(ab, ad);
    dot_a_n          = dot(simplex.pts[0], n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[3]) >= 0 ? -n : n));
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.pts[0], simplex.pts[1], simplex.pts[3], n});
    faces[n_faces++] = {simplex.pts[0], simplex.pts[1], simplex.pts[3], n, dist};

    n                = cross(ac, ad);
    dot_a_n          = dot(simplex.pts[0], n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[3]) >= 0 ? -n : n));
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.pts[0], simplex.pts[2], simplex.pts[3], n});
    faces[n_faces++] = {simplex.pts[0], simplex.pts[2], simplex.pts[3], n, dist};

    n                = cross(bc, bd);
    dot_a_n          = dot(simplex.pts[3], n);
    n                = dot_a_n > 0 ? n : (dot_a_n < 0 ? -n : (dot(n, simplex.pts[2]) >= 0 ? -n : n)); // simplex.pts[3] -> [2]?
    n                = (1 / norm(n)) * n;
    dist             = distmin_origin({simplex.pts[1], simplex.pts[2], simplex.pts[3], n});
    faces[n_faces++] = {simplex.pts[1], simplex.pts[2], simplex.pts[3], n, dist};
  }
  real min_dist;
  int  idx;

  for (int j = 0; j < MAX_ITERATIONS; j++) {
    // Find closest face
    min_dist = INF_REAL;
    for (int k = 0; k < n_faces; k++) {
      if (faces[k].dist < min_dist) {
        min_dist = faces[k].dist;
        idx      = k;
      }
    }

    // obtain a new support point in the direction of the edge normal
    Vec3 support1 = GJK_get_support(v1, count1, -faces[idx].n);
    Vec3 support2 = GJK_get_support(v2, count2, faces[idx].n);
    Vec3 p        = support2 - support1;

    // If the vertex does not expand the polytope in the direction of the normal, the minimum distance
    // is with the closest face (unchanged). Compute and return.
    dist = dot(p - faces[idx].p1, faces[idx].n);
    if (dist < BLAST_GJK_EPSILON) {
      break;
    }

    // Get all the faces which are "seen" by the new point

    // Create new convex hull by determining which line segments are contained in the faces "seen" by new points.
    // These segments must be deleted and the ones which are only contained in one face will be added
    // as parts of the new faces (with the new point forming the remaining two line segments).
    TwoPts current_edge[3];
    n_edges = 0;

    for (int i = 0; i < n_faces;) {
      if (dot(faces[i].n, p - faces[i].p1) > 0) {
        int indx = i;

        current_edge[0] = {faces[indx].p1, faces[indx].p2};
        current_edge[1] = {faces[indx].p1, faces[indx].p3};
        current_edge[2] = {faces[indx].p2, faces[indx].p3};

        faces[indx] = std::move(faces[--n_faces]);

        int  loose_edge_idx = 0;
        bool found_edge     = false;
        for (const auto& edge: current_edge) {
          found_edge = false;
          for (int k = 0; k < n_edges; k++) { // Is the current edge already in loose_edges ?
            if ((edge.p1 == loose_edges[k].p1 && edge.p2 == loose_edges[k].p2) ||
                (edge.p2 == loose_edges[k].p1 && edge.p1 == loose_edges[k].p2)) {
              // edge is already in list
              found_edge     = true;
              loose_edge_idx = k;
              break;
            }
          }

          if (found_edge == false)
            loose_edges[n_edges++] = edge;
          else {
            loose_edges[loose_edge_idx] = std::move(loose_edges[--n_edges]);
          }
        }
      } else {
        i++;
      }
    }

    // rebuild simplex with new faces
    // PolytopeFace new_face;
    Vec3 n;
    real dot_p1_n;
    for (int i = 0; i < n_edges; i++) {
      n        = cross(p - loose_edges[i].p1, p - loose_edges[i].p2);
      dot_p1_n = dot(loose_edges[i].p1, n);
      n        = dot_p1_n > 0 ? n : -n;
      n        = (1 / norm(n)) * n;

      dist             = distmin_origin({loose_edges[i].p1, loose_edges[i].p2, p, n});
      faces[n_faces++] = {loose_edges[i].p1, loose_edges[i].p2, p, n, dist};
      Assert(n_faces < 64);
      if (n_faces == 64) {
        std::cout << "Max capacity faces" << std::endl;
        break;
      }
    }
    Assert(n_faces > 0);
    if (j == MAX_ITERATIONS - 1) {
      std::cout << "Max iterations reached - messy EPA" << std::endl;
    }
  }
  return -std::abs(min_dist);
}

inline host_fn real general_GJK(Vec3* set1, size_t count1, Vec3* set2, size_t count2) {
  // Implemented from the basic algorithm described in Collision Detection manual by Ericson.

  // 1. Initializing simplex to a point from a random direction
  Vec3 direction = set2[0] - set1[0];

  Vec3 a1 = GJK_get_support(set1, count1, -direction);
  Vec3 a2 = GJK_get_support(set2, count2, direction);
  Vec3 V  = a2 - a1;

  Simplex simplex;
  simplex.pts[0] = V;
  simplex.size   = 1;
  real dot_P;

  for (int i = 0; i < MAX_ITERATIONS; i++) { // max iteration count
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
    if (dot_P < BLAST_GJK_EPSILON_SQ && simplex.size >= 3) {
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
    Vec3 direction   = -simplex.P;
    direction        = (1 / norm(direction)) * direction; // necessary in CI-float build
    real old_support = dot(simplex.P, direction);
    real new_support = dot(V, direction);
    if (new_support - old_support <= BLAST_GJK_EPSILON)
      break;

    // 7. Add V to simplex and go to 2.
    Assert(simplex.size <= 3);

    simplex.pts[simplex.size] = V;
    simplex.size++;

    if (i == MAX_ITERATIONS - 1) {
      std::cout << "Max iterations reached - ericson GJK" << std::endl;
    }
  }

  return std::sqrt(dot_P);
}

inline host_fn real distance_GJK(const Capsule& caps, const Box& box) {
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
