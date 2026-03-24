#pragma once

#include <blast>
#include <cmath>

namespace blast {

struct Segment {
  Vec3 p1;
  Vec3 p2;
};

struct Surface {
  Vec3 p;
  Vec3 d1;
  Vec3 d2;
};

struct Plane {
  Vec3 p;
  Vec3 n;
};

struct TwoPts {
  Vec3 p1;
  Vec3 p2;
};

struct Circle {
  Vec3 p;
  Vec3 n;
  real r = 0;
};

struct Triangle {
  Vec3 p1;
  Vec3 p2;
  Vec3 p3;
  real distance = 0;
};

struct OriginVertices {
  Vec3 vert[8];
};

struct Directions {
  Vec3 dir[3];
};

struct Faces {
  Surface surface[12];
};

struct ActiveEdges {
  Segment segments[24];
};

struct SegFace {
  Segment segments[4];
};

constexpr real COLLISION_EPSILON = 1e-9;

// Internal functions not directly exposed to the API
inline blast_fn bool point_in_triangle(Vec3 triangle_v1, Vec3 triangle_v2, Vec3 triangle_origin, Vec3 point);
inline blast_fn bool point_in_surface(Vec3 surf_v1, Vec3 surf_v2, Vec3 surf_origin, Vec3 point);

inline blast_fn Vec3   closest_point(Segment segment, Vec3 point);
inline blast_fn Vec3   closest_point_origin(Segment segment);
inline blast_fn Vec3   closest_point_plane(Vec3 q, Plane p);
inline blast_fn TwoPts closest_points(Segment segment1, Segment segment2);

inline blast_fn real distance(Segment segment, Vec3 point);
inline blast_fn real distance(Surface surface, Vec3 point);
inline blast_fn real distance(Segment segment1, Segment segment2);


inline blast_fn bool point_in_triangle(Vec3 triangle_v1, Vec3 triangle_v2, Vec3 triangle_origin, Vec3 point) {
  // u=P2−P1
  Vec3 u = triangle_v1;
  // v=P3−P1
  Vec3 v = triangle_v2;
  // n=u×v
  Vec3 n = cross(u, v);
  // w=P−P1
  Vec3 w = point - triangle_origin;
  // Barycentric coordinates of the projection point′of point onto T:
  // γ=[(u×w)⋅n]/n²
  real gamma = dot(cross(u, w), n) / dot(n, n);
  // β=[(w×v)⋅n]/n²
  real beta  = dot(cross(w, v), n) / dot(n, n);
  real alpha = 1 - gamma - beta;
  // The point P′ lies inside T if:
  return ((0 <= alpha) && (alpha <= 1) && (0 <= beta) && (beta <= 1) && (0 <= gamma) && (gamma <= 1));
}

inline blast_fn bool point_in_surface(Vec3 surf_v1, Vec3 surf_v2, Vec3 surf_origin, Vec3 point) {
  bool tri1 = point_in_triangle(surf_v1, surf_v2, surf_origin, point);
  bool tri2 = point_in_triangle(-surf_v1, -surf_v2, surf_origin + surf_v1 + surf_v2, point);
  return (tri1 || tri2);
}

inline blast_fn Vec3 closest_point(Segment segment, Vec3 point) {
  Vec3 ab = segment.p2 - segment.p1;
  real t  = dot(point - segment.p1, ab) / dot(ab, ab);
  t       = clamp(t, 0, 1);
  Vec3 d  = segment.p1 + t * ab;
  return d;
}

inline blast_fn Vec3 closest_point_origin(Segment segment) {
  Vec3 ab = segment.p2 - segment.p1;
  real t  = dot(-segment.p1, ab) / dot(ab, ab);

  t = clamp(t, 0, 1);

  Vec3 d = segment.p1 + t * ab;
  return d;
}

inline blast_fn Vec3 closest_point_plane(Vec3 q, Plane p) {
  return q - (dot(p.n, q - p.p) / dot(p.n, p.n)) * p.n;
  ;
}

inline blast_fn real distance(Surface surface, Vec3 point) {
  const Vec3 n = cross(surface.d1, surface.d2);
  if (dot(n, n) < COLLISION_EPSILON)
    return -INF_REAL; // the plane is non-existant
  const Vec3 n_unit = 1 / norm(n) * n;

  // If the surface is a rectangular shape, it is much easier to treat.
  if (std::abs(dot(surface.d1, surface.d2)) <= COLLISION_EPSILON) {
    // norm^2 d1 & d2
    const real val_d1_sq = dot(surface.d1, surface.d1);
    const real val_d2_sq = dot(surface.d2, surface.d2);
    // direction vector from p to point
    const Vec3 val_direction = point - surface.p;

    // How far along is the point in d1 & d2 direction
    const real val_t1          = dot(val_direction, surface.d1) / val_d1_sq;
    const real val_t2          = dot(val_direction, surface.d2) / val_d2_sq;
    const real val_t1_clamped  = clamp(val_t1, 0, 1);
    const real val_t2_clamped  = clamp(val_t2, 0, 1);
    const real val_normal_dist = dot(val_direction, n_unit);

    const real val_dist1 = (val_t1 - val_t1_clamped) * (val_t1 - val_t1_clamped) * val_d1_sq;
    const real val_dist2 = (val_t2 - val_t2_clamped) * (val_t2 - val_t2_clamped) * val_d2_sq;

    const real val_result = (real) sqrt(val_dist1 + val_dist2 + val_normal_dist * val_normal_dist);

    const real val_distance = (val_t1 >= 0 && val_t1 <= 1 && val_t2 >= 0 && val_t2 <= 1) ? val_normal_dist : val_result;
    return val_distance == val_normal_dist ? val_normal_dist : (val_normal_dist < 0 ? -val_distance : val_distance);
  } else {
    const real normal_distance = dot(point - surface.p, n_unit);
    // Vec3       direction  = point - normal_distance * n_unit - surface.p;

    if (point_in_surface(surface.d1, surface.d2, surface.p, point)) {
      return normal_distance;
    }

    Segment segment[4];
    segment[0].p1 = surface.p;
    segment[0].p2 = surface.p + surface.d1;
    segment[1].p1 = surface.p;
    segment[1].p2 = surface.p + surface.d2;
    segment[2].p1 = surface.p + surface.d1;
    segment[2].p2 = surface.p + surface.d1 + surface.d2;
    segment[3].p1 = surface.p + surface.d2;
    segment[3].p2 = surface.p + surface.d2 + surface.d1;

    real dist_min = INF_REAL;
    int  idx      = 0;
    for (int i = 0; i < 4; i++) {
      auto dist_tmp = distance(segment[i], point);
      dist_min      = dist_tmp < dist_min ? dist_tmp : dist_min;
      idx           = dist_tmp < dist_min ? i : idx;
    }
    real resultant = normal_distance < 0 ? -dist_min : dist_min;
    return resultant;
  }
}

inline blast_fn TwoPts closest_points(Segment segment1, Segment segment2) {
  TwoPts result;
  real   s;
  real   t;
  Vec3   c1;
  Vec3   c2;

  Vec3 d1 = segment1.p2 - segment1.p1; // Direction vector of segment S1
  Vec3 d2 = segment2.p2 - segment2.p1; // Direction vector of segment S2
  Vec3 r  = segment1.p1 - segment2.p1;
  real a  = dot(d1, d1);               // Squared length of segment S1, always nonnegative
  real e  = dot(d2, d2);               // Squared length of segment S2, always nonnegative
  real f  = dot(d2, r);
  // Check if either or both segments degenerate into points
  if (a <= COLLISION_EPSILON && e <= COLLISION_EPSILON) {
    // Both segments degenerate into points
    c1        = segment1.p2;
    c2        = segment2.p2;
    result.p1 = c1;
    result.p2 = c2;
    return result;
  }
  if (a <= COLLISION_EPSILON) {
    // First segment degenerates into a point
    s = 0.0f;
    t = f / e; // s = 0 => t = (b*s + f) / e = f / e
    t = clamp(t, 0.0f, 1.0f);
  } else {
    real c = dot(d1, r);
    if (e <= COLLISION_EPSILON) {
      // Second segment degenerates into a point
      t = 0.0f;
      s = clamp(-c / a, 0.0f, 1.0f); // t = 0 => s = (b*t - c) / a = -c / a
    } else {
      // The general nondegenerate case starts here
      real b = dot(d1, d2);
      // If segments not parallel, compute the closest point on L1 to L2 and
      // clamp to segment S1. Else pick arbitrary s (here 0)
      if (real denominator = a * e - b * b; denominator != 0.0f) {
        s = clamp((b * f - c * e) / denominator, 0.0f, 1.0f);
      } else
        s = 0.0f;
      // Compute point on L2 closest to S1(s) using
      // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
      t = (b * s + f) / e;
      // If t in [0,1] done. Else clamp t, recompute s for the new value
      // of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a
      // and clamp s to [0, 1]
      if (t < 0.0f) {
        t = 0.0f;
        s = clamp(-c / a, 0.0f, 1.0f);
      } else if (t > 1.0f) {
        t = 1.0f;
        s = clamp((b - c) / a, 0.0f, 1.0f);
      }
    }
  }
  c1        = segment1.p1 + d1 * s;
  c2        = segment2.p1 + d2 * t;
  result.p1 = c1;
  result.p2 = c2;
  return result;
}

inline blast_fn real distance(Segment segment1, Segment segment2) {
  auto [p1, p2] = closest_points(segment1, segment2);
  return norm(p1 - p2);
}

inline blast_fn real distance(Segment segment, Vec3 point) {
  Vec3 d = closest_point(segment, point);
  return norm(d - point);
}

inline blast_fn real distance(const Capsule& capsule, const Box& box) {
  // Transferring points to OBB frame
  Mat3 R_trans = transpose(box.rotation);

  Vec3 p1 = R_trans * (capsule.p1 - box.center);
  Vec3 p2 = R_trans * (capsule.p2 - box.center);

  Vec3 vec   = p2 - p1;
  Vec3 point = vec.z < 0 ? p1 : p2; // take point of highest z value
  vec        = vec.z < 0 ? -vec : vec;

  // A
  // real x_min = - box.extents.x;
  // real x_max =   box.extents.x;
  // real y_min = - box.extents.y;
  // real y_max =   box.extents.y;
  // real z_min = - box.extents.z;
  // real z_max =   box.extents.z;

  // Creating the eight original vertices
  OriginVertices origin_vertices;
  origin_vertices.vert[0] = {-box.extents.x, -box.extents.y, -box.extents.z};
  origin_vertices.vert[1] = {-box.extents.x, box.extents.y, -box.extents.z};
  origin_vertices.vert[2] = {box.extents.x, -box.extents.y, -box.extents.z};
  origin_vertices.vert[3] = {box.extents.x, box.extents.y, -box.extents.z};
  origin_vertices.vert[4] = {-box.extents.x, -box.extents.y, box.extents.z};
  origin_vertices.vert[5] = {-box.extents.x, box.extents.y, box.extents.z};
  origin_vertices.vert[6] = {box.extents.x, -box.extents.y, box.extents.z};
  origin_vertices.vert[7] = {box.extents.x, box.extents.y, box.extents.z};

  // Thus, we get the three main directions
  Directions directions;
  directions.dir[0] = {2 * box.extents.x, 0, 0};
  directions.dir[1] = {0, 2 * box.extents.y, 0};
  directions.dir[2] = {0, 0, 2 * box.extents.z};

  // some faces depend only on the direction of vec in x
  Faces face;
  // the top four vertices will always be extended, meaning that the top and bottom faces (faces 0-1) will always be of the same format:
  face.surface[0].p  = origin_vertices.vert[4] + vec;
  face.surface[0].d1 = directions.dir[0];
  face.surface[0].d2 = directions.dir[1];

  face.surface[1].p  = origin_vertices.vert[0];
  face.surface[1].d1 = directions.dir[1];
  face.surface[1].d2 = directions.dir[0];

  // The other faces depend on the orientation of the vector either in x (faces 2-5),
  // in y (faces 6-9), or a combination of both (faces 10-11)
  bool positive_x = vec.x >= 0;
  bool positive_y = vec.y >= 0;
  // bool positive_z = vec.z >= 0;

  /*face.surface[2].p  = positive_x ? origin_vertices.vert[5] : origin_vertices.vert[6];
  face.surface[2].d1 = positive_x ? -directions.dir[1] : directions.dir[1];
  face.surface[2].d2 = vec;

  face.surface[3].p  = positive_x ? origin_vertices.vert[1] : origin_vertices.vert[2];
  face.surface[3].d1 = positive_x ? -directions.dir[1] : directions.dir[1];
  face.surface[3].d2 = directions.dir[2];

  face.surface[4].p  = positive_x ? origin_vertices.vert[2] + vec : origin_vertices.vert[1] + vec;
  face.surface[4].d1 = positive_x ? directions.dir[1] : -directions.dir[1];
  face.surface[4].d2 = directions.dir[2];

  face.surface[5].p  = positive_x ? origin_vertices.vert[2] : origin_vertices.vert[1];
  face.surface[5].d1 = positive_x ? directions.dir[1] : -directions.dir[1];
  face.surface[5].d2 = vec;

  face.surface[6].p  = positive_y ? origin_vertices.vert[4] : origin_vertices.vert[7];
  face.surface[6].d1 = positive_y ? directions.dir[0] : -directions.dir[0];
  face.surface[6].d2 = vec;

  face.surface[7].p  = positive_y ? origin_vertices.vert[0] : origin_vertices.vert[3];
  face.surface[7].d1 = positive_y ? directions.dir[0] : -directions.dir[0];
  face.surface[7].d2 = directions.dir[2];

  face.surface[8].p  = positive_y ? origin_vertices.vert[3] + vec : origin_vertices.vert[0] + vec;
  face.surface[8].d1 = positive_y ? -directions.dir[0] : directions.dir[0];
  face.surface[8].d2 = directions.dir[2];

  face.surface[9].p  = positive_y ? origin_vertices.vert[3] : origin_vertices.vert[0];
  face.surface[9].d1 = positive_y ? -directions.dir[0] : directions.dir[0];
  face.surface[9].d2 = vec;

  face.surface[10].p  = positive_x ? (positive_y ? origin_vertices.vert[1] : origin_vertices.vert[4]) : (positive_y ? origin_vertices.vert[7] : origin_vertices.vert[2]);
  face.surface[10].d1 = positive_x ? (positive_y ? directions.dir[2] : -directions.dir[2]) : (positive_y ? -directions.dir[2] : directions.dir[2]);
  face.surface[10].d2 = vec;

  face.surface[11].p  = positive_x ? (positive_y ? origin_vertices.vert[6] : origin_vertices.vert[3]) : (positive_y ? origin_vertices.vert[0] : origin_vertices.vert[5]);
  face.surface[11].d1 = positive_x ? (positive_y ? -directions.dir[2] : directions.dir[2]) : (positive_y ? directions.dir[2] : -directions.dir[2]);
  face.surface[11].d2 = vec;
  */

  // Faces 2-5 depend on the x coordinate
  if (positive_x) {
    face.surface[2].p  = origin_vertices.vert[5];
    face.surface[2].d1 = -directions.dir[1];
    face.surface[2].d2 = vec;

    face.surface[3].p  = origin_vertices.vert[1];
    face.surface[3].d1 = -directions.dir[1];
    face.surface[3].d2 = directions.dir[2];

    face.surface[4].p  = origin_vertices.vert[2] + vec;
    face.surface[4].d1 = directions.dir[1];
    face.surface[4].d2 = directions.dir[2];

    face.surface[5].p  = origin_vertices.vert[2];
    face.surface[5].d1 = directions.dir[1];
    face.surface[5].d2 = vec;
  } else {
    face.surface[2].p  = origin_vertices.vert[6];
    face.surface[2].d1 = directions.dir[1];
    face.surface[2].d2 = vec;

    face.surface[3].p  = origin_vertices.vert[2];
    face.surface[3].d1 = directions.dir[1];
    face.surface[3].d2 = directions.dir[2];

    face.surface[4].p  = origin_vertices.vert[1] + vec;
    face.surface[4].d1 = -directions.dir[1];
    face.surface[4].d2 = directions.dir[2];

    face.surface[5].p  = origin_vertices.vert[1];
    face.surface[5].d1 = -directions.dir[1];
    face.surface[5].d2 = vec;
  }

  // Faces 6-9 depend on the y coordinate
  if (positive_y) {
    face.surface[6].p  = origin_vertices.vert[4];
    face.surface[6].d1 = directions.dir[0];
    face.surface[6].d2 = vec;

    face.surface[7].p  = origin_vertices.vert[0];
    face.surface[7].d1 = directions.dir[0];
    face.surface[7].d2 = directions.dir[2];

    face.surface[8].p  = origin_vertices.vert[3] + vec;
    face.surface[8].d1 = -directions.dir[0];
    face.surface[8].d2 = directions.dir[2];

    face.surface[9].p  = origin_vertices.vert[3];
    face.surface[9].d1 = -directions.dir[0];
    face.surface[9].d2 = vec;
  } else {
    face.surface[6].p  = origin_vertices.vert[7];
    face.surface[6].d1 = -directions.dir[0];
    face.surface[6].d2 = vec;

    face.surface[7].p  = origin_vertices.vert[3];
    face.surface[7].d1 = -directions.dir[0];
    face.surface[7].d2 = directions.dir[2];

    face.surface[8].p  = origin_vertices.vert[0] + vec;
    face.surface[8].d1 = directions.dir[0];
    face.surface[8].d2 = directions.dir[2];

    face.surface[9].p  = origin_vertices.vert[0];
    face.surface[9].d1 = directions.dir[0];
    face.surface[9].d2 = vec;
  }

  // Faces 10-11 depend on the x and y coordinate
  if (positive_x && positive_y) {
    face.surface[10].p  = origin_vertices.vert[1];
    face.surface[10].d1 = directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[6];
    face.surface[11].d1 = -directions.dir[2];
    face.surface[11].d2 = vec;
  } else if (positive_x && !positive_y) {
    face.surface[10].p  = origin_vertices.vert[4];
    face.surface[10].d1 = -directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[3];
    face.surface[11].d1 = directions.dir[2];
    face.surface[11].d2 = vec;
  } else if (!positive_x && positive_y) {
    face.surface[10].p  = origin_vertices.vert[7];
    face.surface[10].d1 = -directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[0];
    face.surface[11].d1 = directions.dir[2];
    face.surface[11].d2 = vec;
  } else {
    face.surface[10].p  = origin_vertices.vert[2];
    face.surface[10].d1 = directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[5];
    face.surface[11].d1 = -directions.dir[2];
    face.surface[11].d2 = vec;
  }

  ActiveEdges active_edges;
  int         n_active_edges = 0;
  real        normal_distance[12];
  real        max_normal_dist = -INF_REAL;

  SegFace seg_face;

  for (int i = 0; i < 12; i++) {
    Vec3 n_current     = cross(face.surface[i].d1, face.surface[i].d2);
    n_current          = 1 / norm(n_current) * n_current;
    normal_distance[i] = dot(n_current, point - face.surface[i].p);
    max_normal_dist    = normal_distance[i] > max_normal_dist ? normal_distance[i] : max_normal_dist;

    if (normal_distance[i] >= 0) {
      if (point_in_surface(face.surface[i].d1, face.surface[i].d2, face.surface[i].p, point))
        return normal_distance[i] - capsule.radius;

      seg_face.segments[0].p1 = face.surface[i].p;
      seg_face.segments[0].p2 = face.surface[i].p + face.surface[i].d1;
      seg_face.segments[1].p1 = face.surface[i].p;
      seg_face.segments[1].p2 = face.surface[i].p + face.surface[i].d2;
      seg_face.segments[2].p1 = face.surface[i].p + face.surface[i].d1;
      seg_face.segments[2].p2 = face.surface[i].p + face.surface[i].d1 + face.surface[i].d2;
      seg_face.segments[3].p1 = face.surface[i].p + face.surface[i].d2;
      seg_face.segments[3].p2 = face.surface[i].p + face.surface[i].d1 + face.surface[i].d2;

      // active_faces[n_active_faces++] = face[i];
      for (const auto& current_seg: seg_face.segments) {
        bool is_in_list = false;
        for (int k = 0; k < n_active_edges; k++) {
          if (current_seg.p1 == active_edges.segments[k].p1 && current_seg.p2 == active_edges.segments[k].p2 ||
              current_seg.p1 == active_edges.segments[k].p2 && current_seg.p2 == active_edges.segments[k].p1) {
            is_in_list = true;
            break;
          }
        }
        if (!is_in_list)
          active_edges.segments[n_active_edges++] = current_seg;
      }
    }
  }

  // If point is inside
  if (n_active_edges == 0) {
    return max_normal_dist - capsule.radius;
  }

  real dist_min = INF_REAL;
  for (int i = 0; i < n_active_edges; i++) {
    const real current_distance = distance(active_edges.segments[i], point);
    dist_min                    = current_distance < dist_min ? current_distance : dist_min;
  }
  return dist_min - capsule.radius;
}

inline blast_fn real distance(const Capsule& capsule, const Sphere& sphere) {
  Segment segment;
  segment.p1       = capsule.p1;
  segment.p2       = capsule.p2;
  real dist_seg_pt = distance(segment, sphere.center);
  return dist_seg_pt - capsule.radius - sphere.radius;
}

inline blast_fn real distance(const Capsule& capsule1, const Capsule& capsule2) {
  Segment segment1{capsule1.p1, capsule1.p2};
  Segment segment2{capsule2.p1, capsule2.p2};
  real    dist_seg_seg = distance(segment1, segment2);
  return dist_seg_seg - capsule1.r - capsule2.r;
}

inline blast_fn Array test_collisions(const ObjMatrix<Capsule>& robot_capsules, const World* world, u32 n_lowest_distances, real start_time, real end_time) {
  Array dist_min(n_lowest_distances, INF_REAL);
  real  dist;

  u32  n_pts     = (u32) (robot_capsules).cols;
  real time_step = n_pts == 1 ? 0 : (end_time - start_time) / (n_pts - 1); // if only one point, p stays 0 and time_step has no meaning anyways

  {
    blast_time_block("External Collision Constraints : Static");
    for (int p = 0; p < n_pts; p++) {
      for (int c = 0; c < robot_capsules.rows; c++) {
        // --- Static tests
        for (const auto& box: world->boxes) {
          dist = distance(robot_capsules(c, p), box);
          for (u32 j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (u32 k = n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }

        for (auto caps: world->capsules) {
          dist = distance(robot_capsules(c, p), caps);
          for (u32 j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (u32 k = n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }

        for (auto sphere: world->spheres) {
          dist = distance(robot_capsules(c, p), sphere);
          for (u32 j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (u32 k = n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }
    }
  }

  {
    blast_time_block("External Collision Constraints : Dynamic");
    for (int p = 0; p < n_pts; p++) {
      // --- Dynamic tests ---
      for (const auto& box: world->dynamic_boxes) {
        auto current_box = box.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_box);
          for (u32 j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (u32 k = n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }

      for (const auto& caps: world->dynamic_capsules) {
        auto current_caps = caps.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_caps);
          for (u32 j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (u32 k = n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }

      for (const auto& sphere: world->dynamic_spheres) {
        auto current_sphere = sphere.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_sphere);
          for (int j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (int k = (int) n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }

      for (auto door: world->dynamic_doors) {
        auto current_box = door.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_box);
          for (int j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (int k = (int) n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }
    }
  }

  return dist_min;
}

// ---------------- Acceleration functions -----------------------------------------

struct CapsuleData {
  Vec3 direction;
  Vec3 point;
  u32  point_idx = 0;
};

inline blast_fn real distance_per_point(const Capsule& capsule, CapsuleData& capsule_data, const Box& box) {
  // Transferring points to OBB frame
  Mat3 R_trans = transpose(box.rotation);

  Vec3 p1 = R_trans * (capsule.p1 - box.center);
  Vec3 p2 = R_trans * (capsule.p2 - box.center);

  Vec3 vec   = p2 - p1;
  Vec3 point = vec.z < 0 ? p1 : p2; // take point of highest z value
  vec        = vec.z < 0 ? -vec : vec;

  // A
  // real x_min = - box.extents.x;
  // real x_max =   box.extents.x;
  // real y_min = - box.extents.y;
  // real y_max =   box.extents.y;
  // real z_min = - box.extents.z;
  // real z_max =   box.extents.z;

  // Creating the eight original vertices
  OriginVertices origin_vertices;
  origin_vertices.vert[0] = {-box.extents.x, -box.extents.y, -box.extents.z};
  origin_vertices.vert[1] = {-box.extents.x, box.extents.y, -box.extents.z};
  origin_vertices.vert[2] = {box.extents.x, -box.extents.y, -box.extents.z};
  origin_vertices.vert[3] = {box.extents.x, box.extents.y, -box.extents.z};
  origin_vertices.vert[4] = {-box.extents.x, -box.extents.y, box.extents.z};
  origin_vertices.vert[5] = {-box.extents.x, box.extents.y, box.extents.z};
  origin_vertices.vert[6] = {box.extents.x, -box.extents.y, box.extents.z};
  origin_vertices.vert[7] = {box.extents.x, box.extents.y, box.extents.z};

  // Thus, we get the three main directions
  Directions directions;
  directions.dir[0] = {2 * box.extents.x, 0, 0};
  directions.dir[1] = {0, 2 * box.extents.y, 0};
  directions.dir[2] = {0, 0, 2 * box.extents.z};

  // some faces depend only on the direction of vec in x
  Faces face;
  // the top four vertices will always be extended, meaning that the top and bottom faces (faces 0-1) will always be of the same format:
  face.surface[0].p  = origin_vertices.vert[4] + vec;
  face.surface[0].d1 = directions.dir[0];
  face.surface[0].d2 = directions.dir[1];

  face.surface[1].p  = origin_vertices.vert[0];
  face.surface[1].d1 = directions.dir[1];
  face.surface[1].d2 = directions.dir[0];

  // The other faces depend on the orientation of the vector either in x (faces 2-5),
  // in y (faces 6-9), or a combination of both (faces 10-11)
  bool positive_x = vec.x >= 0;
  bool positive_y = vec.y >= 0;
  // bool positive_z = vec.z >= 0;

  /*face.surface[2].p  = positive_x ? origin_vertices.vert[5] : origin_vertices.vert[6];
  face.surface[2].d1 = positive_x ? -directions.dir[1] : directions.dir[1];
  face.surface[2].d2 = vec;

  face.surface[3].p  = positive_x ? origin_vertices.vert[1] : origin_vertices.vert[2];
  face.surface[3].d1 = positive_x ? -directions.dir[1] : directions.dir[1];
  face.surface[3].d2 = directions.dir[2];

  face.surface[4].p  = positive_x ? origin_vertices.vert[2] + vec : origin_vertices.vert[1] + vec;
  face.surface[4].d1 = positive_x ? directions.dir[1] : -directions.dir[1];
  face.surface[4].d2 = directions.dir[2];

  face.surface[5].p  = positive_x ? origin_vertices.vert[2] : origin_vertices.vert[1];
  face.surface[5].d1 = positive_x ? directions.dir[1] : -directions.dir[1];
  face.surface[5].d2 = vec;

  face.surface[6].p  = positive_y ? origin_vertices.vert[4] : origin_vertices.vert[7];
  face.surface[6].d1 = positive_y ? directions.dir[0] : -directions.dir[0];
  face.surface[6].d2 = vec;

  face.surface[7].p  = positive_y ? origin_vertices.vert[0] : origin_vertices.vert[3];
  face.surface[7].d1 = positive_y ? directions.dir[0] : -directions.dir[0];
  face.surface[7].d2 = directions.dir[2];

  face.surface[8].p  = positive_y ? origin_vertices.vert[3] + vec : origin_vertices.vert[0] + vec;
  face.surface[8].d1 = positive_y ? -directions.dir[0] : directions.dir[0];
  face.surface[8].d2 = directions.dir[2];

  face.surface[9].p  = positive_y ? origin_vertices.vert[3] : origin_vertices.vert[0];
  face.surface[9].d1 = positive_y ? -directions.dir[0] : directions.dir[0];
  face.surface[9].d2 = vec;

  face.surface[10].p  = positive_x ? (positive_y ? origin_vertices.vert[1] : origin_vertices.vert[4]) : (positive_y ? origin_vertices.vert[7] : origin_vertices.vert[2]);
  face.surface[10].d1 = positive_x ? (positive_y ? directions.dir[2] : -directions.dir[2]) : (positive_y ? -directions.dir[2] : directions.dir[2]);
  face.surface[10].d2 = vec;

  face.surface[11].p  = positive_x ? (positive_y ? origin_vertices.vert[6] : origin_vertices.vert[3]) : (positive_y ? origin_vertices.vert[0] : origin_vertices.vert[5]);
  face.surface[11].d1 = positive_x ? (positive_y ? -directions.dir[2] : directions.dir[2]) : (positive_y ? directions.dir[2] : -directions.dir[2]);
  face.surface[11].d2 = vec;
  */

  // Faces 2-5 depend on the x coordinate
  if (positive_x) {
    face.surface[2].p  = origin_vertices.vert[5];
    face.surface[2].d1 = -directions.dir[1];
    face.surface[2].d2 = vec;

    face.surface[3].p  = origin_vertices.vert[1];
    face.surface[3].d1 = -directions.dir[1];
    face.surface[3].d2 = directions.dir[2];

    face.surface[4].p  = origin_vertices.vert[2] + vec;
    face.surface[4].d1 = directions.dir[1];
    face.surface[4].d2 = directions.dir[2];

    face.surface[5].p  = origin_vertices.vert[2];
    face.surface[5].d1 = directions.dir[1];
    face.surface[5].d2 = vec;
  } else {
    face.surface[2].p  = origin_vertices.vert[6];
    face.surface[2].d1 = directions.dir[1];
    face.surface[2].d2 = vec;

    face.surface[3].p  = origin_vertices.vert[2];
    face.surface[3].d1 = directions.dir[1];
    face.surface[3].d2 = directions.dir[2];

    face.surface[4].p  = origin_vertices.vert[1] + vec;
    face.surface[4].d1 = -directions.dir[1];
    face.surface[4].d2 = directions.dir[2];

    face.surface[5].p  = origin_vertices.vert[1];
    face.surface[5].d1 = -directions.dir[1];
    face.surface[5].d2 = vec;
  }

  // Faces 6-9 depend on the y coordinate
  if (positive_y) {
    face.surface[6].p  = origin_vertices.vert[4];
    face.surface[6].d1 = directions.dir[0];
    face.surface[6].d2 = vec;

    face.surface[7].p  = origin_vertices.vert[0];
    face.surface[7].d1 = directions.dir[0];
    face.surface[7].d2 = directions.dir[2];

    face.surface[8].p  = origin_vertices.vert[3] + vec;
    face.surface[8].d1 = -directions.dir[0];
    face.surface[8].d2 = directions.dir[2];

    face.surface[9].p  = origin_vertices.vert[3];
    face.surface[9].d1 = -directions.dir[0];
    face.surface[9].d2 = vec;
  } else {
    face.surface[6].p  = origin_vertices.vert[7];
    face.surface[6].d1 = -directions.dir[0];
    face.surface[6].d2 = vec;

    face.surface[7].p  = origin_vertices.vert[3];
    face.surface[7].d1 = -directions.dir[0];
    face.surface[7].d2 = directions.dir[2];

    face.surface[8].p  = origin_vertices.vert[0] + vec;
    face.surface[8].d1 = directions.dir[0];
    face.surface[8].d2 = directions.dir[2];

    face.surface[9].p  = origin_vertices.vert[0];
    face.surface[9].d1 = directions.dir[0];
    face.surface[9].d2 = vec;
  }

  // Faces 10-11 depend on the x and y coordinate
  if (positive_x && positive_y) {
    face.surface[10].p  = origin_vertices.vert[1];
    face.surface[10].d1 = directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[6];
    face.surface[11].d1 = -directions.dir[2];
    face.surface[11].d2 = vec;
  } else if (positive_x && !positive_y) {
    face.surface[10].p  = origin_vertices.vert[4];
    face.surface[10].d1 = -directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[3];
    face.surface[11].d1 = directions.dir[2];
    face.surface[11].d2 = vec;
  } else if (!positive_x && positive_y) {
    face.surface[10].p  = origin_vertices.vert[7];
    face.surface[10].d1 = -directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[0];
    face.surface[11].d1 = directions.dir[2];
    face.surface[11].d2 = vec;
  } else {
    face.surface[10].p  = origin_vertices.vert[2];
    face.surface[10].d1 = directions.dir[2];
    face.surface[10].d2 = vec;

    face.surface[11].p  = origin_vertices.vert[5];
    face.surface[11].d1 = -directions.dir[2];
    face.surface[11].d2 = vec;
  }

  ActiveEdges active_edges;
  int         n_active_edges = 0;
  real        normal_distance[12];
  real        max_normal_dist = -INF_REAL;

  SegFace seg_face;

  Vec3 n_max_normal_dist;

  for (int i = 0; i < 12; i++) {
    Vec3 n_current     = cross(face.surface[i].d1, face.surface[i].d2);
    n_current          = 1 / norm(n_current) * n_current;
    normal_distance[i] = dot(n_current, point - face.surface[i].p);
    n_max_normal_dist  = normal_distance[i] > max_normal_dist ? n_current : n_max_normal_dist;
    max_normal_dist    = normal_distance[i] > max_normal_dist ? normal_distance[i] : max_normal_dist;

    if (normal_distance[i] >= 0) {
      if (point_in_surface(face.surface[i].d1, face.surface[i].d2, face.surface[i].p, point)) {
        capsule_data.direction = n_current;
        return normal_distance[i] - capsule.radius;
      }

      seg_face.segments[0].p1 = face.surface[i].p;
      seg_face.segments[0].p2 = face.surface[i].p + face.surface[i].d1;
      seg_face.segments[1].p1 = face.surface[i].p;
      seg_face.segments[1].p2 = face.surface[i].p + face.surface[i].d2;
      seg_face.segments[2].p1 = face.surface[i].p + face.surface[i].d1;
      seg_face.segments[2].p2 = face.surface[i].p + face.surface[i].d1 + face.surface[i].d2;
      seg_face.segments[3].p1 = face.surface[i].p + face.surface[i].d2;
      seg_face.segments[3].p2 = face.surface[i].p + face.surface[i].d1 + face.surface[i].d2;

      // active_faces[n_active_faces++] = face[i];
      for (const auto& current_seg: seg_face.segments) {
        bool is_in_list = false;
        for (int k = 0; k < n_active_edges; k++) {
          if (current_seg.p1 == active_edges.segments[k].p1 && current_seg.p2 == active_edges.segments[k].p2 ||
              current_seg.p1 == active_edges.segments[k].p2 && current_seg.p2 == active_edges.segments[k].p1) {
            is_in_list = true;
            break;
          }
        }
        if (!is_in_list)
          active_edges.segments[n_active_edges++] = current_seg;
      }
    }
  }

  // If point is inside
  if (n_active_edges == 0) {
    capsule_data.direction = -n_max_normal_dist;
    return max_normal_dist - capsule.radius;
  }

  real dist_min = INF_REAL;
  Vec3 current_direction;
  for (int i = 0; i < n_active_edges; i++) {

    Vec3 d                = closest_point(active_edges.segments[i], point);
    real current_distance = dot(d - point, d - point);
    current_direction     = current_distance < dist_min ? current_direction : d - point; // todo: check point - d ?
    dist_min              = current_distance < dist_min ? current_distance : dist_min;
  }
  capsule_data.direction = current_direction;
  return sqrt(dist_min) - capsule.radius;
}

inline blast_fn real test_collisions_per_point(const std::array<Capsule, MAX_CAPSULES>& robot_capsules, const World* world) { // todo: FIX THIS... when MAX_CAPSULE is higher than number of capsules, this returns 0 breaking the optimization
  real dist_min = INF_REAL;

  for (const auto& capsule: robot_capsules) {
    // --- Static tests
    for (const auto& box: world->boxes) {
      if (const auto dist = distance(capsule, box); dist < dist_min)
        dist_min = dist;
    }

    for (auto caps: world->capsules) {
      if (const auto dist = distance(capsule, caps); dist < dist_min)
        dist_min = dist;
    }

    for (auto sphere: world->spheres) {
      if (const auto dist = distance(capsule, sphere); dist < dist_min)
        dist_min = dist;
    }
  }

  /*{ todo: no dynamic collisions yet
    blast_time_block("External Collision Constraints : Dynamic");
    for (int p = 0; p < n_pts; p++) {
      // --- Dynamic tests ---
      for (const auto& box: world->dynamic_boxes) {
        auto current_box = box.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_box);
          for (u32 j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (u32 k = n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }

      for (const auto& caps: world->dynamic_capsules) {
        auto current_caps = caps.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_caps);
          for (u32 j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (u32 k = n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }

      for (const auto& sphere: world->dynamic_spheres) {
        auto current_sphere = sphere.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_sphere);
          for (int j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (int k = (int) n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }

      for (auto door: world->dynamic_doors) {
        auto current_box = door.lookup(start_time + time_step * p);
        for (int c = 0; c < robot_capsules.rows; c++) {
          dist = distance(robot_capsules(c, p), current_box);
          for (int j = 0; j < n_lowest_distances; j++) {
            if (dist < dist_min[j]) {
              for (int k = (int) n_lowest_distances - 1; k > j; k--) {
                dist_min[k] = dist_min[k - 1];
              }
              dist_min[j] = dist;
              break;
            }
          }
        }
      }
    }
  }*/

  return dist_min;
}

} // namespace blast
