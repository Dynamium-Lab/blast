#pragma once
#include <cmath>
#include <algorithm>

#include "blast"

namespace blast {

struct Segment {
    Vec3 p1;
    Vec3 p2;
};

struct surface {
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
    real r;
};

struct Triangle {
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;
    real distance;
};

const real COLLISION_EPSILON = (real)1e-9;

// --- Basic functions ---

inline bool point_in_triangle(Vec3 triangle_v1, Vec3 triangle_v2, Vec3 triangle_origin, Vec3 point) {
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
    real beta = dot(cross(w, v), n) / dot(n, n);
    real alpha = 1 - gamma - beta;
    // The point P′ lies inside T if:
    return ((0 <= alpha) && (alpha <= 1) && (0 <= beta) && (beta <= 1) && (0 <= gamma) && (gamma <= 1));
}

inline bool point_in_surface(Vec3 surf_v1, Vec3 surf_v2, Vec3 surf_origin, Vec3 point) {
    bool tri1 = point_in_triangle(surf_v1, surf_v2, surf_origin, point);
    bool tri2 = point_in_triangle(-surf_v1, -surf_v2, surf_origin+surf_v1+surf_v2, point);
    return (tri1 || tri2);
}

inline Vec3 closest_point(Segment segment, Vec3 point) {
    Vec3 ab = segment.p2 - segment.p1;
    real t = dot(point - segment.p1, ab) / dot(ab, ab);
    t = clamp(t, 0, 1);
    Vec3 d = segment.p1 + t * ab;
    return d;
}

inline real distance(Segment segment, Vec3 point) {
    Vec3 d = closest_point(segment, point);
    return norm(d - point);
}

inline Vec3 closest_point_origin(Segment segment) {
    Vec3 ab = segment.p2 - segment.p1;
    real t = dot(- segment.p1, ab) / dot(ab, ab);

    t = clamp(t, 0, 1);

    Vec3 d = segment.p1 + t * ab;
    return d;
}

inline Vec3 closest_point_plane(Vec3 q, Plane p) {
    return q - (dot(p.n, q - p.p) / dot(p.n, p.n)) * p.n;;
}

inline real distance(surface surface, Vec3 point) {
    Vec3 n = cross(surface.d1, surface.d2);
    if (dot(n, n) < COLLISION_EPSILON)
        return -INF_REAL; // the plane is non-existant
    Vec3 n_unit = 1/norm(n)*n;

    // If the surface is a rectangular shape, it is much easier to treat.
    if (std::abs(dot(surface.d1, surface.d2)) <= COLLISION_EPSILON) {
        // norm^2 d1 & d2
        real val_d1_sq = dot(surface.d1, surface.d1);
        real val_d2_sq = dot(surface.d2, surface.d2);
        // direction vector from p to point
        Vec3 val_direction = point - surface.p;

        // How far along is the point in d1 & d2 direction
        real val_t1 = dot(val_direction, surface.d1) / val_d1_sq;
        real val_t2 = dot(val_direction, surface.d2) / val_d2_sq;
        real val_t1_clamped = clamp(val_t1, 0, 1);
        real val_t2_clamped = clamp(val_t2, 0, 1);
        real val_normaldist = dot(val_direction, n_unit);

        real val_dist1 = (val_t1 - val_t1_clamped)*(val_t1 - val_t1_clamped)*val_d1_sq;
        real val_dist2 = (val_t2 - val_t2_clamped)*(val_t2 - val_t2_clamped)*val_d2_sq;

        real val_result = (real)sqrt(val_dist1 + val_dist2 + val_normaldist*val_normaldist);

        real val_distance = (val_t1 >= 0 && val_t1 <= 1 && val_t2 >= 0 && val_t2 <= 1) ? val_normaldist : val_result;
        return val_distance == val_normaldist ? val_normaldist : (val_normaldist < 0 ? -val_distance : val_distance);
    }
    else {
        real normaldist = dot(point - surface.p, n_unit);
        Vec3 direction = point - normaldist*n_unit - surface.p;

        bool is_inside = point_in_surface(surface.d1, surface.d2, surface.p, point);
        if (is_inside) {
            return normaldist;
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
        int idx;
        for (u32 i = 0; i < 4; i++) {
            auto dist_tmp = distance(segment[i], point);
            dist_min = dist_tmp < dist_min ? dist_tmp : dist_min;
            idx = dist_tmp < dist_min ? i : idx;
        }
        real resultant = normaldist < 0 ? -dist_min : dist_min;
        return resultant;
    }
}

inline TwoPts closest_points(Segment segment1, Segment segment2) {
    TwoPts result;
    real s;
    real t;
    Vec3 c1;
    Vec3 c2;

    Vec3 d1 = segment1.p2 - segment1.p1; // Direction vector of segment S1
    Vec3 d2 = segment2.p2 - segment2.p1; // Direction vector of segment S2
    Vec3 r = segment1.p1 - segment2.p1;
    real a = dot(d1, d1); // Squared length of segment S1, always nonnegative
    real e = dot(d2, d2); // Squared length of segment S2, always nonnegative
    real f = dot(d2, r);
    // Check if either or both segments degenerate into points
    if (a <= COLLISION_EPSILON && e <= COLLISION_EPSILON) {
        // Both segments degenerate into points
        s = t = 0.0f;
        c1 = segment1.p2;
        c2 = segment2.p2;
        result.p1 = c1;
        result.p2 = c2;
        return result;
    }
    if (a <= COLLISION_EPSILON) {
        // First segment degenerates into a point
        s = 0.0f;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = clamp(t, 0.0f, 1.0f);
    }
    else {
        real c = dot(d1, r);
        if (e <= COLLISION_EPSILON) {
            // Second segment degenerates into a point
            t = 0.0f;
            s = clamp(-c / a, 0.0f, 1.0f); // t = 0 => s = (b*t - c) / a = -c / a
        }
        else {
            // The general nondegenerate case starts here
            real b = dot(d1, d2);
            real denom = a * e - b * b; // Always nonnegative
            // If segments not parallel, compute closest point on L1 to L2 and
            // clamp to segment S1. Else pick arbitrary s (here 0)
            if (denom != 0.0f) {
                s = clamp((b * f - c * e) / denom, 0.0f, 1.0f);
            }
            else s = 0.0f;
            // Compute point on L2 closest to S1(s) using
            // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
            t = (b * s + f) / e;
            // If t in [0,1] done. Else clamp t, recompute s for the new value
            // of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a
            // and clamp s to [0, 1]
            if (t < 0.0f) {
                t = 0.0f;
                s = clamp(-c / a, 0.0f, 1.0f);
            }
            else if (t > 1.0f) {
                t = 1.0f;
                s = clamp((b - c) / a, 0.0f, 1.0f);
            }
        }
    }
    c1 = segment1.p1 + d1 * s;
    c2 = segment2.p1 + d2 * t;
    result.p1 = c1;
    result.p2 = c2;
    return result;
}

inline real distance(Segment segment1, Segment segment2) {
    TwoPts close_pt = closest_points(segment1, segment2);
    return norm(close_pt.p1 - close_pt.p2);
}

// --- Primitive tests ---

inline real distance(const Capsule& capsule, const Box& box) {
    // Tranferring points to OBB frame
    Mat3 Rtrans = transpose(box.R);

    Vec3 p1 = Rtrans * (capsule.p1 - box.c);
    Vec3 p2 = Rtrans * (capsule.p2 - box.c);

    Vec3 vec = p2 - p1;
    Vec3 point = vec.z < 0 ? p1 : p2; // take point of highest z value
    vec = vec.z < 0 ? -vec : vec;

    // A
    // real xmin = - box.e.x;
    // real xmax =   box.e.x;
    // real ymin = - box.e.y;
    // real ymax =   box.e.y;
    // real zmin = - box.e.z;
    // real zmax =   box.e.z;

    // Creating the eight original vertices
    Vec3 orgvert[8];
    orgvert[0] = { - box.e.x, - box.e.y, - box.e.z };
    orgvert[1] = { - box.e.x,   box.e.y, - box.e.z };
    orgvert[2] = {   box.e.x, - box.e.y, - box.e.z };
    orgvert[3] = {   box.e.x,   box.e.y, - box.e.z };
    orgvert[4] = { - box.e.x, - box.e.y,   box.e.z };
    orgvert[5] = { - box.e.x,   box.e.y,   box.e.z };
    orgvert[6] = {   box.e.x, - box.e.y,   box.e.z };
    orgvert[7] = {   box.e.x,   box.e.y,   box.e.z };

    // Thus, we get the three main directions
    Vec3 dir[3];
    dir[0] = {2*box.e.x, 0, 0};
    dir[1] = {0, 2*box.e.y, 0};
    dir[2] = {0, 0, 2*box.e.z};

    // some faces depend only on the direction of vec in x
    surface face[12];
    // the top four vertices will always be extended, meaning that the top and bottom faces (faces 0-1) will always be of the same format:
    face[0].p = orgvert[4] + vec;
    face[0].d1 = dir[0];
    face[0].d2 = dir[1];

    face[1].p = orgvert[0];
    face[1].d1 = dir[1];
    face[1].d2 = dir[0];

    // The other faces depend on the orientation of the vector either in x (faces 2-5),
    // in y (faces 6-9), or a combination of both (faces 10-11)
    bool vecx = vec.x >= 0;
    bool vecy = vec.y >= 0;
    bool vecz = vec.z >= 0;

    // Faces 2-5 depend on the x coordinate
    face[2].p = vecx ? orgvert[5] : orgvert[6];
    face[2].d1 = vecx ? - dir[1] : dir[1];
    face[2].d2 = vec;

    face[3].p = vecx ? orgvert[1] : orgvert[2];
    face[3].d1 = vecx ? - dir[1] : dir[1];
    face[3].d2 = dir[2];

    face[4].p = vecx ? orgvert[2] + vec : orgvert[1] + vec;
    face[4].d1 = vecx ? dir[1] : -dir[1];
    face[4].d2 = dir[2];

    face[5].p = vecx ? orgvert[2] : orgvert[1];
    face[5].d1 = vecx ? dir[1] : -dir[1];
    face[5].d2 = vec;

    // Faces 6-9 depend on the y coordinate
    face[6].p = vecy ? orgvert[4] : orgvert[7];
    face[6].d1 = vecy ? dir[0] : -dir[0];
    face[6].d2 = vec;

    face[7].p = vecy ? orgvert[0] : orgvert[3];
    face[7].d1 = vecy ? dir[0] : -dir[0];
    face[7].d2 = dir[2];

    face[8].p = vecy ? orgvert[3] + vec : orgvert[0] + vec;
    face[8].d1 = vecy ? -dir[0] : dir[0];
    face[8].d2 = dir[2];

    face[9].p = vecy ? orgvert[3] : orgvert[0];
    face[9].d1 = vecy ? -dir[0] : dir[0];
    face[9].d2 = vec;

    // Faces 10-11 depend on the x and y coordinate
    face[10].p = vecx ? (vecy ? orgvert[1] : orgvert[4]) : (vecy ? orgvert[7] : orgvert[2]);
    face[10].d1 = vecx ? (vecy ? dir[2] : -dir[2]) : (vecy ? -dir[2] : dir[2]);
    face[10].d2 = vec;

    face[11].p = vecx ? (vecy ? orgvert[6] : orgvert[3]) : (vecy ? orgvert[0] : orgvert[5]);
    face[11].d1 = vecx ? (vecy ? -dir[2] : dir[2]) : (vecy ? dir[2] : -dir[2]);
    face[11].d2 = vec;

    Segment active_edges[24];
    int n_active_edges = 0;
    Vec3 n_current;
    real normaldist[12];
    real max_normal_dist = -INF_REAL;
    for (int i = 0; i < 12; i++) {
        n_current = cross(face[i].d1, face[i].d2);
        n_current = 1 / norm(n_current) * n_current;
        normaldist[i] = dot(n_current, point - face[i].p);
        max_normal_dist = normaldist[i] > max_normal_dist ? normaldist[i] : max_normal_dist;

        if (normaldist[i] >= 0) {
            bool is_inside = point_in_surface(face[i].d1, face[i].d2, face[i].p, point);
            if (is_inside)
                return normaldist[i] - capsule.r;

            Segment seg_face[4];
            seg_face[0].p1 = face[i].p;
            seg_face[0].p2 = face[i].p + face[i].d1;
            seg_face[1].p1 = face[i].p;
            seg_face[1].p2 = face[i].p + face[i].d2;
            seg_face[2].p1 = face[i].p + face[i].d1;
            seg_face[2].p2 = face[i].p + face[i].d1 + face[i].d2;
            seg_face[3].p1 = face[i].p + face[i].d2;
            seg_face[3].p2 = face[i].p + face[i].d1 + face[i].d2;

            // active_faces[n_active_faces++] = face[i];
            for (int j = 0; j < 4; j++) {
                bool is_in_list = false;
                for (int k = 0; k < n_active_edges; k++) {
                    if (seg_face[j].p1 == active_edges[k].p1 && seg_face[j].p2 == active_edges[k].p2 ||
                            seg_face[j].p1 == active_edges[k].p2 && seg_face[j].p2 == active_edges[k].p1) {
                        is_in_list = true;
                        break;
                    }
                }
                if (!is_in_list)
                    active_edges[n_active_edges++] = seg_face[j];
            }
        }
    }

    // If point is inside
    if (n_active_edges == 0) {
        return max_normal_dist - capsule.r;
    }

    real dist_min = INF_REAL;
    real distcurrent;
    for (int i = 0; i < n_active_edges; i++) {
        distcurrent = distance(active_edges[i], point);
        dist_min = distcurrent < dist_min ? distcurrent : dist_min;
    }
    return dist_min - capsule.r;
}

inline real distance(const Capsule& capsule, const Sphere& sphere) {
    Segment segment;
    segment.p1 = capsule.p1;
    segment.p2 = capsule.p2;
    real dist_seg_pt = distance(segment, sphere.c);
    return dist_seg_pt - capsule.r - sphere.r;
}

inline real distance(const Capsule& capsule1, const Capsule& capsule2) {
    Segment segment1 {capsule1.p1, capsule1.p2};
    Segment segment2 {capsule2.p1, capsule2.p2};
    real dist_seg_seg = distance(segment1, segment2);
    return dist_seg_seg - capsule1.r - capsule2.r;
}

inline Array test_collisions(const ObjMatrix<Capsule>& robot_capsules, const World* world, const u32 n_lowest_distances, const real start_time, const real end_time) {
    Array dist_min(n_lowest_distances, INF_REAL);
    real dist;

    u32 n_pts = (u32)(robot_capsules).cols;
    real time_step = n_pts == 1 ? 1 : (end_time - start_time)/(n_pts-1);

    {
        blast_time_block("External Collision Constraints : Static");
        for (u32 p = 0; p < n_pts; p++) {
            for (u32 c = 0; c < robot_capsules.rows; c++) {
                // --- Static tests
                for (auto box : world->boxes) {
                    dist = distance(robot_capsules(c, p), box);
                    for (u32 j = 0; j < n_lowest_distances; j++) {
                        if (dist < dist_min[j]) {
                            for (u32 k = n_lowest_distances - 1; k > j; k--) {
                                dist_min[k] = dist_min[k-1];
                            }
                            dist_min[j] = dist;
                            break;
                        }
                    }
                }

                for (auto caps : world->capsules) {
                    dist = distance(robot_capsules(c, p), caps);
                    for (u32 j = 0; j < n_lowest_distances; j++) {
                        if (dist < dist_min[j]) {
                            for (u32 k = n_lowest_distances - 1; k > j; k--) {
                                dist_min[k] = dist_min[k-1];
                            }
                            dist_min[j] = dist;
                            break;
                        }
                    }
                }

                for (auto sphere : world->spheres) {
                    dist = distance(robot_capsules(c, p), sphere);
                    for (u32 j = 0; j < n_lowest_distances; j++) {
                        if (dist < dist_min[j]) {
                            for (u32 k = n_lowest_distances - 1; k > j; k--) {
                                dist_min[k] = dist_min[k-1];
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
        for (u32 p = 0; p < n_pts; p++) {
            // --- Dynamic tests ---
            for (auto box : world->dynamic_boxes) {
                auto current_box = box.lookup(start_time + time_step*p);
                for (u32 c = 0; c < robot_capsules.rows; c++) {
                    dist = distance(robot_capsules(c, p), current_box);
                    for (u32 j = 0; j < n_lowest_distances; j++) {
                        if (dist < dist_min[j]) {
                            for (u32 k = n_lowest_distances - 1; k > j; k--) {
                                dist_min[k] = dist_min[k-1];
                            }
                            dist_min[j] = dist;
                            break;
                        }
                    }
                }
            }

            for (auto caps : world->dynamic_capsules) {
                auto current_caps = caps.lookup(start_time + time_step*p);
                for (u32 c = 0; c < robot_capsules.rows; c++) {
                    dist = distance(robot_capsules(c, p), current_caps);
                    for (u32 j = 0; j < n_lowest_distances; j++) {
                        if (dist < dist_min[j]) {
                            for (u32 k = n_lowest_distances - 1; k > j; k--) {
                                dist_min[k] = dist_min[k-1];
                            }
                            dist_min[j] = dist;
                            break;
                        }
                    }
                }
            }

            for (auto sphere : world->dynamic_spheres) {
                auto current_sphere = sphere.lookup(start_time + time_step*p);
                for (u32 c = 0; c < robot_capsules.rows; c++) {
                    dist = distance(robot_capsules(c, p), current_sphere);
                    for (u32 j = 0; j < n_lowest_distances; j++) {
                        if (dist < dist_min[j]) {
                            for (u32 k = n_lowest_distances - 1; k > j; k--) {
                                dist_min[k] = dist_min[k-1];
                            }
                            dist_min[j] = dist;
                            break;
                        }
                    }
                }
            }

            for (auto door : world->dynamic_doors) {
                auto current_box = door.get_object_at_time(start_time + time_step*p);
                for (u32 c = 0; c < robot_capsules.rows; c++) {
                    dist = distance(robot_capsules(c, p), current_box);
                    for (u32 j = 0; j < n_lowest_distances; j++) {
                        if (dist < dist_min[j]) {
                            for (u32 k = n_lowest_distances - 1; k > j; k--) {
                                dist_min[k] = dist_min[k-1];
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

} // namespace blast