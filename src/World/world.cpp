#include "blast.h"
#include <cmath>
#include <algorithm>
#include "blast_error.h"

using std::vector;

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

// ======================================
//            Basic functions
// ======================================

bool point_in_triangle(Vec3 triangle_v1, Vec3 triangle_v2, Vec3 triangle_origin, Vec3 point) {
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

bool point_in_surface(Vec3 surf_v1, Vec3 surf_v2, Vec3 surf_origin, Vec3 point) {
    bool tri1 = point_in_triangle(surf_v1, surf_v2, surf_origin, point);
    bool tri2 = point_in_triangle(-surf_v1, -surf_v2, surf_origin+surf_v1+surf_v2, point);
    return (tri1 || tri2);
}

Vec3 closest_point(Segment segment, Vec3 point) {
    Vec3 ab = segment.p2 - segment.p1;
    real t = dot(point - segment.p1, ab) / dot(ab, ab);
    t = clamp(t, 0, 1);
    Vec3 d = segment.p1 + t * ab;
    return d;
}

real distance(Segment segment, Vec3 point) {
    Vec3 d = closest_point(segment, point);
    return norm(d - point);
}

Vec3 closest_point_origin(Segment segment) {
    Vec3 ab = segment.p2 - segment.p1;
    real t = dot(- segment.p1, ab) / dot(ab, ab);

    t = clamp(t, 0, 1);

    Vec3 d = segment.p1 + t * ab;
    return d;
}

Vec3 closest_point_plane(Vec3 q, Plane p) {
    return q - (dot(p.n, q - p.p) / dot(p.n, p.n)) * p.n;;
}

real distance(surface surface, Vec3 point) {
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

// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t.
// From : Real-Time Detection Collision (Christer Ericson)
TwoPts closest_points(Segment segment1, Segment segment2) {
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

real distance(Segment segment1, Segment segment2) {
    TwoPts close_pt = closest_points(segment1, segment2);
    return norm(close_pt.p1 - close_pt.p2);
}

TwoPts projection(Circle circle, Segment segment) {
    Vec3 p1 = segment.p1 - circle.p;
    Vec3 p2 = segment.p2 - circle.p;

    Vec3 d1 = p1;
    Vec3 d2 = cross(circle.n, d1);

    Vec3 d1_unit = (1/norm(d1))*d1;
    Vec3 d2_unit = (1/norm(d2))*d2;

    real x1 = norm(d1);
    real y1 = 0;
    real x2 = dot(p2, d1_unit);
    real y2 = dot(p2, d2_unit);

    real dx = x2 - x1;
    real dy = y2 - y1;
    real dr_sq = dx * dx + dy * dy;
    real D = x1 * y2 - x2 * y1;

    real det = (circle.r * circle.r * dr_sq - D * D);

    if (det < 0) {
        Vec3 point = closest_point(segment, circle.p);
        Vec3 pointcirc = circle.p + (circle.r / norm(point - circle.p)) * (point - circle.p);
        return { pointcirc, pointcirc };
    }

    det = (real)sqrt(det); // note: only sqrt if det is non-negative
    real sign_dy = dy > 0 ? (real)1 : (real)-1;
    real x_1 = (1 / dr_sq) * (D * dy + sign_dy * dx * det);
    real x_2 = (1 / dr_sq) * (D * dy - sign_dy * dx * det);

    real y_1 = (1 / dr_sq) * (-D * dx + std::abs(dy) * det);
    real y_2 = (1 / dr_sq) * (-D * dx - std::abs(dy) * det);

    Vec3 p_1 = circle.p + x_1 * d1_unit + y_1 * d2_unit;
    Vec3 p_2 = circle.p + x_2 * d1_unit + y_2 * d2_unit;
    return { p_1, p_2 };
}

// ======================================
//            Primitive tests
// ======================================

real distance(Capsule capsule, Cylinder cylinder) {
    Segment segment1;
    Segment segment2;
    segment1.p1 = capsule.p1;
    segment1.p2 = capsule.p2;
    segment2.p1 = cylinder.p1;
    segment2.p2 = cylinder.p2;

    TwoPts points = closest_points(segment1, segment2);
    real cond1 = dot(points.p2 - cylinder.p1, points.p2 - cylinder.p1);
    real cond2 = dot(points.p2 - cylinder.p2, points.p2 - cylinder.p2);

    // Depending on which point of the cylinder is closest, we will change the face which we check
    if (cond1 < COLLISION_EPSILON || cond2 < COLLISION_EPSILON) {
        Vec3 cent;
        Vec3 other;

        bool corner = cond2 < cond1;
        cent = corner ? cylinder.p2 : cylinder.p1;
        other = corner ? cylinder.p1 : cylinder.p2;

        // Check if both points project on the circle plane
        Vec3 n = cent - other;
        Vec3 n_unit = (1 / norm(n)) * n;

        Plane face1;
        face1.n = n_unit;
        face1.p = cent;

        Vec3 proj1 = closest_point_plane(capsule.p1, face1);
        Vec3 proj2 = closest_point_plane(capsule.p2, face1);

        real rad_sq1 = dot(proj1 - cent, proj1 - cent);
        real rad_sq2 = dot(proj2 - cent, proj2 - cent);

        real dist_norm_sq1 = dot(capsule.p1 - proj1, capsule.p1 - proj1);
        real dist_norm_sq2 = dot(capsule.p2 - proj2, capsule.p2 - proj2);

        real dist_norm_sq_min = (dist_norm_sq2 < dist_norm_sq1) ? dist_norm_sq2 : dist_norm_sq1;

        // if both points project on the circle plane, then the distance will be the minimum of the two normal distances calculated
        if (rad_sq1 <= cylinder.r * cylinder.r && rad_sq2 <= cylinder.r * cylinder.r)
            return sqrt(dist_norm_sq_min) - capsule.r;

        // if one point projects on the circle plane, it is necessary to check the normal distance as well
        real testnormal = (rad_sq1 <= cylinder.r * cylinder.r) ? dist_norm_sq1 : (rad_sq2 <= cylinder.r * cylinder.r) ? dist_norm_sq2 : INF_REAL;

        // project the points on the plan
        Segment proj;
        proj.p1 = proj1;
        proj.p2 = proj2;
        Circle circle;
        circle.p = cent;
        circle.r = cylinder.r;
        circle.n = n_unit;

        TwoPts pts = projection(circle, proj);

        real dist1 = distance(segment1, pts.p1);
        real dist2 = distance(segment1, pts.p2);

        if (dist2 <= dist1 && dist2 * dist2 < testnormal)
            return dist2 - capsule.r;

        if (dist1 < dist2 && dist1 * dist1 < testnormal)
            return dist1 - capsule.r;

        return sqrt(testnormal) - capsule.r;
    }
    else
        return norm(points.p1 - points.p2) - capsule.r - cylinder.r;
}

real distance(Capsule capsule, Sphere sphere) {
    Segment segment;
    segment.p1 = capsule.p1;
    segment.p2 = capsule.p2;
    real dist_seg_pt = distance(segment, sphere.c);
    return dist_seg_pt - capsule.r - sphere.r;
}

real distance(Capsule capsule1, Capsule capsule2) {
    Segment segment1 {capsule1.p1, capsule1.p2};
    Segment segment2 {capsule2.p1, capsule2.p2};
    real dist_seg_seg = distance(segment1, segment2);
    return dist_seg_seg - capsule1.r - capsule2.r;
}

// Box capsule tests basic
real distance_old(Box box, Capsule capsule) {
    Segment segment;
    segment.p1 = capsule.p1;
    segment.p2 = capsule.p2;

    Mat3 Rtrans = transpose(box.R);

    Vec3 p1 = Rtrans * (segment.p1 - box.c);
    Vec3 p2 = Rtrans * (segment.p2 - box.c);

    Vec3 vec = p2 - p1;
    Vec3 point = p2;

    if (vec.z < 0) {
        vec = -vec;
        point = p1;
    }

    real xmin = - box.e[0];
    real xmax = + box.e[0];
    real ymin = - box.e[1];
    real ymax = + box.e[1];
    real zmin = - box.e[2];
    real zmax = + box.e[2];

    // Creating the eight original vertices
    Vec3 orgvert[8];

    orgvert[0] = { xmin, ymin, zmin };
    orgvert[1] = { xmin, ymax, zmin };
    orgvert[2] = { xmax, ymin, zmin };
    orgvert[3] = { xmax, ymax, zmin };
    orgvert[4] = { xmin, ymin, zmax };
    orgvert[5] = { xmin, ymax, zmax };
    orgvert[6] = { xmax, ymin, zmax };
    orgvert[7] = { xmax, ymax, zmax };

    // Thus, we get the three main directions
    Vec3 dir[3];

    dir[0] = orgvert[2] - orgvert[0];
    dir[1] = orgvert[1] - orgvert[0];
    dir[2] = orgvert[4] - orgvert[0];

    // Creating the added vertices
    Vec3 vert[8];
    for (int i = 0; i < 8; i++) {
        vert[i] = orgvert[i] + vec;
    }

    // some faces depend only on the direction of vec in x
    surface face[12];

    // the top four vertices will always be extended, meaning that the top and bottom faces (faces 0-1) will always be of the same format:
    face[0].p = vert[4];
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

    face[4].p = vecx ? vert[2] : vert[1];
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

    face[8].p = vecy ? vert[3] : vert[0];
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

    real dist[12];
    dist[0] = distance(face[0], point);

    real distcurrent = dist[0];

    for (int i = 1; i < 12; i++) {
        dist[i] = distance(face[i], point);

        bool dist_over = dist[i] >= 0;
        bool dist_under = dist[i] < 0;
        // When dist_over , distcurrent == the closest positive value to 0
        // When dist_under, distcurrent == the closest negative value to 0
        distcurrent = dist_over && ((dist[i] < distcurrent) || (distcurrent < 0)) ? dist[i] :
                      dist_under && (dist[i] > distcurrent) ? dist[i] : distcurrent;
    }
    return distcurrent - capsule.r;
}

// Box capsule tests accelerated
real distance(Capsule capsule, Box box) {
    Segment segment;
    segment.p1 = capsule.p1;
    segment.p2 = capsule.p2;

    Mat3 Rtrans = transpose(box.R);

    Vec3 p1 = Rtrans * (segment.p1 - box.c);
    Vec3 p2 = Rtrans * (segment.p2 - box.c);

    Vec3 vec = p2 - p1;
    Vec3 point = p2;

    if (vec.z < 0) {
        vec = -vec;
        point = p1;
    }

    real xmin = - box.e[0];
    real xmax = + box.e[0];
    real ymin = - box.e[1];
    real ymax = + box.e[1];
    real zmin = - box.e[2];
    real zmax = + box.e[2];

    // Creating the eight original vertices
    Vec3 orgvert[8];
    orgvert[0] = { xmin, ymin, zmin };
    orgvert[1] = { xmin, ymax, zmin };
    orgvert[2] = { xmax, ymin, zmin };
    orgvert[3] = { xmax, ymax, zmin };
    orgvert[4] = { xmin, ymin, zmax };
    orgvert[5] = { xmin, ymax, zmax };
    orgvert[6] = { xmax, ymin, zmax };
    orgvert[7] = { xmax, ymax, zmax };

    // Thus, we get the three main directions
    Vec3 dir[3];
    dir[0] = orgvert[2] - orgvert[0];
    dir[1] = orgvert[1] - orgvert[0];
    dir[2] = orgvert[4] - orgvert[0];

    // Creating the added vertices
    Vec3 vert[8];
    for (int i = 0; i < 8; i++) {
        vert[i] = orgvert[i] + vec;
    }

    // some faces depend only on the direction of vec in x
    surface face[12];
    // the top four vertices will always be extended, meaning that the top and bottom faces (faces 0-1) will always be of the same format:
    face[0].p = vert[4];
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

    face[4].p = vecx ? vert[2] : vert[1];
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

    face[8].p = vecy ? vert[3] : vert[0];
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

void add_box(Vec3 center_point, Vec3 half_width, Mat3 rotation_matrix, World* world) {
    Box new_box;
    new_box.c = center_point;
    new_box.e = half_width;
    new_box.R = rotation_matrix;
    world->boxes.push_back(new_box);
}

void add_sphere(Vec3 center_point, real radius, World* world) {
    Sphere new_sphere;
    new_sphere.c = center_point;
    new_sphere.r = radius;
    world->spheres.push_back(new_sphere);
}

void add_cylinder(Vec3 point1, Vec3 point2, real radius, World* world) {
    Cylinder new_cylinder;
    new_cylinder.p1 = point1;
    new_cylinder.p1 = point2;
    new_cylinder.r = radius;
    world->cylinders.push_back(new_cylinder);
}

void add_capsule(Vec3 point1, Vec3 point2, real radius, World* world) {
    Capsule new_capsule;
    new_capsule.p1 = point1;
    new_capsule.p1 = point2;
    new_capsule.r = radius;
    world->capsules.push_back(new_capsule);
}

Array test_collision(std::vector<Capsule>* robot_capsule_list, World* world, u32 n_lowest_distances) {
    Array dist_min(n_lowest_distances, INF_REAL);
    real dist;

    for (u32 c = 0; c < (*robot_capsule_list).size(); c++) {
        // Box collisions
        for (u32 i = 0; i < world->boxes.size(); i++) {
            dist = distance((*robot_capsule_list)[c], world->boxes[i]);
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

        // Capsule collisions
        for (u32 i = 0; i < world->capsules.size(); i++) {
            dist = distance((*robot_capsule_list)[c], world->capsules[i]);
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

        // Cylinder collisions
        for (u32 i = 0; i < world->cylinders.size(); i++) {
            dist = distance((*robot_capsule_list)[c], world->cylinders[i]);
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

        // Sphere collisions
        for (u32 i = 0; i < world->spheres.size(); i++) {
            dist = distance((*robot_capsule_list)[c], world->spheres[i]);
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

    return dist_min;
}

// Returns distance between an box and a point
real distance(const Box &box, const Vec3 &point) {
    Mat3 Rtrans = transpose(box.R);

    Vec3 point_box = Rtrans * (point - box.c);

    Vec3 proj = {clamp(point_box.x, -box.e.x, box.e.x), clamp(point_box.y, -box.e.y, box.e.y), clamp(point_box.z, -box.e.z, box.e.z)};
    Array dist_in(3);
    dist_in[0] = std::abs(point.x) - box.e.x;
    dist_in[1] = std::abs(point.y) - box.e.y;
    dist_in[2] = std::abs(point.z) - box.e.z;
    real result_if_inside = max(dist_in);
    return result_if_inside > 0 ? norm(proj - point_box) : result_if_inside;
}

// Returns distance between a capsule and a point
real distance(const Capsule &capsule, const Vec3 &point) {
    Segment segment;
    segment.p1 = capsule.p1;
    segment.p2 = capsule.p2;
    return distance(segment, point) - capsule.r;
}

// Returns distance between a sphere and a point
real distance(const Sphere &sph_test, const Vec3 &point) {
    return norm(point - sph_test.c) - sph_test.r;
}

// Gets the point by linear interpolation from caps_list according to values of x
Vec3 get_point(const Array& x, const Matrix &capsule_list) {
    real t = x[0]*(capsule_list.cols-1);
    int t_step = (int)floor(t);
    int t_step_plus1 = (t_step == (capsule_list).cols-1) ? t_step : t_step + 1;
    real s = x[1]*(capsule_list.rows/3-1);
    int s_step = (int)floor(s);
    int s_step_plus1 = (s_step == (capsule_list).rows/3-1) ? s_step : s_step + 1;

    Vec3 p1_1 = {capsule_list(3*s_step, t_step),             capsule_list(3*s_step + 1, t_step),             capsule_list(3*s_step + 2, t_step)};
    Vec3 p1_2 = {capsule_list(3*s_step, t_step_plus1),       capsule_list(3*s_step + 1, t_step_plus1),       capsule_list(3*s_step + 2, t_step_plus1)};
    Vec3 p2_1 = {capsule_list(3*s_step_plus1, t_step),       capsule_list(3*s_step_plus1 + 1, t_step),       capsule_list(3*s_step_plus1 + 2, t_step)};
    Vec3 p2_2 = {capsule_list(3*s_step_plus1, t_step_plus1), capsule_list(3*s_step_plus1 + 1, t_step_plus1), capsule_list(3*s_step_plus1 + 2, t_step_plus1)};

    Vec3 p1 = (1 - (t - t_step)) * p1_1 + (t - t_step) * p1_2;
    Vec3 p2 = (1 - (t - t_step)) * p2_1 + (t - t_step) * p2_2;

    Vec3 p = (1 - (s - s_step)) * p1 + (s - s_step) * p2;
    return p;
}

// Calls get_point and tests this point against the full world
real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const World* world) {
    Vec3 p = get_point(x, robot_cartesian_positions);

    real distmin = INF_REAL;
    real current_dist = 0;
    for (auto box: (*world).boxes) {
        current_dist = distance(box, p);
        distmin = current_dist < distmin ? current_dist : distmin;
    }
    for (auto caps: (*world).capsules) {
        current_dist = distance(caps, p);
        distmin = current_dist < distmin ? current_dist : distmin;
    }
    for (auto sph: (*world).spheres) {
        current_dist = distance(sph, p);
        distmin = current_dist < distmin ? current_dist : distmin;
    }

    return distmin;
}

// Calls get_point and tests this point against the box
real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const Box &box) {
    Vec3 p = get_point(x, robot_cartesian_positions);
    return distance(box, p);
}

// Gets the point by linear interpolation from caps_list according to values of x
Vec3 get_point_with_directions(const Array& x, const Matrix &capsule_list, Array *directions) {
    Assert((*directions).size == 6);

    real t = x[0]*(capsule_list.cols-1);
    int t_step = (int)floor(t);
    int t_step_plus1 = (t_step == (capsule_list).cols-1) ? t_step : t_step + 1;
    real s = x[1]*(capsule_list.rows/3-1);
    int s_step = (int)floor(s);
    int s_step_plus1 = (s_step == (capsule_list).rows/3-1) ? s_step : s_step + 1;

    // t = 0, s = 0
    Vec3 p1_1 = {capsule_list(3*s_step, t_step),             capsule_list(3*s_step + 1, t_step),             capsule_list(3*s_step + 2, t_step)};
    // t = 1, s = 0
    Vec3 p1_2 = {capsule_list(3*s_step, t_step_plus1),       capsule_list(3*s_step + 1, t_step_plus1),       capsule_list(3*s_step + 2, t_step_plus1)};
    // t = 0, s = 1
    Vec3 p2_1 = {capsule_list(3*s_step_plus1, t_step),       capsule_list(3*s_step_plus1 + 1, t_step),       capsule_list(3*s_step_plus1 + 2, t_step)};
    // t = 1, s = 1
    Vec3 p2_2 = {capsule_list(3*s_step_plus1, t_step_plus1), capsule_list(3*s_step_plus1 + 1, t_step_plus1), capsule_list(3*s_step_plus1 + 2, t_step_plus1)};

    Vec3 p1 = (1 - (t - t_step)) * p1_1 + (t - t_step) * p1_2;
    Vec3 p2 = (1 - (t - t_step)) * p2_1 + (t - t_step) * p2_2;

    Vec3 p = (1 - (s - s_step)) * p1 + (s - s_step) * p2;

    // Directions for gradient information
    Vec3 dy = p2 - p1;

    Vec3 p1x = (1 - (s - s_step)) * p1_1 + (s - s_step) * p2_1;
    Vec3 p2x = (1 - (s - s_step)) * p1_2 + (s - s_step) * p2_2;
    Vec3 dx = p2x - p1x;

    Vec3 dx_normalized = dx / (capsule_list.cols-1);
    Vec3 dy_normalized = dy / (capsule_list.rows/3-1);

    Array dir(6);
    dir = { dx_normalized[0], dx_normalized[1], dx_normalized[2], dy_normalized[0], dy_normalized[1], dy_normalized[2] };
    (*directions) = dir;
    
    return p;
}

// Returns distance between an box and a point and displacement vector
real distance_with_directions(const Box &box, const Vec3 &point, Vec3* d) {
    Mat3 Rtrans = transpose(box.R);

    Vec3 point_box = Rtrans * (point - box.c);

    Vec3 proj = {clamp(point_box.x, -box.e.x, box.e.x), clamp(point_box.y, -box.e.y, box.e.y), clamp(point_box.z, -box.e.z, box.e.z)};
    Array dist_in(3);
    dist_in[0] = std::abs(point.x) - box.e.x;
    dist_in[1] = std::abs(point.y) - box.e.y;
    dist_in[2] = std::abs(point.z) - box.e.z;
    real result_if_inside = max(dist_in);
    u32 idx = argmax(dist_in);
    Vec3 dir_inside;
    dir_inside[idx] = 1;

    // Condition ? outside : inside
    (*d) = result_if_inside > 0 ? (proj - point_box) : dir_inside;
    
    // Condition ? outside : inside
    return result_if_inside > 0 ? norm(proj - point_box) : result_if_inside;
}

// Returns distance between a capsule and a point and displacement vector
real distance_with_directions(const Capsule &capsule, const Vec3 &point, Vec3* d) {
    Vec3 ab = capsule.p2 - capsule.p1;
    real t = dot(point - capsule.p1, ab) / dot(ab, ab);
    t = clamp(t, 0, 1);
    Vec3 p_capsule = capsule.p1 + t * ab;

    (*d) = p_capsule - point;

    return norm((*d)) - capsule.r;
}

// Returns distance between a sphere and a point and displacement vector
real distance_with_directions(const Sphere &sph_test, const Vec3 &point, Vec3* d) {
    (*d) = sph_test.c - point;
    return norm((*d)) - sph_test.r;
}

// Calls get_point and tests this point against the full world
real obj_function_gradient(const Array& x, const Matrix &robot_cartesian_positions, const World* world, Array* gradient) {
    Array directions(6);
    Vec3 p = get_point_with_directions(x, robot_cartesian_positions, &directions);
    Vec3 d;

    real distmin = INF_REAL;
    Vec3 dirmin;
    real current_dist = 0;
    for (auto box: (*world).boxes) {
        current_dist = distance_with_directions(box, p, &d);
        dirmin = current_dist < distmin ? d : dirmin;
        distmin = current_dist < distmin ? current_dist : distmin;
    }
    for (auto caps: (*world).capsules) {
        current_dist = distance_with_directions(caps, p, &d);
        dirmin = current_dist < distmin ? d : dirmin;
        distmin = current_dist < distmin ? current_dist : distmin;
    }
    for (auto sph: (*world).spheres) {
        current_dist = distance_with_directions(sph, p, &d);
        dirmin = current_dist < distmin ? d : dirmin;
        distmin = current_dist < distmin ? current_dist : distmin;
    }

    Vec3 direction_x0 = { directions[0], directions[1], directions[2] };
    Vec3 direction_x1 = { directions[3], directions[4], directions[5] };

    (*gradient)[0] = dot(direction_x0, dirmin);
    (*gradient)[1] = dot(direction_x1, dirmin);

    return distmin;
}

// Calls get_point and tests this point against the box
real obj_function_gradient(const Array& x, const Matrix &robot_cartesian_positions, const Box &box, Array* gradient) {
    Array directions(6);
    Vec3 p = get_point_with_directions(x, robot_cartesian_positions, &directions);
    Vec3 dirmin;
    real res = distance_with_directions(box, p, &dirmin);

    Vec3 direction_x0 = { directions[0], directions[1], directions[2] };
    Vec3 direction_x1 = { directions[3], directions[4], directions[5] };

    (*gradient)[0] = dot(direction_x0, dirmin);
    (*gradient)[1] = dot(direction_x1, dirmin);
    return res;
}

// ======================================
//            GJK algorithm
// ======================================

struct Simplex {
    Vec3 P;
    Vec3 a;
    Vec3 b;
    Vec3 c;
    Vec3 d;
    int size;
    int max_size = 1;
};

struct EPAHull {
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;
    Vec3 n;
};

Vec3 gjk_get_support(std::vector<Vec3> vertices, Vec3 direction) {
    real largest_dot = dot(vertices[0], direction);
    Vec3 largest_vertex = vertices[0];
    for (auto& vertex : vertices) {
        const real current_dot = dot(vertex, direction);
        largest_vertex  = current_dot > largest_dot ? vertex : largest_vertex;
        largest_dot     = current_dot > largest_dot ? current_dot : largest_dot;
    }
    Assert(isinf(largest_dot) == false);
    return largest_vertex;
}

// From Ericson's manual
void gjk_solve_simplex2(Simplex* simplex) {
    Vec3 ab = (*simplex).b - (*simplex).a;
    real t = dot(- (*simplex).a, ab) / dot(ab, ab);

    t = clamp(t, 0, 1);

    (*simplex).P = (*simplex).a + t * ab;

    // Vec3 a_temp = (*simplex).a;
    (*simplex).a = t == 1 ? (*simplex).b : (*simplex).a;
    (*simplex).b = t == 1 ? (*simplex).a : (*simplex).b;
    (*simplex).size = (t == 0 || t == 1) ? 1 : 2;
    return;
}

// From Ericson's manual
void gjk_solve_simplex3(Simplex* simplex) {
    // Check if P in vertex region outside A
    Vec3 a = (*simplex).a;
    Vec3 b = (*simplex).b;
    Vec3 c = (*simplex).c;

    Vec3 ac = c - a;
    Vec3 ab = b - a;
    real d1 = dot(ab, - a);
    real d2 = dot(ac, - a);

    // Check if P in vertex region outside B
    real d3 = dot(ab, - b);
    real d4 = dot(ac, - b);

    // // Check if P in edge region of AB, if so return projection of P onto AB
    real vc = d1*d4 - d3*d2;

    // Check if P in vertex region outside C
    real d5 = dot(ab, - c);
    real d6 = dot(ac, - c);
    if (d6 >= 0.0f && d5 <= d6) {
        (*simplex).a = c; // barycentric coordinates (0,0,1)
        (*simplex).c = a;
        (*simplex).size = 1;
        (*simplex).P = c;
        return;
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    real va = d3*d6 - d5*d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
        real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        (*simplex).size = 2;
        (*simplex).a = b;
        (*simplex).b = c;
        (*simplex).c = a;
        (*simplex).P = b + w * (c - b); // barycentric coordinates (0,1-w,w)
        return;
    }

    // Check if P in edge region of AC, if so return projection of P onto AC
    real vb = d5*d2 - d1*d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        real w = d2 / (d2 - d6);
        (*simplex).size = 2;
        (*simplex).b = c;
        (*simplex).c = a;
        (*simplex).P = a + w * ac; // barycentric coordinates (1-w,0,w)
        return;
    }

    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    real denom = 1.0f / (va + vb + vc);
    real v = vb * denom;
    real w = vc * denom;
    (*simplex).P = a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
    return;
}

int point_outside_of_plane(Vec3 point, Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
    // The last input (d) is the one that will be tested
    real signp = dot(point - a, cross(b - a, c - a)); // [AP AB AC]
    real signd = dot(d - a, cross(b - a, c - a)); // [AD AB AC]
    // Points on opposite sides if expression signs are opposite
    return signp * signd < 0.0f;
}

// From Ericson's manual
void gjk_solve_simplex4(Simplex* simplex) {
    Vec3 a = (*simplex).a;
    Vec3 b = (*simplex).b;
    Vec3 c = (*simplex).c;
    Vec3 d = (*simplex).d;
    Vec3 p = {0, 0, 0};
    (*simplex).P = p;

    Simplex simplex_temp = *simplex;

    // Start out assuming point inside all halfspaces, so closest to itself
    Vec3 closestPt = p;
    real bestSqDist = INF_REAL;

    Vec3 a_temp;
    Vec3 b_temp;
    Vec3 c_temp;
    Vec3 d_temp;
    // If point outside face bdc then compute closest point on bcd
    if (point_outside_of_plane(p, b, c, d, a)) {
        gjk_solve_simplex3(&simplex_temp);
        real sqDist = dot(simplex_temp.P, simplex_temp.P);
        if (sqDist < bestSqDist) {
            a_temp = b;
            b_temp = c;
            c_temp = d;
            d_temp = a;

            bestSqDist = sqDist;
            closestPt = simplex_temp.P;
        }
    }
    // Repeat test for face acd
    if (point_outside_of_plane(p, a, c, d, b)) {
        gjk_solve_simplex3(&simplex_temp);
        real sqDist = dot(simplex_temp.P, simplex_temp.P);
        if (sqDist < bestSqDist) {
            a_temp = a;
            b_temp = c;
            c_temp = d;
            d_temp = b;

            bestSqDist = sqDist;
            closestPt = simplex_temp.P;
        }
        // Simplex simplex_temp = *simplex;
    }
    // Repeat test for face adb
    if (point_outside_of_plane(p, a, b, d, c)) {
        gjk_solve_simplex3(&simplex_temp);
        real sqDist = dot(simplex_temp.P, simplex_temp.P);
        if (sqDist < bestSqDist) {
            a_temp = a;
            b_temp = b;
            c_temp = d;
            d_temp = c;

            bestSqDist = sqDist;
            closestPt = simplex_temp.P;
        }

        // Simplex simplex_temp = *simplex;
    }
    // Repeat test for face abc
    if (point_outside_of_plane(p, a, b, c, d)) {
        gjk_solve_simplex3(&simplex_temp);
        real sqDist = dot(simplex_temp.P, simplex_temp.P);
        // Update best closest point if (squared) distance is less than current best
        if (sqDist < bestSqDist) {
            a_temp = a;
            b_temp = b;
            c_temp = c;
            d_temp = d;

            bestSqDist = sqDist;
            closestPt = simplex_temp.P;
        }
        // Simplex simplex_temp = *simplex;
    }

    if (dot(closestPt, closestPt) > COLLISION_EPSILON) {
        (*simplex).size = 3;
        (*simplex).a = a_temp;
        (*simplex).b = b_temp;
        (*simplex).c = c_temp;
        (*simplex).d = d_temp;
        (*simplex).P = closestPt;

        gjk_solve_simplex3(simplex);
    }
    return;
}

real distmin_origin(EPAHull face) {
    Vec3 a = face.p1;
    Vec3 b = face.p2;
    Vec3 c = face.p3;

    Vec3 ab = b - a;
    Vec3 ac = c - a;
    Vec3 bc = c - b;

    // Compute parametric position s for projection P’ of P on AB,
    // P’ = A + s*AB, s = snom/(snom+sdenom)
    real snom = dot(- a, ab), sdenom = dot(- b, a - b);

    // Compute parametric position t for projection P’ of P on AC,
    // P’ = A + t*AC, s = tnom/(tnom+tdenom)
    real tnom = dot(- a, ac), tdenom = dot(- c, a - c);

    if (snom <= 0.0f && tnom <= 0.0f)
        return norm(a); // Vertex region early out

    // Compute parametric position u for projection P’ of P on BC,
    // P’ = B + u*BC, u = unom/(unom+udenom)
    real unom = dot(- b, bc), udenom = dot(- c, b - c);
    if (sdenom <= 0.0f && unom <= 0.0f)
        return norm(b); // Vertex region early out

    if (tdenom <= 0.0f && udenom <= 0.0f)
        return norm(c); // Vertex region early out

    // P is outside (or on) AB if the triple scalar product [N PA PB] <= 0
    // Vec3 n = face.n;
    Vec3 n = cross(b - a, c - a);
    real vc = dot(n, cross(a, b));

    // If P outside AB and within feature region of AB,
    // return projection of P onto AB
    if (vc <= 0.0f && snom >= 0.0f && sdenom >= 0.0f)
        return norm(a + snom / (snom + sdenom) * ab);

    // P is outside (or on) BC if the triple scalar product [N PB PC] <= 0
    real va = dot(n, cross(b, c));

    // If P outside BC and within feature region of BC,
    // return projection of P onto BC
    if (va <= 0.0f && unom >= 0.0f && udenom >= 0.0f)
        return norm(b + unom / (unom + udenom) * bc);
    // P is outside (or on) CA if the triple scalar product [N PC PA] <= 0
    real vb = dot(n, cross(c, a));

    // If P outside CA and within feature region of CA,
    // return projection of P onto CA
    if (vb <= 0.0f && tnom >= 0.0f && tdenom >= 0.0f)
        return norm(a + tnom / (tnom + tdenom) * ac);
    // P must project inside face region. Compute Q using barycentric coordinates
    real u = va / (va + vb + vc);
    real v = vb / (va + vb + vc);
    real w = 1.0f - u - v; // = vc / (va + vb + vc)
    return norm(u * a + v * b + w * c);
}

// (Function is BROKEN - infinite loop)
real solve_epa_algorithm(Simplex simplex, std::vector<Vec3> v1, std::vector<Vec3> v2) {
    // create a face vector that has three points and a normal
    std::vector<EPAHull> faces;
    faces.reserve(12);

    // const int max_size = simplex.max_size;
    // Vec3 simplex_tmp[3 - (max_size - 1)];
    // int add_idx = 0;
    // while (simplex.max_size + add_idx != 4) {
    //     Vec3 sup_tmp1 = GJK_get_support(v1, -faces[idx].n);
    //     Vec3 sup_tmp2 = GJK_get_support(v2, faces[idx].n);
    //     Vec3 p_tmp = sup_tmp2 - sup_tmp1;

    //     simplex_tmp[max_size - 1 + add_idx] = dot(simplex.b, simplex.b) <= 1e-2 ? p_tmp : simplex.b;
    //     add_idx++;
    // }

    Vec3 ab = simplex.b - simplex.a;
    Vec3 ac = simplex.c - simplex.a;
    Vec3 ad = simplex.d - simplex.a;
    Vec3 bc = simplex.c - simplex.b;
    Vec3 bd = simplex.d - simplex.b;

    Vec3 n1 = cross(ab, ac);
    Vec3 n2 = cross(ab, ad);
    Vec3 n3 = cross(ac, ad);
    Vec3 n4 = cross(bc, bd);

    real dot_a_n1 = dot(simplex.a, n1);
    real dot_a_n2 = dot(simplex.a, n2);
    real dot_a_n3 = dot(simplex.a, n3);
    real dot_d_n4 = dot(simplex.d, n4);

    n1 = dot_a_n1 > 0 ? n1 : (dot_a_n1 < 0 ? - n1 : (dot(n1, simplex.d) >= 0 ? -n1 : n1));
    n2 = dot_a_n2 > 0 ? n2 : (dot_a_n2 < 0 ? - n2 : (dot(n2, simplex.d) >= 0 ? -n2 : n2));
    n3 = dot_a_n3 > 0 ? n3 : (dot_a_n3 < 0 ? - n3 : (dot(n3, simplex.d) >= 0 ? -n3 : n3));
    n4 = dot_d_n4 > 0 ? n4 : (dot_d_n4 < 0 ? - n4 : (dot(n4, simplex.a) >= 0 ? -n4 : n4));

    n1 = (1 / norm(n1)) * n1;
    n2 = (1 / norm(n2)) * n2;
    n3 = (1 / norm(n3)) * n3;
    n4 = (1 / norm(n4)) * n4;

    faces.push_back({simplex.a, simplex.b, simplex.c, n1});
    faces.push_back({simplex.a, simplex.b, simplex.d, n2});
    faces.push_back({simplex.a, simplex.c, simplex.d, n3});
    faces.push_back({simplex.b, simplex.c, simplex.d, n4});

    real min_dist;
    int idx;
    real dist;

    while (true) {
        // Find closest face

        min_dist = INF_REAL;
        for (int i = 0; i < size(faces); i++) {
            real current_dist = distmin_origin(faces[i]);
            if (current_dist < min_dist) {
                min_dist = current_dist;
                idx = i;
            }
        }

        // obtain a new support point in the direction of the edge normal
        Vec3 support1 = gjk_get_support(v1, -faces[idx].n);
        Vec3 support2 = gjk_get_support(v2, faces[idx].n);
        Vec3 p = support2 - support1;

        // If the vertex does not expand the polytope in the direction of the normal, the minimum distance
        // is with the closest face (unchanged). Compute and return.
        dist = dot(p - faces[idx].p1, faces[idx].n);
        if (dist < 1e-3) {
            break;
        }

        // Get all the faces which are "seen" by the new point, add them to deleted_faces and delete them.
        std::vector<EPAHull> deleted_faces;
        deleted_faces.reserve(size(faces)-1);
        int kept_face_idx = 0;
        // bool found_kept_face = false;
        bool condition;
        for (int i = 0; i < size(faces); i++) {
            condition = dot(faces[i].n, p - faces[i].p1) > 0;
            if (dot(faces[i].n, p - faces[i].p1) > 0) {
                deleted_faces.push_back(faces[i]);
                faces.erase(faces.begin() + i);
            }
            kept_face_idx = (condition == false && kept_face_idx == 0) ? i : kept_face_idx;
        }

        // Create new convex hull by determining which line segments are contained in the deleted faces.
        // These segments must be deleted and the ones which are only contained in one face will be added
        // as parts of the new faces (with the new point forming the remaining two line segments).
        TwoPts current_edge[3];
        std::vector<TwoPts> loose_edges;
        loose_edges.reserve(size(deleted_faces)+1);
        for (int i = 0; i < size(deleted_faces); i++) { // for all deleted faces found
            current_edge[0] = {deleted_faces[i].p1, deleted_faces[i].p2};
            current_edge[1] = {deleted_faces[i].p1, deleted_faces[i].p3};
            current_edge[2] = {deleted_faces[i].p2, deleted_faces[i].p3};

            int loose_edge_idx = 0;
            bool found_edge = false;
            for (int j = 0; j < 3; j++) { // for all three edges of each face
                found_edge = false;
                for (int k = 0; k < size(loose_edges); k++) { // Is the current edge already in loose_edges ?
                    if ((current_edge[j].p1 == loose_edges[k].p1 && current_edge[j].p2 == loose_edges[k].p2) ||
                            (current_edge[j].p2 == loose_edges[k].p1 && current_edge[j].p1 == loose_edges[k].p2)) {
                        // edge is already in list
                        found_edge = true;
                        loose_edge_idx = k;
                        break;
                    }
                }
                if (found_edge == false)
                    loose_edges.push_back(current_edge[j]);
                else {
                    loose_edges.erase(loose_edges.begin() + loose_edge_idx);
                }
            }
        }

        // rebuild simplex with new faces
        EPAHull new_face;
        Vec3 n;
        real dot_p1_n;
        for (int i = 0; i < size(loose_edges); i++) {
            n = cross(p - loose_edges[i].p1, p - loose_edges[i].p2);
            dot_p1_n = dot(loose_edges[i].p1, n);
            n = dot_p1_n > 0 ? n : (dot_p1_n < 0 ? - n : (dot(n, faces[kept_face_idx].p1) >= 0 ? -n : n));
            // if (dot_p1_n == 0.0) {
            //     int test = 0;
            // }
            n = (1 / norm(n)) * n;

            real check1 = dot(n, loose_edges[i].p1);
            real check2 = dot(n, loose_edges[i].p2);
            real check3 = dot(n, p);
            new_face = {loose_edges[i].p1, loose_edges[i].p2, p, n};
            faces.push_back(new_face);
        }
    }
    return -std::abs(min_dist);
}

// (Function is BROKEN - infinite loop in EPA) Tests 2 sets of points using GJK
real general_gjk(std::vector<Vec3> point_set1, std::vector<Vec3> point_set2) {
    // This version of the GJK algorithm is implemented from the basic algorithm described in Collision Detection
    // manual by Ericson.

    // 1. Initializing simplex to a point from a random direction
    Vec3 V = point_set1[0] - point_set2[0];
    Simplex simplex;
    simplex.a = V;
    simplex.size = 1;

    while (true) {
        // 2. Computing the point P of minimum norm in CH(Q)
        switch (simplex.size) {
        case 1:
            simplex.P = simplex.a;
            break;
        case 2:
            gjk_solve_simplex2(&simplex);
            break;
        case 3:
            gjk_solve_simplex3(&simplex);
            break;
        case 4:
            gjk_solve_simplex4(&simplex);
            break;
        default:
            Assert(false);
        }

        // 3. If P is the origin itself, the origin is clearly contained in the Minkowski difference of A and B.
        // Stop and return A and B as intersecting.
        if (dot(simplex.P, simplex.P) < 1e-6) {
            real dist = solve_epa_algorithm(simplex, point_set1, point_set2);
            return dist;
        }

        // 4. Reduce Q to the smallest subset Q' of Q such that P is still in CH(Q). That is, remove any points
        // from Q not determining the subsimplex of Q in which P lies.
        // (This is done automatically in GJK_solve_simplex functions)

        // 5. Find the next supporting point in direction -P
        Vec3 a1 = gjk_get_support(point_set1, simplex.P);
        Vec3 a2 = gjk_get_support(point_set2, -simplex.P);
        V = a2 - a1;

        // 6. If V is no more exremal in direction -P than P itself, stop and return A and B as not intersecting.
        // The length of the vector from the origin to P is the separation distance of A and B.

        real ans1 = dot(V, -simplex.P)/dot(simplex.P, simplex.P);

        if (ans1 + 1 <= COLLISION_EPSILON)   // No more progress is being made (<5 %). If we do not do this it can cause problems
            break;              // in the evaluation of PointOutsidePlane in evaluation of closest point on simplex4

        // 7. Add V to Q and go to 2.
        // todo : optimize this part of code using std::vector instead of Vec3 for a, b, c, d
        Assert(simplex.size <=3);

        if (simplex.size == 1)
            simplex.b = V;
        if (simplex.size == 2)
            simplex.c = V;
        if (simplex.size == 3)
            simplex.d = V;
        simplex.size += 1;
        simplex.max_size = simplex.size > simplex.max_size ? simplex.size : simplex.max_size;
    }
    return norm(simplex.P);
}

// (Function is BROKEN - infinite loop in GJK)
real gjk_box_capsule(Capsule capsule, Box box) {
    // Initialization of the eight Box points
    Vec3 size_x_org = {box.e.x, 0, 0};
    Vec3 size_y_org = {0, box.e.y, 0};
    Vec3 size_z_org = {0, 0, box.e.z};
    Vec3 size_x = box.R*size_x_org;
    Vec3 size_y = box.R*size_y_org;
    Vec3 size_z = box.R*size_z_org;

    std::vector<Vec3> v1(2);
    v1[0] = capsule.p1;
    v1[1] = capsule.p2;

    std::vector<Vec3> v2(8);
    v2[0] = box.c + size_x + size_y + size_z;
    v2[1] = box.c + size_x + size_y - size_z;
    v2[2] = box.c + size_x - size_y + size_z;
    v2[3] = box.c + size_x - size_y - size_z;
    v2[4] = box.c - size_x + size_y + size_z;
    v2[5] = box.c - size_x + size_y - size_z;
    v2[6] = box.c - size_x - size_y + size_z;
    v2[7] = box.c - size_x - size_y - size_z;

    real dist = general_gjk(v1, v2) - capsule.r;
    return dist;
}

#ifdef BLAST_ENABLE_TESTS

// capsule - Sphere collision tests
struct collision_test_sph {
    Sphere sphere;
    Capsule capsule;
    real expected_dist;
};
collision_test_sph test[] = {
    { { { 4.27, 1.73, 3.74 }, 1 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, 0.42 },
    { { { 2.75, 1.48, 3.44}, 1 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, -0.25 },
    { { { 5.72, 6.47, 3.78 }, 1 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, 1.4 },
    { { { 4.74, 5.26, 4.52 }, 1 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, -0.18 },
};

// capsule - capsule collision tests
struct collision_test_caps {
    Capsule capsule1;
    Capsule capsule2;
    real expected_dist;
};
collision_test_caps test_caps[] = {
    { { { 0.5, 0.95, 6.73 }, { 4, 4.9, 3.51 }, 0.5 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, -0.54 },
    { { { 5.16, 4.9, 4.37 }, { 3.55, 0.95, 8.84 }, 0.5 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, 0.16 },
    { { { 5.16, 4.9, 4.37 }, { 3.55, 0.95, 8.84 }, 0.5 }, { { 5.45, 4.78, 4.57 }, { 6.57, 1, -0.2 }, 0.5 }, -0.66 },
    { { { 2.59, 1.88, 4.22 }, { 18.85, 10.07, 6.1 }, 1 }, { { 0.55, 21.32, 3.08 }, { 5.07, 3.62, 4.19 }, 1 }, -1.44 },
    { { { -9.26, 12.81, 7.98 }, { 6.63, 4.62, 4.04 }, 1 }, { { 0.55, 21.32, 3.08 }, { 5.07, 3.62, 4.19 }, 1 }, -1.52 },
};

// capsule - cylinder collision tests
struct collision_test_cyl {
    Cylinder cylinder;
    Capsule capsule;
    real expected_dist;
};
collision_test_cyl test_cyl[] = {
    /*Test 1*/ { { { -13.46, -3.61, 189 }, { -12.46, -5.11, 190 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, -0.538346 },
    /*Test 2*/ { { { -13.41, -2.79, 189 }, { -15.16, -3.2, 190 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, 1.32728712 },
    /*Test 3*/ { { { -16.3, -3.97, 200.87 }, { -17.3, -2.47, 199.87 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, 1.53080644 },
    /*Test 4*/ { { { -7.25, -10.11, 174.29 }, { -7.86, -9.44, 176.14 }, 2 }, { { -14.66, -13.24, 198.94 }, { -11.16, -7.94, 181.77 }, 1 }, 5.52119334 },
    /*Test 5*/ { { { -14.88, -2.74, 198.21 }, { -15.49, -2.08, 200.06 }, 2 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, -0.44703889 },
};

// capsule - Box collision tests
struct collision_test_box {
    Box box;
    Capsule capsule;
    real expected_dist;
};
collision_test_box test_obb[] = {
    /*Test1*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { -1.21, -5.18, 18.05 }, { -3.89, 8.59, 6.3 }, 1 }, 1.63659624 },
    /*Test2*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 10.46, 0.97, 3.7 }, { 7.79, 14.74, -8.04 }, 1 }, 1.50169942 },
    /*Test3*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 10.05, 1.11, 12.87 }, { 7.37, 14.89, 1.13 }, 1 }, 0.33397901 },
    /*Test4*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.82, 6.9, 12.84 }, { 4.7, -7.14, 24.59 }, 1 }, 4.00838054 },
    /*Test5*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { -3.06, 5.73, 1.67 }, { -0.39, -8.05, 13.41 }, 1 }, 0.89513819 },
    /*Test6*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.97, 2.79, 12.23 }, { 2.29, 16.57, 0.48 }, 1 }, 1.69981481 },
    /*Test7*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 11.49, -0.06, 8.32 }, { 14.17, -13.84, 20.07 }, 1 }, 0.49000000 },
    /*Test8*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 6.76, -1.55, 8 }, { 9.44, -15.33, 19.74 }, 1 }, 0.45000000 },
    /*Test9*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { -1.15, 13.92, -7.48 }, { 1.53, 0.14, 4.27 }, 1 }, 0.73046237 },
    /*Test10*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.54, 0.1, 10.59 }, { 1.86, 13.87, -1.15 }, 1 }, -1.00000000 },
    /*Test11*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 10.76, 0.28, 8.08 }, { 13.44, -13.5, 19.83 }, 1 }, -0.21961454 },
    /*Test12*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.96, 0.43, 8.3 }, { 7.64, -13.35, 20.05 }, 1 }, -1.53000000 },
    /*Test13*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.48, -4.07, 0.76 }, { 8.64, -4.41, 18.58 }, 1 }, 3.06923695 },
    /*Test14*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 17.48, 2.95, 13.77 }, { -0.82, 3.11, 13.77 }, 1 }, 2.41054264 },
    /*Test15*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 11.18, 8.56, 4.82 }, { 11.04, -6.44, 15.29 }, 1 }, 0.09912546 },
    // Box.R and capsule.r changed
    /*Test16*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 2.94, -6.06, 5.74 }, { 4.01, -9.86, 0.98 }, 0.5 }, 1.00474115 },
    /*Test17*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 1.64, 9.52, 11.7 }, { 2.71, 5.72, 6.94 }, 0.5 }, 0.22669533 },
    /*Test18*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 3, 5.69, 6.2 }, { 4.07, 1.89, 1.44 }, 0.5 }, 0.42102531 },
    /*Test19*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 3.1, 4.23, 6.22 }, { 4.17, 0.43, 1.46 }, 0.5 },  0.10695205 },
    /*Test20*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 6.35, 8.57, 12.85 }, { 7.42, 4.77, 8.09 }, 0.5 }, 0.85902522 },
    /*Test21*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 1.88, -2.47, 7.07 }, { 2.95, -6.27, 2.31 }, 0.5 }, 0.37892454 },
    /*Test22*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 1.8, 3.83, 14.57 }, { 2.87, 0.03, 9.81 }, 0.5 }, 1.47889394 },
    /*Test23*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 4.28, 6.32, 8.33 }, { 3.21, 10.12, 13.09 }, 0.5 }, 0.82000000 },
    /*Test24*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 2.36, 2.3, 6.31 }, { 3.43, -1.5, 1.55 }, 0.5 }, 0.26887914 },
    /*Test25*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 2.93, 7.6, 12.29 }, { 4, 3.8, 7.53 }, 0.5 }, -0.93234019 },
    /*Test26*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 7.55, 0.92, 10.24 }, { 6.48, 4.72, 15 }, 0.5 }, -0.32852683 },
    /*Test27*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 6.79, 5, 13.29 }, { 7.86, 1.2, 8.52 }, 0.5 }, -0.40212621 },
    /*Test28*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 5.66, -0.92, 7.44 }, { 8.04, 4.16, 9.96 }, 0.5 }, 0.87078210 },
    /*Test29*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 8.94, 3.63, 11.01 }, { 8.94, -2.55, 11.01 }, 0.5 }, 1.24844065 },
    /*Test30*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 5.91, 5.9, 7.39 }, { 4.14, 5.9, 13.32 }, 0.5 }, 0.40000000 },

    /*Test31*/ { { { 0, 0, 0 }, { 1.421459, 0.796365, 0.574752 }, { 1, 0, 0, 0, 1, 0, 0, 0, 1 } }, { { 2.447825, 6.30643376, -5.224298 }, { -5.4486154, -5.149872, 5.1825604 }, 0 }, -0.05877295 },
    /*Test32*/ { { { 0, 0, 0 }, { 1.7236252662212059, 1.6611691154387149, 0.20442662553572166 }, { 1, 0, 0, 0, 1, 0, 0, 0, 1 } }, { { -4.6785661670641119, 1.5869561202501099, -4.7958108710825522 }, { 6.9379851583657839, -2.4864207406657761, 4.5569750611707658 }, 0 }, -0.438675 },
    /*Test33*/ { { { 0, 0, 0 }, { 1.7865788377999734, 1.4615853859974308, 1.7957946652540584 }, {0.16899262598095025, 0.69208966716438747, -0.70175022976009782, -0.73366055668912455, 0.56377691511991157, 0.37933860538637498, 0.65816690886329632, 0.45074103716231956, 0.60303303184416857 } }, { { -0.050034305603475548, 0.40002630593895194, 2.3444724393109695 }, { 1.3614165219982919, 0.62636518062708535, -1.2242779633088015 }, 0 }, -1.3557069221709168 },
};

TEST_CASE("Collisions", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

    std::vector<Vec3> v1 = {    { -4.6607382574035014, -3.4686124971115451,  0.89349474531514073 },
        {  6.4517093458638977,  6.9310459212756506, -0.24683715666348305 }
    };
    std::vector<Vec3> v2 = {    {  2.4373332611708229,  2.5910471347625212, -0.04489407117738407 },
        {  2.5547766616001675,  4.0525983620611852,  0.19409068785959360 },
        {  1.8335195262912714,  2.7351427068670993, -0.62940465909289878 },
        {  1.9509629267206159,  4.1966939341657632, -0.39041990005592114 },
        {  1.6790895426140411,  2.5264996514086859,  0.72247776841327949 },
        {  1.7965329430433856,  3.9880508787073494,  0.96146252745025707 },
        {  1.0752758077344895,  2.6705952235132639,  0.13796718049776469 },
        {  1.1927192081638340,  4.1321464508119279,  0.37695193953474238 }
    };
    // real dist = general_gjk(v1, v2);

    for (auto t : test) {
        real dist = distance(t.capsule, t.sphere);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_caps) {
        real dist = distance(t.capsule1, t.capsule2);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_cyl) {
        real dist = distance(t.capsule, t.cylinder);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_obb) {
        real dist = distance_old(t.box, t.capsule);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_obb) {
        real dist = distance(t.box, t.capsule);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }
}

TEST_CASE("Collision method benchmarks", "[World]") {
    using namespace blast;

    // real TESTCOLL_EPSILON = 1e-2;

    BENCHMARK("box - capsule Box test with vectors") {
        real dist;
        for (auto t : test_obb) {
            dist = distance_old(t.box, t.capsule);
            // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
        }
        return dist;
    };

    BENCHMARK("box - capsule Box test with vectors accelerated") {
        real dist;
        for (auto t : test_obb) {
            dist = distance(t.box, t.capsule);
            // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
        }
        return dist;
    };

    // BENCHMARK("box - capsule Box test with GJK") {
    //     real dist;
    //     for (auto t : test_obb) {
    //         dist = gjk_box_capsule(t.capsule, t.box);
    //         // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
    //     }
    //     return dist;
    // };
}

TEST_CASE("Collision method comparison exhaustive (Box-cpasules)", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

    vector<Box> obb_list;
    vector<capsule> caps_list;
    int n = 1000;
    // bool error_info = false;

    Array s(3);
    Array c(3);
    Array angles(3);

    Array obb_e(3);
    Array obb_c(3);

    Array seg_p1(3);
    Array seg_p2(3);

    for(u32 i = 0; i < n; i++) {
        fill_random(angles, PI);
        blast::sincos(angles, s, c);

        Mat3 Rz = {c[0], -s[0], 0, s[0], c[0], 0, 0, 0, 1};
        Mat3 Ry = {c[1], 0, s[1], 0, 1, 0, -s[1], 0, c[1]};
        Mat3 Rx = {1, 0, 0, 0, c[2], -s[2], 0, s[2], c[2]};
        auto R = Rz*Ry*Rx;
        CHECK(abs(1.0 - sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2])) < TESTCOLL_EPSILON);
        CHECK(abs(1.0 - sqrt(R[3]*R[3] + R[4]*R[4] + R[5]*R[5])) < TESTCOLL_EPSILON);
        CHECK(abs(1.0 - sqrt(R[6]*R[6] + R[7]*R[7] + R[8]*R[8])) < TESTCOLL_EPSILON);

        // Generate n random Box
        fill_random(obb_e, 2);
        obb_e = array_abs(obb_e);
        fill_random(obb_c, 5);
        obb_c = array_abs(obb_c);
        obb_list.push_back({{obb_c[0], obb_c[1], obb_c[2]}, {obb_e[0], obb_e[1], obb_e[2]}, R});

        // Generate n random capsule
        fill_random(seg_p1, 7);
        fill_random(seg_p2, 7);
        auto r = 0.0;
        caps_list.push_back({{seg_p1[0], seg_p1[1], seg_p1[2]}, {seg_p2[0], seg_p2[1], seg_p2[2]}, r});
    }

    vector<capsule> caps_failed;
    vector<Box> obb_failed;
    for (int cap = 0; cap < caps_list.size(); cap++) {
        // Box collisions
        for (int i = 0; i < obb_list.size(); i++) {
            auto dist_min_vec = distance_old(obb_list[i], caps_list[cap]);
            auto dist_min_vector_acc = distance(obb_list[i], caps_list[cap]);
            // auto dist_min_new = distmin_hull(obb_list[i], caps_list[cap]);
            // auto dist_min_gjk = gjk_box_capsule(caps_list[cap], obb_list[i]);

            CHECK(abs(dist_min_vec - dist_min_vector_acc) < TESTCOLL_EPSILON);
            if (abs(dist_min_vec - dist_min_vector_acc) > TESTCOLL_EPSILON) {
                // save and test capsule and obb in future
                caps_failed.push_back(caps_list[cap]);
                obb_failed.push_back(obb_list[i]);
            }
        }
    }

    // real total_error = 0;
    // real max_error = 0;
    // bool error_when_neg = false;
    // bool error_when_pos = false;
    // std::vector<Vec3> error_distmin_distminnew(caps_failed.size());
    // for (u32 i = 0; i < caps_failed.size(); i++) {
    //     auto dist_min_vec = distance_old(obb_failed[i], caps_failed[i]);
    //     auto dist_min_new = distmin_hull(obb_failed[i], caps_failed[i]);
    //     auto dist_min_gjk = gjk_box_capsule(caps_failed[i], obb_failed[i]);
    //     auto dist_min_pierce = distmin_pierce(obb_failed[i], caps_failed[i]);
    //     real error = abs(dist_min_new - dist_min_vec) * 100 / dist_min;
    //     total_error += error;
    //     int dist = 0;
    //     max_error = error > max_error ? error : max_error;
    //     if (error > 100)
    //         real found_error = 0;
    //     if (dist_min_vec >=0)
    //         error_when_pos = true;
    //     if (dist_min_vec < 0) {
    //         error_when_neg = true;
    //         real check = dist_min_new + dist_min_vec;
    //     }
    //     if (error_info) {
    //         Vec3 error_point = { error, dist_min_vec, dist_min_new };
    //         error_distmin_distminnew.push_back(error_point);
    //     }
    // }
    // real avg_error = total_error / caps_failed.size();
    // real avg_error_over_all = total_error / (n*n);
    // real percent_error = caps_failed.size() / (n*n);
    // std::cout << caps_failed.size() << " failed out of " << n*n << " tests ( " << percent_error << " % ) \n";
    // std::cout << "Average error of " << avg_error << "% and a max error of " << max_error << "%. \n";
    // std::cout << "Counting the tests which passed, we find an average error of : " << avg_error_over_all << " %.  \n";
    // std::cout << "Error when positive dist (true or false): " << error_when_pos << " \n";
    // std::cout << "Error when negative dist (true or false): " << error_when_neg << " \n";
    // if (error_info) {
    //     #include <fstream>
    //     // Open a CSV file for writing
    //     std::ofstream csvFile("C:\\Users\\thoma\\Desktop\\example.csv");
    //     // Check if the file is open
    //     if (!csvFile.is_open()) {
    //         std::cerr << "Error opening file!" << std::endl;
    //     }
    //     // Write headers to the CSV file
    //     csvFile << "error (%), dist_min_vec, dist_min_new" << std::endl;
    //     // Write data
    //     for (int i = 0; i < error_distmin_distminnew.size(); i++) {
    //         Vec3 error_pt = error_distmin_distminnew[i];
    //         csvFile << (error_pt).x << ";" << (error_pt).y << ";" << (error_pt).z << std::endl;
    //     }
    //     // Close the file
    //     csvFile.close();
    //     std::cout << "CSV file created successfully." << std::endl;
    // }
}

#endif

} // namespace blast
