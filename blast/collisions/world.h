#pragma once
#include <iostream>
#include <vector>
#include "blast_math.hpp"

using std::vector;

// todo: cleanup required
// todo: make this fast

namespace blast {

#define COLLISION_EPSILON 1e-9

struct segment {
    Vec3 p1;
    Vec3 p2;
};

struct OBB {
    Vec3 c; // OBB center point
    Vec3 e; // Positive halfwidth extents of OBB along each axis
    Mat3 R; // Local x-, y-, and z-axes (Rotation matrix)
};

struct surf {
    Vec3 p;
    Vec3 d1;
    Vec3 d2;
};

struct capsule {
    Vec3 p1;
    Vec3 p2;
    real r;
};

struct sphere {
    Vec3 c;
    real r;
};

struct cylinder {
    Vec3 p1;
    Vec3 p2;
    real r;
};

struct plane {
    Vec3 p;
    Vec3 n;
};

struct two_pts {
    Vec3 p1;
    Vec3 p2;
};

struct circle {
    Vec3 p;
    Vec3 n;
    real r;
};

struct triangle {
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;
    real distance;
};

struct objlist {
    std::vector<OBB> OBBlist;
    std::vector<cylinder> cyllist;
    std::vector<sphere> sphlist;
    std::vector<capsule> capslist;
};

struct capslist {
    std::vector<capsule> caps; // list of capsules
};

// ======================================
//            Basic functions
// ======================================

host_fn bool pointInTriangle(Vec3 V1, Vec3 V2, Vec3 o, Vec3 P) {
    // u=P2−P1
    Vec3 u = V1;
    // v=P3−P1
    Vec3 v = V2;
    // n=u×v
    Vec3 n = cross(u, v);
    // w=P−P1
    Vec3 w = P - o;
    // Barycentric coordinates of the projection P′of P onto T:
    // γ=[(u×w)⋅n]/n²
    real gamma = dot(cross(u, w), n) / dot(n, n);
    // β=[(w×v)⋅n]/n²
    real beta = dot(cross(w, v), n) / dot(n, n);
    real alpha = 1 - gamma - beta;
    // The point P′ lies inside T if:
    return ((0 <= alpha) && (alpha <= 1) && (0 <= beta) && (beta <= 1) && (0 <= gamma) && (gamma <= 1));
}

host_fn bool pointInSurf(Vec3 V1, Vec3 V2, Vec3 o, Vec3 P) {
    bool tri1 = pointInTriangle(V1, V2, o, P);
    bool tri2 = pointInTriangle(-V1, -V2, o+V1+V2, P);
    return (tri1 || tri2);
}

// todo: change name or add description
host_fn Vec3 ptint(segment seg, Vec3 point) {
    Vec3 ab = seg.p2 - seg.p1;
    real t = dot(point - seg.p1, ab) / dot(ab, ab);
    t = clamp(t, 0, 1);
    Vec3 d = seg.p1 + t * ab;
    return d;
}

host_fn real distmin(segment seg, Vec3 point) {
    Vec3 d = ptint(seg, point);
    return norm(d - point);
}

host_fn Vec3 closept_origin(segment seg) {
    Vec3 ab = seg.p2 - seg.p1;
    real t = dot(- seg.p1, ab) / dot(ab, ab);

    t = clamp(t, 0, 1);

    Vec3 d = seg.p1 + t * ab;
    return d;
}

//todo: change to use the snake case function name convention
host_fn Vec3 ClosestPtPointPlane(Vec3 q, plane p) {
    return q - (dot(p.n, q - p.p) / dot(p.n, p.n)) * p.n;;
}

blast_fn real distmin(surf surf, Vec3 point) {

    Vec3 n = cross(surf.d1, surf.d2);
    if (dot(n, n) < COLLISION_EPSILON)
        return -INF_REAL; // the plane is non-existant
    Vec3 n_unit = 1/norm(n)*n;

    // If the surface is a rectangular shape, it is much easier to treat.
    if (abs(dot(surf.d1, surf.d2)) <= COLLISION_EPSILON) {
        // norm^2 d1 & d2
        real val_d1_sq = dot(surf.d1, surf.d1);
        real val_d2_sq = dot(surf.d2, surf.d2);
        // direction vector from p to point
        Vec3 val_direction = point - surf.p;

        // How far along is the point in d1 & d2 direction
        real val_t1 = dot(val_direction, surf.d1) / val_d1_sq;
        real val_t2 = dot(val_direction, surf.d2) / val_d2_sq;
        real val_t1_clamped = clamp(val_t1, 0, 1);
        real val_t2_clamped = clamp(val_t2, 0, 1);
        real val_normaldist = dot(val_direction, n_unit);

        real val_dist1 = (val_t1 - val_t1_clamped)*(val_t1 - val_t1_clamped)*val_d1_sq;
        real val_dist2 = (val_t2 - val_t2_clamped)*(val_t2 - val_t2_clamped)*val_d2_sq;

        real val_result = sqrt(val_dist1 + val_dist2 + val_normaldist*val_normaldist);

        real val_distance = (val_t1 >= 0 && val_t1 <= 1 && val_t2 >= 0 && val_t2 <= 1) ? val_normaldist : val_result;
        return val_distance == val_normaldist ? val_normaldist : (val_normaldist < 0 ? -val_distance : val_distance);
    }

    //todo: remove??
    // else {
    //     real normaldist = dot(point - surf.p, n_unit);
    //     Vec3 direction = point - normaldist*n_unit - surf.p;
    //     bool is_inside = pointInSurf(surf.d1, surf.d2, surf.p, point);
    //     if (is_inside) {
    //         return normaldist;
    //     }
    //     segment seg[4];
    //     seg[0].p1 = surf.p;
    //     seg[0].p2 = surf.p + surf.d1;
    //     seg[1].p1 = surf.p;
    //     seg[1].p2 = surf.p + surf.d2;
    //     seg[2].p1 = surf.p + surf.d1;
    //     seg[2].p2 = surf.p + surf.d1 + surf.d2;
    //     seg[3].p1 = surf.p + surf.d2;
    //     seg[3].p2 = surf.p + surf.d2 + surf.d1;
    //     real dist_min = INF_REAL;
    //     for (u32 i = 0; i < 4; i++) {
    //         auto dist_tmp = distmin(seg[i], point);
    //         dist_min = dist_tmp < dist_min ? dist_tmp : dist_min;
    //     }
    //     real resultant = normaldist < 0 ? -dist_min : dist_min;
    //     return resultant;
    // }
    // else {
    //     real normaldist = dot(point - surf.p, n_unit);
    //     Vec3 direction = (point - normaldist*n_unit) - surf.p;
    //     Vec3 dy = cross(n_unit, surf.d1);
    //     Vec3 dy_unit = 1 / norm(dy) * dy;
    //     real t2 = dot(dy, direction) / dot(dy, surf.d2);
    //     real t1 = (dot(surf.d1, direction) - t2*dot(surf.d1, surf.d2)) / dot(surf.d1, surf.d1);
    //     // If the point projects on the face
    //     if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
    //         return normaldist;
    //     }
    //     Vec3 empty3 = {0, 0, 0};
    //     // If the point is closest to a corner of the shape
    //     if ((t1 <= 0 || t1 >= 1) && (t2 <= 0 || t2 >= 1)) {
    //         Vec3 testpoint = t1 <= 0 ? (t2 <= 0 ? empty3 :  surf.d2) : (t2 <= 0 ? surf.d1 : surf.d1 + surf.d2);
    //         real distance = sqrt(dot(direction - testpoint, direction - testpoint) + normaldist*normaldist);
    //         return normaldist < 0 ? -distance : distance;
    //     }
    //     real t1_clamped = clamp(t1, 0, 1);
    //     real t2_clamped = clamp(t2, 0, 1);
    //     Vec3 Projection_vector = (t1 < 1 && t1 > 0) ? dy_unit : cross(1/norm(surf.d2)*surf.d2, n_unit);
    //     Projection_vector = dot(Projection_vector, direction) < 0 ? -Projection_vector : Projection_vector;
    //     Vec3 testpoint = (t1 < 0 || t2 < 0) ? empty3 : (surf.d1 + surf.d2);
    //     real dist_plan = dot(Projection_vector, direction - testpoint);
    //     Vec3 projected_point = direction - dist_plan*Projection_vector;
    //     Vec3 clamping_direction = (t1 < 0 || t2 < 0) ? ((t1 < 1 && t1 > 0) ? surf.d1 : surf.d2) : ((t1 < 1 && t1 > 0) ? -surf.d1 : -surf.d2);
    //     real t_current = dot(direction - testpoint, clamping_direction) / dot(clamping_direction, clamping_direction);
    //     real t_current_clamped = clamp(t_current, 0, 1);
    //     real difference_t_current = (t_current - t_current_clamped);
    //     real dist_current = difference_t_current*difference_t_current*dot(clamping_direction, clamping_direction);
    //     real result = sqrt(dist_current + dist_plan*dist_plan + normaldist*normaldist);
    //     real distance = (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) ? normaldist : result;
    //     return distance == normaldist ? normaldist : (normaldist < 0 ? -distance : distance);
    // }
    else {
        real normaldist = dot(point - surf.p, n_unit);
        Vec3 direction = point - normaldist*n_unit - surf.p;

        bool is_inside = pointInSurf(surf.d1, surf.d2, surf.p, point);
        if (is_inside) {
            return normaldist;
        }

        segment seg[4];
        seg[0].p1 = surf.p;
        seg[0].p2 = surf.p + surf.d1;
        seg[1].p1 = surf.p;
        seg[1].p2 = surf.p + surf.d2;
        seg[2].p1 = surf.p + surf.d1;
        seg[2].p2 = surf.p + surf.d1 + surf.d2;
        seg[3].p1 = surf.p + surf.d2;
        seg[3].p2 = surf.p + surf.d2 + surf.d1;

        real dist_min = INF_REAL;
        int idx;
        for (u32 i = 0; i < 4; i++) {
            auto dist_tmp = distmin(seg[i], point);
            dist_min = dist_tmp < dist_min ? dist_tmp : dist_min;
            idx = dist_tmp < dist_min ? i : idx;
        }
        real resultant = normaldist < 0 ? -dist_min : dist_min;

        //todo: remove??
        // real dist;
        // // real normaldist = dot(point - surf.p, n_unit);
        // Vec3 proj = (point - normaldist*n_unit) - surf.p;
        // Vec3 dy = cross(n_unit, surf.d1);
        // real y_max = dot(dy, surf.d2);
        // real y_pt = dot(dy, proj);
        // real x_min =    y_pt < 0 ? 0 :
        //                 (y_pt > 1 ? dot(surf.d2, surf.d1) :
        //                 y_pt/y_max * dot(surf.d2, surf.d1));
        // real x_max = x_min + dot(surf.d1, surf.d1);
        // real x_pt = dot(surf.d1, proj);
        // if (x_pt >= x_min && x_pt <= x_max) {
        //     real diff_y = clamp(y_pt, 0, y_max);
        //     real dist_plan = (y_pt - diff_y) / norm(dy);
        //     dist = sqrt(dist_plan*dist_plan + normaldist*normaldist);
        //     dist = normaldist < 0 ? -dist : dist;
        // } else {
        //     segment seg[2];
        //     seg[0].p1 = surf.p;
        //     seg[0].p2 = surf.p + surf.d2;
        //     seg[1].p1 = surf.p + surf.d1;
        //     seg[1].p2 = surf.p + surf.d1 + surf.d2;
        //     int idx = x_pt < x_min ? 0 : 1;
        //     dist = distmin(seg[idx], point);
        //     dist = normaldist < 0 ? - dist : dist;
        //     // Vec3 temp = x_pt > x_max ? proj - surf.d1 : proj;
        //     // real t = dot(temp, surf.d2) / dot(surf.d2, surf.d2);
        //     // t = clamp(t, 0, 1);
        //     // Vec3 test_pt = x_pt > x_max ? t*surf.d2 + surf.d1 : t*surf.d2;
        //     // dist = normaldist < 0 ? -norm(test_pt - point) : norm(test_pt - point);
        // }
        return resultant;
    }
}

// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t.
// From : Real-Time Detection Collision (Christer Ericson)
host_fn two_pts closept(segment seg1, segment seg2) {
    two_pts result;
    real s;
    real t;
    Vec3 c1;
    Vec3 c2;

    Vec3 d1 = seg1.p2 - seg1.p1; // Direction vector of segment S1
    Vec3 d2 = seg2.p2 - seg2.p1; // Direction vector of segment S2
    Vec3 r = seg1.p1 - seg2.p1;
    real a = dot(d1, d1); // Squared length of segment S1, always nonnegative
    real e = dot(d2, d2); // Squared length of segment S2, always nonnegative
    real f = dot(d2, r);
    // Check if either or both segments degenerate into points
    if (a <= COLLISION_EPSILON && e <= COLLISION_EPSILON) {
        // Both segments degenerate into points
        s = t = 0.0f;
        c1 = seg1.p2;
        c2 = seg2.p2;
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
    c1 = seg1.p1 + d1 * s;
    c2 = seg2.p1 + d2 * t;
    result.p1 = c1;
    result.p2 = c2;
    return result;
}

// returns squared distance
host_fn real distmin_sq(segment seg1, segment seg2) {
    two_pts close_pt = closept(seg1, seg2);
    return dot(close_pt.p1 - close_pt.p2, close_pt.p1 - close_pt.p2);
}

host_fn two_pts intersection(circle circ, segment seg) {
    Vec3 p1 = seg.p1 - circ.p;
    Vec3 p2 = seg.p2 - circ.p;

    Vec3 d1 = p1;
    Vec3 d2 = cross(circ.n, d1);

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

    real det = (circ.r * circ.r * dr_sq - D * D);

    if (det < 0) {
        Vec3 point = ptint(seg, circ.p);
        Vec3 pointcirc = circ.p + (circ.r / norm(point - circ.p)) * (point - circ.p);
        return { pointcirc, pointcirc };
    }

    det = sqrt(det); // note: only sqrt if det is non-negative
    real sign_dy = dy > 0 ? 1 : -1;
    real x_1 = (1 / dr_sq) * (D * dy + sign_dy * dx * det);
    real x_2 = (1 / dr_sq) * (D * dy - sign_dy * dx * det);

    real y_1 = (1 / dr_sq) * (-D * dx + abs(dy) * det);
    real y_2 = (1 / dr_sq) * (-D * dx - abs(dy) * det);

    Vec3 p_1 = circ.p + x_1 * d1_unit + y_1 * d2_unit;
    Vec3 p_2 = circ.p + x_2 * d1_unit + y_2 * d2_unit;
    return { p_1, p_2 };
}

// ======================================
//            Primitive tests
// ======================================

host_fn real distmin(capsule caps, cylinder cyl) {
    segment seg1;
    segment seg2;
    seg1.p1 = caps.p1;
    seg1.p2 = caps.p2;
    seg2.p1 = cyl.p1;
    seg2.p2 = cyl.p2;

    two_pts points = closept(seg1, seg2);
    real cond1 = dot(points.p2 - cyl.p1, points.p2 - cyl.p1);
    real cond2 = dot(points.p2 - cyl.p2, points.p2 - cyl.p2);

    // Depending on which point of the cylinder is closest, we will change the face which we check
    if (cond1 < COLLISION_EPSILON || cond2 < COLLISION_EPSILON) {
        Vec3 cent;
        Vec3 other;

        bool corner = cond2 < cond1;
        cent = corner ? cyl.p2 : cyl.p1;
        other = corner ? cyl.p1 : cyl.p2;

        // Check if both points project on the circle plane
        Vec3 n = cent - other;
        Vec3 n_unit = (1 / norm(n)) * n;

        plane face1;
        face1.n = n_unit;
        face1.p = cent;

        Vec3 proj1 = ClosestPtPointPlane(caps.p1, face1);
        Vec3 proj2 = ClosestPtPointPlane(caps.p2, face1);

        real rad_sq1 = dot(proj1 - cent, proj1 - cent);
        real rad_sq2 = dot(proj2 - cent, proj2 - cent);

        real dist_norm_sq1 = dot(caps.p1 - proj1, caps.p1 - proj1);
        real dist_norm_sq2 = dot(caps.p2 - proj2, caps.p2 - proj2);

        real dist_norm_sq_min = (dist_norm_sq2 < dist_norm_sq1) ? dist_norm_sq2 : dist_norm_sq1;

        // if both points project on the circle plane, then the distance will be the minimum of the two normal distances calculated
        if (rad_sq1 <= cyl.r * cyl.r && rad_sq2 <= cyl.r * cyl.r)
            return sqrt(dist_norm_sq_min) - caps.r;

        // if one point projects on the circle plane, it is necessary to check the normal distance as well
        real testnormal = (rad_sq1 <= cyl.r * cyl.r) ? dist_norm_sq1 : (rad_sq2 <= cyl.r * cyl.r) ? dist_norm_sq2 : INF_REAL;

        // project the points on the plan
        segment proj;
        proj.p1 = proj1;
        proj.p2 = proj2;
        circle circ;
        circ.p = cent;
        circ.r = cyl.r;
        circ.n = n_unit;

        two_pts pts = intersection(circ, proj);

        real dist1 = distmin(seg1, pts.p1);
        real dist2 = distmin(seg1, pts.p2);

        if (dist2 <= dist1 && dist2 * dist2 < testnormal)
            return dist2 - caps.r;

        if (dist1 < dist2 && dist1 * dist1 < testnormal)
            return dist1 - caps.r;

        return sqrt(testnormal) - caps.r;
    }
    else
        return norm(points.p1 - points.p2) - caps.r - cyl.r;
}

host_fn real distmin(capsule caps, sphere sph) {
    segment seg;
    seg.p1 = caps.p1;
    seg.p2 = caps.p2;
    real dist_seg_pt = distmin(seg, sph.c);
    return dist_seg_pt - caps.r - sph.r;
}

host_fn real distmin(capsule caps1, capsule caps2) {
    //todo: why does p2 become p1? Comment if intentionnal.
    segment seg1 {caps1.p2, caps1.p1};
    segment seg2 {caps2.p2, caps2.p1};
    real dist_seg_seg = /*two_segment_distance_sqr(caps1.p1, caps1.p2, caps2.p1, caps2.p2);*/ sqrt(distmin_sq(seg1, seg2));
    return dist_seg_seg - caps1.r - caps2.r;
}

// OBB caps tests basic
host_fn real distmin_vec(OBB OBB, capsule caps) {
    segment seg;
    seg.p1 = caps.p1;
    seg.p2 = caps.p2;

    Mat3 Rtrans = transpose(OBB.R);
    // Mat3 Rtrans = OBB.R;

    Vec3 p1 = Rtrans * (seg.p1 - OBB.c);
    Vec3 p2 = Rtrans * (seg.p2 - OBB.c);

    Vec3 vec = p2 - p1;
    Vec3 point = p2;

    if (vec.z < 0) {
        vec = -vec;
        point = p1;
    }

    real xmin = - OBB.e[0];
    real xmax = + OBB.e[0];
    real ymin = - OBB.e[1];
    real ymax = + OBB.e[1];
    real zmin = - OBB.e[2];
    real zmax = + OBB.e[2];

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
    surf face[12];

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
    dist[0] = distmin(face[0], point);

    real distcurrent = dist[0];

    for (int i = 1; i < 12; i++) {
        dist[i] = distmin(face[i], point);

        bool dist_over = dist[i] >= 0;
        bool dist_under = dist[i] < 0;
        // When dist_over , distcurrent == the closest positive value to 0
        // When dist_under, distcurrent == the closest negative value to 0
        distcurrent = dist_over && ((dist[i] < distcurrent) || (distcurrent < 0)) ? dist[i] :
                      dist_under && (dist[i] > distcurrent) ? dist[i] : distcurrent;
    }
    return distcurrent - caps.r;
}

// OBB caps tests accelerated
host_fn real distmin_vec_acc(OBB OBB, capsule caps) {
    segment seg;
    seg.p1 = caps.p1;
    seg.p2 = caps.p2;

    Mat3 Rtrans = transpose(OBB.R);

    Vec3 p1 = Rtrans * (seg.p1 - OBB.c);
    Vec3 p2 = Rtrans * (seg.p2 - OBB.c);

    Vec3 vec = p2 - p1;
    Vec3 point = p2;

    if (vec.z < 0) {
        vec = -vec;
        point = p1;
    }

    real xmin = - OBB.e[0];
    real xmax = + OBB.e[0];
    real ymin = - OBB.e[1];
    real ymax = + OBB.e[1];
    real zmin = - OBB.e[2];
    real zmax = + OBB.e[2];

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
    surf face[12];
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

    segment active_edges[24];
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
            bool is_inside = pointInSurf(face[i].d1, face[i].d2, face[i].p, point);
            if (is_inside)
                return normaldist[i] - caps.r;

            segment seg_face[4];
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
        return max_normal_dist - caps.r;
    }

    real dist_min = INF_REAL;
    real distcurrent;
    for (int i = 0; i < n_active_edges; i++) {
        distcurrent = distmin(active_edges[i], point);
        dist_min = distcurrent < dist_min ? distcurrent : dist_min;
    }
    return dist_min - caps.r;
}

// OBB caps tests using hulls (NOT WORKING)
host_fn bool comparePoints(const Vec3& a, const Vec3& b) {
    return (a.x < b.x) || (a.x == b.x && a.y < b.y) || (a.x == b.x && a.y == b.y && a.z < b.z);
}

host_fn bool sort_Vec3(const Vec3& a, const Vec3&b) {
    return (dot(a, a) > dot(b, b));
}

host_fn bool isLeftTurn(const Vec3& p, const Vec3& q, const Vec3& r) {
    return (q.x - p.x) * (r.y - p.y) - (r.x - p.x) * (q.y - p.y) > 0;
}

host_fn real distmin_hull(OBB OBB, capsule caps) {
    Mat3 Rtrans = transpose(OBB.R);

    segment seg;
    seg.p1 = Rtrans * (caps.p1 - OBB.c);
    seg.p2 = Rtrans * (caps.p2 - OBB.c);

    real xmin = - OBB.e[0];
    real xmax = + OBB.e[0];
    real ymin = - OBB.e[1];
    real ymax = + OBB.e[1];
    real zmin = - OBB.e[2];
    real zmax = + OBB.e[2];

    auto point = ptint(seg, {0, 0, 0});

    Vec3 signe;
    signe.x = point.x > 0 ? 1.0 : -1.0;
    signe.y = point.y > 0 ? 1.0 : -1.0;
    signe.z = point.z > 0 ? 1.0 : -1.0;

    //todo: use copysign function
    Vec3 p[5];
    Vec3 p1 = {signe.x*xmax, signe.y*ymax, signe.z*zmax};

    real t[6];
    t[0] = (seg.p1.x - p1.x) / (seg.p1.x - seg.p2.x); // x+ plane
    t[1] = (seg.p1.x + p1.x) / (seg.p1.x - seg.p2.x); // x- plane
    t[2] = (seg.p1.y - p1.y) / (seg.p1.y - seg.p2.y); // y+ plane
    t[3] = (seg.p1.y + p1.y) / (seg.p1.y - seg.p2.y); // y- plane
    t[4] = (seg.p1.z - p1.z) / (seg.p1.z - seg.p2.z); // z+ plane
    t[5] = (seg.p1.z + p1.z) / (seg.p1.z - seg.p2.z); // z- plane

    Vec3 ab = seg.p2 - seg.p1;
    Vec3 inter_x1 = seg.p1 + t[0] * ab;
    Vec3 inter_x2 = seg.p1 + t[1] * ab;
    Vec3 inter_y1 = seg.p1 + t[2] * ab;
    Vec3 inter_y2 = seg.p1 + t[3] * ab;
    Vec3 inter_z1 = seg.p1 + t[4] * ab;
    Vec3 inter_z2 = seg.p1 + t[5] * ab;

    bool pierce_x1 = (t[0] < 1 && t[0] > 0) && (inter_x1.y >= ymin && inter_x1.y <= ymax) && (inter_x1.z >= zmin && inter_x1.z <= zmax);
    bool pierce_x2 = (t[1] < 1 && t[1] > 0) && (inter_x2.y >= ymin && inter_x2.y <= ymax) && (inter_x2.z >= zmin && inter_x2.z <= zmax);
    bool pierce_y1 = (t[2] < 1 && t[2] > 0) && (inter_y1.x >= xmin && inter_y1.x <= xmax) && (inter_y1.z >= zmin && inter_y1.z <= zmax);
    bool pierce_y2 = (t[3] < 1 && t[3] > 0) && (inter_y2.x >= xmin && inter_y2.x <= xmax) && (inter_y2.z >= zmin && inter_y2.z <= zmax);
    bool pierce_z1 = (t[4] < 1 && t[4] > 0) && (inter_z1.y >= ymin && inter_z1.y <= ymax) && (inter_z1.x >= xmin && inter_z1.x <= xmax);
    bool pierce_z2 = (t[5] < 1 && t[5] > 0) && (inter_z2.y >= ymin && inter_z2.y <= ymax) && (inter_z2.x >= xmin && inter_z2.x <= xmax);

    Vec3 p1_OBB = {clamp(seg.p1.x, xmin, xmax), clamp(seg.p1.y, ymin, ymax), clamp(seg.p1.z, zmin, zmax)};
    Vec3 p2_OBB = {clamp(seg.p2.x, xmin, xmax), clamp(seg.p2.y, ymin, ymax), clamp(seg.p2.z, zmin, zmax)};
    bool p1_is_inside = (p1_OBB == seg.p1);
    bool p2_is_inside = (p2_OBB == seg.p2);

    bool pierce = pierce_x1 || pierce_x2 || pierce_y1 || pierce_y2 || pierce_z1 ||pierce_z2 || p1_is_inside || p2_is_inside;

    if (pierce) {
        // Define the OBB points in the plan where its normal is the capsule inner segment
        Vec3 OBB_point[8];
        OBB_point[0] = {xmin, ymin, zmin};
        OBB_point[1] = {xmin, ymin, zmax};
        OBB_point[2] = {xmin, ymax, zmin};
        OBB_point[3] = {xmin, ymax, zmax};
        OBB_point[4] = {xmax, ymin, zmin};
        OBB_point[5] = {xmax, ymin, zmax};
        OBB_point[6] = {xmax, ymax, zmin};
        OBB_point[7] = {xmax, ymax, zmax};

        // std::vector<Vec3> points; // Alternate implementation
        std::vector<Vec3> hull_points;
        hull_points.reserve(8);
        real t_OBB;
        t_OBB = dot(ab, - seg.p1) / dot(ab, ab);
        Vec3 OBB_center = - t_OBB*ab;
        for (int i = 0; i < sizeof(OBB_point)/sizeof(Vec3); i++) {
            t_OBB = dot(ab, OBB_point[i] - seg.p1) / dot(ab, ab);
            // points.push_back(OBB_point[i] - t_OBB*ab - OBB_center); // Alternate implementation
            hull_points.push_back(OBB_point[i] - t_OBB*ab - OBB_center);
        }

        // Using Andrew's algorithm for computing the convex hull
        // Sort points lexicographically
        // We can find which two points should be deleted
        std::sort(hull_points.begin(), hull_points.end(), sort_Vec3);
        hull_points.pop_back();
        hull_points.pop_back();

        for (int i = 0; i < hull_points.size(); i++) {
            hull_points[i] = hull_points[i] + OBB_center - seg.p1; // Set the origin as seg.p1
        }

        std::sort(hull_points.begin(), hull_points.end(), comparePoints);
        // Lower hull
        std::vector<Vec3> lower_hull;
        lower_hull.reserve(6);
        for (int i = 0; i < hull_points.size(); ++i) {
            while (lower_hull.size() >= 2 && !isLeftTurn(lower_hull[lower_hull.size() - 2], lower_hull.back(), hull_points[i]))
                lower_hull.pop_back();
            lower_hull.push_back(hull_points[i]);
        }
        // Upper hull
        std::vector<Vec3> upper_hull;
        upper_hull.reserve(6);
        for (int i = hull_points.size() - 1; i >= 0; --i) {
            while (upper_hull.size() >= 2 && !isLeftTurn(upper_hull[upper_hull.size() - 2], upper_hull.back(), hull_points[i]))
                upper_hull.pop_back();
            upper_hull.push_back(hull_points[i]);
        }
        // Concatenate lower and upper hulls (excluding duplicate endpoints)
        lower_hull.pop_back();  // Remove last point of lower hull as it's the same as the first point of upper hull
        lower_hull.insert(lower_hull.end(), upper_hull.begin(), upper_hull.end() - 1);

        // Alternate implementation
        // std::vector<Vec3> hull;
        // hull.reserve(6);
        // // std::vector<Vec3> temp.insert(temp.begin(), points.begin(), points.end());
        // real min_dot = INF_REAL;
        // int idx;
        // Vec3 current_point = points.back();   // Initialize at last point
        // hull.push_back(current_point);
        // points.pop_back();
        // // For every point, we will add the one which is closest to it next in (hull) and delete it from temp.
        // for (int i = 5; i > 1; i--) {
        //     min_dot = INF_REAL;
        //     for (int j = points.size() - 1; j > 0; j--) {
        //         real current_dot = dot(current_point - points[j], current_point - points[j]);
        //         if (current_dot < min_dot) {
        //             min_dot = current_dot;
        //             idx = j;
        //         }
        //     }
        //     hull.push_back(points[idx]);
        //     points.erase(points.begin() + idx);
        //     current_point = points.back();
        // }
        // hull.push_back(current_point);
        // sort(points.begin(), points.end(), comparePoints);

        // real current_dist;
        real current_dist_hull;
        real min_dist = -INF_REAL;
        real min_dist_hull = -INF_REAL;
        // Vec3 closest_point;
        Vec3 closest_point_hull;
        // segment seg_OBB;

        for (int i = 0; i < lower_hull.size(); i++) {
            closest_point_hull = closept_origin({lower_hull[i], lower_hull[(i + 1) % 6]});
            current_dist_hull = -dot(closest_point_hull, closest_point_hull);
            // Alternate implementation
            // closest_point = closept_origin({hull[i], hull[(i + 1) % 6]});
            // current_dist = - dot(closest_point, closest_point);
            // seg_OBB = {points[i], points[(i + 1) % 6]};
            // current_dist = - distmin(seg_OBB, seg.p1);
            // seg_OBB = {lower_hull[i], lower_hull[(i+1) %6]};
            // current_dist = - distmin(seg_OBB, seg.p1);

            // Compare distance with min_dist
            // min_dist = current_dist > min_dist ? current_dist : min_dist;
            min_dist_hull = std::max(current_dist_hull, min_dist_hull);
        }
        min_dist = - sqrt(-min_dist_hull);

        // The two endpoints could also get out of a plane.
        real t_current_plus;
        real t_current_minus;
        Vec3 test_point_plus;
        Vec3 test_point_minus;
        real dist_points;
        real dist_points_opp;

        for (int i = 0; i < 3; i++) {
            std::vector<real> t_active = {t[(2*i+2) % 6], t[(2*i+3) % 6], t[(2*i+4) % 6], t[(2*i+5) % 6]};

            sort(t_active.begin(), t_active.end());

            t_current_plus = clamp(t_active[1], 0, 1);
            t_current_minus = clamp(t_active[2], 0, 1);

            test_point_plus = (1 - t_current_plus) * seg.p1 + t_current_plus*seg.p2;
            test_point_minus = (1 - t_current_minus) * seg.p1 + t_current_minus*seg.p2;

            auto dist_plus = test_point_plus[i] + OBB_point[0][i];
            auto dist_minus = test_point_minus[i] + OBB_point[0][i];
            dist_points = dist_plus < dist_minus ? dist_plus : dist_minus;
            min_dist = dist_points < 0 && dist_points > min_dist ? dist_points : min_dist;

            auto dist_plus_opp = OBB_point[0][i] - test_point_plus[i];
            auto dist_minus_opp = OBB_point[0][i] - test_point_minus[i];
            dist_points_opp = dist_plus_opp < dist_minus_opp ? dist_plus_opp : dist_minus_opp;
            min_dist = dist_points_opp <= 0 && dist_points_opp > min_dist ? dist_points_opp : min_dist;
        }

        return min_dist - caps.r;
    }
    else {
        p[0] = p1;
        p[1] = {-signe.x*xmax, signe.y*ymax, signe.z*zmax};
        p[2] = {signe.x*xmax, -signe.y*ymax, signe.z*zmax};
        p[3] = {signe.x*xmax, signe.y*ymax, -signe.z*zmax};
        segment segOBB_01 = {p[0], p[1]};
        segment segOBB_02 = {p[0], p[2]};
        segment segOBB_03 = {p[0], p[3]};
        two_pts two_point_01 = closept(seg, segOBB_01);
        two_pts two_point_02 = closept(seg, segOBB_02);
        two_pts two_point_03 = closept(seg, segOBB_03);

        vector<two_pts> collision_points(5);
        collision_points = {two_point_01, two_point_02, two_point_03, {seg.p1, p1_OBB}, {seg.p2, p2_OBB}};

        real dist_min = INF_REAL;
        real dist;
        for (auto &two_point:collision_points) {
            dist = dot(two_point.p1 - two_point.p2, two_point.p1 - two_point.p2);
            dist_min = dist < dist_min ? dist : dist_min;
        }

        real dist_min_sqrt = sqrt(dist_min);
        return dist_min_sqrt - caps.r;
    }
}

// Objective function for optimization-based approaches
host_fn real dist_OBB_point(OBB OBB_test, Vec3 point) {
    Mat3 Rtrans = transpose(OBB_test.R);
 
    Vec3 point_OBB = Rtrans * (point - OBB_test.c);
 
    Vec3 proj = {clamp(point_OBB.x, -OBB_test.e.x, OBB_test.e.x), clamp(point_OBB.y, -OBB_test.e.y, OBB_test.e.y), clamp(point_OBB.z, -OBB_test.e.z, OBB_test.e.z)};
    Array dist_in = {abs(point.x) - OBB_test.e.x, abs(point.y) - OBB_test.e.y, abs(point.z) - OBB_test.e.z};
    real result_if_inside = array_max(dist_in);
    return result_if_inside > 0 ? norm(proj - point_OBB) : result_if_inside;
}

host_fn real dist_caps_point(capsule caps_test, Vec3 point) {
    sphere sph;
    sph.c = point;
    sph.r = 0;
    return distmin(caps_test, sph);
}

host_fn real dist_sph_point(sphere sph_test, Vec3 point) {
    return norm(point - sph_test.c) - sph_test.r;
}

host_fn Vec3 get_point(const Array& x, const Matrix &caps_list) {
    real t = x[0]*(caps_list.cols-1);
    int t_step = t;
    int t_step_plus1 = (t_step == (caps_list).cols-1) ? t_step : t_step + 1;
    real s = x[1]*(caps_list.rows/3-1);
    int s_step = s;
    int s_step_plus1 = (s_step == (caps_list).rows/3-1) ? s_step : s_step + 1;
 
    Vec3 p1_1 = {caps_list(3*s_step, t_step),             caps_list(3*s_step + 1, t_step),             caps_list(3*s_step + 2, t_step)};
    Vec3 p1_2 = {caps_list(3*s_step, t_step_plus1),       caps_list(3*s_step + 1, t_step_plus1),       caps_list(3*s_step + 2, t_step_plus1)};
    Vec3 p2_1 = {caps_list(3*s_step_plus1, t_step),       caps_list(3*s_step_plus1 + 1, t_step),       caps_list(3*s_step_plus1 + 2, t_step)};
    Vec3 p2_2 = {caps_list(3*s_step_plus1, t_step_plus1), caps_list(3*s_step_plus1 + 1, t_step_plus1), caps_list(3*s_step_plus1 + 2, t_step_plus1)};
 
    Vec3 p1 = (1 - (t - t_step)) * p1_1 + (t - t_step) * p1_2;
    Vec3 p2 = (1 - (t - t_step)) * p2_1 + (t - t_step) * p2_2;
 
    Vec3 p = (1 - (s - s_step)) * p1 + (s - s_step) * p2;
    return p;
}

host_fn real OBJ_function(const Array& x, const Matrix &caps_list, const objlist* world) {
    Vec3 p = get_point(x, caps_list);

    real distmin = INF_REAL;
    real current_dist = 0;
    for (auto OBB: (*world).OBBlist) {
        current_dist = dist_OBB_point(OBB, p);
        distmin = (distmin < 0) ? (current_dist > distmin && current_dist < 0 ? current_dist : distmin) :
                                  (current_dist < distmin ? current_dist : distmin);
    }
    for (auto caps: (*world).capslist) {
        current_dist = dist_caps_point(caps, p);
        distmin = (distmin < 0) ? (current_dist > distmin ? current_dist : distmin) :
                                  (current_dist < distmin ? current_dist : distmin);
    }
    for (auto sph: (*world).sphlist) {
        current_dist = dist_sph_point(sph, p);
        distmin = (distmin < 0) ? (current_dist > distmin ? current_dist : distmin) :
                                  (current_dist < distmin ? current_dist : distmin);
    }
 
    return distmin;
}

host_fn real OBJ_function(const Array& x, const Matrix &caps_list, const OBB &OBB_test) {
    Vec3 p = get_point(x, caps_list);
    return dist_OBB_point(OBB_test, p);
}

// ======================================
//            Add primitives
// ======================================

host_fn void add_OBB(Vec3 c, Vec3 e, Mat3 R, objlist* world) {
    OBB new_OBB;
    new_OBB.c = c;
    new_OBB.e = e;
    new_OBB.R = R;
    world->OBBlist.push_back(new_OBB);
}

host_fn void add_sphere(Vec3 c, real r, objlist* world) {
    sphere new_sph;
    new_sph.c = c;
    new_sph.r = r;
    world->sphlist.push_back(new_sph);
}

host_fn void add_cylinder(Vec3 p1, Vec3 p2, real r, objlist* world) {
    cylinder new_cyl;
    new_cyl.p1 = p1;
    new_cyl.p1 = p2;
    new_cyl.r = r;
    world->cyllist.push_back(new_cyl);
}

host_fn void add_capsule(Vec3 p1, Vec3 p2, real r, objlist* world) {
    capsule new_caps;
    new_caps.p1 = p1;
    new_caps.p1 = p2;
    new_caps.r = r;
    world->capslist.push_back(new_caps);
}

host_fn std::vector<real> test_collision(capslist* caps_list, objlist* world, int n_var) {
    std::vector<real> dist_min(n_var, INF_REAL);
    real dist;

    for (int c = 0; c < caps_list->caps.size(); c++) {
        // OBB collisions
        for (int i = 0; i < world->OBBlist.size(); i++) {
            dist = distmin_vec_acc(world->OBBlist[i], caps_list->caps[c]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }

        // Capsule collisions
        for (int i = 0; i < world->capslist.size(); i++) {
            dist = distmin(caps_list->caps[c], world->capslist[i]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }

        // Cylinder collisions
        for (int i = 0; i < world->cyllist.size(); i++) {
            dist = distmin(caps_list->caps[c], world->cyllist[i]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }

        // Sphere collisions
        for (int i = 0; i < world->sphlist.size(); i++) {
            dist = distmin(caps_list->caps[c], world->sphlist[i]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }
    }

    // // Self-collision
    // // note : this does not take into account the collision between two subsequent links of a robot,
    // // which means that angle constraints must be put on the articulations.
    // real n_pts = size(robot.pts);
    // real n_link = n_pts - 1;
    // capsule link[n_link];

    // for (int i = 0; i < n_link; i++) {
    //     link[i].p1 = robot.pts[i];
    //     link[i].p2 = robot.pts[i+1];
    //     link[i].r = robot.r[i];
    // }

    // for (int i = 0; i < n_link; i++) {
    //     for (int j = 2 + i; j < n_link; j++) {
    //         dist = distmin(link[i], link[i+j]);
    //         for (int k = 0; k < n_var; k++){
    //             if (dist < dist_min[j]) {
    //                 for (int l = n_var; l > k; l--) {
    //                     dist_min[l] = dist_min[l-1];
    //                 }
    //                 dist_min[k] = dist;
    //                 break;
    //             }
    //         }
    //     }
    // }
    return dist_min;
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

struct EPA_hull {
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;
    Vec3 n;
};

host_fn Vec3 GJK_get_support(std::vector<Vec3> vertices, Vec3 direction) {
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

host_fn void GJK_solve_simplex2_Ericson(Simplex* simplex) {
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

host_fn void GJK_solve_simplex3_Ericson(Simplex* simplex) {
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

host_fn int PointOutsideOfPlane(Vec3 p, Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
    // The last input (d) is the one that will be tested
    real signp = dot(p - a, cross(b - a, c - a)); // [AP AB AC]
    real signd = dot(d - a, cross(b - a, c - a)); // [AD AB AC]
    // Points on opposite sides if expression signs are opposite
    return signp * signd < 0.0f;
}

host_fn void GJK_solve_simplex4_Ericson(Simplex* simplex) {
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
    if (PointOutsideOfPlane(p, b, c, d, a)) {
        GJK_solve_simplex3_Ericson(&simplex_temp);
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
    if (PointOutsideOfPlane(p, a, c, d, b)) {
        GJK_solve_simplex3_Ericson(&simplex_temp);
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
    if (PointOutsideOfPlane(p, a, b, d, c)) {
        GJK_solve_simplex3_Ericson(&simplex_temp);
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
    if (PointOutsideOfPlane(p, a, b, c, d)) {
        GJK_solve_simplex3_Ericson(&simplex_temp);
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

        GJK_solve_simplex3_Ericson(simplex);
    }
    return;
}

host_fn bool pointInFace(EPA_hull face) {
    // u=P2−P1
    Vec3 u = face.p2 - face.p1;
    // v=P3−P1
    Vec3 v = face.p3 - face.p1;
    // w=−P1
    Vec3 w = - face.p1;
    // Barycentric coordinates of the projection P′of P onto T:
    // γ=[(u×w)⋅n]/n²
    real gamma = dot(cross(u, w), face.n);
    // β=[(w×v)⋅n]/n²
    real beta = dot(cross(w, v), face.n);
    real alpha = 1 - gamma - beta;
    // The point P′ lies inside T if:
    return ((0 <= alpha) && (alpha <= 1) && (0 <= beta) && (beta <= 1) && (0 <= gamma) && (gamma <= 1));
}

host_fn real distmin_origin(EPA_hull face) {
    Vec3 a = face.p1;
    Vec3 b = face.p2;
    Vec3 c = face.p3;

    Vec3 ab = b - a;
    Vec3 ac = c - a;
    Vec3 bc = c - b;

    // Compute parametric position s for projection P’ of P on AB,
    // P’ = A + s*AB, s = snom/(snom+sdenom)
    float snom = dot(- a, ab), sdenom = dot(- b, a - b);

    // Compute parametric position t for projection P’ of P on AC,
    // P’ = A + t*AC, s = tnom/(tnom+tdenom)
    float tnom = dot(- a, ac), tdenom = dot(- c, a - c);

    if (snom <= 0.0f && tnom <= 0.0f)
        return norm(a); // Vertex region early out

    // Compute parametric position u for projection P’ of P on BC,
    // P’ = B + u*BC, u = unom/(unom+udenom)
    float unom = dot(- b, bc), udenom = dot(- c, b - c);
    if (sdenom <= 0.0f && unom <= 0.0f)
        return norm(b); // Vertex region early out

    if (tdenom <= 0.0f && udenom <= 0.0f)
        return norm(c); // Vertex region early out

    // P is outside (or on) AB if the triple scalar product [N PA PB] <= 0
    // Vec3 n = face.n;
    Vec3 n = cross(b - a, c - a);
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

// (Function is BROKEN - infinite loop)
host_fn real solve_EPA_algorithm(Simplex simplex, std::vector<Vec3> v1, std::vector<Vec3> v2) {
    // create a face vector that has three points and a normal
    std::vector<EPA_hull> faces;
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
        Vec3 support1 = GJK_get_support(v1, -faces[idx].n);
        Vec3 support2 = GJK_get_support(v2, faces[idx].n);
        Vec3 p = support2 - support1;

        // If the vertex does not expand the polytope in the direction of the normal, the minimum distance
        // is with the closest face (unchanged). Compute and return.
        dist = dot(p - faces[idx].p1, faces[idx].n);
        if (dist < 1e-3) {
            break;
        }

        // Get all the faces which are "seen" by the new point, add them to deleted_faces and delete them.
        std::vector<EPA_hull> deleted_faces;
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
        two_pts current_edge[3];
        std::vector<two_pts> loose_edges;
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
        EPA_hull new_face;
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
            // if ((abs(check1 - check2) >= 1e-2 || abs(check1 - check3) >= 1e-2))
            //     real test = 1;
            new_face = {loose_edges[i].p1, loose_edges[i].p2, p, n};
            faces.push_back(new_face);
        }
    }
    return -abs(min_dist);
}

// (Function is BROKEN - infinite loop in EPA) Tests 2 sets of points using GJK
host_fn real general_GJK(std::vector<Vec3> Set_1, std::vector<Vec3> Set_2) {
    // This version of the GJK algorithm is implemented from the basic algorithm described in Collision Detection
    // manual by Ericson.

    // 1. Initializing simplex to a point from a random direction
    Vec3 V = Set_1[0] - Set_2[0];
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

        // 3. If P is the origin itself, the origin is clearly contained in the Minkowski difference of A and B.
        // Stop and return A and B as intersecting.
        if (dot(simplex.P, simplex.P) < 1e-6) {
            real dist = solve_EPA_algorithm(simplex, Set_1, Set_2);
            return dist;
        }

        // 4. Reduce Q to the smallest subset Q' of Q such that P is still in CH(Q). That is, remove any points
        // from Q not determining the subsimplex of Q in which P lies.
        // (This is done automatically in GJK_solve_simplex functions)

        // 5. Find the next supporting point in direction -P
        Vec3 a1 = GJK_get_support(Set_1, simplex.P);
        Vec3 a2 = GJK_get_support(Set_2, -simplex.P);
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
host_fn real GJK_OBB_caps(capsule caps, OBB box) {
    // Initialization of the eight OBB points
    Vec3 size_x_org = {box.e.x, 0, 0};
    Vec3 size_y_org = {0, box.e.y, 0};
    Vec3 size_z_org = {0, 0, box.e.z};
    Vec3 size_x = box.R*size_x_org;
    Vec3 size_y = box.R*size_y_org;
    Vec3 size_z = box.R*size_z_org;

    std::vector<Vec3> v1(2);
    v1[0] = caps.p1;
    v1[1] = caps.p2;

    std::vector<Vec3> v2(8);
    v2[0] = box.c + size_x + size_y + size_z;
    v2[1] = box.c + size_x + size_y - size_z;
    v2[2] = box.c + size_x - size_y + size_z;
    v2[3] = box.c + size_x - size_y - size_z;
    v2[4] = box.c - size_x + size_y + size_z;
    v2[5] = box.c - size_x + size_y - size_z;
    v2[6] = box.c - size_x - size_y + size_z;
    v2[7] = box.c - size_x - size_y - size_z;

    real dist = general_GJK(v1, v2) - caps.r;
    return dist;
}

host_fn std::vector<real> test_collision_GJK(capslist* caps_list, objlist* world, int n_var) {
    std::vector<real> dist_min(n_var, INF_REAL);
    real dist;

    for (int c = 0; c < caps_list->caps.size(); c++) {
        // OBB collisions
        for (int i = 0; i < world->OBBlist.size(); i++) {
            dist = GJK_OBB_caps(caps_list->caps[c], world->OBBlist[i]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }

        // Capsule collisions
        for (int i = 0; i < world->capslist.size(); i++) {
            dist = distmin(caps_list->caps[c], world->capslist[i]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }

        // Cylinder collisions
        for (int i = 0; i < world->cyllist.size(); i++) {
            dist = distmin(caps_list->caps[c], world->cyllist[i]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }

        // Sphere collisions
        for (int i = 0; i < world->sphlist.size(); i++) {
            dist = distmin(caps_list->caps[c], world->sphlist[i]);
            for (int j = 0; j < n_var; j++) {
                if (dist < dist_min[j]) {
                    for (int k = n_var - 1; k > j; k--) {
                        dist_min[k] = dist_min[k-1];
                    }
                    dist_min[j] = dist;
                    break;
                }
            }
        }
    }

    // // Self-collision
    // // note : this does not take into account the collision between two subsequent links of a robot,
    // // which means that angle constraints must be put on the articulations.
    // real n_pts = size(robot.pts);
    // real n_link = n_pts - 1;
    // capsule link[n_link];

    // for (int i = 0; i < n_link; i++) {
    //     link[i].p1 = robot.pts[i];
    //     link[i].p2 = robot.pts[i+1];
    //     link[i].r = robot.r[i];
    // }

    // for (int i = 0; i < n_link; i++) {
    //     for (int j = 2 + i; j < n_link; j++) {
    //         dist = distmin(link[i], link[i+j]);
    //         for (int k = 0; k < n_var; k++){
    //             if (dist < dist_min[j]) {
    //                 for (int l = n_var; l > k; l--) {
    //                     dist_min[l] = dist_min[l-1];
    //                 }
    //                 dist_min[k] = dist;
    //                 break;
    //             }
    //         }
    //     }
    // }
    return dist_min;
}

#ifdef BLAST_ENABLE_TESTS

// capsule - sphere collision tests
struct collision_test_sph {
    sphere sph;
    capsule caps;
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
    capsule caps1;
    capsule caps2;
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
    cylinder cyl;
    capsule caps;
    real expected_dist;
};
collision_test_cyl test_cyl[] = {
    /*Test 1*/ { { { -13.46, -3.61, 189 }, { -12.46, -5.11, 190 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, -0.538346 },
    /*Test 2*/ { { { -13.41, -2.79, 189 }, { -15.16, -3.2, 190 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, 1.32728712 },
    /*Test 3*/ { { { -16.3, -3.97, 200.87 }, { -17.3, -2.47, 199.87 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, 1.53080644 },
    /*Test 4*/ { { { -7.25, -10.11, 174.29 }, { -7.86, -9.44, 176.14 }, 2 }, { { -14.66, -13.24, 198.94 }, { -11.16, -7.94, 181.77 }, 1 }, 5.52119334 },
    /*Test 5*/ { { { -14.88, -2.74, 198.21 }, { -15.49, -2.08, 200.06 }, 2 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, -0.44703889 },
};

// capsule - OBB collision tests
struct collision_test_box {
    OBB box;
    capsule caps;
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
    // OBB.R and caps.r changed
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

    // The following test does not pass for dist_min_new
    // /*Test34*/ { { { 0, 0, 0 }, { 1.8543682561036101, 0.23108269265489501, 0.11287955803364236 } , { 1, 0, 0, 0, 1, 0, 0, 0, 1 } }, { { 4.4589046322448933, 6.1459446090413632, -5.0723802406288092 }, { -4.5833194312111516, -3.5836095282891645, 1.8703880978400422 }, 0 }, 0.33293883592789086 },
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
    // real dist = general_GJK(v1, v2);

    for (auto t : test) {
        real dist = distmin(t.caps, t.sph);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_caps) {
        real dist = distmin(t.caps1, t.caps2);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_cyl) {
        real dist = distmin(t.caps, t.cyl);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_obb) {
        real dist = distmin_vec(t.box, t.caps);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    for (auto t : test_obb) {
        real dist = distmin_vec_acc(t.box, t.caps);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }

    // for (auto t : test_obb) {
    //     real dist = distmin_hull(t.box, t.caps);
    //     // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
    //     CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    // }

    // GJK Algorithm (new/Ericson version)
    // for (auto t : test_obb) {
    //     real dist = GJK_OBB_caps(t.caps, t.box);
    //     CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    // }
}

TEST_CASE("Collision method benchmarks", "[World]") {
    using namespace blast;

    // real TESTCOLL_EPSILON = 1e-2;

    BENCHMARK("box - capsule OBB test with vectors") {
        real dist;
        for (auto t : test_obb) {
            dist = distmin_vec(t.box, t.caps);
            // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
        }
        return dist;
    };

    BENCHMARK("box - capsule OBB test with vectors accelerated") {
        real dist;
        for (auto t : test_obb) {
            dist = distmin_vec_acc(t.box, t.caps);
            // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
        }
        return dist;
    };

    // BENCHMARK("box - capsule OBB test with hulls") {
    //     real dist;
    //     for (auto t : test_obb) {
    //         dist = distmin_hull(t.box, t.caps);
    //         // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
    //     }
    //     return dist;
    // };

    // BENCHMARK("box - capsule OBB test with GJK") {
    //     real dist;
    //     for (auto t : test_obb) {
    //         dist = GJK_OBB_caps(t.caps, t.box);
    //         // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
    //     }
    //     return dist;
    // };
}

TEST_CASE("Collision method comparison exhaustive (OBB-cpasules)", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

    vector<OBB> obb_list;
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

        // Generate n random OBB
        fill_random(obb_e, 2);
        obb_e = array_abs(obb_e);
        fill_random(obb_c, 5);
        obb_c = array_abs(obb_c);
        obb_list.push_back({{obb_c[0], obb_c[1], obb_c[2]}, {obb_e[0], obb_e[1], obb_e[2]}, R});

        // Generate n random caps
        fill_random(seg_p1, 7);
        fill_random(seg_p2, 7);
        auto r = 0.0;
        caps_list.push_back({{seg_p1[0], seg_p1[1], seg_p1[2]}, {seg_p2[0], seg_p2[1], seg_p2[2]}, r});
    }

    vector<capsule> caps_failed;
    vector<OBB> obb_failed;
    for (int cap = 0; cap < caps_list.size(); cap++) {
        // OBB collisions
        for (int i = 0; i < obb_list.size(); i++) {
            auto dist_min_vec = distmin_vec(obb_list[i], caps_list[cap]);
            auto dist_min_vector_acc = distmin_vec_acc(obb_list[i], caps_list[cap]);
            // auto dist_min_new = distmin_hull(obb_list[i], caps_list[cap]);
            // auto dist_min_gjk = GJK_OBB_caps(caps_list[cap], obb_list[i]);

            CHECK(abs(dist_min_vec - dist_min_vector_acc) < TESTCOLL_EPSILON);
            if (abs(dist_min_vec - dist_min_vector_acc) > TESTCOLL_EPSILON) {
                // save and test caps and obb in future
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
    //     auto dist_min_vec = distmin_vec(obb_failed[i], caps_failed[i]);
    //     auto dist_min_new = distmin_hull(obb_failed[i], caps_failed[i]);
    //     auto dist_min_gjk = GJK_OBB_caps(caps_failed[i], obb_failed[i]);
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
