#pragma once
#include <iostream>
#include <vector>
#include "blast.hpp"

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

struct EPA_dist_norm {
    Vec3 normal;
    real distance;
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

struct EPA_simplex {
    std::vector<triangle> triangles;
};

struct edge {
    real distance;
    Vec3 normal;
    int index;
};

struct closeface {
    real distance;
    Vec3 normal;
    int index;
};

struct three_pts {
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;
};

struct bool_simplex {
    Vec3 a;
    Vec3 b;
    Vec3 c;
    Vec3 d;
    int count;
};

struct bool_Vec3 {
    bool intersection;
    Vec3 direction;
};

struct bool_real {
    bool intersection;
    real distance;
};

struct bool_simplex_2stdvec {
    bool intersection;
    bool_simplex simplex;
    std::vector<Vec3> v1;
    std::vector<Vec3> v2;
};

host_fn Vec3 normalize(Vec3 a) {
    return (1/norm(a))*a;
}

// A paper in 2005 suggested an elegant solution to the problem of the point in the triangle. This is this solution (adapted from https://gamedev.stackexchange.com/questions/28781/easy-way-to-project-point-onto-triangle-or-plane)
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

host_fn EPA_dist_norm distmin_origin(triangle tri) {
    EPA_dist_norm result;

    Vec3 o = tri.p1;
    Vec3 v1 = tri.p2 - tri.p1;
    Vec3 v2 = tri.p3 - tri.p1;
    Vec3 n = cross(v1, v2);
    if (dot(n, o) < 0)
        n = -n;

    real distmin = INF_REAL;
    if (pointInTriangle(v1, v2, o, {0, 0, 0})) {
        Vec3 n_unit = normalize(n);
        real d = dot(o, n_unit);
        result.normal = n_unit;
        result.distance = d;
        return result;
    }

    segment seg[3];

    seg[0].p1 = tri.p1;
    seg[0].p2 = tri.p2;
    seg[1].p1 = tri.p2;
    seg[1].p2 = tri.p3;
    seg[2].p1 = tri.p3;
    seg[2].p2 = tri.p1;

    for (int i = 0; i < 3; i++) {
        Vec3 d = closept_origin(seg[i]);
        if (dot(d, d) < distmin*distmin) {
            distmin = norm(d);
            result.normal = d;
        }
    }
    result.distance = distmin;
    return result;
}

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
        real val_d1_sq = dot(surf.d1, surf.d1);
        real val_d2_sq = dot(surf.d2, surf.d2);
        Vec3 val_direction = point - surf.p;

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
    else {
        real normaldist = dot(point - surf.p, n_unit);
        Vec3 direction = (point - normaldist*n_unit) - surf.p;

        Vec3 dy = cross(n_unit, surf.d1);
        Vec3 dy_unit = 1 / norm(dy) * dy;

        real t2 = dot(dy, direction) / dot(dy, surf.d2);
        real t1 = (dot(surf.d1, direction) - t2*dot(surf.d1, surf.d2)) / dot(surf.d1, surf.d1);

        // If the point projects on the face
        if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
            return normaldist;
        }

        Vec3 empty3 = {0, 0, 0};

        // If the point is closest to a corner of the shape
        if ((t1 <= 0 || t1 >= 1) && (t2 <= 0 || t2 >= 1)) {
            Vec3 testpoint = t1 <= 0 ? (t2 <= 0 ? empty3 :  surf.d2) : (t2 <= 0 ? surf.d1 : surf.d1 + surf.d2);
            real distance = sqrt(dot(direction - testpoint, direction - testpoint) + normaldist*normaldist);
            return normaldist < 0 ? -distance : distance;
        }

        real t1_clamped = clamp(t1, 0, 1);
        real t2_clamped = clamp(t2, 0, 1);

        Vec3 Projection_vector = (t1 < 1 && t1 > 0) ? dy_unit : cross(1/norm(surf.d2)*surf.d2, n_unit);
        Projection_vector = dot(Projection_vector, direction) < 0 ? -Projection_vector : Projection_vector;
        Vec3 testpoint = (t1 < 0 || t2 < 0) ? empty3 : (surf.d1 + surf.d2);

        real dist_plan = dot(Projection_vector, direction - testpoint);
        Vec3 projected_point = direction - dist_plan*Projection_vector;

        Vec3 clamping_direction = (t1 < 0 || t2 < 0) ? ((t1 < 1 && t1 > 0) ? surf.d1 : surf.d2) : ((t1 < 1 && t1 > 0) ? -surf.d1 : -surf.d2);
        real t_current = dot(direction - testpoint, clamping_direction) / dot(clamping_direction, clamping_direction);
        real t_current_clamped = clamp(t_current, 0, 1);

        real difference_t_current = (t_current - t_current_clamped);
        real dist_current = difference_t_current*difference_t_current*dot(clamping_direction, clamping_direction);
        real result = sqrt(dist_current + dist_plan*dist_plan + normaldist*normaldist);

        real distance = (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) ? normaldist : result;
        return distance == normaldist ? normaldist : (normaldist < 0 ? -distance : distance);
    }
}

host_fn two_pts closept(segment seg1, segment seg2) {
    // From : Real-Time Detection Collision (Christer Ericson)

    // Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
    // S2(t)=P2+t*(Q2-P2), returning s and t.

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
host_fn real distmin(segment seg1, segment seg2) {
    two_pts close_pt = closept(seg1, seg2);
    return dot(close_pt.p1 - close_pt.p2, close_pt.p1 - close_pt.p2); // ATTENTION! CECI DONNE LE RESULTAT AU CARRE
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

        // // To find the point on the side of the circle : (This was an attempt at fixing the current mistake,
        // // which leads to only a very small error. Thus it was deemed ok and this attempt scrapped)
        // Vec3 caps_segment = caps.p2 - caps.p1;
        // Vec3 caps_segment_unit = 1 / norm(caps_segment) * caps_segment;
        // Vec3 cross_prod = cross(n, caps_segment);
        // //Vec3 cross_prod = cross(n_unit, caps_segment_unit);
        // Vec3 cross_prod_unit = 1 / norm(cross_prod) * cross_prod;
        // Vec3 pt_direction = cross(caps_segment, cross_prod);
        // // Vec3 pt_direction = cross(caps_segment_unit, cross_prod_unit);
        // Vec3 pt_direction_unit = 1 / norm(pt_direction) * pt_direction;
        // pt_direction = (dot(pt_direction, points.p1 - points.p2) > 0) ? pt_direction : -pt_direction;
        // pt_direction_unit = 1 / norm(pt_direction) * pt_direction;
        // Vec3 projected_direction = ClosestPtPointPlane(pt_direction_unit, face1);
        // Vec3 projected_direction_unit = 1 / norm(projected_direction - cent) * (projected_direction - cent);
        // real verification_direction = dot(projected_direction_unit, n);
        // Vec3 testpoint = cent + (projected_direction_unit)*cyl.r;
        // Vec3 point_caps = ptint(seg1, testpoint);
        // real dist_point = norm(testpoint - point_caps);
        // // real dist_point = distmin(seg1, testpoint);
        // // dist_point = (dot(point_caps - cent, testpoint - cent) < dot(testpoint - cent, testpoint - cent)) ? -dist_point : dist_point;
        // real distance = dist_point < testnormal ? dist_point - caps.r : testnormal - caps.r;
        // // return distance;

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
    segment seg1;
    segment seg2;
    seg1.p1 = caps1.p2;
    seg1.p2 = caps1.p1;
    seg2.p1 = caps2.p2;
    seg2.p2 = caps2.p1;
    real dist_seg_seg = /*two_segment_distance_sqr(caps1.p1, caps1.p2, caps2.p1, caps2.p2);*/ sqrt(distmin(seg1, seg2));
    return dist_seg_seg - caps1.r - caps2.r;
}

host_fn real distmin(OBB OBB, capsule caps) {
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
        bool cond1 = dist[i] >= 0;
        bool cond21 = (dist[i] - distcurrent < COLLISION_EPSILON);
        bool cond22 = (distcurrent < 0);
        bool cond2 = (cond21 || cond22);
        if (cond1 && cond2) {
            distcurrent = dist[i];
        }
        else if (dist[i] < 0 && dist[i] > distcurrent) {
            distcurrent = dist[i];
        }
    }
    return distcurrent - caps.r;
}

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
            dist = distmin(world->OBBlist[i], caps_list->caps[c]);
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
    int count;
    std::vector<Vec3> fixed;
};

struct gjkresult {
    bool intersection;
    ComplexSimplex final_simplex;
    Vec3 A_closept;
    Vec3 B_closept;
    real minimal_distance;
};

host_fn bool GJK_same_direction(Vec3 a, Vec3 b) {
    return dot(a, b) > 0;
}

host_fn real GJK_triangle_area_2d(real x1, real y1, real x2, real y2, real x3, real y3) {
    return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

host_fn Vec3 GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 c, Vec3 p) {
    Vec3 ba = b - a;
    Vec3 ca = c - a;
    Vec3 abc = cross(ba, ca);
    real nu;
    real nv;
    real ood;

    real x = abs(abc.x);
    real y = abs(abc.y);
    real z = abs(abc.z);

    if ((x >= y) && (x >= z)) {
        // x is the largest so project onto the yz plane
        nu = GJK_triangle_area_2d(p.y, p.z, b.y, b.z, c.y, c.z);
        nv = GJK_triangle_area_2d(p.y, p.z, c.y, c.z, a.y, a.z);
        ood = 1 / abc.x;
    }
    else if ((y >= x) && (y >= z)) {
        // y is the largest so project onto the xz plane
        nu = GJK_triangle_area_2d(p.x, p.z, b.x, b.z, c.x, c.z);
        nv = GJK_triangle_area_2d(p.x, p.z, c.x, c.z, a.x, a.z);
        ood = 1 / -abc.y;
    }
    else {
        // z is the largest so project onto the xy plane
        nu = GJK_triangle_area_2d(p.x, p.y, b.x, b.y, c.x, c.y);
        nv = GJK_triangle_area_2d(p.x, p.y, c.x, c.y, a.x, a.y);
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

host_fn two_pts GJK_get_local_points(ComplexSimplex simplex, Vec3 p) {
    two_pts loc;
    real barycentric;
    Vec3 barycentric3;

    switch (simplex.count) {
    case 1:
        loc.p1 = simplex.a1;
        loc.p2 = simplex.a2;
        break;

    case 2:
        barycentric = GJK_convert_barycentric(simplex.a, simplex.b, p);
        loc.p1 = GJK_convert_cartesian(simplex.a1, simplex.b1, barycentric);
        loc.p2 = GJK_convert_cartesian(simplex.a2, simplex.b2, barycentric);
        break;

    case 3:
        barycentric3 = GJK_convert_barycentric(simplex.a, simplex.b, simplex.c, p);
        loc.p1 = GJK_convert_cartesian(simplex.a1, simplex.b1, simplex.c1, barycentric3);
        loc.p2 = GJK_convert_cartesian(simplex.a2, simplex.b2, simplex.c2, barycentric3);
        Assert(isnan(loc.p1.x) == false);
        break;

    default:
        Assert (false/*"invalid simplex dimension for nearest feature resolution."*/);
    }
    Assert(isnan(loc.p1.x) == false);
    return loc;
}

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

host_fn two_pts GJK_solve_simplex2(ComplexSimplex& simplex) {
    Vec3 ab = simplex.b - simplex.a;
    Vec3 ao = -simplex.a;

    const auto d = dot(ab, ao);
    simplex.count = d > 0 ? 2 : 1;
    ao = d > 0 ? cross(cross(ab, ao), ab) : ao;
    ab = d > 0 ? simplex.a + d / dot(ab, ab) * ab : simplex.a;

    return {ab, ao};
}

host_fn two_pts GJK_solve_simplex3(ComplexSimplex& simplex) {
    Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
    Vec3 ac = simplex.c - simplex.a;
    Vec3 ao = -simplex.a;

    if (GJK_same_direction(cross(abc, ac), ao)) {
        // the origin is in the direction of the triangle normal
        if (GJK_same_direction(ac, ao)) {
            // the origin is nearest to the line ac
            // simplex c remains as c
            // simplex.b_all = simplex.c_all;
            simplex.b = simplex.c;
            simplex.b1 = simplex.c1;
            simplex.b2 = simplex.c2;
            simplex.count = 2;
            real t = dot(ao, ac) / dot(ac, ac);
            return {simplex.a + t * ac, cross(cross(ac, ao), ac)};
        }
        else {
            Vec3 ab = simplex.b - simplex.a;
            if (GJK_same_direction(ab, ao)) {
                // origin is closest to the line ab
                simplex.count = 2;
                real t = dot(ao, ab) / dot(ab, ab);
                return {simplex.a + t * ab, cross(cross(ab, ao), ab)};
            }
            else {
                // the origin is nearest to the point a
                simplex.count = 1;
                return {simplex.a, ao};
            }
        }
    }
    else {
        Vec3 ab = simplex.b - simplex.a;
        if (GJK_same_direction(cross(ab, abc), ao)) {
            if (GJK_same_direction(ab, ao)) {
                // the origin is closest to line ab
                simplex.count = 2;
                real t = dot(ao, ab) / dot(ab, ab);
                return {simplex.a + t * ab, cross(cross(ab, ao), ab)};
            }
            else {
                // the origin is nearest to the point a
                simplex.count = 1;
                return {simplex.a, ao};
            }
        }
        else {
            if (GJK_same_direction(abc, ao)) {
                // the origin is closest to the triangle abc
                Vec3 bo = -simplex.b;
                Vec3 co = -simplex.c;
                real d1 = dot(ab, ao);
                real d2 = dot(ac, ao);
                real d3 = dot(ab, bo);
                real d4 = dot(ac, bo);
                real d5 = dot(ab, co);
                real d6 = dot(ac, co);
                real va = d3 * d6 - d5 * d4;
                real vb = d5 * d2 - d1 * d6;
                real vc = d1 * d4 - d3 * d2;
                real denom = 1 / (va + vb + vc);
                real v = vb * denom;
                real w = vc * denom;
                return {simplex.a + ab * v + ac * w, abc};
            }
            else {
                // the origin is closest to the triangle acb (other side of the triangle)
                // simplex.b_all, simplex.c_all = swap(simplex.b_all, simplex.c_all);
                Vec3 sb = simplex.c;
                Vec3 sb1 = simplex.c1;
                Vec3 sb2 = simplex.c2;

                simplex.c = simplex.b;
                simplex.c1 = simplex.b1;
                simplex.c2 = simplex.b2;

                simplex.b = sb;
                simplex.b1 = sb1;
                simplex.b2 = sb2;

                Vec3 bo = -simplex.b;
                Vec3 co = -simplex.c;
                real d1 = dot(ac, ao);
                real d2 = dot(ab, ao);
                real d3 = dot(ac, co);
                real d4 = dot(ab, co);
                real d5 = dot(ac, bo);
                real d6 = dot(ab, bo);
                real va = d3 * d6 - d5 * d4;
                real vb = d5 * d2 - d1 * d6;
                real vc = d1 * d4 - d3 * d2;
                real denom = 1 / (va + vb + vc);
                real v = vb * denom;
                real w = vc * denom;
                return {simplex.a + ab * v + ac * w, -abc};
            }
        }
    }
}

host_fn two_pts GJK_solve_simplex4(ComplexSimplex& simplex) {
    Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
    Vec3 acd = cross(simplex.c - simplex.a, simplex.d - simplex.a);
    Vec3 adb = cross(simplex.d - simplex.a, simplex.b - simplex.a);
    Vec3 ao = -simplex.a;

    real abc_dir = GJK_same_direction(abc, ao);
    real acd_dir = GJK_same_direction(acd, ao);
    real adb_dir = GJK_same_direction(adb, ao);

    if (!abc_dir && !acd_dir && !adb_dir) {
        // the origin is inside the simplex
        return {{ 0, 0, 0 }, {0, 0, 0}};
    }
    if (abc_dir && !acd_dir && !adb_dir) {
        // the origin is near abc
        simplex.count = 3;
        return GJK_solve_simplex3(simplex);
    }
    if (!abc_dir && acd_dir && !adb_dir) {
        // the origin is near acd
        // simplex.b_all = simplex.c_all;
        simplex.b = simplex.c;
        simplex.b1 = simplex.c1;
        simplex.b2 = simplex.c2;
        // simplex.c_all = simplex.d_all;
        simplex.c = simplex.d;
        simplex.c1 = simplex.d1;
        simplex.c2 = simplex.d2;

        simplex.count = 3;
        return GJK_solve_simplex3(simplex);
    }
    if (!abc_dir && !acd_dir && adb_dir) {
        // the origin is near adb
        // simplex.c_all = simplex.b_all;
        simplex.c = simplex.b;
        simplex.c1 = simplex.b1;
        simplex.c2 = simplex.b2;

        // simplex.b_all = simplex.d_all;
        simplex.b = simplex.d;
        simplex.b1 = simplex.d1;
        simplex.b2 = simplex.d2;

        simplex.count = 3;
        return GJK_solve_simplex3(simplex);
    }

    // the origin potentially falls on multiple triangles
    ComplexSimplex simplex_abc = simplex;
    simplex_abc.count = 3;

    // There may be a better way to do this.

    ComplexSimplex simplex_acd = simplex;
    //simplex_acd.b_all = simplex_acd.c_all;
    simplex_acd.b = simplex.c;
    simplex_acd.b1 = simplex.c1;
    simplex_acd.b2 = simplex.c2;
    //simplex_acd.c_all = simplex_acd.d_all;
    simplex_acd.c = simplex.d;
    simplex_acd.c1 = simplex.d1;
    simplex_acd.c2 = simplex.d2;

    simplex_acd.count = 3;

    ComplexSimplex simplex_adb = simplex;
    // simplex_adb.c_all = simplex_adb.b_all;
    simplex_adb.c = simplex.b;
    simplex_adb.c1 = simplex.b1;
    simplex_adb.c2 = simplex.b2;
    // simplex_adb.c_all = simplex_adb.d_all;
    simplex_adb.b = simplex.d;
    simplex_adb.b1 = simplex.d1;
    simplex_adb.b2 = simplex.d2;

    simplex_adb.count = 3;

    two_pts solved;
    solved = GJK_solve_simplex3(simplex_abc);
    Vec3 p_abc = solved.p1;
    Vec3 dir_abc = solved.p2;
    solved = GJK_solve_simplex3(simplex_adb);
    Vec3 p_adb = solved.p1;
    Vec3 dir_adb = solved.p2;
    solved = GJK_solve_simplex3(simplex_acd);
    Vec3 p_acd = solved.p1;
    Vec3 dir_acd = solved.p2;

    real abc_d2 = dot(p_abc, p_abc);
    real acd_d2 = dot(p_acd, p_acd);
    real adb_d2 = dot(p_adb, p_adb);

    if ((abc_d2 <= acd_d2) && (abc_d2 <= adb_d2)) {
        simplex = simplex_abc;
        return {p_abc, dir_abc};
    }
    else if ((acd_d2 <= abc_d2) && (acd_d2 <= adb_d2)) {
        simplex = simplex_acd;
        return {p_acd, dir_acd};
    }
    else if ((adb_d2 <= abc_d2) && (adb_d2 <= acd_d2)) {
        simplex = simplex_adb;
        return {p_adb, dir_adb};
    }

    Assert(false); // This souldn't be reachable

    // the origin isn't outside of any plane, so it's inside the tetrahedron
    return {{ 0, 0, 0 }, {0, 0, 0}};
}

// ======================================
//            EPA algorithm
// ======================================
// This version of the EPA algorithm is adapted from : https://dyn4j.org/2010/05/epa-expanding-polytope-algorithm/

host_fn closeface findClosestFace(std::vector<triangle> simplex) {
    closeface closest;
    closest.distance = INF_REAL;
    for (int i = 0; i < size(simplex); i++) {
        triangle current_triangle;
        current_triangle.p1 = simplex[i].p1;
        current_triangle.p2 = simplex[i].p2;
        current_triangle.p3 = simplex[i].p3;

        Vec3 vec1 = current_triangle.p2-current_triangle.p1;
        Vec3 vec2 = current_triangle.p3-current_triangle.p1;

        Vec3 origin = {0, 0, 0};
        Vec3 n = cross(vec1, vec2);
        if (norm(n) < COLLISION_EPSILON)
            break;
        if (dot(n, current_triangle.p1) < 0)
            n = - n;
        Vec3 n_unit = (1/norm(n))*n;

        EPA_dist_norm distance = distmin_origin(current_triangle);
        if (distance.distance < closest.distance) {
            closest.distance = distance.distance;
            closest.normal = distance.normal;
            closest.index = i;
        }
    }
    return closest;
}

host_fn bool check_same_point(Vec3 p1, Vec3 p2) {
    return ((p1.x - p2.x)*(p1.x - p2.x) < COLLISION_EPSILON && (p1.y - p2.y)*(p1.y - p2.y) < COLLISION_EPSILON && (p1.z - p2.z)*(p1.z - p2.z) < COLLISION_EPSILON);
}

host_fn bool check_same_triangle(triangle tri1, triangle tri2) {
    if (check_same_point(tri1.p1, tri2.p1) == 0 && check_same_point(tri1.p1, tri2.p2) == 0 && check_same_point(tri1.p1, tri2.p3) == 0)
        return 0;
    if (check_same_point(tri1.p2, tri2.p1) == 0 && check_same_point(tri1.p2, tri2.p2) == 0 && check_same_point(tri1.p2, tri2.p3) == 0)
        return 0;
    if (check_same_point(tri1.p3, tri2.p1) == 0 && check_same_point(tri1.p3, tri2.p2) == 0 && check_same_point(tri1.p3, tri2.p3) == 0)
        return 0;
    return 1;
}

host_fn real solve_EPA_algorithm(ComplexSimplex simplex, std::vector<Vec3> v1, std::vector<Vec3> v2) {
    // create a simplex that can take as many points as necessary (EPA_simplex)
    std::vector<triangle> s;
    s.push_back({simplex.a, simplex.b, simplex.c});
    s.push_back({simplex.a, simplex.b, simplex.d});
    s.push_back({simplex.a, simplex.c, simplex.d});
    s.push_back({simplex.b, simplex.c, simplex.d});

    // loop to find the collision information
    while (true) {
        // obtain the feature (edge for 2D) closest to the
        // origin on the Minkowski Difference
        closeface e = findClosestFace(s);
        // obtain a new support point in the direction of the edge normal
        Vec3 support1 = GJK_get_support(v1, -e.normal);
        Vec3 support2 = GJK_get_support(v2, e.normal);
        Vec3 p = support2 - support1;
        // check the distance from the origin to the edge against the
        // distance p is along e.normal
        real d = dot(p, e.normal);
        if (abs(d) - abs(e.distance) < 1e-2) {
            // the tolerance should be something positive close to zero (ex. 0.00001)

            // if the difference is less than the tolerance then we can
            // assume that we cannot expand the simplex any further and
            // we have our solution
            return -abs(d);
        }
        else {
            // we haven't reached the edge of the Minkowski Difference
            // so continue expanding by adding the new triangle to the simplex
            // in between the points that made the closest triangle. To do this,
            // we must delete the triangle which is currently the closest and create three
            // triangles that are inserted in its place.
            Vec3 p1 = s[e.index].p1;
            Vec3 p2 = s[e.index].p2;
            Vec3 p3 = s[e.index].p3;

            triangle new_tri1 = {p1, p2, p};
            triangle new_tri2 = {p1, p3, p};
            triangle new_tri3 = {p2, p3, p};

            // we delete the face which was closest before. it is now inside the polytope.
            s.erase(s.begin() + e.index);

            // in the case where two iterations return the same support point, it can be the case that two
            // faces are created on top of one another. In this case, the face is always inside the
            // polytope and it should therefore be deleted.
            for (int i = (int)size(s) - 1; i >= 0; i--) {
                if (check_same_triangle(s[i], new_tri1)) {
                    s.erase(s.begin() + i); // if the triangle was already in the list, then we must delete it as this face is now inside the simplex
                    break;
                }
                else if (i == 0) {
                    s.push_back(new_tri1); // if it is a new triangle, add it to the list
                    break;
                }
            }

            for (int i = (int)size(s) - 1; i >= 0; i--) {
                if (check_same_triangle(s[i], new_tri2)) {
                    s.erase(s.begin() + i);
                    break;
                }
                else if (i == 0) {
                    s.push_back(new_tri2); // if it is a new triangle, add it to the list
                    break;
                }
            }

            for (int i = (int)size(s) - 1; i >= 0; i--) {
                if (check_same_triangle(s[i], new_tri3)) {
                    s.erase(s.begin() + i);
                    break;
                }
                else if (i == 0) {
                    s.push_back(new_tri3); // if it is a new triangle, add it to the list
                    break;
                }
            }
        }
    }
}

host_fn real solve_EPA_algorithm_bool(bool_simplex simplex, std::vector<Vec3> v1, std::vector<Vec3> v2) {
    // create a simplex that can take as many points as necessary (EPA_simplex)
    std::vector<triangle> s;
    s.push_back({simplex.a, simplex.b, simplex.c});
    s.push_back({simplex.a, simplex.b, simplex.d});
    s.push_back({simplex.a, simplex.c, simplex.d});
    s.push_back({simplex.b, simplex.c, simplex.d});

    // loop to find the collision information
    while (true) {
        // obtain the feature (edge for 2D) closest to the
        // origin on the Minkowski Difference
        closeface e = findClosestFace(s);
        // obtain a new support point in the direction of the edge normal
        Vec3 support1 = GJK_get_support(v1, -e.normal);
        Vec3 support2 = GJK_get_support(v2, e.normal);
        Vec3 p = support2 - support1;
        // check the distance from the origin to the edge against the
        // distance p is along e.normal
        real d = dot(p, e.normal);
        if (abs(d) - abs(e.distance) < 1e-2) {
            // the tolerance should be something positive close to zero (ex. 0.00001)

            // if the difference is less than the tolerance then we can
            // assume that we cannot expand the simplex any further and
            // we have our solution
            return -abs(d);
        }
        else {
            // we haven't reached the edge of the Minkowski Difference
            // so continue expanding by adding the new triangle to the simplex
            // in between the points that made the closest triangle. To do this,
            // we must delete the triangle which is currently the closest and create three
            // triangles that are inserted in its place.
            Vec3 p1 = s[e.index].p1;
            Vec3 p2 = s[e.index].p2;
            Vec3 p3 = s[e.index].p3;

            triangle new_tri1 = {p1, p2, p};
            triangle new_tri2 = {p1, p3, p};
            triangle new_tri3 = {p2, p3, p};

            // we delete the face which was closest before. it is now inside the polytope.
            s.erase(s.begin() + e.index);

            // in the case where two iterations return the same support point, it can be the case that two
            // faces are created on top of one another. In this case, the face is always inside the
            // polytope and it should therefore be deleted.
            for (int i = (int)size(s) - 1; i >= 0; i--) {
                if (check_same_triangle(s[i], new_tri1)) {
                    s.erase(s.begin() + i); // if the triangle was already in the list, then we must delete it as this face is now inside the simplex
                    break;
                }
                else if (i == 0) {
                    s.push_back(new_tri1); // if it is a new triangle, add it to the list
                    break;
                }
            }

            for (int i = (int)size(s) - 1; i >= 0; i--) {
                if (check_same_triangle(s[i], new_tri2)) {
                    s.erase(s.begin() + i);
                    break;
                }
                else if (i == 0) {
                    s.push_back(new_tri2); // if it is a new triangle, add it to the list
                    break;
                }
            }

            for (int i = (int)size(s) - 1; i >= 0; i--) {
                if (check_same_triangle(s[i], new_tri3)) {
                    s.erase(s.begin() + i);
                    break;
                }
                else if (i == 0) {
                    s.push_back(new_tri3); // if it is a new triangle, add it to the list
                    break;
                }
            }
        }
    }
}

host_fn gjkresult GJK_solve_gjk_simple(capsule caps, OBB box) {
    ComplexSimplex simplex;
    gjkresult results;

    // Initialization of the eight OBB points
    Vec3 size_x_org = {box.e.x, 0, 0};
    Vec3 size_y_org = {0, box.e.y, 0};
    Vec3 size_z_org = {0, 0, box.e.z};
    Vec3 size_x = box.R*size_x_org;
    Vec3 size_y = box.R*size_y_org;
    Vec3 size_z = box.R*size_z_org;

    std::vector<Vec3> v1(8);
    v1[0] = box.c + size_x + size_y + size_z;
    v1[1] = box.c + size_x + size_y - size_z;
    v1[2] = box.c + size_x - size_y + size_z;
    v1[3] = box.c + size_x - size_y - size_z;
    v1[4] = box.c - size_x + size_y + size_z;
    v1[5] = box.c - size_x + size_y - size_z;
    v1[6] = box.c - size_x - size_y + size_z;
    v1[7] = box.c - size_x - size_y - size_z;

    std::vector<Vec3> v2(2);
    v2[0] = caps.p1;
    v2[1] = caps.p2;

    // The general code starts here
    Vec3 direction = v1[0] - v2[0];
    simplex.a1 = GJK_get_support(v1, -direction);
    simplex.a2 = GJK_get_support(v2, direction);
    simplex.a = simplex.a2 - simplex.a1; // todo: a1 and a2 necessary? Same for b1, b2, ...

    simplex.count = 1;

    two_pts solved;
    Vec3 p;
    real old_simplex_count;

    while (true) {
        old_simplex_count = simplex.count;

        if (simplex.count == 1) {
            p = simplex.a;
            direction = -simplex.a;
        }
        else {
            solved = (simplex.count == 2) ? GJK_solve_simplex2(simplex) : ((simplex.count == 3) ? GJK_solve_simplex3(simplex) : GJK_solve_simplex4(simplex));
            p = solved.p1;
            direction = solved.p2;
        }
        // switch (simplex.count) {
        // case 1:
        //     p = simplex.a;
        //     direction = -simplex.a;
        //     break;
        // case 2:
        //     solved = GJK_solve_simplex2(simplex);
        //     p = solved.p1;
        //     direction = solved.p2;
        //     break;
        // case 3:
        //     // if colinear --> find closest to origin and return
        //     solved = GJK_solve_simplex3(simplex);
        //     p = solved.p1;
        //     direction = solved.p2;
        //     break;
        // case 4:
        //     // if 4 points in same plane --> find closest to origin and return
        //     solved = GJK_solve_simplex4(simplex);
        //     p = solved.p1;
        //     direction = solved.p2;
        //     break;
        // default:
        //     // This should never be called
        //     Assert(false);
        // }

        if (dot(p, p) < COLLISION_EPSILON && old_simplex_count == 4) {
            results.intersection = true;
            results.final_simplex = simplex;
            results.A_closept = {};
            results.B_closept = {};
            results.minimal_distance = solve_EPA_algorithm(simplex, v1, v2) - caps.r;
            Assert(isnan(results.minimal_distance) == false);
            break;
        }

        Vec3 support1 = GJK_get_support(v1, -direction);
        Vec3 support2 = GJK_get_support(v2, direction);

        Vec3 support = support2 - support1;
        real goal_post = dot(p, direction);
        real current_highest = dot(support, direction);

        if (current_highest - goal_post <=  COLLISION_EPSILON || abs(dot(simplex.a, simplex.a) - dot(support, simplex.a)) <= COLLISION_EPSILON || (simplex.count >= 2 && dot(simplex.b, simplex.b) - dot(support, simplex.b) <= COLLISION_EPSILON) || (simplex.count >= 3 && abs(dot(simplex.c, simplex.c) - dot(support, simplex.c)) <= COLLISION_EPSILON) || (simplex.count == 4 && abs(dot(simplex.d, simplex.d) - dot(support, simplex.d)) <= COLLISION_EPSILON)) {
            if (simplex.count == 3 && norm(cross(simplex.b - simplex.a, simplex.c - simplex.a)) < 1e-9) {
                segment seg_test;
                // Calculate all points in minkowski difference
                // Make a list of all the points that have the same dot product as our point p (closest point)
                std::vector<three_pts> full_minkowski_list;
                std::vector<Vec3> pts_list;
                for (auto& vertex1 : v1) {
                    for (auto& vertex2 : v2) {
                        Vec3 new_pt = vertex2 - vertex1;
                        full_minkowski_list.push_back({vertex1, vertex2, new_pt});
                        if (abs(goal_post - dot(new_pt, direction)) < 1e-9)
                            pts_list.push_back(new_pt);
                    }
                }

                // Within this list, find the two points which are the most extreme
                // Create a segment with these two points and test segment-point with the origin
                Vec3 new_dir = simplex.b - simplex.a;
                real max_dot = -INF_REAL;
                real min_dot = INF_REAL;
                for (auto& pts : pts_list) {
                    real current_dot = dot(pts, new_dir);
                    if (current_dot > max_dot) {
                        max_dot = current_dot;
                        seg_test.p1 = pts;
                    }
                    else if (current_dot < min_dot) {
                        min_dot = current_dot;
                        seg_test.p2 = pts;
                    }
                }
                simplex.a = seg_test.p1;
                simplex.a1 = GJK_get_support(v1, -simplex.a);
                simplex.a2 = GJK_get_support(v2, simplex.a);

                simplex.b = seg_test.p2;
                simplex.b1 = GJK_get_support(v1, -simplex.b);
                simplex.b2 = GJK_get_support(v2, simplex.b);

                simplex.count = 2;
            }

            two_pts local = GJK_get_local_points(simplex, p);

            results.intersection = false;
            results.A_closept = local.p1;
            results.B_closept = local.p2;
            results.final_simplex = simplex;
            results.minimal_distance = norm(results.A_closept - results.B_closept) - caps.r;
            Assert(isnan(results.minimal_distance) == false);
            break;
        }

        Assert(simplex.count < 4 );/*"You cannot have a simplex that's a tetrahedron at this point."*/

        simplex.d = simplex.c;
        simplex.d1 = simplex.c1;
        simplex.d2 = simplex.c2;

        simplex.c = simplex.b;
        simplex.c1 = simplex.b1;
        simplex.c2 = simplex.b2;

        simplex.b = simplex.a;
        simplex.b1 = simplex.a1;
        simplex.b2 = simplex.a2;

        simplex.count += 1;

        simplex.a1 = support1;
        simplex.a2 = support2;
        simplex.a = support;
    }

    return results;
}

// ======================================
//      Boolean GJK algorithm
// ======================================

host_fn Vec3 GJK_find_direction2(bool_simplex& simplex) {
    Vec3 ab = simplex.b - simplex.a;
    Vec3 dir = -simplex.a;

    const auto d = dot(ab, dir);
    simplex.count = d > 0 ? 2 : 1;
    dir = d > 0 ? cross(cross(ab, dir), ab) : dir;

    return dir;
}

host_fn Vec3 GJK_find_direction3(bool_simplex& simplex) {
    Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
    Vec3 ac = simplex.c - simplex.a;
    Vec3 ao = -simplex.a;

    if (GJK_same_direction(cross(abc, ac), ao)) {
        // the origin is in the direction of the triangle normal
        if (GJK_same_direction(ac, ao)) {
            // the origin is nearest to the line ac
            simplex.b = simplex.c;
            simplex.count = 2;
            return cross(cross(ac, ao), ac);
        }
        else {
            Vec3 ab = simplex.b - simplex.a;
            if (GJK_same_direction(ab, ao)) {
                // origin is closest to the line ab
                simplex.count = 2;
                return cross(cross(ab, ao), ab);
            }
            else {
                // the origin is nearest to the point a
                simplex.count = 1;
                return ao;
            }
        }
    }
    else {
        Vec3 ab = simplex.b - simplex.a;
        if (GJK_same_direction(cross(ab, abc), ao)) {
            if (GJK_same_direction(ab, ao)) {
                // the origin is closest to line ab
                simplex.count = 2;
                return cross(cross(ab, ao), ab);
            }
            else {
                // the origin is nearest to the point a
                simplex.count = 1;
                return ao;
            }
        }
        else {
            if (GJK_same_direction(abc, ao)) {
                // the origin is closest to the triangle abc
                return abc;
            }
            else {
                // the origin is closest to the triangle acb (other side of the triangle)
                return -abc;
            }
        }
    }
}

host_fn bool_Vec3 GJK_find_direction4(bool_simplex& simplex) {
    Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
    Vec3 ao = -simplex.a;

    if (GJK_same_direction(abc, ao)) {
        // the origin is nearest to the triangle abc
        simplex.count = 3;
        return {false, GJK_find_direction3(simplex)};
    }

    Vec3 acd = cross(simplex.c - simplex.a, simplex.d - simplex.a);

    if (GJK_same_direction(acd, ao)) {
        // the origin is nearest to the triangle acd
        simplex.b = simplex.c;

        simplex.c = simplex.d;
        simplex.count = 3;
        return {false, GJK_find_direction3(simplex)};
    }

    Vec3 adb = cross(simplex.d - simplex.a, simplex.b - simplex.a);

    if (GJK_same_direction(adb, ao)) {
        // the origin is nearest to the triangle adb
        simplex.c = simplex.b;
        simplex.b = simplex.d;

        simplex.count = 3;
        return {false, GJK_find_direction3(simplex)};
    }

    Vec3 bdc = cross(simplex.d - simplex.b, simplex.c - simplex.b);
    Vec3 bo = -simplex.b;

    // we normally wouldn't have to check bdc, but since we're allowing cached simplexes we do
    if (GJK_same_direction(bdc, bo)) {
        simplex.a = simplex.b;
        simplex.b = simplex.d;
        simplex.count = 3;
        return {false, GJK_find_direction3(simplex)};
    }

    return {true, {0, 0, 0}};
}

host_fn Vec3 GJK_bool_get_support(std::vector<Vec3> Minkowski_list, Vec3 direction) {
    real largest_dot = -INF_REAL;
    real current_dot;
    Vec3 largest_vertex;

    for (auto& vertex : Minkowski_list) {
        current_dot = dot(vertex, direction);
        largest_vertex  = current_dot > largest_dot ? vertex : largest_vertex;
        largest_dot     = current_dot > largest_dot ? current_dot : largest_dot;
    }
    return largest_vertex;
}

host_fn real GJK_bool_simplexdistance(bool_simplex& simplex) {
    real distance;
    switch (simplex.count) {
    case 1:
        distance = norm(simplex.a);
        break;
    case 2: {
        segment seg = {simplex.a, simplex.b};
        distance = distmin(seg, {0, 0, 0});
        break;
    }
    case 3: {
        plane p;
        p.n = cross(simplex.c - simplex.a, simplex.b - simplex.a);
        p.p = simplex.a;
        distance = norm(ClosestPtPointPlane({0, 0, 0}, p));
        break;
    }
    default:
        // this should not be called
        Assert(false);
    }
    return distance;
}

// Adapted from jai version here : https://github.com/kujukuju/KodaPhysics/blob/master/src/Gjk.jai
host_fn bool GJK_bool_gen(std::vector<Vec3> v1, std::vector<Vec3> v2, real radius) {
    bool_simplex simplex;
    bool intersection;

    int full_size = (int)(size(v1) * size(v2));
    std::vector<Vec3> full_minkowski_list(full_size);

    // if (full_size > 32) {
    //     // The size does not warrant for the full computation of the minkowski difference
    //     simplex.a = v1[0] - v2[0];
    // }
    // else {
    // Calculate all points in minkowski difference and pick the point of lowest norm to start
    real low_norm = INF_REAL;
    for (auto& vertex1 : v1) {
        for (auto& vertex2 : v2) {
            Vec3 new_pt = vertex2 - vertex1;
            full_minkowski_list.push_back(new_pt);
            real current_norm_sq = dot(new_pt, new_pt);
            if (current_norm_sq < low_norm*low_norm) {
                low_norm = sqrt(current_norm_sq);
                simplex.a = new_pt;
            }
        }
    }
    // }

    Vec3 direction;
    simplex.count = 1;

    while (true) {
        switch (simplex.count) {
        case 1:
            direction = -simplex.a;
            break;
        case 2:
            direction = GJK_find_direction2(simplex);
            break;
        case 3:
            direction = GJK_find_direction3(simplex);
            break;
        case 4: {
            bool_Vec3 sol = GJK_find_direction4(simplex);
            if (sol.intersection == true) {
                intersection = true;
                goto end;
            }
            direction = sol.direction;
            break;
        }
        default:
            // This should never be called
            Assert(false);
        }

        Assert(simplex.count < 4 ); /*"You cannot have a simplex that's a tetrahedron at this point."*/

        if (GJK_bool_simplexdistance(simplex) - radius <= COLLISION_EPSILON) {
            intersection = true;
            break;
        }

        // Vec3 support = full_size < 32 ? GJK_get_support(v2, direction) - GJK_get_support(v1, -direction) : GJK_bool_get_support(full_minkowski_list, direction);
        // Vec3 support1 = GJK_get_support(v1, -direction);
        // Vec3 support2 = GJK_get_support(v2, direction);

        Vec3 support = GJK_bool_get_support(full_minkowski_list, direction); //support2 - support1;

        if (dot(support, direction) <= 0) {
            // That means that we did not cross the origin, therefore it is impossible that the
            // origin is inside the minkowski difference.
            intersection = false;
            break;
        }

        simplex.d = simplex.c;
        simplex.c = simplex.b;
        simplex.b = simplex.a;
        simplex.a = support;

        simplex.count += 1;
    }
end :
    return intersection;
}

host_fn bool GJK_bool(capsule caps, OBB box) {
    // Initialization of the eight OBB points
    Vec3 size_x_org = {box.e.x, 0, 0};
    Vec3 size_y_org = {0, box.e.y, 0};
    Vec3 size_z_org = {0, 0, box.e.z};
    Vec3 size_x = box.R*size_x_org;
    Vec3 size_y = box.R*size_y_org;
    Vec3 size_z = box.R*size_z_org;

    std::vector<Vec3> v1(8);
    v1[0] = box.c + size_x + size_y + size_z;
    v1[1] = box.c + size_x + size_y - size_z;
    v1[2] = box.c + size_x - size_y + size_z;
    v1[3] = box.c + size_x - size_y - size_z;
    v1[4] = box.c - size_x + size_y + size_z;
    v1[5] = box.c - size_x + size_y - size_z;
    v1[6] = box.c - size_x - size_y + size_z;
    v1[7] = box.c - size_x - size_y - size_z;

    std::vector<Vec3> v2(2);
    v2[0] = caps.p1;
    v2[1] = caps.p2;

    return GJK_bool_gen(v1, v2, caps.r);
}

// Boolean GJK with depth
// host_fn bool_simplex_2stdvec GJK_bool_gen_informations(std::vector<Vec3> v1, std::vector<Vec3> v2, real radius){
//     bool_simplex simplex;
//     bool_simplex_2stdvec results;
//     int full_size = size(v1) * size(v2);
//     std::vector<Vec3> full_minkowski_list(full_size);
//     real low_norm = INF_REAL;
//     for (auto& vertex1 : v1) {
//         for (auto& vertex2 : v2) {
//             Vec3 new_pt = vertex2 - vertex1;
//             full_minkowski_list.push_back(new_pt);
//             real current_norm_sq = dot(new_pt, new_pt);
//             if (current_norm_sq < low_norm*low_norm) {
//                 low_norm = sqrt(current_norm_sq);
//                 simplex.a = new_pt;
//             }
//         }
//     }
//     Vec3 direction;
//     simplex.count = 1;
//     while (true) {
//         if (simplex.count == 4) {
//             bool_Vec3 sol = GJK_find_direction4(simplex);
//             if (sol.intersection == true) {
//                 results.intersection = true;
//                 results.simplex = simplex;
//                 results.v1 = v1;
//                 results.v2 = v2;
//                 goto end;
//             }
//             direction = sol.direction;
//         }
//         else {
//             direction = simplex.count == 1 ? -simplex.a : (simplex.count == 2 ? GJK_find_direction2(simplex) : (GJK_find_direction3(simplex)));
//         }
//         Assert(simplex.count < 4 ); /*"You cannot have a simplex that's a tetrahedron at this point."*/
//         if (GJK_bool_simplexdistance(simplex) - radius <= COLLISION_EPSILON) {
//             results.intersection = true;
//             results.simplex = simplex;
//             results.v1 = v1;
//             results.v2 = v2;
//             break;
//         }
//         // Vec3 support = full_size < 32 ? GJK_get_support(v2, direction) - GJK_get_support(v1, -direction) : GJK_bool_get_support(full_minkowski_list, direction);
//         // Vec3 support1 = GJK_get_support(v1, -direction);
//         // Vec3 support2 = GJK_get_support(v2, direction);
//         Vec3 support = GJK_bool_get_support(full_minkowski_list, direction); //support2 - support1;
//         if (dot(support, direction) <= 0) {
//             // That means that we did not cross the origin, therefore it is impossible that the
//             // origin is inside the minkowski difference.
//             results.intersection = false;
//             results.simplex = simplex;
//             results.v1 = v1;
//             results.v2 = v2;
//             break;
//         }
//         simplex.d = simplex.c;
//         simplex.c = simplex.b;
//         simplex.b = simplex.a;
//         simplex.a = support;
//         simplex.count += 1;
//     }
//     end :
//     return results;
// }
// host_fn bool_simplex_2stdvec GJK_bool_informations(capsule caps, OBB box) {
//     // Initialization of the eight OBB points
//     Vec3 size_x_org = {box.e.x, 0, 0};
//     Vec3 size_y_org = {0, box.e.y, 0};
//     Vec3 size_z_org = {0, 0, box.e.z};
//     Vec3 size_x = box.R*size_x_org;
//     Vec3 size_y = box.R*size_y_org;
//     Vec3 size_z = box.R*size_z_org;
//     std::vector<Vec3> v1(8);
//     v1[0] = box.c + size_x + size_y + size_z;
//     v1[1] = box.c + size_x + size_y - size_z;
//     v1[2] = box.c + size_x - size_y + size_z;
//     v1[3] = box.c + size_x - size_y - size_z;
//     v1[4] = box.c - size_x + size_y + size_z;
//     v1[5] = box.c - size_x + size_y - size_z;
//     v1[6] = box.c - size_x - size_y + size_z;
//     v1[7] = box.c - size_x - size_y - size_z;
//     std::vector<Vec3> v2(2);
//     v2[0] = caps.p1;
//     v2[1] = caps.p2;
//     return GJK_bool_gen_informations(v1, v2, caps.r);
// }
// host_fn bool_real GJK_bool_dist_pen(capsule caps, OBB box) {
//     bool_simplex_2stdvec results = GJK_bool_informations(caps, box);
//     if (results.intersection) {
//         return {true, solve_EPA_algorithm_bool(results.simplex, results.v1, results.v2) - caps.r};
//     }
//     return {false, INF_REAL};
// }
// // ======================================
// //          méthode des points
// // ======================================
// // in progress
// host_fn real pts_distmin(OBB OBBtest, capsule caps) {
//     segment seg;
//     seg.p1 = caps.p1;
//     seg.p2 = caps.p2;
//     Mat3 Rtrans = transpose(OBBtest.R);
//     // Mat3 Rtrans = OBB.R;
//     Vec3 p1 = Rtrans * (seg.p1 - OBBtest.c);
//     Vec3 p2 = Rtrans * (seg.p2 - OBBtest.c);
//     Vec3 vec = p2 - p1;
//     // This ensures that p1 is always below p2 in z coordinates
//     if (vec.z < 0) {
//         vec = -vec;
//         Vec3 point = p1;
//         p1 = p2;
//         p2 = point;
//     }
//     segment seg_test;
//     seg_test = { p1, p2 };
//     // We must find the points which need to be tested later.
//     Vec3 p[2];
//     // The point p1 clamped on the OBB
//     p[0].x = clamp(p1.x, -OBBtest.e.x, OBBtest.e.x);
//     p[0].y = clamp(p1.y, -OBBtest.e.y, OBBtest.e.y);
//     p[0].z = clamp(p1.z, -OBBtest.e.z, OBBtest.e.z);
//     // The point p2 clamped on the OBB
//     p[1].x = clamp(p2.x, -OBBtest.e.x, OBBtest.e.x);
//     p[1].y = clamp(p2.y, -OBBtest.e.y, OBBtest.e.y);
//     p[1].z = clamp(p2.z, -OBBtest.e.z, OBBtest.e.z);
//     // The point on the closest OBB line
//     Vec3 vec_unit = (1 / norm(vec)) * vec;
//     Vec3 direction = dot(-p1, p2 - p1) * (vec_unit) + p1;
//     int quadrant = direction.x > 0 ? (direction.z > 0 ? 1 : 4) : (direction.z > 0 ? 2 : 3);
//     real xmin = - OBBtest.e[0];
//     real xmax = + OBBtest.e[0];
//     real ymin = - OBBtest.e[1];
//     real ymax = + OBBtest.e[1];
//     real zmin = - OBBtest.e[2];
//     real zmax = + OBBtest.e[2];
//     // Creating the eight original vertices
//     Vec3 orgvert[8];
//     orgvert[0] = { xmax, ymin, zmax };
//     orgvert[1] = { xmax, ymax, zmax };
//     orgvert[2] = { xmin, ymin, zmax };
//     orgvert[3] = { xmin, ymax, zmax };
//     orgvert[4] = { xmin, ymin, zmin };
//     orgvert[5] = { xmin, ymax, zmin };
//     orgvert[6] = { xmax, ymin, zmin };
//     orgvert[7] = { xmax, ymax, zmin };
//     segment seg_OBB;
//     seg_OBB.p1 = orgvert[2 * quadrant - 2];
//     seg_OBB.p2 = orgvert[2 * quadrant - 1];
//     segment corner_p1;
//     corner_p1.p1 = seg_OBB.p1;
//     corner_p1.p2 = p1;
//     Vec3 cross_product = cross(vec, corner_p1.p1 - p1);
//     cross_product.y = quadrant == 2 || quadrant == 3 ? -cross_product.y : cross_product.y;
//     if (cross_product.y > 0) {
//     }
//     real dist_seg = sqrt(distmin(seg, seg_OBB));
//     real dist_p3 = distmin(seg, p[0]);
//     real dist_p4 = distmin(seg, p[1]);
//     real dist_min = (dist_seg < dist_p3 && dist_seg < dist_p4) ? dist_seg : (dist_p3 < dist_p4) ? dist_p3 : dist_p4;
//     return dist_min - caps.r;
// }
// // Points method
// host_fn real dist_OBB_caps(OBB OBBtest, capsule caps) {
//     segment seg;
//     seg.p1 = caps.p1;
//     seg.p2 = caps.p2;
//     Mat3 Rtrans = transpose(OBBtest.R);
//     // Mat3 Rtrans = OBB.R;
//     Vec3 p1 = Rtrans * (seg.p1 - OBBtest.c);
//     Vec3 p2 = Rtrans * (seg.p2 - OBBtest.c);
//     Vec3 vec = p2 - p1;
//     Vec3 point = p2;
//     if (vec.z < 0) {
//         vec = -vec;
//         point = p1;
//     }
//     p2 = p1 + vec;
//     segment seg_test;
//     seg_test = {p1,p2};
//     // We must find the points which need to be tested later.
//     Vec3 p[8];
//     // The point p1 clamped on the OBB
//     p[0].x = clamp(p1.x, -OBBtest.e.x, OBBtest.e.x);
//     p[0].y = clamp(p1.y, -OBBtest.e.y, OBBtest.e.y);
//     p[0].z = clamp(p1.z, -OBBtest.e.z, OBBtest.e.z);
//     // The point p2 clamped on the OBB
//     p[1].x = clamp(p2.x, -OBBtest.e.x, OBBtest.e.x);
//     p[1].y = clamp(p2.y, -OBBtest.e.y, OBBtest.e.y);
//     p[1].z = clamp(p2.z, -OBBtest.e.z, OBBtest.e.z);
//     // Early return : if both points project on the same face, then the smallest of the two distances
//     // will be returned (there is no need to test anything else).
//     // if ((p[1].x > -OBBtest.e.x && p[1].x < OBBtest.e.x && p[2].x > -OBBtest.e.x && p[2].x < OBBtest.e.x) &&
//     // ((p[1].y > -OBBtest.e.y && p[1].y < OBBtest.e.y && p[2].y > -OBBtest.e.y && p[2].y < OBBtest.e.y) || (p[1].z > -OBBtest.e.z && p[1].z < OBBtest.e.z && p[2].z > -OBBtest.e.z && p[2].z < OBBtest.e.z)) ||
//     // ((p[1].y > -OBBtest.e.y && p[1].y < OBBtest.e.y && p[2].y > -OBBtest.e.y && p[2].y < OBBtest.e.y) || (p[1].z > -OBBtest.e.z && p[1].z < OBBtest.e.z && p[2].z > -OBBtest.e.z && p[2].z < OBBtest.e.z)))
//     // {
//     //     real dist1 = distmin(seg, p[1]);
//     //     real dist2 = distmin(seg, p[2]);
//     //     return dist1 < dist2 ? dist1 - caps.r : dist2 - caps.r;
//     // }
//     // We find the face which will be tested : this face is the face in which the vector from the center
//     // of the OBB to the closest point on the line to the center crosses first.
//     // Vec3 closept = closept_origin(seg);
//     // Vec3 t;
//     // t.x = closept.x / OBBtest.e.x;
//     // t.y = closept.y / OBBtest.e.y;
//     // t.z = closept.z / OBBtest.e.z;
//     // Vec3 t_abs = { abs(t.x), abs(t.y), abs(t.z) };
//     // int face = (t_abs.x > t_abs.y && t_abs.x > t_abs.z) ? 0 : (t_abs.y > t_abs.z) ? 1 : 2;
//     // int sign = (face == 1 && t.x < 0) || (face == 2 && t.y < 0) || (face == 3 && t.z < 0) ? -1 : 1;
//     // int face_idx = face*sign;
//     // std::vector<real> dist;
//     // std::vector<real> t(6);
//     // t[0] = (OBBtest.e.x - p1.x) / vec.x;
//     // t[1] = (-OBBtest.e.x - p1.x) / vec.x;
//     // t[2] = (OBBtest.e.y - p1.y) / vec.y;
//     // t[3] = (-OBBtest.e.y - p1.y) / vec.y;
//     // t[4] = (OBBtest.e.z - p1.z) / vec.z;
//     // t[5] = (-OBBtest.e.z - p1.z) / vec.z;
//     // int dim1 = (face == 0) ? 1 : 0;
//     // int dim2 = (face == 2) ? 1 : 2;
//     real s[6];
//     // todo : fix division by 0
//     s[0] = (OBBtest.e.x - p1.x) / vec.x;
//     s[1] = (-OBBtest.e.x - p1.x) / vec.x;
//     s[2] = (OBBtest.e.y - p1.y) / vec.y;
//     s[3] = (-OBBtest.e.y - p1.y) / vec.y;
//     s[4] = (OBBtest.e.z - p1.z) / vec.z;
//     s[5] = (-OBBtest.e.z - p1.z) / vec.z;
//     real s_clamped[6];
//     s_clamped[0] = clamp(s[0], 0, 1);
//     s_clamped[1] = clamp(s[1], 0, 1);
//     s_clamped[2] = clamp(s[2], 0, 1);
//     s_clamped[3] = clamp(s[3], 0, 1);
//     s_clamped[4] = clamp(s[4], 0, 1);
//     s_clamped[5] = clamp(s[5], 0, 1);
//     // We find the points on the segment which will clamp to the closest point on the OBB (with t)
//     Vec3 pt_segment[6];
//     // !!! Is it possible to skip the first step by creating points with logic here?
//     pt_segment[0] = p1 + s_clamped[0]*vec;
//     pt_segment[1] = p1 + s_clamped[1]*vec;
//     pt_segment[2] = p1 + s_clamped[2]*vec;
//     pt_segment[3] = p1 + s_clamped[3]*vec;
//     pt_segment[4] = p1 + s_clamped[4]*vec;
//     pt_segment[5] = p1 + s_clamped[5]*vec;
//     for (int i = 0; i < 6; i++) {
//         p[i + 2].x = clamp(pt_segment[i].x, -OBBtest.e.x, OBBtest.e.x);
//         p[i + 2].y = clamp(pt_segment[i].y, -OBBtest.e.y, OBBtest.e.y);
//         p[i + 2].z = clamp(pt_segment[i].z, -OBBtest.e.z, OBBtest.e.z);
//     }
//     // for debugging : p_general
//     Vec3 p_general[8];
//     Vec3 p_segment_general[8];
//     for (int i = 0; i < 8; i++) {
//         p_general[i] = OBBtest.R*p[i+2] + OBBtest.c;
//         p_segment_general[i] = OBBtest.R*pt_segment[i+2] + OBBtest.c;
//     }
//     real dist[8];
//     // Initializing at dist[0]
//     Vec3 pt_on_segment = ptint(seg_test, p[0]);
//     dist[0] = norm(pt_on_segment - p[0]);
//     dist[0] = (dot(p[0], p[0]) >= dot(pt_on_segment, p[0])) ? -dist[0] : dist[0];
//     real distcurrent = dist[0];
//     for (int i = 1; i < 8; i++) {
//         pt_on_segment = ptint(seg_test, p[i]);
//         dist[i] = norm(pt_on_segment - p[i]);
//         dist[i] = (dot(p[i], p[i]) >= dot(pt_on_segment, pt_on_segment)) ? -dist[i] : dist[i];
//         if (dist[i] >= 0 && (dist[i] < distcurrent || distcurrent < 0)) {
//             distcurrent = abs(dist[i]) < COLLISION_EPSILON ? distcurrent : dist[i];
//         }
//         else if (dist[i] < 0 && dist[i] > distcurrent) {
//             distcurrent = abs(dist[i]) < COLLISION_EPSILON ? distcurrent : dist[i];
//         }
//     }
//     return distcurrent - caps.r;
//     // bool dim1_intersection[2] = { s[dim1][0] >= 0 && s[dim1][0] <= 1, s[dim1][1] >= 0 && s[dim1][1] <= 1};
//     // bool dim2_intersection[2] = { s[dim2][0] >= 0 && s[dim2][0] <= 1, s[dim2][1] >= 0 && s[dim2][1] <= 1};
//     // // how many times does our segment intersect its first dimension axis?
//     // int dim1_number =   dim1_intersection[0] && dim1_intersection[1] ? 2 :
//     //                     dim1_intersection[0] || dim1_intersection[1] ? 1 : 0;
//     // int dim2_number =   dim2_intersection[0] && dim2_intersection[1] ? 2 :
//     //                     dim2_intersection[0] || dim2_intersection[1] ? 1 : 0;
//     // real coord1 = sign*OBBtest.e[face - 1];
//     // int idx = 2;
//     // if (face == 1) {
//     //     dist[idx++] = (dim1_intersection[0] == true) ? distmin(seg, { sign * OBBtest.e.x, p1.y + s[1][0]*vec.y, OBBtest.e.z }) : INF_REAL;
//     //     dist[idx++] = (dim1_intersection[1] == true) ? distmin(seg, { sign * OBBtest.e.x, p1.y + s[1][1]*vec.y, -OBBtest.e.z }) : INF_REAL;
//     //     dist[idx++] = (dim2_intersection[0] == true) ? distmin(seg, { sign * OBBtest.e.x, OBBtest.e.y, p1.z + s[2][0]*vec.z }) : INF_REAL;
//     //     dist[idx++] = (dim2_intersection[1] == true) ? distmin(seg, { sign * OBBtest.e.x, -OBBtest.e.y, p1.z + s[2][1]*vec.z }) : INF_REAL;
//     // }
//     // if (face == 2) {
//     //     // todo : face 2 and 3.
//     //     dist[idx++] = (dim1_intersection[0] == true) ? distmin(seg, { sign * OBBtest.e.x, p1.y + s[1][0]*vec.y, OBBtest.e.z }) : INF_REAL;
//     //     dist[idx++] = (dim1_intersection[1] == true) ? distmin(seg, { sign * OBBtest.e.x, p1.y + s[1][1]*vec.y, -OBBtest.e.z }) : INF_REAL;
//     //     dist[idx++] = (dim2_intersection[0] == true) ? distmin(seg, { sign * OBBtest.e.x, OBBtest.e.y, p1.z + s[2][0]*vec.z }) : INF_REAL;
//     //     dist[idx++] = (dim2_intersection[1] == true) ? distmin(seg, { sign * OBBtest.e.x, -OBBtest.e.y, p1.z + s[2][1]*vec.z }) : INF_REAL;
//     // }
//     // Vec3 pt1_1;
//     // Vec3 pt1_2;
//     // Vec3 pt2_1;
//     // Vec3 pt2_2;
//     // pt1_1.x = face == 1 ? sign * OBBtest.e.x : (s[0][0] >= 0 && s[0][0] <= 1 ? s[0][0] : (s[0][1] >= 0 && s[0][1] <= 1 ? s[0][1] : 0));
//     // pt1_1.y = face == 2 ? sign * OBBtest.e.y : (s[dim1][0] >= 0 && s[dim1][0] <= 1 ? s[dim1][0] : (s[dim1][1] >= 0 && s[dim1][1] <= 1 ? s[dim1][1] : 0));
//     // pt1_1.z = face == 3 ? sign * OBBtest.e.z : (s[dim1][0] >= 0 && s[dim1][0] <= 1 ? s[dim1][0] : (s[dim1][1] >= 0 && s[dim1][1] <= 1 ? s[dim1][1] : 0));
//     // if (s[0] >= 0 && s[0] <= 1) {
//     //     dist[idx++] = distmin(seg, { sign * OBBtest.e.x, p1.y + s[0]*vec.y, OBBtest.e.z });
//     // }
//     // if (s[1] >= 0 && s[1] <= 1) {
//     //     dist[idx++] = distmin(seg, { sign * OBBtest.e.x, p1.y + s[1]*vec.y, -OBBtest.e.z });
//     // }
//     // if (s[2] >= 0 && s[2] <= 1) {
//     //     dist[idx++] = distmin(seg, { sign * OBBtest.e.x, OBBtest.e.y, p1.z + s[2]*vec.z });
//     // }
//     // if (s[3] >= 0 && s[3] <= 1) {
//     //     dist[idx++] = distmin(seg, { sign * OBBtest.e.x, -OBBtest.e.y, p1.z + s[3]*vec.z });
//     // }
//     // If there is no intersection, then we should take the corner closest to the segment.
//     // if (size(dist) == 0) {
//     //     return distmin(seg, )
//     // }
//     // real distcurrent = dist[0];
//     // for (auto& d : dist) {
//     //     if (d >= 0 && (d < distcurrent || distcurrent < 0)) {
//     //         distcurrent = d;
//     //     }
//     //     else if (d < 0 && d > distcurrent) {
//     //         distcurrent = d;
//     //     }
//     // }
// }


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
};


TEST_CASE("Collisions", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

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

    // Minkowski sum method
    for (auto t : test_obb) {
        real dist = distmin(t.box, t.caps);
        // std::cout << "The distance difference is " << abs(dist - t.expected_dist) << ", or " << abs(dist - t.expected_dist) * 100 / abs(t.expected_dist) << " %." << std::endl;
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }
    // GJK Algorithm
    for (auto t : test_obb) {
        gjkresult res = GJK_solve_gjk_simple(t.caps, t.box);
        CHECK(abs(res.minimal_distance - t.expected_dist) < TESTCOLL_EPSILON);
    }
    // boolean GJK algorithm
    for (auto t : test_obb) {
        bool res = GJK_bool(t.caps, t.box);
        bool expected = (t.expected_dist < 0);
        CHECK(res == expected);
    }

}

#endif

} // namespace blast
