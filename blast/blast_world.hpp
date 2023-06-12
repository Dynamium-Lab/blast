#pragma once
#include <iostream>
#include "blast.hpp"

namespace blast {

const real COLLISION_EPSILON = 1e-9;

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

struct objlist {
    std::vector<OBB> OBBlist;
    std::vector<cylinder> cyllist;
    std::vector<sphere> sphlist;
    std::vector<capsule> capslist;
};
objlist world;

struct robot {
    std::vector<Vec3> pts; // pts of the center of every articulation in the given position
    std::vector<real> r; // radii of the respective link
};


// A paper in 2005 suggested an elegant solution to the problem of the point in the triangle. This is this solution (adapted from https://gamedev.stackexchange.com/questions/28781/easy-way-to-project-point-onto-triangle-or-plane)
bool pointInTriangle(Vec3 V1, Vec3 V2, Vec3 o, Vec3 P) {
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

bool pointInSurf(Vec3 V1, Vec3 V2, Vec3 o, Vec3 P){
    bool tri1 = pointInTriangle(V1,V2,o,P);
    bool tri2 = pointInTriangle(-V1,-V2,o+V1+V2,P);
    return (tri1 || tri2);
}

Vec3 ptint(segment seg, Vec3 point) {
    Vec3 ab = seg.p2 - seg.p1;
    real t = dot(point - seg.p1, ab) / dot(ab, ab);

    t = clamp(t, 0, 1);

    Vec3 d = seg.p1 + t * ab;
    return d;
}

real distmin(segment seg, Vec3 point) {
    Vec3 d = ptint(seg, point);
    return norm(d - point);
}

real distmin(surf surf, Vec3 point) {
    Vec3 n = cross(surf.d1, surf.d2);
    Vec3 n_unit = (1/norm(n)*n);

    if (norm(n) < COLLISION_EPSILON) return -INF_REAL;

    bool pt_in_surf = pointInSurf(surf.d1, surf.d2, surf.p, point);

    real normaldist = dot(point - surf.p, n_unit);

    if (pt_in_surf == true) return normaldist;
    else {
        segment seg1;
        segment seg2;
        segment seg3;
        segment seg4;

        seg1.p1 = surf.p;
        seg1.p2 = seg1.p1 + surf.d1;
        seg2.p1 = seg1.p2;
        seg2.p2 = seg2.p1 + surf.d2;
        seg3.p1 = seg2.p2;
        seg3.p2 = seg3.p1 - surf.d1;
        seg4.p1 = seg3.p2;
        seg4.p2 = seg4.p1 - surf.d2;

        real dist1 = distmin(seg1, point);
        real dist2 = distmin(seg2, point);
        real dist3 = distmin(seg3, point);
        real dist4 = distmin(seg4, point);

        real distmin = dist1;
        if (dist2 < distmin) distmin = dist2;
        if (dist3 < distmin) distmin = dist3;
        if (dist4 < distmin) distmin = dist4;

        if (normaldist < 0) distmin = -distmin;

        return distmin;
    }
}

two_pts closept(segment seg1, segment seg2) {
    // From : Real-Time Detection Collision (Christer Ericson)

    // Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
    // S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
    // distance between between S1(s) and S2(t)

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

real distmin(segment seg1, segment seg2) {
    two_pts close_pt = closept(seg1,seg2);
    return dot(close_pt.p1 - close_pt.p2, close_pt.p1 - close_pt.p2); // ATTENTION! CECI DONNE LE RESULTAT AU CARRE
}

Vec3 ClosestPtPointPlane(Vec3 q, plane p) {
    return q - (dot(p.n, q - p.p) / dot(p.n, p.n)) * p.n;;
}

two_pts intersection(circle circ, segment seg) {
    Vec3 p1 = seg.p1 - circ.p;
    Vec3 p2 = seg.p2 - circ.p;
    
    Vec3 d1 = p1;
    Vec3 d2 = cross(circ.n,d1);

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

real distmin(capsule caps, cylinder cyl) {
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

        cent = cyl.p1;
        other = cyl.p2;

        if (cond2 < cond1) {
            cent = cyl.p2;
            other = cyl.p1;
        }

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

        real dist_norm_sq_min = dist_norm_sq1;
        if  (dist_norm_sq2 < dist_norm_sq1) {
            dist_norm_sq_min = dist_norm_sq2;
        }

        // if both points project on the circle plane, then the distance will be the minimum of the two normal distances calculated
        if (rad_sq1 <= cyl.r * cyl.r && rad_sq2 <= cyl.r * cyl.r)
            return sqrt(dist_norm_sq_min) - caps.r;

        // if one point projects on the circle plane, it is necessary to check the normal distance as well
        real testnormal = INF_REAL;
        if (rad_sq1 <= cyl.r * cyl.r){
            testnormal = dist_norm_sq1;
        }
        else if (rad_sq2 <= cyl.r * cyl.r) {
            testnormal = dist_norm_sq2;
        }

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

real distmin(capsule caps, sphere sph) {
    segment seg;
    seg.p1 = caps.p1;
    seg.p2 = caps.p2;
    real dist_seg_pt = distmin(seg, sph.c);
    return dist_seg_pt - caps.r - sph.r;
}

real distmin(capsule caps1, capsule caps2) {
    segment seg1;
    segment seg2;
    seg1.p1 = caps1.p2;
    seg1.p2 = caps1.p1;
    seg2.p1 = caps2.p2;
    seg2.p2 = caps2.p1;
    real dist_seg_seg = /*two_segment_distance_sqr(caps1.p1, caps1.p2, caps2.p1, caps2.p2);*/ sqrt(distmin(seg1, seg2));
    return dist_seg_seg - caps1.r - caps2.r;
}

real distmin(OBB OBB, capsule caps) {
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
    for (int i = 0; i < 8; i++){
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

    // The other faces depend on the orientation of the vector either in x (faces 2-5), in y (faces 6-9), or a combination of both (faces 10-11)

    // Faces 2-5
    if (vec.x >= 0) {
        face[2].p = orgvert[5];
        face[2].d1 = - dir[1];
        face[2].d2 = vec;

        face[3].p = orgvert[1];
        face[3].d1 = - dir[1];
        face[3].d2 = dir[2];

        face[4].p = vert[2];
        face[4].d1 = dir[1];
        face[4].d2 = dir[2];

        face[5].p = orgvert[2];
        face[5].d1 = dir[1];
        face[5].d2 = vec;
    }

    if (vec.x < 0) {
        face[2].p = orgvert[6];
        face[2].d1 = dir[1];
        face[2].d2 = vec;

        face[3].p = orgvert[2];
        face[3].d1 = dir[1];
        face[3].d2 = dir[2];

        face[4].p = vert[1];
        face[4].d1 = -dir[1];
        face[4].d2 = dir[2];

        face[5].p = orgvert[1];
        face[5].d1 = -dir[1];
        face[5].d2 = vec;
    }

    // Faces 6-9
    if (vec.y >= 0) {
        face[6].p = orgvert[4];
        face[6].d1 = dir[0];
        face[6].d2 = vec;

        face[7].p = orgvert[0];
        face[7].d1 = dir[0];
        face[7].d2 = dir[2];

        face[8].p = vert[3];
        face[8].d1 = -dir[0];
        face[8].d2 = dir[2];

        face[9].p = orgvert[3];
        face[9].d1 = -dir[0];
        face[9].d2 = vec;
    }

    if (vec.y < 0) {
        face[6].p = orgvert[7];
        face[6].d1 = -dir[0];
        face[6].d2 = vec;

        face[7].p = orgvert[3];
        face[7].d1 = -dir[0];
        face[7].d2 = dir[2];

        face[8].p = vert[0];
        face[8].d1 = dir[0];
        face[8].d2 = dir[2];

        face[9].p = orgvert[0];
        face[9].d1 = dir[0];
        face[9].d2 = vec;
    }

    // Faces 10-11
    if (vec.x >= 0 && vec.y >= 0) {
        face[10].p = orgvert[1];
        face[10].d1 = dir[2];
        face[10].d2 = vec;

        face[11].p = orgvert[6];
        face[11].d1 = -dir[2];
        face[11].d2 = vec;
    }

    if (vec.x >= 0 && vec.y < 0) {
        face[10].p = orgvert[4];
        face[10].d1 = -dir[2];
        face[10].d2 = vec;

        face[11].p = orgvert[3];
        face[11].d1 = dir[2];
        face[11].d2 = vec;
    }

    if (vec.x < 0 && vec.y >= 0) {
        face[10].p = orgvert[7];
        face[10].d1 = -dir[2];
        face[10].d2 = vec;

        face[11].p = orgvert[0];
        face[11].d1 = dir[2];
        face[11].d2 = vec;
    }

    if (vec.x < 0 && vec.y < 0) {
        face[10].p = orgvert[2];
        face[10].d1 = dir[2];
        face[10].d2 = vec;
        
        face[11].p = orgvert[5];
        face[11].d1 = -dir[2];
        face[11].d2 = vec;        
    }

    real dist[12];
    dist[0] = distmin(face[0], point);

    real distcurrent = dist[0];
    
    for (int i = 1; i < 12; i++) {
        dist[i] = distmin(face[i], point);
        if (dist[i] >= 0 && (dist[i] < distcurrent || distcurrent < 0)) {
            distcurrent = dist[i];
        }
        else if (dist[i] < 0 && dist[i] > distcurrent) {
            distcurrent = dist[i];
        }
    }
    return distcurrent - caps.r;
}

void add_OBB(Vec3 c, Vec3 e, Mat3 R) {
    OBB new_OBB;
    new_OBB.c = c;
    new_OBB.e = e;
    new_OBB.R = R;
    world.OBBlist.push_back(new_OBB);
}

void add_sph(Vec3 c, real r) {
    sphere new_sph;
    new_sph.c = c;
    new_sph.r = r;
    world.sphlist.push_back(new_sph);
}

void add_cyl(Vec3 p1, Vec3 p2, real r) {
    cylinder new_cyl;
    new_cyl.p1 = p1;
    new_cyl.p1 = p2;
    new_cyl.r = r;
    world.cyllist.push_back(new_cyl);
}

void add_caps(Vec3 p1, Vec3 p2, real r) {
    capsule new_caps;
    new_caps.p1 = p1;
    new_caps.p1 = p2;
    new_caps.r = r;
    world.capslist.push_back(new_caps);
}
 
// real test_collision(/*robot, */objlist world, int n_var){
//     real dist_min[n_var];
//     int idx = 0;

//     // OBB collisions
//     for (int i = 0; i < size(world.OBBlist); i++) {
//         real dist = distmin(t.caps, world.OBBlist[i]);
//         for (int j = 0; j < n_var; j++){
//                 if (dist < dist_min[j]) {
//                     for (int k = n_var; k > j; k--) {
//                         dist_min[k] = dist_min[k-1];
//                     }
//                     dist_min[j] = dist;
//                     break
//                 }
//         }
//     }

//     // Capsule collisions
//     for (int i = 0; i < size(world.capslist); i++) {
//         real dist = distmin(t.caps, world.capslist[i]);
//         for (int j = 0; j < n_var; j++){
//                 if (dist < dist_min[j]) {
//                     for (int k = n_var; k > j; k--) {
//                         dist_min[k] = dist_min[k-1];
//                     }
//                     dist_min[j] = dist;
//                     break
//                 }
//         }
//     }

//     // Cylinder collisions
//     for (int i = 0; i < size(world.cyllist); i++) {
//         real dist = distmin(t.caps, world.cyllist[i]);
//         for (int j = 0; j < n_var; j++){
//                 if (dist < dist_min[j]) {
//                     for (int k = n_var; k > j; k--) {
//                         dist_min[k] = dist_min[k-1];
//                     }
//                     dist_min[j] = dist;
//                     break
//                 }
//         }
//     }

//     // Sphere collisions
//     for (int i = 0; i < size(world.sphlist); i++) {
//         real dist = distmin(t.caps, world.sphlist[i]);
//         for (int j = 0; j < n_var; j++){
//                 if (dist < dist_min[j]) {
//                     for (int k = n_var; k > j; k--) {
//                         dist_min[k] = dist_min[k-1];
//                     }
//                     dist_min[j] = dist;
//                     break
//                 }
//         }
//     }

//     // Self-collision
//     // note : this does not take into account the collision between two subsequent links of a robot, 
//     // which means that angle constraints must be put on the articulations.
//     real n_pts = size(robot.pts);
//     real n_link = n_pts - 1;
//     capsule link[n_link];

//     for (int i = 0; i < n_pts - 1; i++) {
//         link[i].p1 = robot.pts[i];
//         link[i].p2 = robot.pts[i+1];
//         link[i].r = robot.r[i];
//     }

//     for (int i = 0; i < n_link; i++) {
//         for (int j = 2 + i; j < n_link; j++) {
//             dist = distmin(link[i], link[i+j]);
//             for (int k = 0; k < n_var; k++){
//                 if (dist < dist_min[j]) {
//                     for (int l = n_var; l > k; l--) {
//                         dist_min[l] = dist_min[l-1];
//                     }
//                     dist_min[k] = dist;
//                     break
//                 }
//             }
//         }
//     }

//     return dist_min;
// }

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


real GJK_triangle_area_2d(real x1, real y1, real x2, real y2, real x3, real y3) {
    return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

Vec3 GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 c, Vec3 p) {
    Vec3 abc = cross(b - a, c - a);
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

    Vec3 result;
    result.x = nu * ood;
    result.y = nv * ood;
    result.z = 1 - result.x - result.y;
    return result;
}

real GJK_convert_barycentric(Vec3 a, Vec3 b, Vec3 p) {
    Vec3 ab = b - a;
    Vec3 ap = p - a;

    real t = dot(ap, ab) / dot(ab, ab);
    return t;
}

Vec3 GJK_convert_cartesian(Vec3 a, Vec3 b, real barycentric) {
    Vec3 result;
    result.x = a.x + barycentric * (b.x - a.x);
    result.y = a.y + barycentric * (b.y - a.y);
    result.z = a.z + barycentric * (b.z - a.z);
    return result;

}

Vec3 GJK_convert_cartesian(Vec3 a, Vec3 b, Vec3 c, Vec3 barycentric) {
    Vec3 result;
    result.x = a.x * barycentric.x + b.x * barycentric.y + c.x * barycentric.z;
    result.y = a.y * barycentric.x + b.y * barycentric.y + c.y * barycentric.z;
    result.z = a.z * barycentric.x + b.z * barycentric.y + c.z * barycentric.z;
    return result;
}

two_pts GJK_get_local_points(ComplexSimplex simplex, Vec3 p) {
    two_pts loc;
    real barycentric;
    Vec3 barycentric3;

    switch (simplex.count) {
    case '1':
        loc.p1 = simplex.a1;
        loc.p2 = simplex.a2;
        return loc;

    case '2':
        barycentric = GJK_convert_barycentric(simplex.a, simplex.b, p);
        loc.p1 = GJK_convert_cartesian(simplex.a1, simplex.b1, barycentric);
        loc.p2 = GJK_convert_cartesian(simplex.a2, simplex.b2, barycentric);
        return loc;

    case '3':
        barycentric3 = GJK_convert_barycentric(simplex.a, simplex.b, simplex.c, p);
        loc.p1 = GJK_convert_cartesian(simplex.a1, simplex.b1, simplex.c1, barycentric3);
        loc.p2 = GJK_convert_cartesian(simplex.a2, simplex.b2, simplex.c2, barycentric3);
        return loc;

    default:
        Assert (false/*"invalid simplex dimension for nearest feature resolution."*/);
        return loc;
    }
}

Vec3 GJK_get_support(KodaVertices vertices, Vec3 direction) {
    real largest_dot = dot(vertices.fixed[0], direction);
    Vec3 largest_vertex = vertices.fixed[0];
    for (int i = 0; i < vertices.count; i++) {
        Vec3 vertex = vertices.fixed[i];

        real current_dot = dot(vertex, direction);
        if (current_dot > largest_dot) {
            largest_dot = current_dot;
            largest_vertex = vertex;
        }
    }
    return largest_vertex;
}

Vec3 GJK_get_support(sphere sph, Vec3 direction) {
    return sph.c + sph.r * direction;
}

Vec3 GJK_get_support(capsule caps, Vec3 direction) {
    Vec3 ab = caps.p2 - caps.p1;
    if (norm(cross(ab,direction)) < COLLISION_EPSILON)
        return caps.p1 + caps.r*direction; // In the case where the direction is perpendicular to the capsule.
    else if (GJK_same_direction(ab, direction))
        return caps.p2 + caps.r*direction;
    else
        return caps.p1 + caps.r*direction;
}

Vec3 GJK_get_support(cylinder cyl, Vec3 direction) {
    Vec3 ab = cyl.p2 - cyl.p1;
    if (norm(cross(ab,direction)) < COLLISION_EPSILON)
        return cyl.p1 + cyl.r*direction; // In the case where the direction is perpendicular to the cylinder.
    plane circ;
    circ.n = ab;
    circ.p = cyl.p2;
    Vec3 direction_unit = (1/norm(direction))*direction; // maybe already done before?
    Vec3 proj = ClosestPtPointPlane(direction_unit, circ);
    if (dot(proj - cyl.p2,proj - cyl.p2) - cyl.r*cyl.r < COLLISION_EPSILON) // the direction is similar to the cylinder and the projected point should be on the face
        if (GJK_same_direction(ab, direction))
            return proj;
        else
            return proj - ab;
    else {
        if (GJK_same_direction(ab, direction)){
            Vec3 dir_proj = proj - cyl.p2;
            Vec3 dir_proj_unit = (1/norm(dir_proj))*dir_proj;
            return cyl.p2 + cyl.r*dir_proj_unit;
        }
        else {
            Vec3 dir_proj = proj - cyl.p1;
            Vec3 dir_proj_unit = (1/norm(dir_proj))*dir_proj;
            return cyl.p1 + cyl.r*dir_proj_unit;
        }
    }
}

bool GJK_same_direction(Vec3 a, Vec3 b) {
    return dot(a,b) > 0;
}

Vec3 GJK_solve_simplex2(ComplexSimplex simplex) {
    Vec3 ab = simplex.b - simplex.a;
    Vec3 ao = -simplex.a;

    if (GJK_same_direction(ab, ao)) {
        // the origin falls on the line, project the origin onto the line and return
        real t = dot(ao, ab) / dot(ab, ab);
        return simplex.a + t * ab;
    }
    else {
        // the origin is closer to a
        simplex.count = 1;
        return simplex.a;
    }
}

Vec3 GJK_solve_simplex3(ComplexSimplex simplex) {
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
            return simplex.a + t * ac;
        }
        else {
            Vec3 ab = simplex.b - simplex.a;
            if (GJK_same_direction(ab, ao)) {
                // origin is closest to the line ab
                simplex.count = 2;
                real t = dot(ao, ab) / dot(ab, ab);
                return simplex.a + t * ab;
            }
            else {
                // the origin is nearest to the point a
                simplex.count = 1;
                return simplex.a;
            }
        }
    }
    else {
        Vec3 ab = simplex.b - simplex.a;
        if (GJK_same_direction(cross(ab, abc), ao)) {
            if (GJK_same_direction(ab, ao)) {
                // the origin is closest to line ab
                simplex.count = 1;
                real t = dot(ao, ab) / dot(ab, ab);
                return simplex.a + t * ab;
             }
            else {
            // the origin is nearest to the point a
            simplex.count = 1;
            return simplex.a;
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
                return simplex.a + ab * v + ac * w;
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
                return simplex.a + ab * v + ac * w;
            }
        }
    }
}

Vec3 GJK_solve_simplex4(ComplexSimplex simplex) {
    Vec3 abc = cross(simplex.b - simplex.a, simplex.c - simplex.a);
    Vec3 acd = cross(simplex.c - simplex.a, simplex.d - simplex.a);
    Vec3 adb = cross(simplex.d - simplex.a, simplex.b - simplex.a);
    Vec3 ao = -simplex.a;

    real abc_dir = GJK_same_direction(abc, ao);
    real acd_dir = GJK_same_direction(acd, ao);
    real adb_dir = GJK_same_direction(adb, ao);

    if (!abc_dir && !acd_dir && !adb_dir) {
        // the origin is inside the simplex
        return { 0, 0, 0 };
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

    ComplexSimplex simplex_acd;
    //simplex_acd.b_all = simplex_acd.c_all;
    simplex_acd.b = simplex.c;
    simplex_acd.b1 = simplex.c1;
    simplex_acd.b2 = simplex.c2;
    //simplex_acd.c_all = simplex_acd.d_all;
    simplex_acd.c = simplex.d;
    simplex_acd.c1 = simplex.d1;
    simplex_acd.c2 = simplex.d2;

    simplex_acd.count = 3;

    ComplexSimplex simplex_adb;
    // simplex_adb.b_all = simplex_adb.b_all;
    simplex_adb.b = simplex.c;
    simplex_adb.b1 = simplex.c1;
    simplex_adb.b2 = simplex.c2;
    // simplex_adb.c_all = simplex_adb.d_all;
    simplex_adb.c = simplex.d;
    simplex_adb.c1 = simplex.d1;
    simplex_adb.c2 = simplex.d2;

    simplex_adb.count = 3;

    Vec3 p_abc = GJK_solve_simplex3(simplex_abc);
    Vec3 p_acd = GJK_solve_simplex3(simplex_acd);
    Vec3 p_adb = GJK_solve_simplex3(simplex_adb);

    real abc_d2 = dot(p_abc, p_abc);
    real acd_d2 = dot(p_acd, p_acd);
    real adb_d2 = dot(p_adb, p_adb);

    if ((abc_d2 <= acd_d2) && (abc_d2 <= adb_d2)){
        simplex = simplex_abc;
        return p_abc;
    }
    else if ((acd_d2 <= abc_d2) && (acd_d2 <= adb_d2)) {
        simplex = simplex_acd;
        return p_acd;
    }
    else if ((adb_d2 <= abc_d2) && (adb_d2 <= acd_d2)) {
        simplex = simplex_adb;
        return p_adb;
    }

    Assert (false);

    // the origin isn't outside of any plane, so it's inside the tetrahedron
    return { 0,0,0 };
}

gjkresult GJK_solve_gjk_simple(KodaVertices v1, KodaVertices v2) {
    ComplexSimplex simplex;
    gjkresult results;

    Vec3 direction = v1.fixed[0] - v2.fixed[0];
    simplex.a1 = GJK_get_support(v1, -direction);
    simplex.a2 = GJK_get_support(v2, direction);
    simplex.a = simplex.a2 - simplex.a1;
    simplex.count = 1;

    while (true) {
        Vec3 p;
        switch (simplex.count) {
        case '1':
            p = simplex.a;
        case '2':
            p = GJK_solve_simplex2(simplex);
        case '3':
            p = GJK_solve_simplex3(simplex);
        case '4':
            p = GJK_solve_simplex4(simplex);
        default:
            std::cout << "error: the simplex size is not valid!" << std::endl;
        }
        if (dot(p, p) < COLLISION_EPSILON) {
            results.intersection = true;
            results.final_simplex = simplex;
            results.minimal_distance = 0;
            return results;
        }

        direction = -p;

        Vec3 support1 = GJK_get_support(v1, -direction);
        Vec3 support2 = GJK_get_support(v2, direction);
        
        Vec3 support = support2 - support1;

        if (dot(p, direction) >= dot(support, direction)) {
            two_pts local = GJK_get_local_points(simplex, p);

            results.intersection = false;
            results.A_closept = local.p1;
            results.B_closept = local.p2;
            results.final_simplex = simplex;
            results.minimal_distance = norm(results.A_closept - results.B_closept);
            return results;
        }
        Assert(simplex.count < 4 /*"You cannot have a simplex that's a tetrahedron at this point."*/);

        simplex.d = simplex.c;
        simplex.d1 = simplex.c1;
        simplex.d2 = simplex.c2;

        simplex.c = simplex.b;
        simplex.c1 = simplex.b1;
        simplex.c2 = simplex.b2;

        simplex.b = simplex.c;
        simplex.b1 = simplex.c1;
        simplex.b2 = simplex.c2;

        simplex.count += 1;

        simplex.a1 = support1;
        simplex.a2 = support2;
        simplex.a = support;
    }

    /*results.intersection = false;
    results.final_simplex = simplex;
    results.A_closept = {};
    results.B_closept = {};
    results.minimal_distance = norm(results.A_closept - results.B_closept);
    return results;*/
}


int main()
{
    std::cout << "Hello World!\n";
}




} // namespace blast

#ifdef BLAST_ENABLE_TESTS

TEST_CASE("Sphere collision", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

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

    for (auto t : test) {
        real dist = distmin(t.caps, t.sph);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }
}

TEST_CASE("Capsule collision", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

    struct collision_test_caps {
        capsule caps1;
        capsule caps2;
        real expected_dist;
    };

    collision_test_caps test[] = {
        { { { 0.5, 0.95, 6.73 }, { 4, 4.9, 3.51 }, 0.5 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, -0.54 },
        { { { 5.16, 4.9, 4.37 }, { 3.55, 0.95, 8.84 }, 0.5 }, { { 4, 4.2, 4.3 }, { 0.5, 0.4, 0.9 }, 0.5 }, 0.16 },
        { { { 5.16, 4.9, 4.37 }, { 3.55, 0.95, 8.84 }, 0.5 }, { { 5.45, 4.78, 4.57 }, { 6.57, 1, -0.2 }, 0.5 }, -0.66 },
        { { { 2.59, 1.88, 4.22 }, { 18.85, 10.07, 6.1 }, 1 }, { { 0.55, 21.32, 3.08 }, { 5.07, 3.62, 4.19 }, 1 }, -1.44 },
        { { { -9.26, 12.81, 7.98 }, { 6.63, 4.62, 4.04 }, 1 }, { { 0.55, 21.32, 3.08 }, { 5.07, 3.62, 4.19 }, 1 }, -1.52 },
    };

    for (auto t : test) {
        real dist = distmin(t.caps1, t.caps2);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }
}

TEST_CASE("Cylinder collision", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

    struct collision_test_cyl {
        cylinder cyl;
        capsule caps;
        real expected_dist;
    };

    collision_test_cyl test[] = {
        /*Test 1*/ { { { -13.46, -3.61, 189 }, { -12.46, -5.11, 190 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, -0.54 },
        /*Test 2*/ { { { -13.41, -2.79, 189 }, { -15.16, -3.2, 190 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, 1.33 },
        /*Test 3*/ { { { -16.3, -3.97, 200.87 }, { -17.3, -2.47, 199.87 }, 1 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, 1.53 },
        /*Test 4*/ { { { -7.25, -10.11, 174.29 }, { -7.86, -9.44, 176.14 }, 2 }, { { -14.66, -13.24, 198.94 }, { -11.16, -7.94, 181.77 }, 1 }, 5.53 },
        /*Test 5*/ { { { -14.88, -2.74, 198.21 }, { -15.49, -2.08, 200.06 }, 2 }, { { -10.7, -8.99, 180.98 }, { -14.2, -3.19, 197.98 }, 1 }, -0.45 },
    };

    for (auto t : test) {
        real dist = distmin(t.caps, t.cyl);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }
}

TEST_CASE("OBB collision", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 6e-4;

    struct collision_test_box {
        OBB box;
        capsule caps;
        real expected_dist;
    };

    collision_test_box test[] = {
        /*Test1*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { -1.21, -5.18, 18.05 }, { -3.89, 8.59, 6.3 }, 1 }, 1.63659624 },
        /*Test2*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 10.46, 0.97, 3.7 }, { 7.79, 14.74, -8.04 }, 1 }, 1.50169942 },
        /*Test3*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, 1, 0, -1, 0, 0, 0, 0, 1 } }, { { 10.05, 1.11, 12.87 }, { 7.37, 14.89, 1.13 }, 1 }, 0.33397901 },
        /*Test5*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.82, 6.9, 12.84 }, { 4.7, -7.14, 24.59 }, 1 }, 4.00838054 },
        /*Test6*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { -3.06, 5.73, 1.67 }, { -0.39, -8.05, 13.41 }, 1 }, 0.89513819 },
        /*Test7*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.97, 2.79, 12.23 }, { 2.29, 16.57, 0.48 }, 1 }, 1.69981481 },
        /*Test8*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 11.49, -0.06, 8.32 }, { 14.17, -13.84, 20.07 }, 1 }, 0.49000000 },
        /*Test9*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 6.76, -1.55, 8 }, { 9.44, -15.33, 19.74 }, 1 }, 0.45000000 },
        /*Test10*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { -1.15, 13.92, -7.48 }, { 1.53, 0.14, 4.27 }, 1 }, 0.73046237 },
        /*Test11*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.54, 0.1, 10.59 }, { 1.86, 13.87, -1.15 }, 1 }, -1.00000000 },
        /*Test12*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 10.76, 0.28, 8.08 }, { 13.44, -13.5, 19.83 }, 1 }, -0.21961454 },
        /*Test13*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.96, 0.43, 8.3 }, { 7.64, -13.35, 20.05 }, 1 }, -1.53000000 },
        /*Test14*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 4.48, -4.07, 0.76 }, { 8.64, -4.41, 18.58 }, 1 }, 3.06923695 },
        /*Test15*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 17.48, 2.95, 13.77 }, { -0.82, 3.11, 13.77 }, 1 }, 2.41054264 },
        /*Test16*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0, -1, 0, 1, 0, 0, 0, 0, 1 } }, { { 11.18, 8.56, 4.82 }, { 11.04, -6.44, 15.29 }, 1 }, 0.09912546 },

        // OBB.R and caps.r changed
        /*Test17*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 2.94, -6.06, 5.74 }, { 4.01, -9.86, 0.98 }, 0.5 }, 1.00474115 },
        /*Test18*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 1.64, 9.52, 11.7 }, { 2.71, 5.72, 6.94 }, 0.5 }, 0.22669533 },
        /*Test19*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 3, 5.69, 6.2 }, { 4.07, 1.89, 1.44 }, 0.5 }, 0.42102531 },
        /*Test20*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 3.1, 4.23, 6.22 }, { 4.17, 0.43, 1.46 }, 0.5 },  0.10695205 },
        /*Test21*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 6.35, 8.57, 12.85 }, { 7.42, 4.77, 8.09 }, 0.5 }, 0.85902522 },
        /*Test22*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 1.88, -2.47, 7.07 }, { 2.95, -6.27, 2.31 }, 0.5 }, 0.37892454 },
        /*Test23*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 1.8, 3.83, 14.57 }, { 2.87, 0.03, 9.81 }, 0.5 }, 1.47889394 },
        /*Test24*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 4.28, 6.32, 8.33 }, { 3.21, 10.12, 13.09 }, 0.5 }, 0.82000000 },
        /*Test25*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 2.36, 2.3, 6.31 }, { 3.43, -1.5, 1.55 }, 0.5 }, 0.26887914 },
        /*Test26*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 2.93, 7.6, 12.29 }, { 4, 3.8, 7.53 }, 0.5 }, -0.93234019 },
        /*Test27*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 7.55, 0.92, 10.24 }, { 6.48, 4.72, 15 }, 0.5 }, -0.32852683 },
        /*Test28*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 6.79, 5, 13.29 }, { 7.86, 1.2, 8.52 }, 0.5 }, -0.40212621 },
        /*Test29*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 5.66, -0.92, 7.44 }, { 8.04, 4.16, 9.96 }, 0.5 }, 0.87078210 },
        /*Test30*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 8.94, 3.63, 11.01 }, { 8.94, -2.55, 11.01 }, 0.5 }, 1.24844065 },
        /*Test31*/ { { { 5, 0, 9 }, { 0.1, 5, 3 }, { 0.707, 0, -0.707, 0, 1, 0, 0.707, 0, 0.707 } }, { { 5.91, 5.9, 7.39 }, { 4.14, 5.9, 13.32 }, 0.5 }, 0.40000000 },
    };

    for (auto t : test) {
        real dist = distmin(t.box, t.caps);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }
}

// ======================================
//          GJK algorithm tests
// ======================================

TEST_CASE("Sphere collision", "[World]") {
    using namespace blast;

    real TESTCOLL_EPSILON = 1e-2;

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

    for (auto t : test) {
        real dist = distmin(t.caps, t.sph);
        CHECK(abs(dist - t.expected_dist) < TESTCOLL_EPSILON);
    }
}


#endif