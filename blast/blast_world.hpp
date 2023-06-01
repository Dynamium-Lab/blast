#include <iostream>
#include "blast/blast.hpp"

using namespace blast;
namespace blast {
    struct segment {
    Vec3 p1;
    Vec3 p2;
};

struct OBB {
    Vec3 c; // OBB center point
    Mat3 R; // Local x-, y-, and z-axes (Rotation matrix)
    Vec3 e; // Positive halfwidth extents of OBB along each axis
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
    Vec3 p;
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
    real r;
};

Vec3 point; 
Vec3 Vecprin;
Vec3 vec;
real face_crit;
Vec3 orgvert;
Vec3 vert;
real EPSILON = 1e-9;


// A paper in 2005 suggested an elegant solution to the problem of the point in the triangle. This is this solution (adapted from https://gamedev.stackexchange.com/questions/28781/easy-way-to-project-point-onto-triangle-or-plane)
bool pointInTriangle(Vec3 V1, Vec3 V2, Vec3 o, Vec3 P)
{
    // u=P2−P1
    Vec3 u = V1;
    // v=P3−P1
    Vec3 v = V2;
    // n=u×v
    Vec3 n = cross(u,v);
    // w=P−P1
    Vec3 w = P - o;
    // Barycentric coordinates of the projection P′of P onto T:
    // γ=[(u×w)⋅n]/n²
    real gamma = dot(cross(u, w),n) / dot(n,n);
    // β=[(w×v)⋅n]/n²
    real beta = dot(cross(w, v),n) / dot(n,n);
    real alpha = 1 - gamma - beta;
    // The point P′ lies inside T if:
    return ((0 <= alpha) && (alpha <= 1) && (0 <= beta) && (beta <= 1) && (0 <= gamma) && (gamma <= 1));
}

real distmin(segment obj1, Vec3 point) {
    Vec3 ab = obj1.p2 - obj1.p1;
    real t = dot(point - obj1.p1, ab) / dot(ab, ab);
    
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;

    Vec3 d = obj1.p1 + t * ab;
    return norm(d - point);
}

Vec3 ptint(segment obj1, Vec3 point) {
    Vec3 ab = obj1.p2 - obj1.p1;
    real t = dot(point - obj1.p1, obj1.p2 - obj1.p1) / dot(obj1.p2 - obj1.p1, obj1.p2 - obj1.p1);

    t = clamp(t, 0, 1);

    Vec3 d = obj1.p1 + t * ab;
    return d;
}

real distmin(surf surf, Vec3 point) {
    Vec3 n = cross(surf.d1, surf.d2);

    if (norm(n) < EPSILON) return -INF_REAL;

    bool tri1 = pointInTriangle(surf.d1, surf.d2, surf.p, point);
    bool tri2 = pointInTriangle(-surf.d1, -surf.d2, surf.p + surf.d1 + surf.d2, point);

    real normaldist = dot(point - surf.p, n) / norm(n);

    if (tri1==1 || tri2 == 1) return normaldist;
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

real distmin(segment seg1, segment seg2) {
    // From : Real-Time Detection Collision (Christer Ericson)
    
    // Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
    // S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
    // distance between between S1(s) and S2(t)
    
    real s;
    real t;
    Vec3 c1;
    Vec3 c2;

    Vec3 d1 = seg1.p1 - seg1.p2; // Direction vector of segment S1
    Vec3 d2 = seg2.p1 - seg2.p2; // Direction vector of segment S2
    Vec3 r = seg1.p1 - seg2.p1;
    real a = dot(d1, d1); // Squared length of segment S1, always nonnegative
    real e = dot(d2, d2); // Squared length of segment S2, always nonnegative
    real f = dot(d2, r);
    // Check if either or both segments degenerate into points
    if (a <= EPSILON && e <= EPSILON) {
        // Both segments degenerate into points
        s = t = 0.0f;
        c1 = seg1.p1;
        c2 = seg2.p1;
        return dot(c1 - c2, c1 - c2);
    }
    if (a <= EPSILON) {
        // First segment degenerates into a point
        s = 0.0f;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = clamp(t, 0.0f, 1.0f);
    }
    else {
        real c = dot(d1, r);
        if (e <= EPSILON) {
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
    return dot(c1 - c2, c1 - c2); // ATTENTION! CECI DONNE LE RESULTAT AU CARRE
}
two_pts closept(segment seg1, segment seg2) {
    // From : Real-Time Detection Collision (Christer Ericson)

    // Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
    // S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
    // distance between between S1(s) and S2(t)

    real s;
    real t;
    Vec3 c1;
    Vec3 c2;

    two_pts result;

    Vec3 d1 = seg1.p2 - seg1.p1; // Direction vector of segment S1
    Vec3 d2 = seg2.p2 - seg2.p1; // Direction vector of segment S2
    Vec3 r = seg1.p1 - seg2.p1;
    real a = dot(d1, d1); // Squared length of segment S1, always nonnegative
    real e = dot(d2, d2); // Squared length of segment S2, always nonnegative
    real f = dot(d2, r);
    // Check if either or both segments degenerate into points
    if (a <= EPSILON && e <= EPSILON) {
        // Both segments degenerate into points
        s = t = 0.0f;
        c1 = seg1.p1;
        c2 = seg2.p1;
        result.p1 = c1;
        result.p2 = c2;
        return result;
    }
    if (a <= EPSILON) {
        // First segment degenerates into a point
        s = 0.0f;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = clamp(t, 0.0f, 1.0f);
    }
    else {
        real c = dot(d1, r);
        if (e <= EPSILON) {
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

Vec3 ClosestPtPointPlane(Vec3 q, plane p)
{
    return q - (dot(p.n, q - p.p) / dot(p.n, p.n)) * p.n;;
}

int signx(real x) {
    if (x < 0) {
        return -1;
    }
    return 1;
}

two_pts intersection(circle circ, segment seg) {
    real dx = seg.p2.x - seg.p1.x;
    real dy = seg.p2.y - seg.p1.y;
    real dr_sq = dx * dx + dy * dy;
    real D = seg.p1.x * seg.p2.y - seg.p2.x * seg.p1.y;

    real det = sqrt(circ.r * circ.r * dr_sq - D * D);

    if (det < 0) {
        Vec3 point = ptint(seg, circ.p);
        Vec3 pointcirc = (circ.r / norm(point - circ.p)) * (point - circ.p);
        return { pointcirc, pointcirc };
    }

    real x1 = (1 / dr_sq) * (D * dy + signx(dy) * dx * det);
    real x2 = (1 / dr_sq) * (D * dy - signx(dy) * dx * det);

    real y1 = (1 / dr_sq) * (-D * dx + abs(dy) * det);
    real y2 = (1 / dr_sq) * (-D * dx - abs(dy) * det);

    Vec3 p1 = { x1,y1,0 };
    Vec3 p2 = { x2,y2,0 };
    return { p1, p2 };
}

real distmin(capsule caps, cylinder cyl) {
    segment seg1;
    segment seg2;
    seg1.p1 = caps.p1;
    seg1.p2 = caps.p2;
    seg2.p1 = cyl.p1;
    seg2.p2 = cyl.p2;
    
    two_pts points = closept(seg1, seg2); // !! Erreur dans cette formule !!
    real cond1 = dot(points.p2 - cyl.p1, points.p2 - cyl.p1);
    real cond2 = dot(points.p2 - cyl.p2, points.p2 - cyl.p2);

    // Depending on which point of the cylinder is closest, we will change the face which we check
    if (cond1 < EPSILON || cond2 < EPSILON) {
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

        real dist_sq;
        dist_sq = rad_sq1;
        if (rad_sq2 < rad_sq1) {
            dist_sq = rad_sq2;
        }

        // if both points project on the circle plane, then the distance will be the minimum of the two normal distances calculated
        if (rad_sq1 < cyl.r * cyl.r && rad_sq2 < cyl.r * cyl.r) {
            return sqrt(dist_sq) - caps.r;
        }

        // project the points on the plan
        segment proj;
        proj.p1 = proj1;
        proj.p2 = proj2;
        circle circ;
        circ.p = cent;
        circ.r = cyl.r;

        two_pts pts = intersection(circ, proj);

        real dist1 = distmin(seg1, pts.p1);
        real dist2 = distmin(seg1, pts.p2);

        if (dist2 < dist1 && dist2*dist2 < dist_sq) {
            return dist2 - caps.r;
        }
        if (dist1 < dist2 && dist1 * dist1 < dist_sq) {
            return dist1 - caps.r;
        }
        return sqrt(dist_sq) - caps.r;
    }
    else {
        return norm(points.p1 - points.p2) - caps.r - cyl.r;
    }
}

real distmin(capsule caps, sphere sph) {
    segment seg;
    seg.p1 = caps.p1;
    seg.p2 = caps.p2;
    real dist_seg_pt = distmin(seg, sph.p);
    return dist_seg_pt - caps.r - sph.r;
}

real distmin(capsule caps1, capsule caps2) {
    segment seg1;
    segment seg2;
    seg1.p1 = caps1.p1;
    seg1.p2 = caps1.p2;
    seg2.p1 = caps2.p1;
    seg2.p2 = caps2.p2;
    real dist_seg_seg = distmin(seg1, seg2);
    return dist_seg_seg - caps1.r - caps2.r;
}

real distmin(OBB OBB, segment seg) {
    Vecprin = seg.p2 - seg.p1;

    Mat3 Rtrans = transpose(OBB.R);

    vec = Rtrans * (seg.p2 - seg.p1);
    point = Rtrans * (seg.p2 - OBB.c);

    if (vec.z < 0) {
        vec = -vec;
        point = Rtrans * seg.p1;
    }

    real xmin = - OBB.e[0];
    real xmax = + OBB.e[0];
    real ymin = - OBB.e[1];
    real ymax = + OBB.e[1];
    real zmin = - OBB.e[2];
    real zmax = + OBB.e[2];
    
    // Creating the eight original vertices
    Vec3 orgvert[8];

    orgvert[0] = { xmin,ymin,zmin };
    orgvert[1] = { xmin,ymax,zmin };
    orgvert[2] = { xmax,ymin,zmin };
    orgvert[3] = { xmax,ymax,zmin };
    orgvert[4] = { xmin,ymin,zmax };
    orgvert[5] = { xmin,ymax,zmax };
    orgvert[6] = { xmax,ymin,zmax };
    orgvert[7] = { xmax,ymax,zmax };

    // Thus, we get the three main directions
    Vec3 dir[3];

    dir[0] = orgvert[2] - orgvert[0];
    dir[1] = orgvert[1] - orgvert[0];
    dir[2] = orgvert[4] - orgvert[0];

    // Creating the added vertices
    Vec3 vert[8];

    // bottom four
    vert[0] = orgvert[0] + vec;
    vert[1] = orgvert[1] + vec;
    vert[2] = orgvert[2] + vec;
    vert[3] = orgvert[3] + vec;

    // top four
    vert[4] = orgvert[4] + vec;
    vert[5] = orgvert[5] + vec;
    vert[6] = orgvert[6] + vec;
    vert[7] = orgvert[7] + vec;
    
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
    
    real distmax = -INFINITY;

    for (int i = 0; i < 12; i++) {
        dist[i] = distmin(face[i], point);
        // std::cout << "La distance de la face " << i << " est " << dist[i] << std::endl;
        if (dist[i] > distmax) {
            distmax = dist[i];
            face_crit = i;
        }
    }
    return distmax;
}

}

int main()
{
    // Initialisation de l'OBB test
    OBB OBBTEST;
    OBBTEST.c = { 0, 0, 0 };
    OBBTEST.e = { 1, 2, 3 };
    OBBTEST.R = { 1,0,0, 0,1,0, 0,0,1 };

    // Initialisation du segment test
    segment SEGMENTTEST;
    SEGMENTTEST.p1 = { 3,3,0.5 };
    SEGMENTTEST.p2 = { 2,2,2 };

    real rayon = 0;
    real dist = distmin(OBBTEST, SEGMENTTEST) - rayon;
    std::cout << "Le vecteur est defini dans les axes principaux par un point en " << SEGMENTTEST.p2[0] << ", " << SEGMENTTEST.p2[1] << ", " << SEGMENTTEST.p2[2] << " et un vecteur aux valeurs " << Vecprin[0] << ", " << Vecprin[1] << ", " << Vecprin[2] << std::endl;
    std::cout << "Le vecteur est defini dans les axes de l OBB par un point en " << point[0] << ", " << point[1] << ", " << point[2] << " et un vecteur aux valeurs (signes inverses) " << vec[0] << ", " << vec[1] << ", " << vec[2] << std::endl;
    std::cout << "La distance maximale entre les faces de l OBB et le segment est de " << dist << " et celle-ci est avec la face " << face_crit << std::endl;

    // FORMAT TEST
    //// Initialisation de l'OBB test
    //OBB OBBTEST;
    //OBBTEST.c = { 0, 0, 0 };
    //OBBTEST.e = { 1, 2, 3 };
    //OBBTEST.R = { 1,0,0, 0,1,0, 0,0,1 };
    //
    //// Initialisation du segment test
    // Vec3 size = transpose(OBBTEST.R)*OBBTEST.e;
    // Vec3 q1 = { -1.2, -1.1, -1 };
    // Vec3 q2 = { -1.9, -1.5, -2 };
    // real rayon = 0;
    // 
    // segment SEGMENTTEST;
    // SEGMENTTEST.p1 = OBBTEST.c +q1.x*size.x + q1.y*size.y + q1.z*size.z;
    // SEGMENTTEST.p2 = OBBTEST.c +q2.x*size.x + q2.y*size.y + q2.z*size.z;
    // 
    // real dist = distmin(OBBTEST, SEGMENTTEST) - rayon;

    //// Test du programme (point-plan)
    //surf SURFTEST;
    //SURFTEST.d1 = { 1,0,0 };
    //SURFTEST.d2 = { 0,1,0 };
    //SURFTEST.p = { 0,0,0 };

    //Vec3 POINTTEST = { 2,2,0 };
    //real dist = distmin(SURFTEST, POINTTEST);

    //std::cout << "La distance entre le point et le plan est " << dist << std::endl;

    //// Test de Triangle point
    //Vec3 V1 = { 0,1,1 };
    //Vec3 V2 = { 1,0,1 };
    //Vec3 o = { 0,0,0 };
    //Vec3 P = { 1,1,1 };
    //
    //bool point = pointInTriangle(V1, V2, o, P);

    //std::cout << point;
}