#pragma once
#include <blast>

namespace blast {


// Returns distance between an box and a point
inline real distance(const Box &box, const Vec3 &point) {
    Mat3 Rtrans = transpose(box.rotation);

    Vec3 point_box = Rtrans * (point - box.center);

    Vec3 proj = {clamp(point_box.x, -box.extents.x, box.extents.x), clamp(point_box.y, -box.extents.y, box.extents.y), clamp(point_box.z, -box.extents.z, box.extents.z)};
    Array dist_in(3);
    dist_in[0] = std::abs(point_box.x) - box.extents.x;
    dist_in[1] = std::abs(point_box.y) - box.extents.y;
    dist_in[2] = std::abs(point_box.z) - box.extents.z;
    real result_if_inside = max(dist_in);
    if (max(dist_in) == 0.29048075959276026) {
        int stop = 0;
    }
    return result_if_inside > 0 ? norm(proj - point_box) : result_if_inside;
}

// Returns distance between a capsule and a point
inline real distance(const Capsule &capsule, const Vec3 &point) {
    Segment segment;
    segment.p1 = capsule.p1;
    segment.p2 = capsule.p2;
    return distance(segment, point) - capsule.radius;
}

// Returns distance between a sphere and a point
inline real distance(const Sphere &sphere, const Vec3 &point) {
    return norm(point - sphere.center) - sphere.radius;
}

// Gets the point by linear interpolation from caps_list according to values of x
inline Vec3 get_point(const Array& x, const Matrix &capsule_list) {
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
inline real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const World* world) {
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
inline real obj_function(const Array& x, const Matrix &robot_cartesian_positions, const Box &box) {
    Vec3 p = get_point(x, robot_cartesian_positions);
    return distance(box, p);
}

// Gets the point by linear interpolation from caps_list according to values of x
inline Vec3 get_point_with_directions(const Array& x, const Matrix &capsule_list, Array *directions) {
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
inline real distance_with_directions(const Box &box, const Vec3 &point, Vec3* d) {
    Mat3 Rtrans = transpose(box.rotation);

    Vec3 point_box = Rtrans * (point - box.center);

    Vec3 proj = {clamp(point_box.x, -box.extents.x, box.extents.x), clamp(point_box.y, -box.extents.y, box.extents.y), clamp(point_box.z, -box.extents.z, box.extents.z)};
    Array dist_in(3);
    dist_in[0] = std::abs(point.x) - box.extents.x;
    dist_in[1] = std::abs(point.y) - box.extents.y;
    dist_in[2] = std::abs(point.z) - box.extents.z;
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
inline real distance_with_directions(const Capsule &capsule, const Vec3 &point, Vec3* d) {
    Vec3 ab = capsule.p2 - capsule.p1;
    real t = dot(point - capsule.p1, ab) / dot(ab, ab);
    t = clamp(t, 0, 1);
    Vec3 p_capsule = capsule.p1 + t * ab;

    (*d) = p_capsule - point;

    return norm((*d)) - capsule.radius;
}

// Returns distance between a sphere and a point and displacement vector
inline real distance_with_directions(const Sphere &sph_test, const Vec3 &point, Vec3* d) {
    (*d) = sph_test.center - point;
    return norm((*d)) - sph_test.radius;
}

// Calls get_point and tests this point against the full world
inline real obj_function_gradient(const Array& x, const Matrix &robot_cartesian_positions, const World* world, Array* gradient) {
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
inline real obj_function_gradient(const Array& x, const Matrix &robot_cartesian_positions, const Box &box, Array* gradient) {
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

} // namespace blast