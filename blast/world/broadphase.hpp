#pragma once
#include <blast>

#if ENABLE_TRACY
#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"
#endif

namespace blast {

struct AABBPair {
  int  aabb_obj;
  int  aabb_cap;
  real dist;
};

struct ComparePair {
  inline bool operator()(const AABBPair& aabb1, const AABBPair& aabb2) {
    return aabb1.dist > aabb2.dist;
  }
};

struct QueuePair : public std::priority_queue<AABBPair, std::vector<AABBPair>, ComparePair> {
  QueuePair(size_t reserve_capacity = 0) {
    if (reserve_capacity > 0) {
      this->c.reserve(reserve_capacity);
    }
  }

  void clear_and_reserve(size_t capacity) {
    this->c.clear();
    this->c.reserve(capacity);
  }
};

inline host_fn void create_AABB_from_sphere(const Sphere& sphere, AxisAlignedBoundingBox& aabb, const void* ptr, int point_in_segment = -1) {
  aabb.center           = sphere.center;
  aabb.extents.x        = sphere.radius;
  aabb.extents.y        = sphere.radius;
  aabb.extents.z        = sphere.radius;
  aabb.child_type       = CollisionObjectType::sphere;
  aabb.child_ptr        = ptr;
  aabb.point_in_segment = point_in_segment;
}

inline host_fn void create_AABB_from_box(const Box& box, AxisAlignedBoundingBox& aabb, const void* ptr, int point_in_segment = -1) {
  aabb.center = box.center;

  Vec3 x = {box.extents.x, 0, 0};
  Vec3 y = {0, box.extents.y, 0};
  Vec3 z = {0, 0, box.extents.z};

  Vec3 x0 = box.rotation * x;
  Vec3 y0 = box.rotation * y;
  Vec3 z0 = box.rotation * z;

  aabb.extents.x = std::abs(x0[0]) + std::abs(y0[0]) + std::abs(z0[0]);
  aabb.extents.y = std::abs(x0[1]) + std::abs(y0[1]) + std::abs(z0[1]);
  aabb.extents.z = std::abs(x0[2]) + std::abs(y0[2]) + std::abs(z0[2]);

  aabb.child_type       = CollisionObjectType::box;
  aabb.child_ptr        = ptr;
  aabb.point_in_segment = point_in_segment;
}

inline host_fn void create_AABB_from_door(const Box& box, AxisAlignedBoundingBox& aabb, const void* ptr, int point_in_segment = -1) {
  aabb.center = box.center;

  Vec3 x = {box.extents.x, 0, 0};
  Vec3 y = {0, box.extents.y, 0};
  Vec3 z = {0, 0, box.extents.z};

  Vec3 x0 = box.rotation * x;
  Vec3 y0 = box.rotation * y;
  Vec3 z0 = box.rotation * z;

  aabb.extents.x = std::abs(x0[0]) + std::abs(y0[0]) + std::abs(z0[0]);
  aabb.extents.y = std::abs(x0[1]) + std::abs(y0[1]) + std::abs(z0[1]);
  aabb.extents.z = std::abs(x0[2]) + std::abs(y0[2]) + std::abs(z0[2]);

  aabb.child_type       = CollisionObjectType::door;
  aabb.child_ptr        = ptr;
  aabb.point_in_segment = point_in_segment;
}

inline host_fn void create_AABB_from_capsule(const Capsule& capsule, AxisAlignedBoundingBox& aabb, const void* ptr, int point_in_segment = -1) {
  aabb.center     = (capsule.p1 + capsule.p2) / 2;
  Vec3 halfwidths = (capsule.p2 - capsule.p1) / 2;
  aabb.extents.x  = std::abs(halfwidths.x) + capsule.radius;
  aabb.extents.y  = std::abs(halfwidths.y) + capsule.radius;
  aabb.extents.z  = std::abs(halfwidths.z) + capsule.radius;

  aabb.child_type       = CollisionObjectType::capsule;
  aabb.child_ptr        = ptr;
  aabb.point_in_segment = point_in_segment;
}

inline host_fn void create_AABB_from_AABBs(AxisAlignedBoundingBox& b1, AxisAlignedBoundingBox& b2, AxisAlignedBoundingBox& aabb, int id1, int id2) {
  Vec3 a   = b1.center + b1.extents;
  Vec3 b   = b2.center + b2.extents;
  Vec3 max = {std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2])};

  a        = b1.center - b1.extents;
  b        = b2.center - b2.extents;
  Vec3 min = {std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2])};

  aabb.center  = (max + min) * 0.5;
  aabb.extents = (max - min) * 0.5;

  aabb.children[0] = id1;
  aabb.children[1] = id2;

  aabb.child_type = CollisionObjectType::aabb;
}

inline host_fn void find_boxes_to_merge(const BoundingVolumeHierarchy& BVH, std::vector<int>& boxes, int& box1, int& box2, int num_boxes) {
  real closest_dist, dist; // squared distances
  closest_dist = INF_REAL;
  Vec3 diff;

  for (int i = 0; i < num_boxes - 1; i++) {
    for (int j = i + 1; j < num_boxes; j++) {
      diff = BVH.leaves[boxes[i]].center - BVH.leaves[boxes[j]].center;
      dist = dot(diff, diff);
      if (dist < closest_dist) {
        closest_dist = dist;
        box1         = i;
        box2         = j;
      }
    }
  }
}

// Use if worlds regularly contain more than 50 obstacles
inline host_fn void find_boxes_to_merge_sap(const BoundingVolumeHierarchy& BVH, std::vector<int>& boxes, int& box1, int& box2, int num_boxes) {
  if (num_boxes == 2) {
    box1 = 0;
    box2 = 1;
    return;
  }

  std::vector<int> indices(num_boxes);
  for (int i = 0; i < num_boxes; i++) {
    indices[i] = i;
  }

  std::stable_sort(indices.begin(), indices.end(), [&BVH, &boxes](int a, int b) { return BVH.leaves[boxes[a]].center.x < BVH.leaves[boxes[b]].center.x; });

  Vec3 diff;
  real closest_dist = INF_REAL, dist, dx;

  for (int i = 0; i < num_boxes - 1; i++) {
    int idx1 = indices[i];
    for (int j = i + 1; j < num_boxes; j++) {
      int idx2 = indices[j];
      dx       = BVH.leaves[boxes[idx2]].center.x - BVH.leaves[boxes[idx1]].center.x;
      if (dx * dx >= closest_dist) {
        break;
      }
      diff = BVH.leaves[boxes[idx1]].center - BVH.leaves[boxes[idx2]].center;
      dist = dot(diff, diff);

      if (dist < closest_dist) {
        closest_dist = dist;
        box1         = idx1;
        box2         = idx2;
      }
    }
  }
}

inline host_fn void create_bounding_volume_hierarchy(std::vector<int>& objects, BoundingVolumeHierarchy& BVH, int n_leaves) {
  int num = objects.size();
  int i, j, min, max;
  while (num > 1) {
    find_boxes_to_merge(BVH, objects, i, j, num);
    create_AABB_from_AABBs(BVH.leaves[objects[i]], BVH.leaves[objects[j]], BVH.leaves[n_leaves], objects[i], objects[j]);
    n_leaves++;
    min = i;
    max = j;
    if (max < min)
      max = i, min = j;

    // removes merged aabbs and adds the new one
    objects[min] = n_leaves - 1;
    objects[max] = objects[num - 1];
    num--;
  }

  // Assign remaining aabb to BVH root node
  BVH.root = 2 * BVH.num_objects - 2;
}

inline host_fn void create_bounding_volume_hierarchy_sap(std::vector<int>& objects, BoundingVolumeHierarchy& BVH, int n_leaves) {
  int num = objects.size();
  int i, j, min, max;
  while (num > 1) {
    find_boxes_to_merge_sap(BVH, objects, i, j, num);
    create_AABB_from_AABBs(BVH.leaves[objects[i]], BVH.leaves[objects[j]], BVH.leaves[n_leaves], objects[i], objects[j]);
    n_leaves++;
    min = i;
    max = j;
    if (max < min)
      max = i, min = j;

    // removes merged aabbs and adds the new one
    objects[min] = n_leaves - 1;
    objects[max] = objects[num - 1];
    num--;
  }
  // Assign remaining aabb to BVH parent aabb
  BVH.root = 2 * BVH.num_objects - 2;
}

inline host_fn void create_static_bounding_volume_hierarchy(World& world, BoundingVolumeHierarchy& BVH) {
  // Create leaves
  BVH.num_objects = world.boxes.size() + world.spheres.size() + world.capsules.size();
  BVH.leaves.resize(2 * BVH.num_objects - 1);
  std::vector<int> objects;

  int n_leaves = 0;
  // Add objects
  for (auto& box: world.boxes) {
    create_AABB_from_box(box, BVH.leaves[n_leaves], &box);
    n_leaves++;
  }
  for (auto& sphere: world.spheres) {
    create_AABB_from_sphere(sphere, BVH.leaves[n_leaves], &sphere);
    n_leaves++;
  }
  for (auto& capsule: world.capsules) {
    create_AABB_from_capsule(capsule, BVH.leaves[n_leaves], &capsule);
    n_leaves++;
  }

  objects.resize(BVH.num_objects);
  for (int i = 0; i < BVH.num_objects; i++) {
    objects[i] = i;
  }

  create_bounding_volume_hierarchy(objects, BVH, n_leaves);
}

inline host_fn void create_static_bounding_volume_hierarchy_sap(World& world, BoundingVolumeHierarchy& BVH) {
  // Create leaves
  BVH.num_objects = world.boxes.size() + world.spheres.size() + world.capsules.size();
  BVH.leaves.resize(2 * BVH.num_objects - 1);
  std::vector<int> objects;

  int n_leaves = 0;
  // Add objects
  for (auto& box: world.boxes) {
    create_AABB_from_box(box, BVH.leaves[n_leaves], &box);
    n_leaves++;
  }
  for (auto& sphere: world.spheres) {
    create_AABB_from_sphere(sphere, BVH.leaves[n_leaves], &sphere);
    n_leaves++;
  }
  for (auto& capsule: world.capsules) {
    create_AABB_from_capsule(capsule, BVH.leaves[n_leaves], &capsule);
    n_leaves++;
  }

  objects.resize(BVH.num_objects);
  for (int i = 0; i < BVH.num_objects; i++) {
    objects[i] = i;
  }

  create_bounding_volume_hierarchy_sap(objects, BVH, n_leaves);
}

inline host_fn void create_dynamic_bounding_volume_hierarchy(World& world, BoundingVolumeHierarchy& BVH, real time = 0.0) {
  // Create leaves
  BVH.num_objects = world.dynamic_boxes.size() + world.dynamic_spheres.size() + world.dynamic_capsules.size() + world.dynamic_doors.size();
  BVH.time        = time;
  if (BVH.num_objects == 0)
    return;

  BVH.leaves.resize(2 * BVH.num_objects - 1);
  std::vector<int> objects;

  int n_leaves = 0;
  // Add dynamic objects
  for (auto& box: world.dynamic_boxes) {
    create_AABB_from_box(box.lookup(time), BVH.leaves[n_leaves], &box);
    n_leaves++;
  }
  for (auto& sphere: world.dynamic_spheres) {
    create_AABB_from_sphere(sphere.lookup(time), BVH.leaves[n_leaves], &sphere);
    n_leaves++;
  }
  for (auto& capsule: world.dynamic_capsules) {
    create_AABB_from_capsule(capsule.lookup(time), BVH.leaves[n_leaves], &capsule);
    n_leaves++;
  }
  for (auto& door: world.dynamic_doors) {
    create_AABB_from_box(door.lookup(time), BVH.leaves[n_leaves], &door);
    BVH.leaves[n_leaves].child_type = CollisionObjectType::door;
    n_leaves++;
  }

  objects.resize(BVH.num_objects);
  for (int i = 0; i < BVH.num_objects; i++) {
    objects[i] = i;
  }

  create_bounding_volume_hierarchy(objects, BVH, n_leaves);
}

inline host_fn void create_time_bvh_dynamic_objects(World& world, BoundingVolumeHierarchy& BVH, int start_point_in_segment, int n_points_per_segment,
                                                    real opt_time, int n_segments, int trajectory_start_time) {
  // Create leaves
  BVH.num_objects = (world.dynamic_boxes.size() + world.dynamic_spheres.size() + world.dynamic_capsules.size() + world.dynamic_doors.size()) * n_points_per_segment;

  if (BVH.num_objects == 0)
    return;

  BVH.leaves.resize(2 * BVH.num_objects - 1);
  std::vector<int> objects;
  int              max_point = n_segments * n_points_per_segment - 1;

  int n_leaves = 0;
  // Add objects
  for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
    int  current_point = start_point_in_segment + point_in_segment;
    real current_time  = opt_time * ((real) current_point / (real) max_point) + trajectory_start_time; // trajectory time * progression along trajectory

    for (const auto& box: world.dynamic_boxes) {
      create_AABB_from_box(box.lookup(current_time), BVH.leaves[n_leaves], &box, point_in_segment);
      n_leaves++;
    }
    for (const auto& sphere: world.dynamic_spheres) {
      create_AABB_from_sphere(sphere.lookup(current_time), BVH.leaves[n_leaves], &sphere, point_in_segment);
      n_leaves++;
    }
    for (const auto& capsule: world.dynamic_capsules) {
      create_AABB_from_capsule(capsule.lookup(current_time), BVH.leaves[n_leaves], &capsule, point_in_segment);
      n_leaves++;
    }
    for (const auto& door: world.dynamic_doors) {
      create_AABB_from_door(door.lookup(current_time), BVH.leaves[n_leaves], &door, point_in_segment);
      BVH.leaves[n_leaves].child_type = CollisionObjectType::door;
      n_leaves++;
    }
  }

  objects.resize(BVH.num_objects);
  for (int i = 0; i < BVH.num_objects; i++) {
    objects[i] = i;
  }

  create_bounding_volume_hierarchy_sap(objects, BVH, n_leaves);
}

// One BVH per capsule because one collision constraint needed per capsule
inline host_fn void create_time_bounding_volume_hierarchies(std::array<BoundingVolumeHierarchy, MAX_CAPSULES>& BVH, std::vector<std::array<Capsule, MAX_CAPSULES>>& capsules,
                                                            const int num_capsules, const int num_points_in_segment) {
  std::vector<int> objects(num_points_in_segment);

  for (int caps_id = 0; caps_id < num_capsules; caps_id++) {
    BVH[caps_id].num_objects = num_points_in_segment;
    BVH[caps_id].leaves.resize(2 * num_points_in_segment - 1);

    for (int point = 0; point < num_points_in_segment; point++) {
      create_AABB_from_capsule(capsules[point][caps_id], BVH[caps_id].leaves[point], &capsules[point][caps_id]);
      BVH[caps_id].leaves[point].point_in_segment = point;
      objects[point]                              = point;
    }

    create_bounding_volume_hierarchy(objects, BVH[caps_id], num_points_in_segment);
  }
}

inline host_fn void minimum_distance(const Capsule& capsule, BoundingVolumeHierarchy& BVH, real& dist_min, CollisionEntities& collision_objects, int point_in_segment) {
  BVH.queue.clear_and_reserve(BVH.num_objects);

  real                   dist = INF_REAL;
  AxisAlignedBoundingBox arm;
  create_AABB_from_capsule(capsule, arm, &capsule);

  std::vector<AxisAlignedBoundingBox>& leaves = BVH.leaves;
  AxisAlignedBoundingBox*              object = &leaves[BVH.root];

  leaves[object->children[0]].dist = distance(arm, leaves[object->children[0]]);
  leaves[object->children[1]].dist = distance(arm, leaves[object->children[1]]);
  BVH.queue.push(object->children[0]);
  BVH.queue.push(object->children[1]);

  while (BVH.queue.size() > 0) {
    object = &leaves[BVH.queue.top()];
    BVH.queue.pop();

    // STOPPING CRITERIA: NO FURTHER IMPROVEMENT POSSIBLE
    if (object->dist >= dist_min) {
      return;
    }

    switch (object->child_type) {
      case CollisionObjectType::aabb: {
        leaves[object->children[0]].dist = distance(arm, leaves[object->children[0]]);
        leaves[object->children[1]].dist = distance(arm, leaves[object->children[1]]);
        BVH.queue.push(object->children[0]);
        BVH.queue.push(object->children[1]);
        break;
      }
      case CollisionObjectType::box: {
        Box& box = *((Box*) object->child_ptr);
        dist     = distance(capsule, box);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::box;
          collision_objects.box               = box;
          collision_objects.point_in_segment  = point_in_segment;
        }
        break;
      }
      case CollisionObjectType::sphere: {
        Sphere& sphere = *((Sphere*) object->child_ptr);
        dist           = distance(capsule, sphere); // MAYBE: make function for aabb/sphere dist_min
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::sphere;
          collision_objects.sphere            = sphere;
          collision_objects.point_in_segment  = point_in_segment;
        }
        break;
      }
      case CollisionObjectType::capsule: {
        Capsule& caps = *((Capsule*) object->child_ptr);
        dist          = distance(capsule, caps);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::capsule;
          collision_objects.capsule           = caps;
          collision_objects.point_in_segment  = point_in_segment;
        }
        break;
      }
    }
  }
}

inline host_fn void minimum_distance_dynamic(const Capsule& capsule, BoundingVolumeHierarchy& BVH, real& dist_min, CollisionEntities& collision_objects, int point_in_segment) {
  BVH.queue.clear_and_reserve(BVH.num_objects);

  real                   dist = INF_REAL;
  AxisAlignedBoundingBox arm;
  create_AABB_from_capsule(capsule, arm, &capsule);

  std::vector<AxisAlignedBoundingBox>& leaves = BVH.leaves;
  AxisAlignedBoundingBox*              object = &leaves[BVH.root];

  leaves[object->children[0]].dist = distance(arm, leaves[object->children[0]]);
  leaves[object->children[1]].dist = distance(arm, leaves[object->children[1]]);
  BVH.queue.push(object->children[0]);
  BVH.queue.push(object->children[1]);

  while (BVH.queue.size() > 0) {
    object = &leaves[BVH.queue.top()];
    BVH.queue.pop();

    // STOPPING CRITERIA: NO FURTHER IMPROVEMENT POSSIBLE
    if (object->dist >= dist_min) {
      return;
    }

    switch (object->child_type) {
      case CollisionObjectType::aabb: {
        leaves[object->children[0]].dist = distance(arm, leaves[object->children[0]]);
        leaves[object->children[1]].dist = distance(arm, leaves[object->children[1]]);
        BVH.queue.push(object->children[0]);
        BVH.queue.push(object->children[1]);
        break;
      }
      case CollisionObjectType::box: {
        Box box = ((DynamicBox*) object->child_ptr)->lookup(BVH.time);
        dist    = distance(capsule, box);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::box;
          collision_objects.box               = box;
          collision_objects.point_in_segment  = point_in_segment;
        }
        break;
      }
      case CollisionObjectType::sphere: {
        Sphere sphere = ((DynamicSphere*) object->child_ptr)->lookup(BVH.time);
        dist          = distance(capsule, sphere);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::sphere;
          collision_objects.sphere            = sphere;
          collision_objects.point_in_segment  = point_in_segment;
        }
        break;
      }
      case CollisionObjectType::capsule: {
        Capsule caps = ((DynamicCapsule*) object->child_ptr)->lookup(BVH.time);
        dist         = distance(capsule, caps);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::capsule;
          collision_objects.capsule           = caps;
          collision_objects.point_in_segment  = point_in_segment;
        }
        break;
      }
      case CollisionObjectType::door: {
        Box door = ((DynamicDoor*) object->child_ptr)->lookup(BVH.time);
        dist     = distance(capsule, door);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::box;
          collision_objects.box               = door;
          collision_objects.point_in_segment  = point_in_segment;
        }
      }
    }
  }
}

inline host_fn void minimum_distance_static_objects_time(BoundingVolumeHierarchy& BVH_obj, BoundingVolumeHierarchy& BVH_time, real& dist_min,
                                                         CollisionEntities& collision_objects) {
  QueuePair queue(BVH_obj.num_objects * BVH_time.num_objects * 4);
  dist_min  = INF_REAL;
  real dist = INF_REAL;

  std::vector<AxisAlignedBoundingBox>& leaves_obj  = BVH_obj.leaves;
  std::vector<AxisAlignedBoundingBox>& leaves_time = BVH_time.leaves;

  AxisAlignedBoundingBox* object1 = &leaves_obj[BVH_obj.root];
  AxisAlignedBoundingBox* object2 = &leaves_time[BVH_time.root];

  AABBPair pair;

  // Test BVHs' children against eachother and add to priority queue
  pair.aabb_obj = object1->children[0];
  pair.aabb_cap = object2->children[0];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_cap = object2->children[1];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_obj = object1->children[1];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_cap = object2->children[0];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  CollisionObjectType child_type1;
  CollisionObjectType child_type2;

  while (queue.size() > 0) {
    pair        = queue.top();
    object1     = &leaves_obj[pair.aabb_obj];
    child_type1 = object1->child_type;
    object2     = &leaves_time[pair.aabb_cap];
    child_type2 = object2->child_type;
    queue.pop();

    // STOPPING CRITERIA: NO FURTHER IMPROVEMENT POSSIBLE
    if (pair.dist >= dist_min) {
      return;
    }

    if (child_type1 == CollisionObjectType::aabb && child_type2 == CollisionObjectType::capsule) {
      pair.aabb_obj = object1->children[0];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);

      pair.aabb_obj = object1->children[1];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);
    } else if (child_type1 == CollisionObjectType::box && child_type2 == CollisionObjectType::capsule) {
      Box&           box  = *((Box*) object1->child_ptr);
      const Capsule& caps = *((Capsule*) object2->child_ptr);
      dist                = distance(caps, box);
      if (dist < dist_min) {
        dist_min                            = dist;
        collision_objects.other_object_type = CollisionObjectType::box;
        collision_objects.box               = box;
        collision_objects.point_in_segment  = object2->point_in_segment;
      }
    } else if (child_type1 == CollisionObjectType::aabb && child_type2 == CollisionObjectType::aabb) {
      // Faster to only enter bigger AABB
      if ((object1->extents.x * object1->extents.y * object1->extents.z) > (object2->extents.x * object2->extents.y * object2->extents.z)) {
        pair.aabb_obj = object1->children[0];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);

        pair.aabb_obj = object1->children[1];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);
      } else {
        pair.aabb_cap = object2->children[0];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);

        pair.aabb_cap = object2->children[1];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);
      }
    } else if (child_type1 != CollisionObjectType::aabb && child_type2 == CollisionObjectType::aabb) {
      pair.aabb_cap = object2->children[0];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);

      pair.aabb_cap = object2->children[1];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);
    } else if (child_type1 == CollisionObjectType::capsule && child_type2 == CollisionObjectType::capsule) {
      Capsule& caps1 = *((Capsule*) object1->child_ptr);
      Capsule& caps2 = *((Capsule*) object2->child_ptr);
      dist           = distance(caps1, caps2);
      if (dist < dist_min) {
        dist_min                            = dist;
        collision_objects.other_object_type = CollisionObjectType::capsule;
        collision_objects.capsule           = caps1;
        collision_objects.point_in_segment  = object2->point_in_segment;
      }
    } else if (child_type1 == CollisionObjectType::sphere && child_type2 == CollisionObjectType::capsule) {
      Sphere&  sphere = *((Sphere*) object1->child_ptr);
      Capsule& caps   = *((Capsule*) object2->child_ptr);
      dist            = distance(caps, sphere);
      if (dist < dist_min) {
        dist_min                            = dist;
        collision_objects.other_object_type = CollisionObjectType::sphere;
        collision_objects.sphere            = sphere;
        collision_objects.point_in_segment  = object2->point_in_segment;
      }
    }
  }
}

inline host_fn void minimum_distance_dynamic_objects_time(BoundingVolumeHierarchy& BVH_obj, BoundingVolumeHierarchy& BVH_time, real& dist_min,
                                                          CollisionEntities& collision_objects) {
  QueuePair queue(BVH_obj.num_objects * BVH_time.num_objects * 4);

  real dist = INF_REAL;

  std::vector<AxisAlignedBoundingBox>& leaves_obj  = BVH_obj.leaves;
  std::vector<AxisAlignedBoundingBox>& leaves_time = BVH_time.leaves;

  AxisAlignedBoundingBox* object1 = &leaves_obj[BVH_obj.root];
  AxisAlignedBoundingBox* object2 = &leaves_time[BVH_time.root];

  AABBPair pair;

  // Test BVHs' children against eachother and add to priority queue
  pair.aabb_obj = object1->children[0];
  pair.aabb_cap = object2->children[0];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_cap = object2->children[1];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_obj = object1->children[1];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_cap = object2->children[0];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);
  queue.push(pair);

  CollisionObjectType child_type1;
  CollisionObjectType child_type2;

  while (queue.size() > 0) {
    pair        = queue.top();
    object1     = &leaves_obj[pair.aabb_obj];
    child_type1 = object1->child_type;
    object2     = &leaves_time[pair.aabb_cap];
    child_type2 = object2->child_type;
    queue.pop();

    // STOPPING CRITERIA: NO FURTHER IMPROVEMENT POSSIBLE
    if (pair.dist >= dist_min) {
      return;
    }

    if (child_type1 == CollisionObjectType::aabb && child_type2 == CollisionObjectType::capsule) {
      pair.aabb_obj = object1->children[0];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);

      pair.aabb_obj = object1->children[1];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);
    } else if (child_type1 == CollisionObjectType::box && child_type2 == CollisionObjectType::capsule && object1->point_in_segment == object2->point_in_segment) {
      Box&     box  = *((Box*) object1->child_ptr);
      Capsule& caps = *((Capsule*) object2->child_ptr);
      dist          = distance(caps, box);
      if (dist < dist_min) {
        dist_min                            = dist;
        collision_objects.other_object_type = CollisionObjectType::box;
        collision_objects.box               = box;
        collision_objects.point_in_segment  = object2->point_in_segment;
      }
    } else if (child_type1 == CollisionObjectType::aabb && child_type2 == CollisionObjectType::aabb) {
      if ((object1->extents.x * object1->extents.y * object1->extents.z) > (object2->extents.x * object2->extents.y * object2->extents.z)) {
        pair.aabb_obj = object1->children[0];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);

        pair.aabb_obj = object1->children[1];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);
      } else {
        pair.aabb_cap = object2->children[0];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);

        pair.aabb_cap = object2->children[1];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);
      }
    } else if (child_type1 != CollisionObjectType::aabb && child_type2 == CollisionObjectType::aabb) {
      pair.aabb_cap = object2->children[0];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);

      pair.aabb_cap = object2->children[1];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);
    } else if (child_type1 == CollisionObjectType::capsule && child_type2 == CollisionObjectType::capsule && object1->point_in_segment == object2->point_in_segment) {
      Capsule& caps1 = *((Capsule*) object1->child_ptr);
      Capsule& caps2 = *((Capsule*) object2->child_ptr);
      dist           = distance(caps1, caps2);
      if (dist < dist_min) {
        dist_min                            = dist;
        collision_objects.other_object_type = CollisionObjectType::capsule;
        collision_objects.capsule           = caps1;
        collision_objects.point_in_segment  = object2->point_in_segment;
      }
    } else if (child_type1 == CollisionObjectType::sphere && child_type2 == CollisionObjectType::capsule && object1->point_in_segment == object2->point_in_segment) {
      Sphere&  sphere = *((Sphere*) object1->child_ptr);
      Capsule& caps   = *((Capsule*) object2->child_ptr);
      dist            = distance(caps, sphere);
      if (dist < dist_min) {
        dist_min                            = dist;
        collision_objects.other_object_type = CollisionObjectType::sphere;
        collision_objects.sphere            = sphere;
        collision_objects.point_in_segment  = object2->point_in_segment;
      }
    }
  }
}

inline host_fn void minimum_distance_double_time_bvh(BoundingVolumeHierarchy& BVH_obj, BoundingVolumeHierarchy& BVH_time, real& dist_min, CollisionEntities& collision_objects,
                                                     int start_point_in_segment, int n_points_per_segment, real opt_time, int n_segments, int trajectory_start_time) {

  QueuePair queue(BVH_obj.num_objects * BVH_time.num_objects * 4);

  real dist = INF_REAL;

  std::vector<AxisAlignedBoundingBox>& leaves_obj  = BVH_obj.leaves;
  std::vector<AxisAlignedBoundingBox>& leaves_time = BVH_time.leaves;

  AxisAlignedBoundingBox* object1 = &leaves_obj[BVH_obj.root];
  AxisAlignedBoundingBox* object2 = &leaves_time[BVH_time.root];

  AABBPair pair;

  // Test BVHs' children against eachother and add to priority queue
  pair.aabb_obj = object1->children[0];
  pair.aabb_cap = object2->children[0];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_cap = object2->children[1];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_obj = object1->children[1];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);

  pair.aabb_cap = object2->children[0];
  pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
  queue.push(pair);
  queue.push(pair);

  CollisionObjectType child_type1;
  CollisionObjectType child_type2;

  while (queue.size() > 0) {
    pair        = queue.top();
    object1     = &leaves_obj[pair.aabb_obj];
    child_type1 = object1->child_type;
    object2     = &leaves_time[pair.aabb_cap];
    child_type2 = object2->child_type;
    queue.pop();

    // STOPPING CRITERIA: NO FURTHER IMPROVEMENT POSSIBLE
    if (pair.dist >= dist_min) {
      return;
    }

    if (child_type1 == CollisionObjectType::aabb && child_type2 == CollisionObjectType::capsule) {
      pair.aabb_obj = object1->children[0];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);

      pair.aabb_obj = object1->children[1];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);
    } else if (child_type1 == CollisionObjectType::aabb && child_type2 == CollisionObjectType::aabb) {
      // Faster to only enter bigger AABB
      if ((object1->extents.x * object1->extents.y * object1->extents.z) > (object2->extents.x * object2->extents.y * object2->extents.z)) {
        pair.aabb_obj = object1->children[0];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);

        pair.aabb_obj = object1->children[1];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);
      } else {
        pair.aabb_cap = object2->children[0];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);

        pair.aabb_cap = object2->children[1];
        pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
        queue.push(pair);
      }
      // }
    } else if (child_type1 != CollisionObjectType::aabb && child_type2 == CollisionObjectType::aabb) {
      pair.aabb_cap = object2->children[0];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);

      pair.aabb_cap = object2->children[1];
      pair.dist     = distance(leaves_obj[pair.aabb_obj], leaves_time[pair.aabb_cap]);
      queue.push(pair);
    } else if (leaves_obj[pair.aabb_obj].point_in_segment == leaves_time[pair.aabb_cap].point_in_segment) {
      int  current_point = start_point_in_segment + leaves_obj[pair.aabb_obj].point_in_segment;
      int  max_point     = n_segments * n_points_per_segment - 1;
      real current_time  = opt_time * ((real) current_point / (real) max_point) + trajectory_start_time;

      if (child_type1 == CollisionObjectType::box && child_type2 == CollisionObjectType::capsule) {
        Box            box  = ((DynamicBox*) object1->child_ptr)->lookup(current_time);
        const Capsule& caps = *((Capsule*) object2->child_ptr);
        dist                = distance(caps, box);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::box;
          collision_objects.box               = box;
          collision_objects.point_in_segment  = object2->point_in_segment;
        }
      } else if (child_type1 == CollisionObjectType::sphere && child_type2 == CollisionObjectType::capsule) {
        Sphere   sphere = ((DynamicSphere*) object1->child_ptr)->lookup(current_time);
        Capsule& caps   = *((Capsule*) object2->child_ptr);
        dist            = distance(caps, sphere);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::sphere;
          collision_objects.sphere            = sphere;
          collision_objects.point_in_segment  = object2->point_in_segment;
        }
      } else if (child_type1 == CollisionObjectType::capsule && child_type2 == CollisionObjectType::capsule) {
        Capsule  caps1 = ((DynamicCapsule*) object1->child_ptr)->lookup(current_time);
        Capsule& caps2 = *((Capsule*) object2->child_ptr);
        dist           = distance(caps1, caps2);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::capsule;
          collision_objects.capsule           = caps1;
          collision_objects.point_in_segment  = object2->point_in_segment;
        }
      } else if (child_type1 == CollisionObjectType::door && child_type2 == CollisionObjectType::capsule) {
        Box      door = ((DynamicDoor*) object1->child_ptr)->lookup(current_time);
        Capsule& caps = *((Capsule*) object2->child_ptr);
        dist          = distance(caps, door);
        if (dist < dist_min) {
          dist_min                            = dist;
          collision_objects.other_object_type = CollisionObjectType::capsule;
          collision_objects.box               = door;
          collision_objects.point_in_segment  = object2->point_in_segment;
        }
      }
    }
  }
}
} // namespace blast
