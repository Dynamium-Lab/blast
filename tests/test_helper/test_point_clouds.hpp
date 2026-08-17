#pragma once
#include <blast>

namespace blast {
PointCloud point_cloud_from_box(Box box, int n_extra = 0) {
  PointCloud cloud;

  // Initialization of the eight OBB points
  Vec3 size_x_org = {box.extents.x, 0, 0};
  Vec3 size_y_org = {0, box.extents.y, 0};
  Vec3 size_z_org = {0, 0, box.extents.z};
  Vec3 size_x     = box.rotation * size_x_org;
  Vec3 size_y     = box.rotation * size_y_org;
  Vec3 size_z     = box.rotation * size_z_org;

  std::vector<Vec3> v2(8);
  v2[0] = box.center + size_x + size_y + size_z;
  v2[1] = box.center + size_x + size_y - size_z;
  v2[2] = box.center + size_x - size_y + size_z;
  v2[3] = box.center + size_x - size_y - size_z;
  v2[4] = box.center - size_x + size_y + size_z;
  v2[5] = box.center - size_x + size_y - size_z;
  v2[6] = box.center - size_x - size_y + size_z;
  v2[7] = box.center - size_x - size_y - size_z;


  // Add random points inside box
  cloud.points = v2;
  if (n_extra) {
    cloud.points.reserve(cloud.points.size() + n_extra);
    for (int i = 0; i < n_extra; i++) {
      cloud.points.push_back(box.center + random_real() * size_x + random_real() * size_y + random_real() * size_z);
    }
  }

  return cloud;
}

std::vector<std::pair<Box, PointCloud>> get_boxes_and_point_clouds(World world) {
  std::vector<std::pair<Box, PointCloud>> pairs;

  for (auto box: world.boxes) {
    pairs.push_back({box, point_cloud_from_box(box)});
  }

  return pairs;
}
} // namespace blast
