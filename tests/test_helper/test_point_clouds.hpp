#pragma once
#include <blast>

namespace blast {
PointCloud point_cloud_from_box(Box box) {
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

  cloud.set = v2;
  // // Add box corners
  // for (int i = 0; i < 8; i++) {
  //   Vec3 rotated_corner = box.rotation * corners[i];
  //   cloud.set.push_back(box.center + rotated_corner);
  // }

  // // Add random points inside box
  // for (int i = 0; i < 0; i++) { // can be up to 8 extra points
  //   Vec3 scaled_corner = corners[i];
  //   for (int j = 0; j < 3; j++) {
  //     scaled_corner[j] *= random_real();
  //   }
  //   scaled_corner = box.rotation * scaled_corner;
  //   cloud.set.push_back(box.center + scaled_corner);
  // }

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
