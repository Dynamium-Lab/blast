#pragma once
#include <blast>

namespace blast {

struct AxisAlignedBoundingBox {
    Vec3 center;
    Vec3 dimensions;
    std::vector<AxisAlignedBoundingBox*> children;
};

struct BoundingVolumeHierarchy {
    std::vector<AxisAlignedBoundingBox> boxes;
    int top_box;
};

BoundingVolumeHierarchy create_bounding_volume_hierarchy(World world) {

}

} // namespace blast