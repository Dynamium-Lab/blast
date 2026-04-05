// Minimal consumer — verifies that the installed Blast headers and the
// exported nlopt target are all reachable via find_package(Blast).
#include <blast>
#include <iostream>

int main() {
  // Exercise a type from blast_math.hpp (included transitively via <blast>)
  blast::Vec3 v{1.0, 2.0, 3.0};
  std::cout << "Blast::Vec3 constructed: "
            << v.x << " " << v.y << " " << v.z << "\n";
  return 0;
}
