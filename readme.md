# Blast
A fast robot manipulator motion optimization library that generates smooth, collision-free trajectories that respect kinematic constraints, dynamic constraints, and custom constraints. 

# Download
## Releases
You can find the latest release from the [Releases](https://github.com/Dynamium-Lab/blast/releases) page. Download and extract the latest version inside your project, ex: `extern/blast/`.

## Git submodule
To add `Blast` to a project as a submodule, run the following command from your git repository:
```sh
git submodule add https://github.com/Dynamium-Lab/blast.git <path>/blast/
```
replace `<path>` with the path to where you store dependancies in your project, ex:
```sh
git submodule add https://github.com/Dynamium-Lab/blast.git extern/blast/
```


# Build

## Cmake
The easiest method to build and use `Blast` is to use `cmake`. Add the following to your `CMakeLists.txt` file:
```Cmake
add_subdirectory(<path>/blast)  # replace <path> with the path to where you stored 'Blast'.
```

Then, to link `Blast` to a target:
```Cmake
target_link_libraries(<target> PRIVATE blast)  # replace <target> with the target name given to add_executable()
```

## Manual build (less recommended)
`Blast`'s only dependancy is [NLopt](https://nlopt.readthedocs.io/en/latest/). If `NLopt` is installed on the system, make sure to link to it.

All files are in the `blast` subfolder. All `*.cpp` files must be built. The `blast/CMakeLists.txt` file can help you locate them.

# Usage
Once built and linked, everything is in a namespace called `blast::` in the main `blast.h` header file.
```c++
#include "blast.h"

int main() {
    blast::Vec3 a = {1,2,3};
    blast::Vec3 b = {4,3,2};
    blast::print(a+b);

    return 0;
}
```

# Examples
A set of examples are provided in the `/examples` folder.
