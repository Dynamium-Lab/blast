# Blast

    1: Make it work
    2: Make it faster
    3: repeat step 2

`Blast` is a fast robot manipulator motion optimization library that generates smooth, collision-free trajectories that respect kinematic constraints, dynamic constraints, and custom constraints. 
## Usage
Once built and linked, simply include the main `blast.h` header file. All functions and types are namespaced with `blast::`.

```c++
#include "blast.h"

int main() {
    blast::Vec3 a = {1,2,3};
    blast::Vec3 b = {4,3,2};
    blast::print(a+b);

    return 0;
}
```

## Examples
A set of examples are provided in the `/examples` folder.

## Download
### Latest release
You can find the latest release from the [Releases](https://github.com/Dynamium-Lab/blast/releases) page. Download and extract the latest version inside your project, ex: `extern/blast/`.

### As a git submodule
To add `Blast` to a project as a submodule, run the following command from your git repository:
```sh
git submodule add https://github.com/Dynamium-Lab/blast.git <path>/blast/
```
replace `<path>` with the path to where you store dependancies in your project, ex:
```sh
git submodule add https://github.com/Dynamium-Lab/blast.git extern/blast/
```


## Build

### Cmake (easiest)
The easiest method to build and use `Blast` is to use `cmake`. Add the following to your `CMakeLists.txt` file:
```Cmake
# replace <path> with the path to where you stored 'Blast'
add_subdirectory(<path>/blast)

# replace <target> with appropriate target name
add_executable(<target> example.cpp)

# link blast to the new target
target_link_libraries(<target> PRIVATE blast nlopt)
```

#### Build your application (windows + MSVC):
```sh
cmake -B build
cmake --build build --config Release
```

#### Build your application (almost every other build system):
```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Manual build (more versatile)
`Blast`'s only dependancy is [NLopt](https://nlopt.readthedocs.io/en/latest/). If `NLopt` is installed on the system, make sure to link to it. If not, a working version of `NLopt` ships with the repository in `extern/nlopt`. Compile `.cpp` files in their `src/algs` and `src/api` folders.

All files for `Blast` are in the `blast` subfolder. All `.cpp` files must be built. The `blast/CMakeLists.txt` file can help you locate them.


