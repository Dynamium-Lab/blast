# Blast
A motion optimization header library that generates smooth trajectories that respect kinematic and dynamic constraints. 

## How to use Blast
To use blast, it is not necessary to 'install' anything. The library lives in headers and can be compiled with the build tools of your choice simply by including the main Blast header into your source files:

```c++
#include "blast.hpp"
```

Optionally, some utility functions are provided in the `/optional` folder:

```c++
#include "optional/blast_optional_utilities.hpp"
```

To use in a project, simply copy the `/blast` folder to your project or point your compiler to it.

__Note: There is a bug in the gcc compiler for windows (MinGW-w64) that makes it difficult to use AVX instructions. MinGW is therefore not recommended for Windows.__

## Examples
A set of examples are provided in the `/examples` folder. To build the examples, it is recommended, but not required to use CMake. For convenience, a build script that uses CMake is provided:

```
./build_examples.sh
```

This script builds a debug version of the examples.

## Tests and Benchmarks
Internally, tests and benchmarks are used to ensure performance of the library. These can be found in the `/tests` and `/benchmarks` folders, respectively. To build these, CMake is prefered. For convenience, `build_tests.sh` and `build_benchmarks.sh` scripts can be used.

The __tests__ use the [GoogleTest](https://github.com/google/googletest) framework. If using cmake, the Google Test framework will be downloaded and compiled automatically. Otherwise, it must be installed and configured manually.

The __benchmarks__ use the [Google Benchmark](https://github.com/google/benchmark) framework. As far as we know, this must be installed manually following this [page](https://github.com/google/benchmark#installation). The scripts here presume that this has been done in the `../lib/benchmark` folder. This hardcoded path will be fixed in a future release.

