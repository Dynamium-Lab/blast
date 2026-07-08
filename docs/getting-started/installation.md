# Installation

## Requirements

- CMake 3.21+
- A C++17 compiler (GCC 9+, Clang 10+, MSVC 2019+)
- x86-64 with at least SSE4 (AVX2 recommended for full performance)

No other system dependencies are required — NLopt ships with the repository.

## C++

### Option A — `find_package` (recommended for installed use)

Configure and install Blast once:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build
cmake --install build
```

Then, in your project's `CMakeLists.txt`:

```cmake
find_package(Blast REQUIRED)
target_link_libraries(my_app PRIVATE Blast::blast)
```

### Option B — `add_subdirectory` (submodule / embedded)

```sh
git submodule add https://github.com/Dynamium-Lab/blast.git extern/blast
```

```cmake
add_subdirectory(extern/blast)
target_link_libraries(my_app PRIVATE blast)
```

Developer targets (examples, tests, benchmarks) are automatically excluded when Blast is
consumed this way.

## Python

Requires Python 3.9+, a C++17 compiler, and CMake 3.21+ on your `PATH`.

Create an isolated environment first:

```sh
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
```

Then build and install the bindings:

```sh
pip install .
```

This compiles the native extension and installs the `blast` package into the active
environment.

## Updating

**`find_package` install** — pull the latest, rebuild, and reinstall:

```sh
git pull
cmake --build build
cmake --install build
```

**Submodule** — update to the latest commit on `main`:

```sh
git submodule update --remote extern/blast
git add extern/blast
git commit -m "Update blast"
```
