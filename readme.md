# Blast

    1: Make it work
    2: Make it faster
    3: repeat step 2

**Blast** is a fast trajectory optimization library for robot manipulators. Given a robot model and a motion task (start pose, goal pose, and constraints), Blast finds a smooth, time-optimal, collision-free B-spline trajectory that satisfies joint position, velocity, and acceleration limits.

Blast is the software accompanying the ICRA 2026 paper *"Accelerating Trajectory Optimization by Exploiting B-Spline Gradient Structure"*.

**Full documentation:** [docs.dynamium.ca](https://docs.dynamium.ca)

---

## Requirements

- CMake 3.21+
- A C++17 compiler (GCC 9+, Clang 10+, MSVC 2019+)
- x86-64 with at least SSE4 (AVX2 recommended for full performance)

No other system dependencies are required — NLopt ships with the repository.

---

## Installation

### Option A — find_package (recommended for installed use)

Configure and install Blast once:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build
cmake --install build
```

In your project's `CMakeLists.txt`:

```cmake
find_package(Blast REQUIRED)
target_link_libraries(my_app PRIVATE Blast::blast)
```

### Option B — add_subdirectory (submodule / embedded)

```sh
git submodule add https://github.com/Dynamium-Lab/blast.git extern/blast
```

```cmake
add_subdirectory(extern/blast)
target_link_libraries(my_app PRIVATE blast)
```

Developer targets (examples, tests, benchmarks) are automatically excluded when Blast is consumed this way.

### Updating

**find_package install** — pull the latest, rebuild, and reinstall:

```sh
git pull
cmake --build build
cmake --install build
```

**Submodule** — update the submodule to the latest commit on `main`:

```sh
git submodule update --remote extern/blast
git add extern/blast
git commit -m "Update blast"
```

### Building examples locally

```sh
cmake --preset dev        # configure (RelWithDebInfo + -march=native)
cmake --build --preset dev
./build/dev/examples/example_01_forward_kinematics
```

### Option C — Python bindings

Requires Python 3.9+, a C++17 compiler, and CMake 3.21+ on your `PATH`.

Create an isolated environment first so the build doesn't affect your system packages:

```sh
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
```

Then build and install the bindings:

```sh
pip install .
```

This compiles the native extension and installs the `blast` package into the active environment.

---

## Quick start

```cpp
#define BLAST_USE_NATIVE_SQP  // use the bundled SQP solver
#include <blast>
#include <iostream>

int main() {
    using namespace blast;

    // 1. Build the robot model (UR5e shown; define your own Manipulator similarly)
    Manipulator robot = make_ur5e();

    // 2. Define the motion task: stop-to-stop from start to goal joint angles
    Task task = Task::stop_to_stop(
        {0.0, -1.57, 1.57, -1.57, -1.57, 0.0},  // start (rad)
        {1.57, -1.0,  1.0, -1.57, -1.57, 0.5}    // goal  (rad)
    );

    // 3. Set up the optimization problem
    Optimization opt(robot, task);
    opt.bspline              = Bspline(16, 110, 5, robot.n_joints);
    opt.constraints.position     = true;
    opt.constraints.velocity     = true;
    opt.constraints.acceleration = true;
    opt.success_tolerance    = 0.01;
    opt.guess.type           = Guess::random;
    opt.guess.n_random_shots = 50;

    // 4. Solve
    Result result = optimize(&opt);

    std::cout << (result.success ? "Success" : "Failed") << "\n";
    std::cout << "Duration: " << result.x.back() << " s\n";
    std::cout << "Solve time: " << result.compute_time << " ms\n";
    return 0;
}
```

### Python

The example below uses a UR5e. All numeric values match `examples/ur5e.hpp`.

```python
import numpy as np
import blast

# --- Limits ---
limits = blast.ManipulatorLimits()
limits.position_max     = blast.array([6.2832, 6.2832, 3.1416, 6.2832, 6.2832, 6.2832])
limits.position_min     = -limits.position_max
limits.velocity_max     = blast.array([3.1416] * 6)
limits.acceleration_max = blast.array([13.96]  * 6)
limits.torque_max       = blast.array([150., 150., 150., 28., 28., 28.])
limits.tool_speed_max   = 2.0

# --- Kinematics ---
kin = blast.ManipulatorKinematics()
kin.joint_offsets = [
    blast.Vec3(0, 0, 0),
    blast.Vec3(-0.425, 0, 0),
    blast.Vec3(-0.3922, 0, 0.1333),
    blast.Vec3(0, -0.0997, 0),
    blast.Vec3(0, 0.0996, 0),
    blast.Vec3(0, 0, 0),
]
kin.joint_axes           = [blast.Vec3(0, 0, 1)] * 6
kin.first_joint_position = blast.Vec3(0.0, 0.0, 0.1625)

# Mat3 is column-major: Mat3(col0r0, col0r1, col0r2, col1r0, ...)
rots = [blast.Mat3()] * blast.MAX_JOINTS
rots[0] = blast.Mat3(-1, 0, 0,  0, -1, 0,  0,  0, 1)
rots[1] = blast.Mat3( 1, 0, 0,  0,  0, 1,  0, -1, 0)
rots[2] = blast.Mat3( 1, 0, 0,  0,  1, 0,  0,  0, 1)
rots[3] = blast.Mat3( 1, 0, 0,  0,  1, 0,  0,  0, 1)
rots[4] = blast.Mat3( 1, 0, 0,  0,  0, 1,  0, -1, 0)
rots[5] = blast.Mat3( 1, 0, 0,  0,  0, -1, 0,  1, 0)
kin.static_rotations = rots + [blast.Mat3()]  # MAX_JOINTS + 1 entries

# --- Dynamics ---
dyn = blast.ManipulatorDynamics()
dyn.link_masses = blast.array([3.7, 8.393, 2.275, 1.219, 1.219, 0.1879, 0, 0])
dyn.inertia_tensors = [
    blast.Mat3(0.0103, 0, 0,  0, 0.0103, 0,  0, 0, 0.0067),
    blast.Mat3(0.1339, 0, 0,  0, 0.1339, 0,  0, 0, 0.0151),
    blast.Mat3(0.0312, 0, 0,  0, 0.0312, 0,  0, 0, 0.0041),
    blast.Mat3(0.0026, 0, 0,  0, 0.0026, 0,  0, 0, 0.0022),
    blast.Mat3(0.0026, 0, 0,  0, 0.0026, 0,  0, 0, 0.0022),
    blast.Mat3(0.0001, 0, 0,  0, 0.0001, 0,  0, 0, 0.0001),
    blast.Mat3(), blast.Mat3(),
]
dyn.cog_offsets = [
    blast.Vec3(0, 0, 0),
    blast.Vec3(-0.2125, 0, 0.138),
    blast.Vec3(-0.1961, 0, 0.007),
    blast.Vec3(0, 0, 0),
    blast.Vec3(0, 0, 0),
    blast.Vec3(0, 0, -0.0229),
    blast.Vec3(0, 0, 0),
    blast.Vec3(0, 0, 0),
]

# --- Collision geometry ---
caps_cfg = blast.ManipulatorCapsules()

base_sphere        = blast.Sphere()
base_sphere.center = blast.Vec3(0, 0, 0.0325)
base_sphere.radius = 0.09
caps_cfg.base_sphere    = base_sphere
caps_cfg.collision_base = blast.array([0, 0, 0, 1, 1, 1, 1])

col_mat = np.zeros((7, 7), dtype=np.uint8)
col_mat[3, 0] = col_mat[4, 0] = col_mat[5, 0] = col_mat[6, 0] = 1
col_mat[4, 1] = col_mat[5, 1] = col_mat[6, 1] = 1
col_mat[6, 3] = 1
caps_cfg.collision_matrix = col_mat

def make_cap(frame, p1, p2, r):
    c = blast.CollisionModelCapsule()
    c.joint_frame = frame
    c.p1 = blast.Vec3(*p1)
    c.p2 = blast.Vec3(*p2)
    c.radius = r
    return c

caps_cfg.capsule_list = [
    make_cap(0, (0, 0, 0),                 (0, -0.15, 0),              0.09),
    make_cap(1, (-0.42, 0, 0.1375),        (0, 0, 0.1375),             0.06),
    make_cap(2, (0, 0, 0.02),              (0, 0, 0.18),               0.065),
    make_cap(2, (-0.373156, 0, 0.0085),    (0.000441, 0.000121, 0.0085), 0.05),
    make_cap(3, (0, 0, 0),                 (0, 0, -0.155),             0.0425),
    make_cap(3, (0, 0.03, 0),              (0, -0.1, 0),               0.0425),
    make_cap(5, (0, 0, -0.14),             (0, 0, -0.01),              0.038),
]

# --- Build robot ---
robot = blast.Manipulator(6, limits, kin, dyn, caps_cfg)

# --- Task ---
start = blast.array([0.0, -1.57,  1.57, -1.57, -1.57,  0.0])
goal  = blast.array([1.57, -1.0,  1.0,  -1.57, -1.57,  0.5])
task  = blast.Task.stop_to_stop(start, goal)

# --- Optimize ---
opt = blast.Optimization(robot, task)
opt.bspline                  = blast.Bspline(16, 110, 5, robot.n_joints)
opt.constraints.position     = True
opt.constraints.velocity     = True
opt.constraints.acceleration = True
opt.success_tolerance        = 0.01
opt.guess.type               = blast.GuessType.random
opt.guess.n_random_shots     = 50

result = blast.optimize(opt)

print("Success:", result.success)
print("Duration:", result.compute_time, "ms")
```

> **Array layout** — `Trajectory.pos`, `.vel`, and `.acc` are Fortran-contiguous arrays of shape `(n_joints, n_points)`. Index as `traj.pos[joint_idx, time_step]`. Use `np.ascontiguousarray(traj.pos)` if you need a C-contiguous copy.
>
> **Dtype** — always construct input arrays with `blast.array(values)` (or `np.asarray(values, dtype=blast.REAL_DTYPE)`) to avoid silent float32/float64 mismatches.

---

## Key concepts

### Manipulator

A `Manipulator` holds the robot's kinematic and dynamic model: DH parameters, joint limits, link masses, inertia tensors, and capsule collision geometry.

```cpp
Manipulator robot;
robot.n_joints        = 6;
robot.limits          = my_limits;
robot.kinematics      = my_kinematics;
robot.dynamics        = my_dynamics;
robot.collision_model = my_capsules;
```

See `examples/ur5e.hpp` for a complete UR5e definition.

### Task

A `Task` specifies the boundary conditions of the motion — start and goal joint angles, velocities, and accelerations.

```cpp
// Most common: start and stop at rest
Task task = Task::stop_to_stop(q_start, q_goal);

// Full control over boundary conditions
Task task;
task.start.position     = q_start;
task.start.velocity     = {0, 0, 0, 0, 0, 0};
task.start.acceleration = {0, 0, 0, 0, 0, 0};
task.goal.position      = q_goal;
task.goal.velocity      = {0, 0, 0, 0, 0, 0};
task.goal.acceleration  = {0, 0, 0, 0, 0, 0};
```

### Bspline

The trajectory is represented as a B-spline. The `Bspline` struct configures the resolution and smoothness of that representation.

```cpp
// Bspline(n_control_points, n_eval_points, degree, n_joints)
opt.bspline = Bspline(16, 110, 5, robot.n_joints);
//                    ^^  ^^^  ^
//                    |   |    degree 5 → C4-continuous
//                    |   evaluation density (constraint samples)
//                    optimizer resolution (decision variables per joint)
```

Higher `n_control_points` gives the optimizer more freedom but increases solve time. Higher `n_eval_points` makes constraint evaluation denser (important for collision avoidance).

### Optimization

`Optimization` bundles robot, task, B-spline, constraints, solver limits, and the initial guess strategy.

```cpp
Optimization opt(robot, task);

// Constraints
opt.constraints.position            = true;
opt.constraints.velocity            = true;
opt.constraints.acceleration        = true;
opt.constraints.self_collisions     = true;
opt.constraints.external_collisions = true;

// Solver settings
opt.success_tolerance = 0.01;   // relative constraint violation threshold
opt.max_eval          = 5000;   // max function evaluations

// Initial guess
opt.guess.type           = Guess::random;
opt.guess.n_random_shots = 50;  // number of random starting points
```

### World (collision)

A `World` holds static obstacles. Attach it to an `Optimization` to enable collision avoidance.

```cpp
World world;
world.add_box(
    Vec3{0.4, 0.0, 0.6},             // centre (m)
    Vec3{0.05, 0.3, 0.3},            // half-extents (m)
    Mat3{1,0,0, 0,1,0, 0,0,1}        // rotation
);
// also: world.add_sphere(...), world.add_capsule(...)

opt.world = world;
opt.constraints.external_collisions = true;
```

### Result

```cpp
Result result = optimize(&opt);

result.success               // bool — did the optimizer satisfy all constraints?
result.compute_time          // double — wall-clock solve time (ms)
result.num_eval              // int — number of objective function evaluations
result.max_constraint_value  // double — largest constraint violation (0 = feasible)
result.x                     // Array — optimized B-spline control points;
                             //         result.x.back() is the trajectory duration (s)
```

### Optimization methods

```cpp
// Default: evaluates one constraint per segment (fast, good for most problems)
Result r = optimize(&opt, OptimizationMethod::with_segments);

// Evaluates constraints at all eval points (slower, denser for collision)
Result r = optimize(&opt, OptimizationMethod::standard);
```

---

## Examples

All examples are in `examples/` and use the UR5e defined in `examples/ur5e.hpp`.

| Example | What it shows |
|---|---|
| `example_01_forward_kinematics` | Build a robot model; read joint frame positions after FK |
| `example_02_trajectory_optimization` | Point-to-point motion with kinematic constraints |
| `example_03_collision_avoidance` | Add a box obstacle; retry loop for harder problems |

Build and run a specific example:

```sh
cmake --preset dev
cmake --build --preset dev --target example_02_trajectory_optimization
./build/dev/examples/example_02_trajectory_optimization
```

---

## Citing Blast

If you use Blast in academic work, please cite:

```bibtex
@inproceedings{blast2026,
  title     = {Accelerating Trajectory Optimization by Exploiting B-Spline Gradient Structure},
  author    = {Doiron, Nikos and Duquette, Thomas and {Tchane Djogdom}, {Gilde Vanel} and Gallant, Andr\'{e}},
  booktitle = {IEEE International Conference on Robotics and Automation (ICRA)},
  year      = {2026}
}
```

---

## Contributing

Contributions are welcome. Please follow the coding style enforced by `.clang-format` — the CI will check formatting on every pull request.
