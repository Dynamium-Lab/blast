# Quick start (C++)

A complete point-to-point optimization. The UR5e model is defined in
`examples/ur5e.hpp`; see the [Concepts guide](../guides/concepts.md) for how to build
your own `Manipulator`.

```cpp
#define BLAST_USE_NATIVE_SQP  // use the bundled SQP solver
#include <blast>
#include <iostream>

int main() {
    using namespace blast;

    // 1. Build the robot model (UR5e shown; define your own Manipulator similarly).
    Manipulator robot = make_ur5e();

    // 2. Define the motion task: stop-to-stop from start to goal joint angles.
    Task task = Task::stop_to_stop(
        {0.0, -1.57, 1.57, -1.57, -1.57, 0.0},  // start (rad)
        {1.57, -1.0,  1.0, -1.57, -1.57, 0.5}    // goal  (rad)
    );

    // 3. Set up the optimization problem.
    Optimization opt(robot, task);
    opt.bspline                  = Bspline(16, 110, 5, robot.n_joints);
    opt.constraints.position     = true;
    opt.constraints.velocity     = true;
    opt.constraints.acceleration = true;
    opt.success_tolerance        = 0.01;
    opt.guess.type               = Guess::random;
    opt.guess.n_random_shots     = 50;

    // 4. Solve.
    Result result = optimize(&opt);

    std::cout << (result.success ? "Success" : "Failed") << "\n";
    std::cout << "Duration:   " << result.x.back() << " s\n";
    std::cout << "Solve time: " << result.compute_time << " ms\n";
    return 0;
}
```

## Building and running the examples (Linux or similar build tools)

```sh
cmake --preset dev        # configure (RelWithDebInfo + -march=native)
cmake --build --preset dev --target example_02_trajectory_optimization
./build/dev/examples/example_02_trajectory_optimization
```

## Building and running the examples (MSVC)

```sh
cmake --preset dev-msvc        # configure (RelWithDebInfo + -march=native)
cmake --build --preset dev-msvc --target example_02_trajectory_optimization
./build/dev/examples/example_02_trajectory_optimization
```

| Example                              | What it shows                                            |
|--------------------------------------|----------------------------------------------------------|
| `example_01_forward_kinematics`      | Build a robot model; read joint frame positions after FK |
| `example_02_trajectory_optimization` | Point-to-point motion with kinematic constraints         |
| `example_03_collision_avoidance`     | Add a box obstacle; retry loop for harder problems       |

Next: the [Concepts guide](../guides/concepts.md) explains `Manipulator`, `Task`,
`Bspline`, `Optimization`, `World`, and `Result` in detail.
