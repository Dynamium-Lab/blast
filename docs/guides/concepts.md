# Concepts

Blast is organized around a few core types. This page gives the mental model; exact
fields and methods live in the [C++](../api/cpp/index.md) and
[Python](../api/python/index.md) API references.

## Manipulator

A `Manipulator` holds the robot's kinematic and dynamic model: joint geometry, joint
limits, link masses, inertia tensors, and capsule collision geometry.

```cpp
Manipulator robot;
robot.n_joints        = 6;
robot.limits          = my_limits;
robot.kinematics      = my_kinematics;
robot.dynamics        = my_dynamics;
robot.collision_model = my_capsules;
```

See `examples/ur5e.hpp` for a complete UR5e definition.

## Task

A `Task` specifies the boundary conditions of the motion — start and goal joint angles,
velocities, and accelerations.

```cpp
// Most common: start and stop at rest.
Task task = Task::stop_to_stop(q_start, q_goal);

// Full control over boundary conditions.
Task task;
task.start.position     = q_start;
task.start.velocity     = {0, 0, 0, 0, 0, 0};
task.start.acceleration = {0, 0, 0, 0, 0, 0};
task.goal.position      = q_goal;
```

## Bspline

The trajectory is represented as a B-spline. The `Bspline` struct configures the
resolution and smoothness of that representation.

```cpp
// Bspline(n_control_points, n_eval_points, degree, n_joints)
opt.bspline = Bspline(16, 110, 5, robot.n_joints);
//                    ^^  ^^^  ^
//                    |   |    degree 5 -> C4-continuous
//                    |   evaluation density (constraint samples)
//                    optimizer resolution (decision variables per joint)
```

Higher `n_control_points` gives the optimizer more freedom but increases solve time.
Higher `n_eval_points` makes constraint evaluation denser (important for collision
avoidance).

## Optimization

`Optimization` bundles robot, task, B-spline, constraints, solver limits, and the initial
guess strategy. See [Constraints](constraints.md) and
[Objectives & guesses](objectives-and-guesses.md) for the knobs.

## World (collision)

A `World` holds static obstacles. Attach it to an `Optimization` to enable collision
avoidance.

```cpp
World world;
world.add_box(
    Vec3{0.4, 0.0, 0.6},      // centre (m)
    Vec3{0.05, 0.3, 0.3},     // half-extents (m)
    Mat3{1,0,0, 0,1,0, 0,0,1} // rotation
);
// also: world.add_sphere(...), world.add_capsule(...)

opt.world = world;
opt.constraints.external_collisions = true;
```

## Result

```cpp
Result result = optimize(&opt);

result.success               // did the optimizer satisfy all constraints?
result.compute_time          // wall-clock solve time (ms)
result.num_eval              // number of objective function evaluations
result.max_constraint_value  // largest constraint violation (0 = feasible)
result.x                     // optimized control points; x.back() is duration (s)
```
