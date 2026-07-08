# Constraints

Constraints are toggled through the `ConstraintSelection` member of an `Optimization`
(`opt.constraints`). Only the enabled constraints are evaluated by the solver.

```cpp
Optimization opt(robot, task);

opt.constraints.position            = true;
opt.constraints.velocity            = true;
opt.constraints.acceleration        = true;
opt.constraints.torque              = true;
opt.constraints.tool_speed          = true;
opt.constraints.self_collisions     = true;
opt.constraints.external_collisions = true;  // requires opt.world to be set
```

| Constraint            | Meaning                                                 | Source of limits      |
|-----------------------|---------------------------------------------------------|-----------------------|
| `position`            | Joint angles stay within `[position_min, position_max]` | `ManipulatorLimits`   |
| `velocity`            | Joint velocities stay within `velocity_max`             | `ManipulatorLimits`   |
| `acceleration`        | Joint accelerations stay within `acceleration_max`      | `ManipulatorLimits`   |
| `torque`              | Joint torques stay within `torque_max`                  | `ManipulatorLimits`   |
| `tool_speed`          | The tool absolute speed stays within `tool_speed_max`   | `ManipulatorLimits`   |
| `self_collisions`     | Robot capsules do not collide with each other           | `ManipulatorCapsules` |
| `external_collisions` | Robot capsules do not collide with `World` obstacles    | `World`               |

## Constraint sampling: `OptimizationMethod`

How densely constraints are evaluated along the trajectory is controlled by the method
passed to `optimize`:

```cpp
// Default: one constraint per B-spline segment (fast, good for most problems).
Result r = optimize(&opt, OptimizationMethod::with_segments);

// Point-based with finite-difference gradients (reference baseline).
Result r = optimize(&opt, OptimizationMethod::baseline);

// Point-based with analytical gradients for position/velocity/acceleration.
Result r = optimize(&opt, OptimizationMethod::with_analytical_pva);

// As above, plus analytical torque-dynamics gradients.
Result r = optimize(&opt, OptimizationMethod::with_analytical_dynamics);
```

`success_tolerance` sets the relative constraint-violation threshold below which a
solution is considered feasible:

```cpp
opt.success_tolerance = 0.01;  // 1% relative violation
opt.max_eval          = 5000;  // cap on function evaluations
```

See the [API reference](../api/cpp/index.md) for the full `ConstraintSelection` and
`OptimizationMethod` definitions.
