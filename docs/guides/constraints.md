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

## Feasibility: `success_tolerance` and `collision_buffer`

`success_tolerance` sets the relative constraint-violation threshold below which a
solution is considered feasible:

```cpp
opt.success_tolerance = 0.01;  // 1% relative violation
opt.max_eval          = 5000;  // cap on function evaluations
```

It is the only feasibility bar: the solver's own constraint tolerance, the margin the
limits are tightened by during the solve, and the post-solve accept gate. Because the
limits are tightened by exactly this amount, **a solve reported successful satisfies the
original, untightened limits** -- the tolerance buys the solver an early exit, not a
violation you have to check for afterwards.

`success_tolerance` is dimensionless (a fraction of a limit), so it cannot express a
collision distance. `collision_buffer` is its counterpart in metres:

```cpp
opt.collision_buffer = 0.001;  // default: 1 mm
```

Capsules and obstacles are each grown by `collision_buffer / 2` during the solve, so a
pairwise check carries exactly this margin, and the collision constraint is rescaled so
`success_tolerance` lands on it. An accepted solve is therefore strictly clear of the
true geometry. The buffer is the slack the solver may consume, not clearance left over:
a 1 mm buffer guarantees no penetration, not 1 mm of remaining gap. Larger values
perturb the geometry more during the solve and can make tightly constrained problems
harder to solve. Set it to `0` to fall back to `success_tolerance`.

`opt.tighten_for_tolerance = false` disables the tightening entirely. It exists to
isolate its effect when measuring; with it off, an accepted solve may violate the real
limits by up to `success_tolerance`.

See the [API reference](../api/cpp/index.md) for the full `ConstraintSelection` and
`OptimizationMethod` definitions.
