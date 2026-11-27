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

How densely constraints are evaluated along the trajectory, and how collisions are searched,
is set by `opt.method` before calling `optimize`:

```cpp
// Default: one constraint per B-spline segment (fast, good for most problems).
opt.method = OptimizationMethod::with_segments;

// Point-based with finite-difference gradients (reference baseline).
opt.method = OptimizationMethod::baseline;

// Point-based with analytical gradients for position/velocity/acceleration.
opt.method = OptimizationMethod::with_analytical_pva;

// As above, plus analytical torque-dynamics gradients.
opt.method = OptimizationMethod::with_analytical_dynamics;

// Segment-based, with a bounding volume hierarchy over the obstacles (many obstacles).
opt.method = OptimizationMethod::broadphase;

// As broadphase, plus a hierarchy over each capsule's samples within a segment.
opt.method = OptimizationMethod::double_broadphase;

Result r = optimize(&opt);
```

The two broadphase methods read a hierarchy over the static obstacles that the caller
builds once, after setting the world (dynamic obstacles are handled internally):

```cpp
opt.world = world;
create_static_bounding_volume_hierarchy(opt.world, opt.world.static_bounding_volume_hierarchy);
```

They give the same constraints and gradients as `with_segments`, and follow the same
tightening and success rules below.

## Custom constraints

Besides the built-in constraints, any method accepts user constraints of two kinds. Each
returns an `Array` of values, each satisfied when `<= 0`, and is registered with the number of
values it returns.

**Per point**: evaluated at every B-spline sample.

```cpp
// Keep the tool above z = 0.2 m.
Array tool_above(Optimization* opt, Array* x, ManipulatorTempData* state, u32 point) {
  return Array(1, (real) 0.2 - state->p_j[opt->manip.n_joints].z);
}
opt.constraints.add_point_constraint(tool_above, /*n_values=*/1);
```

`state` holds the sample's forward kinematics (`rotations`, `rotations_mult`, `p_j`), collision
model (`capsule_list`) and joint efforts (`efforts`). The sample's joint positions, velocities,
accelerations and time are `opt->bspline.traj.pos/vel/acc.col(point)` and
`opt->bspline.traj.t[point]`. A point constraint must not recompute the trajectory.

**Global**: evaluated once per trajectory.

```cpp
// Duration at most max_duration, passed through custom_data.
Array duration_cap(Optimization* opt, Array* x) {
  return Array(1, x->back() - *static_cast<real*>(opt->custom_data));
}
opt.constraints.add_global_constraint(duration_cap, /*n_values=*/1);
```

`x` is the decision vector: the free control points, joint by joint, then the duration
`x->back()`. A global constraint may call `opt->bspline.compute_trajectory(*x, opt->task)`.
(`add_constraint` is the former name of `add_global_constraint`.)

Where they go in the constraint vector:

| Method                                                                  | Point constraints                                       | Global constraints |
|-------------------------------------------------------------------------|---------------------------------------------------------|--------------------|
| `baseline`, `with_analytical_pva`, `with_analytical_dynamics`           | at the end of every sample's block                      | after the last sample |
| `with_segments`, `broadphase`, `double_broadphase`                      | at the end of every segment's block, worst sample of the segment, value by value | after the last segment |

Their gradients are forward finite differences in every method: for a point constraint, only the
sample is perturbed (the control points that affect it, and the duration), so it costs a few
evaluations of the function and of the kinematics per sample, or per segment for the
segment-based methods; for a global constraint, one evaluation per decision variable.

Custom constraints follow the tightening and success rules below: during the solve every
custom value is shifted up by `success_tolerance` (`custom_constraint_shift`), so an accepted
solve has every custom value `<= 0`, and the output trajectory is checked with the point
constraints evaluated at every validation sample. Write the values on a scale where
`success_tolerance` is a sensible margin (for example divide by the limit, as the built-in
constraints do).

## Feasibility: `success_tolerance` and `collision_buffer`

`success_tolerance` sets the relative constraint-violation threshold below which a
solution is considered feasible:

```cpp
opt.success_tolerance = 0.01;  // 1% relative violation
opt.max_eval          = 5000;  // cap on function evaluations
```

It is the only feasibility bar: the solver's own constraint tolerance, the margin the
limits are tightened by during the solve, and the post-solve accept gate. The limits are
tightened by exactly this amount, and the output trajectory is validated on the true limits,
so **a solve reported successful satisfies the original, untightened limits** -- the
tolerance buys the solver an early exit, not a violation you have to check for afterwards.

`success_tolerance` is dimensionless (a fraction of a limit), so it cannot express a
collision distance. `collision_buffer` is its counterpart in metres:

```cpp
opt.collision_buffer = 0.001;  // default: 1 mm
```

Capsules and obstacles are each grown by `collision_buffer / 2` during the solve, so a
pairwise check carries exactly `collision_buffer`, and the collision constraint is rescaled so
`success_tolerance` lands on it. An accepted solve is therefore strictly clear of the
true geometry. `collision_buffer` is the slack the solver may consume, not clearance left over:
a 1 mm `collision_buffer` guarantees no penetration, not 1 mm of remaining gap. Larger values
perturb the geometry more during the solve and can make tightly constrained problems
harder to solve. It must be positive and should stay small for tight environments: the
default `1e-3` works well, while smaller values steepen the collision gradients
(`collision_scale = success_tolerance / collision_buffer`) and cost convergence. A value of
`0` or less is replaced by the default `1e-3`, with a printed notice.

Two exceptions keep the tightening from making a task impossible. The start and goal of a
stop-to-stop task are fixed, so constraints there cannot be improved by the solver:

- **Position** limits are tightened like the others. If the start or goal of a joint lies
  outside the tightened bound, that bound is moved out for the whole trajectory, to halfway
  between the endpoint and the true limit, so the fixed endpoint is strictly inside it. Only
  that side moves, and the constraint keeps its slope. The other side still reaches
  `success_tolerance` exactly at its true limit; on the widened side the solver may accept a
  position slightly past the true limit (by less than the tightening), and the check below then
  reports the solve as failed. An endpoint at, past, or within rounding of the true limit is not
  moved: there is no room left to widen into, and no trajectory from an endpoint at the limit can
  pass the check below.
- **Collisions** next to the start or goal (the first and last segment) aim for at most the
  clearance that endpoint itself has, instead of the full `collision_buffer`. A start that
  clears an obstacle by less than `collision_buffer` is still accepted; every other part of the
  trajectory keeps the full `collision_buffer`.

A solve is reported as a success only if the output trajectory (every `output_steps_ms`,
1 ms by default) also meets the true, untightened limits and has no contact with the true
geometry. The solver only sees its constraint points; this check also covers the trajectory
between them.

Velocity, acceleration, torque, and tool speed do not need this exception: a stop-to-stop task
starts and ends at rest, so the endpoint values for those constraints are zero, well within any
tightened limit.

`opt.tighten_for_tolerance = false` disables the tightening entirely. It exists to
isolate its effect when measuring; with it off, an accepted solve may violate the real
limits by up to `success_tolerance`.

## Decision bounds: duration and control points

The solver's variables are the free B-spline control points and the duration `T`. Each has
a box the solver cannot leave, enforced exactly at every step:

```cpp
opt.min_duration = 0.01;               // T in [min_duration, max_duration], seconds
opt.max_duration = 30.0;
opt.control_point_bound_margin = 0.0;  // widen control-point box by this fraction of each joint's range
opt.bound_control_points = true;       // false: duration bounds only, no control-point box
```

With the control points inside the position limits, the whole trajectory is inside them
too (a B-spline stays within the hull of its control points).

The bounds are on by default. Set `bound_control_points = false` to keep only the duration bounds.

**Tasks close to a joint limit**: because the convex-hull property is conservative, control
points capped at the limit produce a trajectory that stays strictly inside it. Raise
`control_point_bound_margin` (e.g. `0.05`) to give the free control points room to create
a smooth curve through a near-limit region.

**Joints without position limits**: the control-point box falls back to `[-2π, 2π]`, which
covers all non-screwing stop-to-stop tasks. Set `bound_control_points = false` to remove
the box entirely, for example for continuous-rotation joints.

An initial guess outside the box is clamped into it.

See the [API reference](../api/cpp/index.md) for the full `ConstraintSelection` and
`OptimizationMethod` definitions.
