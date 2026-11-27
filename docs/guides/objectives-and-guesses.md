# Objectives & initial guesses

## Objective

The `Objective` controls what the optimizer minimizes. By default Blast minimizes
trajectory duration (`time_weight = 1.0`). Custom objective terms can be added:

```cpp
opt.objective.time_weight = 1.0;
opt.objective.add_custom_objective(my_objective_fn, /*weight=*/0.5);
```

A custom objective is a function `real(Optimization*, Array* x)` that returns the value to be
added (scaled by its weight) to the total cost. Its gradient is taken by forward finite
differences, one evaluation per decision variable.

An objective that knows its gradient is registered with `add_custom_objective_with_gradient`. It
has the signature `real(Optimization*, Array* x, Array* grad)`: it returns the value and, when
`grad` is not null, writes `df/dx` into it (already zeroed, one entry per decision variable):

```cpp
opt.objective.add_custom_objective_with_gradient(my_objective_with_gradient, /*weight=*/1.0);
```

`examples/example_04_custom_constraints.cpp` gives the gradients of a sum of absolute velocities
and of absolute accelerations.

## Initial guess

Trajectory optimization is non-convex, so the starting point matters. The `Guess` member
selects the strategy:

```cpp
// Random restarts — try N random initial guesses, keep the best.
opt.guess.type           = Guess::random;
opt.guess.n_random_shots = 50;

// User-supplied starting vector.
opt.guess = Guess(initial_x);          // Guess::custom

// Try each column of a candidate matrix.
opt.guess = Guess(candidate_matrix);   // Guess::from_list

// Shotgun — a fixed number of random shots.
opt.guess = Guess(100u);               // Guess::shotgun

// Straight line start->goal; deterministic, no random draws.
opt.guess.type = Guess::straight_line;
```

| `GuessType` | Behaviour |
|---|---|
| `custom` | Use the user-provided `initial_x` vector |
| `random` | Try `n_random_shots` random guesses, keep the best |
| `shotgun` | Fire a fixed number of random shots |
| `from_list` | Try each column of the `candidates` matrix |
| `straight_line` | Straight line start->goal; time from the manipulator's vel/accel limits |

```{tip}
More random shots improves robustness on hard problems at the cost of solve time. Start
with `random` + 50 shots and increase only if optimization fails to converge.
```

### Starting from another planner's path

A joint-space waypoint path, for example from a sampling-based planner, can seed the
optimizer through `guess_from_path`, which returns a `Guess::custom` starting vector:

```cpp
Matrix path = ...;                             // n_joints x n_waypoints, one column per waypoint, start to goal
opt.guess = Guess(guess_from_path(&opt, path));
```

It fits the B-spline curve to the path by least squares, so that the curve moves along the path
at a steady rate: at 30 % of the duration it should be 30 % of the way along the path. Distance
along the path is the sum of the Euclidean joint-space distances between waypoints, so unevenly
spaced waypoints do not skew the fit. The start and goal control points stay pinned as for any
stop-to-stop task, so the curve starts and ends at rest. The fitted control points sit just
outside the path's corners, so the curve follows the path instead of cutting them.

The duration is the shortest for which the control points of the velocity and acceleration
curves are within the joint limits, a conservative choice for the whole curve. It requires a
stop-to-stop task; a path of zero length gives the straight-line guess, and so does a path
of the wrong shape (not one row per joint, or fewer than two waypoints), with a printed message
in release builds. A debug build stops on it instead.

See the [API reference](../api/cpp/index.md) for the full `Objective` and `Guess`
definitions.
