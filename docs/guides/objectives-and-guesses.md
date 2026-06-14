# Objectives & initial guesses

## Objective

The `Objective` controls what the optimizer minimizes. By default Blast minimizes
trajectory duration (`time_weight = 1.0`). Custom objective terms can be added:

```cpp
opt.objective.time_weight = 1.0;
opt.objective.add_custom_objective(my_objective_fn, /*weight=*/0.5);
```

A custom objective is a function `double(Optimization*)` that returns the value to be
added (scaled by its weight) to the total cost.

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
```

| `GuessType` | Behaviour |
|---|---|
| `custom` | Use the user-provided `initial_x` vector |
| `random` | Try `n_random_shots` random guesses, keep the best |
| `shotgun` | Fire a fixed number of random shots |
| `from_list` | Try each column of the `candidates` matrix |

```{tip}
More random shots improves robustness on hard problems at the cost of solve time. Start
with `random` + 50 shots and increase only if optimization fails to converge.
```

See the [API reference](../api/cpp/index.md) for the full `Objective` and `Guess`
definitions.
