# Quick start (Python)

After `pip install .` (see [Installation](installation.md)), the bindings are importable
as `blast`. The example below uses a UR5e; all numeric values match `examples/ur5e.hpp`.

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

# --- Kinematics, dynamics, collision geometry ---
# See the full UR5e setup in the project README; abbreviated here.
kin  = blast.ManipulatorKinematics()
dyn  = blast.ManipulatorDynamics()
caps = blast.ManipulatorCapsules()
# ... populate kin / dyn / caps ...

# --- Build robot ---
robot = blast.Manipulator(6, limits, kin, dyn, caps)

# --- Task ---
start = blast.array([0.0, -1.57,  1.57, -1.57, -1.57,  0.0])
goal  = blast.array([1.57, -1.0,   1.0, -1.57, -1.57,  0.5])
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

print("Success: ", result.success)
print("Solve time:", result.compute_time, "ms")
```

```{note}
**Array layout** — `Trajectory.pos`, `.vel`, and `.acc` are Fortran-contiguous arrays of
shape `(n_joints, n_points)`. Index as `traj.pos[joint_idx, time_step]`. Use
`np.ascontiguousarray(traj.pos)` if you need a C-contiguous copy.

**Dtype** — always construct input arrays with `blast.array(values)` (or
`np.asarray(values, dtype=blast.REAL_DTYPE)`) to avoid silent float32/float64 mismatches.
```

The full Python API is documented in the [Python API reference](../api/python/index.md).
