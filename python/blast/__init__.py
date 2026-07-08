"""
blast — Python bindings for the Blast trajectory optimization library.

Quick start
-----------
>>> import numpy as np
>>> import blast
>>>
>>> limits = blast.ManipulatorLimits()
>>> limits.position_max = np.deg2rad([360, 360, 180, 360, 360, 360])
>>> limits.position_min = -limits.position_max
>>> limits.velocity_max = np.full(6, np.pi)
>>> limits.acceleration_max = np.full(6, 13.96)
>>> limits.torque_max = np.array([150., 150., 150., 28., 28., 28.])
>>> limits.tool_speed_max = 2.0
>>>
>>> # ... set up kinematics, dynamics, capsules ...
>>> # manip = blast.Manipulator(6, limits, kinematics, dynamics, capsules)
>>> # task  = blast.Task.stop_to_stop(start_q, goal_q)
>>> # opt   = blast.Optimization(manip, task)
>>> # result = blast.optimize(opt)
>>> # print(result.success, result.trajectory.pos.shape)

Note on array layouts
---------------------
``Trajectory.pos``, ``.vel``, ``.acc`` are returned as **Fortran-contiguous**
numpy arrays of shape ``(n_joints, n_points)``.  Index them as
``traj.pos[joint_idx, time_step]``.  Call ``np.ascontiguousarray(traj.pos)``
if you need a C-contiguous copy.

``REAL_DTYPE`` is ``np.float64`` (default build) or ``np.float32``.
Always construct input arrays with ``blast.array(values)`` or
``np.asarray(values, dtype=blast.REAL_DTYPE)`` to avoid silent dtype mismatches.
"""

from ._blast_core import (
    # Math
    Vec3,
    Mat3,
    MAX_JOINTS,
    MAX_CAPSULES,
    REAL_DTYPE,
    REAL_IS_DOUBLE,
    # World
    Box,
    Sphere,
    Capsule,
    DynamicBox,
    DynamicSphere,
    DynamicCapsule,
    World,
    # Manipulator
    CollisionModelCapsule,
    ManipulatorLimits,
    ManipulatorKinematics,
    ManipulatorDynamics,
    ManipulatorCapsules,
    Tool,
    Manipulator,
    # Task
    TaskBoundary,
    Task,
    # Trajectory
    Trajectory,
    Bspline,
    # Optimization
    GuessType,
    OptimizationMethod,
    ConstraintSelection,
    Objective,
    Guess,
    Optimization,
    Result,
    optimize,
)

import numpy as _np


def array(values) -> _np.ndarray:
    """Return a numpy array with the correct dtype for Blast (``REAL_DTYPE``).

    Use this helper whenever you need to pass numeric data to Blast functions
    to avoid silent float32/float64 mismatches.

    Parameters
    ----------
    values : array-like
        Any sequence or array that numpy can convert.

    Returns
    -------
    np.ndarray
        1-D or N-D array with dtype == blast.REAL_DTYPE.
    """
    return _np.asarray(values, dtype=REAL_DTYPE)


__all__ = [
    # Math
    "Vec3", "Mat3", "MAX_JOINTS", "MAX_CAPSULES", "REAL_DTYPE", "REAL_IS_DOUBLE",
    # World
    "Box", "Sphere", "Capsule", "DynamicBox", "DynamicSphere", "DynamicCapsule", "World",
    # Manipulator
    "CollisionModelCapsule", "ManipulatorLimits", "ManipulatorKinematics",
    "ManipulatorDynamics", "ManipulatorCapsules", "Tool", "Manipulator",
    # Task
    "TaskBoundary", "Task",
    # Trajectory
    "Trajectory", "Bspline",
    # Optimization
    "GuessType", "OptimizationMethod", "ConstraintSelection", "Objective",
    "Guess", "Optimization", "Result", "optimize",
    # Helpers
    "array",
]
