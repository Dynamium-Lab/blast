"""Shared test fixtures, mirroring examples/ur5e.hpp in Python."""
import numpy as np
import pytest
import blast


@pytest.fixture(scope="session")
def ur5e() -> blast.Manipulator:
    """Fully configured UR5e manipulator, matching examples/ur5e.hpp."""
    # --- Limits ---
    limits = blast.ManipulatorLimits()
    limits.position_max     = blast.array([6.2832, 6.2832, 3.1416, 6.2832, 6.2832, 6.2832])
    limits.position_min     = blast.array([-6.2832, -6.2832, -3.1416, -6.2832, -6.2832, -6.2832])
    limits.velocity_max     = blast.array([3.1416, 3.1416, 3.1416, 3.1416, 3.1416, 3.1416])
    limits.acceleration_max = blast.array([13.96, 13.96, 13.96, 13.96, 13.96, 13.96])
    limits.torque_max       = blast.array([150.0, 150.0, 150.0, 28.0, 28.0, 28.0])
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
    kin.joint_axes = [blast.Vec3(0, 0, 1)] * 6
    kin.first_joint_position = blast.Vec3(0.0, 0.0, 0.1625)

    # Mat3 is column-by-column: Mat3(col0r0, col0r1, col0r2, col1r0, ...)
    rots = [blast.Mat3()] * blast.MAX_JOINTS  # placeholders; we only set [0..5]
    # Use the same values as ur5e.hpp (column-major: {a,b,c,d,e,f,g,h,i} sets columns)
    rots[0] = blast.Mat3(-1, 0, 0,  0, -1, 0,  0, 0, 1)
    rots[1] = blast.Mat3( 1, 0, 0,  0,  0, 1,  0,-1, 0)
    rots[2] = blast.Mat3( 1, 0, 0,  0,  1, 0,  0, 0, 1)
    rots[3] = blast.Mat3( 1, 0, 0,  0,  1, 0,  0, 0, 1)
    rots[4] = blast.Mat3( 1, 0, 0,  0,  0, 1,  0,-1, 0)
    rots[5] = blast.Mat3( 1, 0, 0,  0,  0,-1,  0, 1, 0)
    # Pad to MAX_JOINTS+1 (static_rotations has MAX_JOINTS+1 entries)
    kin.static_rotations = rots + [blast.Mat3()]

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

    base_sphere = blast.Sphere()
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
        make_cap(0, (0, 0, 0),                   (0, -0.15, 0),               0.09),
        make_cap(1, (-0.42, 0, 0.1375),           (0, 0, 0.1375),              0.06),
        make_cap(2, (0, 0, 0.02),                 (0, 0, 0.18),                0.065),
        make_cap(2, (-0.373156, 0, 0.00850418),   (0.000440611, 0.000121443, 0.00850418), 0.05),
        make_cap(3, (0, 0, 0),                    (0, 0, -0.155),              0.0425),
        make_cap(3, (0, 0.03, 0),                 (0, -0.1, 0),               0.0425),
        make_cap(5, (0, 0, -0.14),                (0, 0, -0.01),               0.038),
    ]

    return blast.Manipulator(6, limits, kin, dyn, caps_cfg)


@pytest.fixture(scope="session")
def ur5e_task() -> blast.Task:
    """Standard stop-to-stop task matching examples/ur5e.hpp."""
    start = blast.array([1.94822,  0.473555, -0.0255247, -0.448375,  0.370356, -3.12883])
    goal  = blast.array([2.5825,   0.0700,   -0.3892,     0.3196,    0.9927,   -3.17328])
    return blast.Task.stop_to_stop(start, goal)
