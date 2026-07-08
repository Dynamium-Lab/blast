"""End-to-end optimization tests using the UR5e fixture."""
import numpy as np
import pytest
import blast


class TestManipulatorConstruction:
    def test_n_joints(self, ur5e):
        assert ur5e.n_joints == 6

    def test_has_no_tool_by_default(self, ur5e):
        assert not ur5e.has_tool


class TestTask:
    def test_stop_to_stop_shape(self, ur5e_task):
        assert ur5e_task.start.position.shape == (6,)
        assert ur5e_task.goal.position.shape  == (6,)

    def test_stop_to_stop_values(self, ur5e_task):
        np.testing.assert_allclose(
            ur5e_task.start.position,
            blast.array([1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883]),
        )


class TestOptimizationConstruction:
    def test_construct(self, ur5e, ur5e_task):
        opt = blast.Optimization(ur5e, ur5e_task)
        assert opt.x_len() > 0

    def test_bspline_default(self, ur5e, ur5e_task):
        opt = blast.Optimization(ur5e, ur5e_task)
        assert opt.bspline.n_joints == 6

    def test_explicit_bspline(self, ur5e, ur5e_task):
        bs = blast.Bspline(16, 110, 5, 6)
        opt = blast.Optimization(ur5e, ur5e_task, bs)
        assert opt.bspline.n_ctrl == 16
        assert opt.bspline.n_points == 110


class TestOptimizationSmokeTest:
    """Run a quick optimization and verify result shape/types.

    We use very permissive settings (large tolerance, few evaluations) so the
    test completes quickly without depending on the optimizer finding a
    feasible solution.
    """

    @pytest.fixture(scope="class")
    def quick_result(self, ur5e, ur5e_task):
        opt = blast.Optimization(ur5e, ur5e_task)
        opt.bspline                 = blast.Bspline(12, 80, 5, 6)
        opt.constraints.position    = True
        opt.constraints.velocity    = True
        opt.constraints.acceleration = True
        opt.success_tolerance       = 0.05
        opt.max_tries               = 1
        opt.max_eval                = 200
        opt.max_time                = 10.0
        opt.guess                   = blast.Guess(10)  # 10 random shots
        return blast.optimize(opt)

    def test_result_type(self, quick_result):
        assert isinstance(quick_result.success, bool)

    def test_compute_time_positive(self, quick_result):
        assert quick_result.compute_time > 0.0

    def test_trajectory_shape(self, quick_result):
        traj = quick_result.trajectory
        pos  = traj.pos
        assert pos.ndim  == 2
        assert pos.shape[0] == 6   # n_joints
        assert pos.shape[1] > 0    # n_points

    def test_trajectory_dtype(self, quick_result):
        assert quick_result.trajectory.pos.dtype == blast.REAL_DTYPE

    def test_trajectory_fortran_contiguous(self, quick_result):
        pos = quick_result.trajectory.pos
        assert np.isfortran(pos), (
            "Trajectory.pos should be Fortran-contiguous (column-major, "
            "matching blast::Matrix storage)"
        )

    def test_time_array_shape(self, quick_result):
        t = quick_result.trajectory.t
        assert t.ndim == 1
        assert t.shape[0] == quick_result.trajectory.pos.shape[1]

    def test_x_shape(self, quick_result):
        x = quick_result.x
        assert x.ndim  == 1
        assert x.shape[0] > 0

    def test_x_dtype(self, quick_result):
        assert quick_result.x.dtype == blast.REAL_DTYPE


class TestWorldObstacles:
    def test_add_box(self, ur5e, ur5e_task):
        opt = blast.Optimization(ur5e, ur5e_task)
        opt.world.add_box(
            blast.Vec3(0.5, 0.0, 0.5),
            blast.Vec3(0.1, 0.1, 0.1),
            blast.Mat3(1, 0, 0,  0, 1, 0,  0, 0, 1),
        )
        opt.constraints.external_collisions = True
        # Just verify it doesn't crash at construction — don't run the full optimizer
        assert opt.x_len() > 0

    def test_add_sphere(self, ur5e, ur5e_task):
        opt = blast.Optimization(ur5e, ur5e_task)
        opt.world.add_sphere(blast.Vec3(0.5, 0.0, 0.5), 0.05)
        assert len(opt.world.spheres) == 1


class TestResultViews:
    """Verify that numpy views returned by Result stay valid."""

    @pytest.fixture(scope="class")
    def quick_result(self, ur5e, ur5e_task):
        opt = blast.Optimization(ur5e, ur5e_task)
        opt.bspline           = blast.Bspline(12, 80, 5, 6)
        opt.max_eval          = 200
        opt.max_time          = 10.0
        opt.guess             = blast.Guess(10)
        opt.success_tolerance = 0.05
        return blast.optimize(opt)

    def test_trajectory_copy_is_independent(self, quick_result):
        import gc
        traj = quick_result.trajectory
        pos  = traj.pos
        del traj
        gc.collect()
        # pos is a copy — it must remain readable after traj is collected
        assert pos.ndim == 2
        _ = pos[0, 0]

    def test_x_view_survives_gc(self, quick_result):
        import gc
        x = quick_result.x
        gc.collect()
        _ = x[0]
