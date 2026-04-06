"""Tests for Vec3, Mat3 construction and numpy round-trips."""
import gc
import numpy as np
import pytest
import blast


# ---------------------------------------------------------------------------
# REAL_DTYPE
# ---------------------------------------------------------------------------

def test_real_dtype_is_float():
    assert blast.REAL_DTYPE in (np.float64, np.float32)


def test_blast_array_helper():
    a = blast.array([1.0, 2.0, 3.0])
    assert a.dtype == blast.REAL_DTYPE
    assert a.shape == (3,)


# ---------------------------------------------------------------------------
# Vec3
# ---------------------------------------------------------------------------

class TestVec3:
    def test_construct_xyz(self):
        v = blast.Vec3(1.0, 2.0, 3.0)
        assert v.x == pytest.approx(1.0)
        assert v.y == pytest.approx(2.0)
        assert v.z == pytest.approx(3.0)

    def test_construct_default(self):
        v = blast.Vec3()
        assert v.x == 0.0 and v.y == 0.0 and v.z == 0.0

    def test_construct_from_numpy(self):
        arr = blast.array([4.0, 5.0, 6.0])
        v = blast.Vec3(arr)
        assert v.x == pytest.approx(4.0)
        assert v.y == pytest.approx(5.0)
        assert v.z == pytest.approx(6.0)

    def test_to_numpy_shape_dtype(self):
        v = blast.Vec3(1.0, 2.0, 3.0)
        a = v.to_numpy()
        assert a.shape == (3,)
        assert a.dtype == blast.REAL_DTYPE

    def test_round_trip(self):
        v = blast.Vec3(1.5, -2.5, 3.5)
        a = v.to_numpy()
        v2 = blast.Vec3.from_numpy(a)
        assert v2.x == pytest.approx(v.x)
        assert v2.y == pytest.approx(v.y)
        assert v2.z == pytest.approx(v.z)

    def test_getitem_setitem(self):
        v = blast.Vec3(1.0, 2.0, 3.0)
        assert v[0] == pytest.approx(1.0)
        assert v[2] == pytest.approx(3.0)
        v[1] = 99.0
        assert v.y == pytest.approx(99.0)

    def test_len(self):
        assert len(blast.Vec3()) == 3

    def test_implicit_conversion_via_constructor(self):
        # Passing numpy array where Vec3 expected (using from_numpy explicitly)
        arr = blast.array([7.0, 8.0, 9.0])
        v = blast.Vec3.from_numpy(arr)
        assert v.x == pytest.approx(7.0)

    def test_wrong_length_raises(self):
        with pytest.raises(Exception):
            blast.Vec3.from_numpy(blast.array([1.0, 2.0]))

    def test_to_numpy_is_copy(self):
        v = blast.Vec3(1.0, 2.0, 3.0)
        a = v.to_numpy()
        a[0] = 999.0
        assert v.x == pytest.approx(1.0)  # original unmodified


# ---------------------------------------------------------------------------
# Mat3
# ---------------------------------------------------------------------------

class TestMat3:
    def test_construct_9_values(self):
        # Column-by-column: Mat3(col0r0, col0r1, col0r2, col1r0, ...)
        # Col0 = (1,2,3), Col1 = (4,5,6), Col2 = (7,8,9)
        m = blast.Mat3(1, 2, 3, 4, 5, 6, 7, 8, 9)
        # m(row, col):
        assert m(0, 0) == pytest.approx(1.0)  # col0, row0
        assert m(1, 0) == pytest.approx(2.0)  # col0, row1
        assert m(2, 0) == pytest.approx(3.0)  # col0, row2
        assert m(0, 1) == pytest.approx(4.0)  # col1, row0

    def test_to_numpy_shape_dtype(self):
        m = blast.Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)  # identity
        a = m.to_numpy()
        assert a.shape == (3, 3)
        assert a.dtype == blast.REAL_DTYPE

    def test_identity_round_trip(self):
        eye_np = np.eye(3, dtype=blast.REAL_DTYPE)
        m = blast.Mat3.from_numpy(eye_np)
        a = m.to_numpy()
        np.testing.assert_allclose(a, eye_np, atol=1e-12)

    def test_column_major_layout(self):
        """Verify that the numpy output matches the logical matrix entries."""
        # Mat3(col0r0=1, col0r1=2, col0r2=3,  col1r0=4, ...)
        m = blast.Mat3(1, 2, 3, 4, 5, 6, 7, 8, 9)
        a = m.to_numpy()
        # a is C-order (3,3): a[r,c] should equal m(r,c)
        assert a[0, 0] == pytest.approx(m(0, 0))
        assert a[2, 1] == pytest.approx(m(2, 1))

    def test_from_numpy_row_major(self):
        a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=blast.REAL_DTYPE)
        m = blast.Mat3.from_numpy(a)
        for r in range(3):
            for c in range(3):
                assert m(r, c) == pytest.approx(a[r, c])

    def test_wrong_shape_raises(self):
        with pytest.raises(Exception):
            blast.Mat3.from_numpy(np.eye(4, dtype=blast.REAL_DTYPE))


# ---------------------------------------------------------------------------
# ManipulatorLimits — Array round-trips
# ---------------------------------------------------------------------------

class TestManipulatorLimits:
    def test_position_max_round_trip(self):
        lim = blast.ManipulatorLimits()
        orig = blast.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        lim.position_max = orig
        back = lim.position_max
        assert back.dtype == blast.REAL_DTYPE
        np.testing.assert_allclose(back[:6], orig)

    def test_setter_accepts_float32_via_conversion(self):
        """nanobind converts float32 → float64 implicitly; verify the data is correct."""
        if not blast.REAL_IS_DOUBLE:
            pytest.skip("Only relevant for double builds")
        lim = blast.ManipulatorLimits()
        lim.position_max = np.ones(6, dtype=np.float32)
        np.testing.assert_allclose(lim.position_max[:6], 1.0)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_max_joints():
    assert blast.MAX_JOINTS == 8

def test_max_capsules():
    assert blast.MAX_CAPSULES == 7
