"""Native solver interface checks independent of whole-model baselines."""

import ctypes

import numpy as np
import pytest

from smew.equations import water_equations


_minpack = pytest.importorskip("smew._native._minpack")


def test_native_water_solver_reports_status_and_final_residual():
    parameters = (0.0, 1e-6, 1e-8, 1e-3, 1e-14)
    initial = np.array([1e-6])
    state, residual, status, evaluations = _minpack.solve_water(initial, parameters)
    assert status == 1
    assert evaluations > 0
    np.testing.assert_array_equal(initial, [1e-6])
    np.testing.assert_allclose(residual, water_equations(state, *parameters), atol=1e-16)


def test_numba_solver_entry_uses_raw_buffers_and_checks_layout():
    solver = ctypes.CDLL(_minpack.__file__).smew_solve
    solver.restype = ctypes.c_int
    solver.argtypes = (
        ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_double, ctypes.c_int,
        ctypes.c_void_p, ctypes.c_int,
    )
    parameters = np.array((0.0, 1e-6, 1e-8, 1e-3, 1e-14))
    state = np.array((1e-6,))
    residual = np.empty(1)
    work = np.empty(8)
    args = (parameters.ctypes.data, 5, state.ctypes.data,
            residual.ctypes.data, 1.4901161193847656e-8, 400,
            work.ctypes.data, work.size)
    assert solver(0, 2, *args) == -2
    assert solver(0, 1, *args) == 1
    np.testing.assert_allclose(residual, water_equations(state, *parameters), atol=1e-16)


def test_workspace_is_reused_and_checked():
    workspace = _minpack.Workspace(1)
    parameters = (0.0, 1e-6, 1e-8, 1e-3, 1e-14)
    _, first, _, _ = _minpack.solve_water([1e-6], parameters, workspace=workspace)
    _, second, _, _ = _minpack.solve_water([2e-6], parameters, workspace=workspace)
    assert first is second
    with pytest.raises(ValueError, match="workspace for 5 unknowns"):
        _minpack.solve_kelland([1.0] * 5, [1.0] * 15, workspace=workspace)
    workspace.work = np.empty(1)
    with pytest.raises(ValueError, match="workspace buffers"):
        _minpack.solve_water([1e-6], parameters, workspace=workspace)


def test_native_residual_failure_is_distinct_from_solver_status():
    with pytest.raises(RuntimeError, match="residual evaluation failed during solve"):
        _minpack.solve_water([0.0], (0.0, 1e-6, 1e-8, 1e-3, 1e-14))


@pytest.mark.parametrize(
    ("initial", "parameters", "message"),
    (
        ([1.0, 2.0], [0.0] * 5, "needs 1 state"),
        ([1.0], [0.0] * 4, "needs 5 parameters"),
    ),
)
def test_native_solver_rejects_incompatible_layouts(
    initial, parameters, message,
):
    with pytest.raises(ValueError, match=message):
        _minpack.solve_water(initial, parameters)
