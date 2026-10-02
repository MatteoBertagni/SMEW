"""Error codes, messages, transport, and composition with compiled callers."""

import pickle

import numpy as np
import pytest
from numba import njit

import smew
from smew.errors import _solver_error


def test_error_codes_are_unique():
    codes = list(smew.ErrorCode)
    assert len({code.value for code in codes}) == len(codes)


@pytest.mark.parametrize("code,base", [
    (smew.ErrorCode.CHEMISTRY_NO_CONVERGENCE, ValueError),
    (smew.ErrorCode.CHEMISTRY_RESIDUAL, RuntimeError),
])
def test_exception_keeps_code_context_and_serializes(code, base):
    message = ("biogeochem.chemistry: failure; step=37, time_days=1.5, "
               "solver_status=4, max_abs_residual=0.25, previous_pH=5.2, trial_H=0.0001")
    with pytest.raises(base) as caught:
        smew.raise_for_error(code, message)
    error = caught.value
    assert isinstance(error, smew.SMEWError)
    assert error.code == code
    assert error.message == message
    text = str(error)
    assert code.name in text
    assert message in text
    restored = pickle.loads(pickle.dumps(error))
    assert type(restored) is type(error)
    assert restored.code == code
    assert restored.message == message
    assert str(restored) == text


def test_success_does_not_raise():
    smew.raise_for_error(0, "")


@njit
def _parent_respiration(values):
    result, code, message = smew.respiration_numba(
        None, 1000., None, 1., "loam", values*.6, values,
        1., .3, -values, 1., 1., 1000.,
    )
    # A parent allocating arrays must propagate the failure without raising.
    if code != smew.ErrorCode.OK:
        return None, code, message
    return result, code, message


def test_compiled_parent_preserves_child_failure():
    result, code, message = _parent_respiration(np.ones(8))
    assert result is None
    assert code == smew.ErrorCode.RESPIRATION_MEAN_ACTIVITY
    assert "step=" not in message  # Aggregate failure, not a timestep failure.
    assert "mean_activity=0" in message
    with pytest.raises(smew.SMEWError, match="mean_activity=0"):
        smew.raise_for_error(code, message)


def test_respiration_reports_initial_state():
    values = np.ones(8)
    with pytest.raises(smew.SMEWError) as caught:
        smew.respiration(None, 1000., .01, 1., "loam", values*.6, values,
                         1., .3, -values, 1., 1.)
    assert caught.value.code == smew.ErrorCode.RESPIRATION_INITIAL_ACTIVITY
    for fragment in ("step=0", "activity=0", "moisture=0.6", "temperature_C=-1", "SOC=1000"):
        assert fragment in caught.value.message


@njit
def _invalid_solver_error():
    # Unwritten scratch memory must not appear in diagnostics.
    return _solver_error(np.int32(-2), np.empty(1), smew.ErrorCode.RAIN_RESIDUAL,
                         smew.ErrorCode.RAIN_SOLVER_INPUT, "biogeochem.rainwater")


def test_invalid_solver_buffers_have_no_residual_measurement():
    result, code, message = _invalid_solver_error()
    assert result is None
    assert code == smew.ErrorCode.RAIN_SOLVER_INPUT
    assert "solver_status=-2" in message
    assert "max_abs_residual=" not in message

