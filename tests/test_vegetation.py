"""Seasonal vegetation caps Euler overshoot and reports numerical failures."""

import numpy as np
import pytest
from numba import config, njit

import smew


def test_zero_capacity_keeps_moisture_finite_after_emergence():
    stages = np.array([1] * 10 + [2] * 10 + [3] * 10 + [4, 0, 0], dtype=np.int8)
    biomass, harvest = smew.veg_seasonal(stages, 0., .05)
    np.testing.assert_array_equal(biomass, 0.)
    np.testing.assert_array_equal(harvest, 0.)
    outputs = smew.moisture_balance(
        np.zeros(stages.size), .28, "loam", np.full(stages.size, .002),
        biomass, 3000., 1, .588, stages.size / 144., 1 / 144.,
    )
    for index in (0, 4, 5, 6):  # Moisture, leakage, transpiration, evaporation.
        assert np.isfinite(outputs[index]).all()


@pytest.mark.parametrize("duration", [6, 12])
@pytest.mark.parametrize("fraction", [.05, .5, .95])
def test_seasonal_growth_matches_classic_vegetation(duration, fraction):
    stages = np.array([2] * duration + [3], dtype=np.int8)
    biomass, harvest = smew.veg_seasonal(stages, 3000., fraction)
    expected = smew.veg(3000. * fraction, float(duration), 3000., 0., np.ones(stages.size), 1.)
    assert biomass.dtype == np.float32
    np.testing.assert_allclose(biomass, expected, rtol=1e-6)
    assert not harvest.any()


@pytest.mark.parametrize("duration", [1, 2, 3, 5])
@pytest.mark.parametrize("fraction", [.05, .5, .95])
def test_short_growth_period_caps_overshoot_and_preserves_harvest(duration, fraction):
    stages = np.array([2] * duration + [3] * 20 + [4, 0, 0], dtype=np.int8)
    biomass, harvest = smew.veg_seasonal(stages, 3000., fraction)
    for values in (biomass, harvest):
        assert values.dtype == np.float32
        assert np.isfinite(values).all()
        assert (values >= 0).all() and (values <= 3000.).all()
    assert biomass[0] == np.float32(3000. * fraction)
    assert (np.diff(biomass[:-1]) >= 0).all()
    assert biomass[-2] == 3000.
    assert biomass[-1] == 0.
    np.testing.assert_array_equal(np.flatnonzero(harvest), [stages.size - 1])
    assert harvest[-1] == biomass[-2]


@pytest.mark.parametrize("capacity", [-1., np.nan, np.inf, -np.inf])
def test_invalid_capacity_is_rejected(capacity):
    with pytest.raises(ValueError, match="k_v"):
        smew.veg_seasonal(np.array([2, 3, 4, 0]), capacity, .05)


@pytest.mark.parametrize("fraction", [-.1, 0., 1., 1.1, np.nan, np.inf])
def test_invalid_initial_fraction_is_rejected(fraction):
    with pytest.raises(ValueError, match="f_v_in"):
        smew.veg_seasonal(np.array([2, 3, 4, 0]), 3000., fraction)


@pytest.mark.parametrize("capacity", [1e-20, .1, np.float32(.1)])
def test_capacity_scale_does_not_make_growth_nonfinite(capacity):
    stages = np.array([2] * 12 + [3] * 20 + [4, 0, 0], dtype=np.int8)
    biomass, harvest = smew.veg_seasonal(stages, capacity, .95)
    for values in (biomass, harvest):
        assert np.isfinite(values).all()
        assert (values >= 0).all() and (values <= np.float32(capacity)).all()


def test_nonfinite_growth_returns_error_instead_of_a_partial_trajectory():
    stages = np.array([2] + [3] * 20 + [4, 0, 0], dtype=np.int8)
    # A representable capacity can still overflow float32 during an Euler update.
    with np.errstate(over="ignore", invalid="ignore"):
        result, code, message = smew.veg_seasonal_numba(stages, 2e38, .5)
    assert result is None
    assert code == smew.ErrorCode.VEGETATION_GROWTH
    assert "nonfinite biomass" in message and "T_v=1" in message


@pytest.mark.parametrize("stages", [[], [0, 0], [3, 4, 0, 1, 4, 0]])
def test_no_growth_period_produces_zero_biomass(stages):
    biomass, harvest = smew.veg_seasonal(np.array(stages, dtype=np.int8), 3000., .05)
    np.testing.assert_array_equal(biomass, np.zeros(len(stages)))
    np.testing.assert_array_equal(harvest, np.zeros(len(stages)))


def test_partial_season_does_not_shift_next_seasons_harvest():
    stages = np.array([3, 4, 0, 1] + [2] * 12 + [3, 4, 0, 0], dtype=np.int8)
    biomass, harvest = smew.veg_seasonal(stages, 10., .1)
    reference, reference_harvest = smew.veg_seasonal(stages[3:], 10., .1)
    np.testing.assert_array_equal(biomass[:3], 0.)
    np.testing.assert_array_equal(biomass[3:], reference)
    np.testing.assert_array_equal(harvest[3:], reference_harvest)
    assert harvest[-1] > 0


def test_normal_growth_preserves_euler_curve_and_harvest_timing():
    stages = np.array([0, 1] + [2] * 12 + [3, 4, 0, 0], dtype=np.int8)
    biomass, harvest = smew.veg_seasonal(stages, 10., .1)
    expected = np.zeros(stages.size, dtype=np.float32)
    expected[2] = 1.
    for i in range(3, 17):
        expected[i] = expected[i-1] + (6 / (12 * 10.)) * expected[i-1] * (10. - expected[i-1])
    np.testing.assert_allclose(biomass, expected, rtol=1e-6)
    np.testing.assert_array_equal(np.flatnonzero(harvest), [17])
    assert harvest[17] == biomass[16]


def test_growth_period_at_end_of_array_has_no_harvest():
    biomass, harvest = smew.veg_seasonal(np.array([0, 1, 2], dtype=np.int8), 10., .1)
    np.testing.assert_array_equal(biomass, [0., 0., 1.])
    assert not harvest.any()


@njit(nogil=True)
def _parent_vegetation(stages, capacity, fraction):
    # The parent owns a Numba allocation when the child reports a failure.
    owned_stages = stages.copy()
    result, code, message = smew.veg_seasonal_numba(owned_stages, capacity, fraction)
    if code != 0:
        return None, code, message
    return result, code, message


@pytest.mark.parametrize("stages,capacity,fraction,expected_code,context", [
    ([2, 3, 4, 0], -1., .05, smew.ErrorCode.VEGETATION_CAPACITY, "k_v=-1"),
    ([2, 3, 4, 0], 10., 2., smew.ErrorCode.VEGETATION_INITIAL_FRACTION, "f_v_in=2"),
])
def test_compiled_parent_propagates_validation_error(stages, capacity, fraction, expected_code, context):
    stages = np.array(stages, dtype=np.int8)
    result, code, message = _parent_vegetation(stages, capacity, fraction)
    assert result is None
    assert code == expected_code
    assert context in message
    # Only the Python boundary raises, preserving the compiled diagnostic.
    with pytest.raises(smew.SMEWError) as caught:
        smew.veg_seasonal(stages, capacity, fraction)
    assert caught.value.code == code
    assert caught.value.message == message


@pytest.mark.parametrize("capacity", [0., 10.])
def test_compiled_parent_returns_same_biomass_as_python_wrapper(capacity):
    stages = np.array([0, 1] + [2] * 12 + [3, 4, 0, 0], dtype=np.int8)
    result, code, message = _parent_vegetation(stages, capacity, .1)
    assert code == 0 and message == ""
    for actual, expected in zip(result, smew.veg_seasonal(stages, capacity, .1)):
        np.testing.assert_array_equal(actual, expected)


@pytest.mark.skipif(config.DISABLE_JIT or not config.NRT_STATS, reason="Requires compiled NRT allocation counters")
@pytest.mark.parametrize("capacity", [-1., 0., 10.])
def test_compiled_parent_releases_arrays_on_success_and_failure(capacity):
    from tests.numba_probe import check_released

    stages = np.array([0, 1, 2, 3, 4, 0, 0], dtype=np.int8)
    check_released(lambda: _parent_vegetation(stages, capacity, .05))


@pytest.mark.skipif(config.DISABLE_JIT or not config.NRT_STATS, reason="Requires compiled NRT allocation counters")
def test_compiled_parent_releases_arrays_after_growth_failure():
    from tests.numba_probe import check_released

    stages = np.array([2] + [3] * 20 + [4, 0, 0], dtype=np.int8)
    result, code, message = _parent_vegetation(stages, 2e38, .5)
    assert result is None and code == smew.ErrorCode.VEGETATION_GROWTH
    check_released(lambda: _parent_vegetation(stages, 2e38, .5))
