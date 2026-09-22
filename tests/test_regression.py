"""Every selected output must match its accepted reference at every timestep."""

import numpy as np
import pytest

from tests.support import CHECKS


def assert_series_matches(actual, expected, time_days, *, label, rtol, atol):
    assert actual.shape == expected.shape, (
        f"{label}: shape changed: {expected.shape} -> {actual.shape}"
    )
    assert actual.ndim >= 1 and actual.shape[-1] == len(time_days), (
        f"{label}: the last axis must correspond to time"
    )
    assert np.isfinite(actual).all(), f"{label}: current output contains NaN or infinity"
    assert np.isfinite(expected).all(), f"{label}: reference contains NaN or infinity"
    difference = np.abs(actual - expected)
    outside = difference > atol + rtol * np.abs(expected)
    if outside.any():
        failing_indices = np.argwhere(outside)
        first = tuple(failing_indices[np.argmin(failing_indices[:, -1])])
        worst = np.unravel_index(np.argmax(difference), difference.shape)
        pytest.fail(
            f"{label}: {outside.sum()}/{actual.size} values exceed "
            f"atol={atol:g}, rtol={rtol:g}. "
            f"First difference at day {time_days[first[-1]]:.8g}, index {first}: "
            f"reference={expected[first]:.12g}, current={actual[first]:.12g}. "
            f"Maximum absolute difference={difference[worst]:.12g} "
            f"at day {time_days[worst[-1]]:.8g}, index {worst}."
        )


def test_time_grid(results, baseline):
    time_days = results["time_days"]
    assert time_days.ndim == 1 and time_days.size > 1
    assert np.isfinite(time_days).all()
    assert np.all(np.diff(time_days) > 0)
    # The time axis must match exactly; no interpolation or resampling.
    np.testing.assert_array_equal(time_days, baseline["time_days"])


def test_reference_variables(baseline):
    assert set(baseline) == {"time_days", *CHECKS["variables"]}, (
        "Reference variable list differs from checks.toml; review before updating it."
    )


@pytest.mark.parametrize("variable", tuple(CHECKS["variables"]))
def test_time_series(scenario_name, results, baseline, variable):
    rule = CHECKS["variables"][variable]
    assert_series_matches(
        results[variable], baseline[variable], results["time_days"],
        label=f"{scenario_name}/{variable} ({rule['units']})",
        rtol=rule.get("rtol", CHECKS["comparison"]["rtol"]), atol=rule["atol"],
    )


def test_notebook_runs_without_plotting(notebook_run):
    assert notebook_run["plot_button"].value is False
    assert "figure" not in notebook_run
    assert "plt" not in notebook_run


def test_full_year_duration(notebook_run):
    assert notebook_run["duration_days"] == 365
    assert notebook_run["timestep_minutes"] == 10
    assert notebook_run["results"]["time_days"].size == 365 * 24 * 6
