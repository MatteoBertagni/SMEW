"""Physical envelopes are independent of agreement with the saved reference."""

import numpy as np
import pytest

from tests.support import CHECKS


def assert_physical_range(values, rule, *, label):
    assert values.size > 0, f"{label}: empty output"
    assert np.isfinite(values).all(), f"{label}: NaN or infinity"
    tolerance = rule.get("bound_atol", 0.0)
    if "minimum" in rule:
        assert np.all(values >= rule["minimum"] - tolerance), (
            f"{label}: minimum {values.min():.12g} is below "
            f"{rule['minimum']} (allowed roundoff {tolerance:g})"
        )
    if "maximum" in rule:
        assert np.all(values <= rule["maximum"] + tolerance), (
            f"{label}: maximum {values.max():.12g} exceeds "
            f"{rule['maximum']} (allowed roundoff {tolerance:g})"
        )


@pytest.mark.parametrize("variable", tuple(CHECKS["variables"]))
def test_physical_ranges(scenario_name, results, variable):
    rule = CHECKS["variables"][variable]
    assert_physical_range(
        results[variable], rule,
        label=f"{scenario_name}/{variable} ({rule['units']})",
    )
