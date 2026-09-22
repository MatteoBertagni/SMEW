"""Check that regressions and invalid values cannot silently pass the harness."""

import numpy as np
import pytest

from tests.test_physical_ranges import assert_physical_range
from tests.test_regression import assert_series_matches


def test_detects_interior_change_with_unchanged_endpoints():
    reference = np.array([1.0, 2.0, 3.0])
    current = np.array([1.0, 2.1, 3.0])
    with pytest.raises(pytest.fail.Exception, match="First difference at day 1"):
        assert_series_matches(
            current, reference, np.arange(3), label="example/pH", rtol=0, atol=1e-7,
        )


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_matching_nonfinite_values_are_rejected(invalid):
    values = np.array([invalid])
    with pytest.raises(AssertionError, match="NaN or infinity"):
        assert_series_matches(values, values, np.array([0]), label="pH", rtol=0, atol=1e-7)
    with pytest.raises(AssertionError, match="NaN or infinity"):
        assert_physical_range(values, {"minimum": 0}, label="s")


def test_out_of_range_value_fails_even_if_it_matches_a_reference():
    with pytest.raises(AssertionError, match="exceeds"):
        assert_physical_range(np.array([1.01]), {"minimum": 0, "maximum": 1}, label="s")
