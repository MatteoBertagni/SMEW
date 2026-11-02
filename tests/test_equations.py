"""Known aqueous balances, checked with and without Numba compilation."""

import numpy as np
import pytest

from smew.equations import water_alk_equations, water_alk_residual


@pytest.mark.parametrize("alkalinity,expected", (
    (0., -0.00101901),
    (0.002, 0.00098099),
    (-0.001, -0.00201901),
    (0.00101901, 0.),
))
def test_water_alk_residual_sign_and_parameter_order(alkalinity, expected):
    # At H=1e-6: HCO3=1e-3, CO3=1e-5, OH=1e-8, hence Alk=0.00101901.
    parameters = (alkalinity, 1e-6, 1e-8, 1e-3, 1e-14)
    assert water_alk_residual(1e-6, *parameters) == pytest.approx(
        expected, rel=1e-13, abs=1e-18,
    )
    assert water_alk_equations(np.array([1e-6]), *parameters) == pytest.approx(
        expected, rel=1e-13, abs=1e-18,
    )
