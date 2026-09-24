"""Check the private Cython equation module when a native wheel is installed."""

import numpy as np
import pytest

from smew import equations as python_equations


native = pytest.importorskip("smew._native._equations")


def test_private_extension_keeps_python_equations_importable():
    assert python_equations.__file__.endswith("equations.py")
    assert set(native.__pyx_capi__) == {
        "water_residual", "h_residual", "biogeochem_residual",
        "cec_calcium_residual", "total_to_cec_residual", "kelland_residual",
    }


@pytest.mark.parametrize(
    ("name", "arguments"),
    (
        ("water_equations", (np.array((1e-6,)), 0.0, 1e-6, 1e-8, 1e-3, 1e-14)),
        ("h_equations", (np.array((1e-6,)), 1e-6, 1e-8, 1e-3, 1e-14, 0.1)),
        ("cec_calcium_equation", (np.array((0.7,)), 1e-6, 1.0, 1e-3,
                                  1e-3, 1e-3, 1e-3, 1e-6, 1.0, 1.0, 1.0,
                                  1.0, 1.0)),
    ),
)
def test_native_equations_match_python(name, arguments):
    np.testing.assert_allclose(
        getattr(native, name)(*arguments),
        getattr(python_equations, name)(*arguments),
        rtol=1e-13,
        atol=1e-13,
    )
