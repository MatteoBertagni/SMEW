"""Numba and Python execution selected once, at process startup."""

import os
import subprocess
import sys

import numpy as np
import pytest

import smew

pytestmark = pytest.mark.xdist_group("numba_modes")


@pytest.fixture(scope="module")
def mode_outputs(tmp_path_factory):
    directory = tmp_path_factory.mktemp("numba-modes")
    outputs = {}
    for mode, disabled in (("python", "1"), ("compiled", "0")):
        path = directory / f"{mode}.npz"
        process = subprocess.run(
            [sys.executable, "-m", "tests.numba_probe", str(path)],
            env={**os.environ, "NUMBA_DISABLE_JIT": disabled},
            capture_output=True, text=True,
        )
        assert process.returncode == 0, process.stdout + process.stderr
        with np.load(path) as archive:
            outputs[mode] = {name: archive[name] for name in archive.files}
    assert outputs["python"].keys() == outputs["compiled"].keys()
    return outputs


@pytest.mark.parametrize("group", (
    "soil", "pores", "mineral", "omega", "_gradient", "CO2_atm", "D_0", "Dw_0",
    "MM", "K_Al", "K_C", "K_GT_CEC", "plant_nutr_f", "carb_weath_const",
    "carb_W", "sil_Wr", "temp", "veg", "season", "boundaries", "up_act",
    "moisture", "respiration", "mov_avg", "psd_evol", "psd_number_from_mass",
    "normalized_cumulative_area", "wet_f_Anand", "wetness", "f_CEC_to_conc",
    "f_CEC_and_conc_to_K", "Amann", "conc_to_f_CEC", "total_to_f_CEC_and_conc",
    "Kelland", "water_equations", "h_equations", "biogeochem_equations",
    "cec_calcium_equation", "total_to_cec_equations", "kelland_equations",
    "water_residual", "h_residual", "biogeochem_residual", "cec_calcium_residual",
    "total_to_cec_residual", "kelland_residual", "simulation", "two_rocks",
))
def test_numerical_outputs_match_between_modes(mode_outputs, group):
    names = [name for name in mode_outputs["python"] if name.split("/")[0] == group]
    assert names, group
    # Seasonal biomass uses float32; solvers accumulate roundoff differently.
    rtol = 1e-6 if group == "season" else 1e-9
    for name in names:
        python = mode_outputs["python"][name]
        compiled = mode_outputs["compiled"][name]
        # min_const uses NaN for unavailable solubility products (its last field).
        if group == "mineral" and name.endswith("/10") and np.isnan(python).all():
            np.testing.assert_array_equal(compiled, python, err_msg=name)
            continue
        assert np.isfinite(python).all(), name
        assert np.isfinite(compiled).all(), name
        np.testing.assert_allclose(compiled, python, rtol=rtol, atol=1e-9, err_msg=name)


def test_vegetation_growth_and_harvest_timing(mode_outputs):
    for outputs in mode_outputs.values():
        growth = outputs["veg"]
        np.testing.assert_array_equal(growth[:4], 0.)
        assert growth[4] == 1.
        assert growth[5] == pytest.approx(1.135)
        v, harvest = outputs["season/2/0"], outputs["season/2/1"]
        np.testing.assert_array_equal(v[:2], 0.)
        assert v[2] == 1.
        np.testing.assert_array_equal(np.flatnonzero(harvest), [17])
        assert harvest[17] == v[16]
        assert v[17] == 0.
        assert not outputs["season/3/1"].any()  # Season continues past the final sample.


def test_pore_gradient_matches_numpy(mode_outputs):
    expected = np.gradient([1., 2., 8., 10.], [0., .2, 1., 2.])
    for outputs in mode_outputs.values():
        np.testing.assert_allclose(outputs["_gradient"], expected, rtol=1e-14)
        for name, values in outputs.items():
            if name.startswith("pores/") and name.endswith("/1"):
                d = outputs[name[:-1] + "0"]
                assert (values >= 0.).all()
                assert np.trapezoid(values, d) == pytest.approx(1.)


def test_compiled_helpers_need_no_native_solver(tmp_path):
    process = subprocess.run(
        [sys.executable, "-m", "tests.numba_probe", str(tmp_path / "helpers.npz"), "--helpers-only"],
        env={**os.environ, "NUMBA_DISABLE_JIT": "0"},
        capture_output=True, text=True,
    )
    assert process.returncode == 0, process.stdout + process.stderr


def test_cumulative_area_preserves_trapezoidal_integral():
    from scipy.integrate import cumulative_trapezoid

    x = np.array((4., 1., 2., 2., 3.))
    y = np.array((2., -1., 1., 7., 3.))
    order = np.argsort(x)
    grid, indices = np.unique(x[order], return_index=True)
    cumulative = cumulative_trapezoid(np.maximum(y[order][indices], 0.), grid, initial=0.)
    actual_grid, actual = smew.normalized_cumulative_area(x, y)
    np.testing.assert_array_equal(actual_grid, grid)
    np.testing.assert_allclose(actual, cumulative / cumulative[-1], rtol=1e-15)
