"""Public Python and compiled solver selection."""

import inspect
import subprocess
import sys

import numpy as np
import pytest

import smew


ENTRY_POINTS = (
    smew.biogeochem_balance,
    smew.biogeochem_balance2psd,
    smew.conc_to_f_CEC,
    smew.total_to_f_CEC_and_conc,
    smew.Kelland,
)


@pytest.mark.parametrize(
    "entry_point",
    (smew.f_CEC_to_conc, smew.f_CEC_and_conc_to_K, smew.Amann),
)
def test_algebraic_initialization_has_no_backend_option(entry_point):
    assert "backend" not in inspect.signature(entry_point).parameters


@pytest.mark.parametrize("entry_point", ENTRY_POINTS)
def test_backend_is_keyword_only_and_defaults_to_python(entry_point):
    parameter = inspect.signature(entry_point).parameters["backend"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default == "python"


@pytest.mark.parametrize("entry_point", ENTRY_POINTS)
def test_unknown_backend_fails_before_model_work(entry_point):
    required = {
        name: None
        for name, parameter in inspect.signature(entry_point).parameters.items()
        if parameter.default is inspect.Parameter.empty
    }
    with pytest.raises(ValueError, match="Unknown backend"):
        entry_point(**required, backend="invalid")


def test_2psd_remains_python_only():
    required = {
        name: None
        for name, parameter in inspect.signature(
            smew.biogeochem_balance2psd
        ).parameters.items()
        if parameter.default is inspect.Parameter.empty
    }
    with pytest.raises(NotImplementedError, match="does not support"):
        smew.biogeochem_balance2psd(**required, backend="compiled")


def test_explicit_python_backend_runs_initialization():
    fractions = np.array((0.6, 0.2, 0.08, 0.05, 0.04, 0.03))
    concentrations, _ = smew.f_CEC_to_conc(
        fractions, 6.0, "loam", 1.0, 1.0,
    )
    default, _ = smew.conc_to_f_CEC(
        concentrations, 6.0, "loam", 1.0, 1.0,
    )
    explicit, _ = smew.conc_to_f_CEC(
        concentrations, 6.0, "loam", 1.0, 1.0, backend="python",
    )
    np.testing.assert_array_equal(explicit, default)


def test_python_simulation_needs_no_native_or_numba_imports():
    script = """
import importlib.abc
import sys

class BlockCompiled(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if (fullname.split('.')[0] in ('numba', 'llvmlite') or
                fullname.startswith(('smew._native', 'smew._simulation_compiled'))):
            raise AssertionError(f'Python mode imported {fullname}')

sys.meta_path.insert(0, BlockCompiled())
import smew
fractions = (0.6, 0.2, 0.08, 0.05, 0.04, 0.03)
concentrations, _ = smew.f_CEC_to_conc(fractions, 6.0, 'loam', 1.0, 1.0)
smew.conc_to_f_CEC(concentrations, 6.0, 'loam', 1.0, 1.0, backend='python')
capacity = 0.1
water_volume = 0.4 * 0.3 * 0.6 * 1000
totals = [concentrations[i] * water_volume + fractions[i] * capacity / charge
          for i, charge in enumerate((2, 2, 1, 1))]
smew.total_to_f_CEC_and_conc(
    totals, 6.0, 0.07, [0.6], 'loam', 0.4, 0.3, capacity, 1.0, 1.0,
    backend='python',
)
totals[0] += 0.01
smew.Kelland(
    totals, 6.0, concentrations, [0.6], 'loam', 0.4, 0.3, capacity, 1.0, 1.0,
    backend='python',
)
from tests.example_marimo_notebook import app
_, definitions = app.run(defs={'duration_days': 2, 'timestep_minutes': 60})
assert definitions['results']['pH'].size == 48
"""
    process = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True,
    )
    assert process.returncode == 0, process.stdout + process.stderr


def test_numba_compatible_cumulative_area_preserves_trapezoidal_integral():
    from scipy.integrate import cumulative_trapezoid
    from smew.weathering import normalized_cumulative_area

    x = np.array((4.0, 1.0, 2.0, 2.0, 3.0))
    y = np.array((2.0, -1.0, 1.0, 7.0, 3.0))
    order = np.argsort(x)
    grid, indices = np.unique(x[order], return_index=True)
    cumulative = cumulative_trapezoid(
        np.maximum(y[order][indices], 0.0), grid, initial=0.0,
    )
    actual_grid, actual = normalized_cumulative_area(x, y)
    np.testing.assert_array_equal(actual_grid, grid)
    np.testing.assert_allclose(actual, cumulative / cumulative[-1], rtol=1e-15)


def _reject_python_solver(*args, **kwargs):
    raise AssertionError("Compiled initialization called a Python solver")


@pytest.mark.parametrize("soil,conv_mol,pH_in", (
    ("sand", 1.0, 6), ("loam", 1.0, 6.0), ("clay", 1e6, 6.0),
))
def test_compiled_and_python_initialization_can_run_successively(
    monkeypatch, soil, conv_mol, pH_in,
):
    native = pytest.importorskip("smew._native._minpack")
    fractions = np.array((0.6, 0.2, 0.08, 0.05, 0.04, 0.03))
    concentrations, _ = smew.f_CEC_to_conc(
        fractions, pH_in, soil, conv_mol, 1.0,
    )
    python, _ = smew.conc_to_f_CEC(
        concentrations, pH_in, soil, conv_mol, 1.0, backend="python",
    )
    with monkeypatch.context() as patch:
        patch.setattr("smew._utils.fsolve", _reject_python_solver)
        patch.setattr(native, "solve_cec_calcium", _reject_python_solver)
        compiled, _ = smew.conc_to_f_CEC(
            concentrations, pH_in, soil, conv_mol, 1.0, backend="compiled",
        )
    np.testing.assert_allclose(compiled, python, rtol=1e-9, atol=1e-12)
    again, _ = smew.conc_to_f_CEC(
        concentrations, pH_in, soil, conv_mol, 1.0, backend="python",
    )
    np.testing.assert_array_equal(again, python)
    from smew._simulation_compiled import compiled_conc_to_f_CEC
    from smew.ic import _conc_to_f_CEC

    assert compiled_conc_to_f_CEC.py_func is _conc_to_f_CEC
    assert compiled_conc_to_f_CEC.nopython_signatures


def test_compiled_initialization_uses_both_vector_systems(monkeypatch):
    native = pytest.importorskip("smew._native._minpack")
    fractions = np.array((0.6, 0.2, 0.08, 0.05, 0.04, 0.03))
    concentrations, _ = smew.f_CEC_to_conc(
        fractions, 6.0, "loam", 1.0, 1.0,
    )
    ca, mg, potassium, na, _ = concentrations
    capacity = 0.1
    water_volume = 0.4 * 0.3 * 0.6 * 1000
    totals = [
        ca * water_volume + fractions[0] / 2 * capacity,
        mg * water_volume + fractions[1] / 2 * capacity,
        potassium * water_volume + fractions[2] * capacity,
        na * water_volume + fractions[3] * capacity,
    ]
    args = (totals, 6.0, fractions[4] + fractions[5], [0.6],
            "loam", 0.4, 0.3, capacity, 1.0, 1.0)
    python = smew.total_to_f_CEC_and_conc(*args, backend="python")
    with monkeypatch.context() as patch:
        patch.setattr("smew._utils.fsolve", _reject_python_solver)
        patch.setattr(native, "solve_total_to_cec", _reject_python_solver)
        compiled = smew.total_to_f_CEC_and_conc(*args, backend="compiled")
    for left, right in zip(python, compiled):
        np.testing.assert_allclose(right, left, rtol=1e-9, atol=1e-12)

    totals[0] += 0.01
    args = (totals, 6.0, concentrations, [0.6], "loam",
            0.4, 0.3, capacity, 1.0, 1.0)
    python = smew.Kelland(*args, backend="python")
    with monkeypatch.context() as patch:
        patch.setattr("smew._utils.fsolve", _reject_python_solver)
        patch.setattr(native, "solve_kelland", _reject_python_solver)
        compiled = smew.Kelland(*args, backend="compiled")
    for left, right in zip(python, compiled):
        np.testing.assert_allclose(right, left, rtol=1e-9, atol=1e-12)

    from smew import _simulation_compiled, ic

    for dispatcher, function in (
        (_simulation_compiled.compiled_total_to_f_CEC_and_conc, ic._total_to_f_CEC_and_conc),
        (_simulation_compiled.compiled_Kelland, ic._Kelland),
    ):
        assert dispatcher.py_func is function
        assert dispatcher.nopython_signatures


def test_compiled_simulation_uses_numba_timestep_loop():
    pytest.importorskip("smew._native._minpack")
    from tests.example_marimo_notebook import app

    _, python_definitions = app.run(defs={"duration_days": 2, "timestep_minutes": 60})

    _, definitions = app.run(defs={
        "duration_days": 2, "timestep_minutes": 60, "backend": "compiled",
    })
    assert definitions["results"]["pH"].size == 48
    assert np.isfinite(definitions["results"]["pH"]).all()
    from smew import _simulation_compiled

    from smew.biogeochem import _biogeochem_balance

    assert _simulation_compiled.compiled_balance.py_func is _biogeochem_balance
    assert _simulation_compiled.compiled_balance.nopython_signatures
    assert definitions["chemistry"].keys() == python_definitions["chemistry"].keys()
    for chemistry in (definitions["chemistry"], python_definitions["chemistry"]):
        assert not {"s", "pH_in", "conc_in", "mineral", "backend", "keyword_ssa"} & chemistry.keys()
    for name in ("pH", "Ca", "Mg", "Alk", "IC_tot", "M_rock", "wet_f"):
        np.testing.assert_allclose(
            definitions["chemistry"][name], python_definitions["chemistry"][name],
            rtol=1e-10, atol=1e-9,
        )


@pytest.mark.parametrize("variant", ("no_rock", "nonlinear_frozen"))
def test_compiled_timestep_branches_match_python(monkeypatch, variant):
    pytest.importorskip("smew._native._minpack")
    from tests.example_marimo_notebook import app

    settings = {"duration_days": 2, "timestep_minutes": 60}
    if variant == "no_rock":
        settings.update(rock_mass_g_m2=0.0, application_day=0.0)
    else:
        settings.update(
            dissolution_factor=1.0,
            mineral_mass_fractions=[1.0], minerals=["forsterite"],
            particle_diameters_um=[50.0, 100.0, 200.0],
            particle_mass_fractions=[0.2, 0.5, 0.3],
            pore_model="ding2016", pore_particle_mixing=1.0,
            rock_density_g_m3=3e6, wet_surface_model="nonlinear",
        )
    original = smew.biogeochem_balance

    def run(backend):
        def balance(**inputs):
            if variant == "nonlinear_frozen":
                inputs["temp_soil"] = inputs["temp_soil"].copy()
                inputs["temp_soil"][24:] = -2.0
            return original(**inputs)

        monkeypatch.setattr(smew, "biogeochem_balance", balance)
        _, definitions = app.run(defs={**settings, "backend": backend})
        return definitions["chemistry"]

    python = run("python")
    compiled = run("compiled")
    if variant == "nonlinear_frozen":
        assert compiled["frozen"].sum() == 24
    for name in ("pH", "Ca", "IC_tot", "SA", "wet_f"):
        np.testing.assert_allclose(
            compiled[name], python[name], rtol=1e-10, atol=1e-9,
        )
