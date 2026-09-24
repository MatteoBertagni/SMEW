"""Public Python and compiled solver selection."""

import inspect

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


def test_compiled_and_python_initialization_can_run_successively():
    pytest.importorskip("smew._native._minpack")
    fractions = np.array((0.6, 0.2, 0.08, 0.05, 0.04, 0.03))
    concentrations, _ = smew.f_CEC_to_conc(
        fractions, 6.0, "loam", 1.0, 1.0,
    )
    python, _ = smew.conc_to_f_CEC(
        concentrations, 6.0, "loam", 1.0, 1.0, backend="python",
    )
    compiled, _ = smew.conc_to_f_CEC(
        concentrations, 6.0, "loam", 1.0, 1.0, backend="compiled",
    )
    np.testing.assert_allclose(compiled, python, rtol=1e-9, atol=1e-12)


def test_compiled_initialization_uses_both_vector_systems():
    pytest.importorskip("smew._native._minpack")
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
    compiled = smew.total_to_f_CEC_and_conc(*args, backend="compiled")
    for left, right in zip(python, compiled):
        np.testing.assert_allclose(right, left, rtol=1e-9, atol=1e-12)

    totals[0] += 0.01
    args = (totals, 6.0, concentrations, [0.6], "loam",
            0.4, 0.3, capacity, 1.0, 1.0)
    python = smew.Kelland(*args, backend="python")
    compiled = smew.Kelland(*args, backend="compiled")
    for left, right in zip(python, compiled):
        np.testing.assert_allclose(right, left, rtol=1e-9, atol=1e-12)


def test_compiled_simulation_runs_with_python_process_loop(monkeypatch):
    pytest.importorskip("smew._native._minpack")
    from tests.example_marimo_notebook import app

    original = smew.biogeochem_balance

    def compiled_balance(**inputs):
        return original(**inputs, backend="compiled")

    monkeypatch.setattr(smew, "biogeochem_balance", compiled_balance)
    _, definitions = app.run(defs={"duration_days": 2, "timestep_minutes": 60})
    assert definitions["results"]["pH"].size == 48
    assert np.isfinite(definitions["results"]["pH"]).all()
