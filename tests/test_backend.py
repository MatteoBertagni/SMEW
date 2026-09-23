"""Public backend selection before the native backend is available."""

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
@pytest.mark.parametrize(
    ("backend", "error", "message"),
    (("compiled", RuntimeError, "not been implemented"),
     ("invalid", ValueError, "Unknown backend")),
)
def test_unavailable_or_unknown_backend_fails_before_model_work(
    entry_point, backend, error, message,
):
    required = {
        name: None
        for name, parameter in inspect.signature(entry_point).parameters.items()
        if parameter.default is inspect.Parameter.empty
    }
    with pytest.raises(error, match=message):
        entry_point(**required, backend=backend)


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
