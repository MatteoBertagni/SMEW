# SPDX-License-Identifier: AGPL-3.0-only
"""Small helpers shared by SMEW modules."""


import warnings

from scipy.optimize import fsolve
from smew.equations import (
    water_equations, h_equations, biogeochem_equations,
    cec_calcium_equation, total_to_cec_equations, kelland_equations,
)


def require_backend(backend):
    """Validate execution choice and check that the native solver is installed."""
    if backend not in ("python", "compiled"):
        raise ValueError(
            f"Unknown backend {backend!r}; choose 'python' or 'compiled'."
        )
    if backend == "compiled":
        try:
            from smew._native import _minpack
        except ImportError as exc:
            raise RuntimeError(
                "The compiled backend is unavailable. Install a native wheel, "
                "build the extensions with 'make install-native', or select "
                "backend='python'."
            ) from exc
        return _minpack
    return None


def _solve_system(system, state, parameters, residual, work, xtol):
    """SciPy solver adapter; Numba replaces this call with the native binding."""
    if system == 0:
        equations = water_equations
    elif system == 1:
        equations = h_equations
    elif system == 2:
        equations = biogeochem_equations
    elif system == 3:
        equations = cec_calcium_equation
    elif system == 4:
        equations = total_to_cec_equations
    elif system == 5:
        equations = kelland_equations
    else:
        raise ValueError("Unknown equation system")
    state[:] = fsolve(equations, state, args=tuple(parameters), xtol=xtol)
    residual[:] = equations(state, *parameters)
    return 1  # fsolve already reports nonconvergence through its own warnings.


def _warn_solver_status(status):
    """Report native nonconvergence outside the compiled calculation."""
    if status != 1:
        warnings.warn(
            f"MINPACK stopped with status {status}", RuntimeWarning, stacklevel=3,
        )
