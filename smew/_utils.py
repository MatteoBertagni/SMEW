# SPDX-License-Identifier: AGPL-3.0-only
"""Small helpers shared by SMEW modules."""


import warnings

from scipy.optimize import fsolve
from numba.extending import overload
from smew.equations import (
    water_equations, h_equations, biogeochem_equations,
    cec_calcium_equation, total_to_cec_equations, kelland_equations,
)


# Integer tags passed to the native dispatcher to select the corresponding C residual
# function pointer inside the compiled cminpack wrapper.
WATER_SYSTEM = 0
HYDROGEN_SYSTEM = 1
BIOGEOCHEM_SYSTEM = 2
CEC_CALCIUM_SYSTEM = 3
TOTAL_TO_CEC_SYSTEM = 4
KELLAND_SYSTEM = 5


# Execution contract: operates as the pure-Python fallback using SciPy's fsolve.
# When running under the compiled backend, this implementation is replaced by Numba's @overload.
def _solve_system(system, state, parameters, residual, work, xtol):
    """SciPy solver adapter; Numba replaces this call with the native binding."""
    if system == WATER_SYSTEM:
        equations = water_equations
    elif system == HYDROGEN_SYSTEM:
        equations = h_equations
    elif system == BIOGEOCHEM_SYSTEM:
        equations = biogeochem_equations
    elif system == CEC_CALCIUM_SYSTEM:
        equations = cec_calcium_equation
    elif system == TOTAL_TO_CEC_SYSTEM:
        equations = total_to_cec_equations
    elif system == KELLAND_SYSTEM:
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


@overload(_solve_system)
def _compiled_solve_system(system, state, parameters, residual, work, xtol):
    # Loaded during JIT compilation; disabled JIT keeps the SciPy implementation.
    from smew._native._numba import _solve_with_cminpack
    return _solve_with_cminpack
