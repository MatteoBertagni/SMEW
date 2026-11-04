# SPDX-License-Identifier: AGPL-3.0-only
"""Small helpers shared by SMEW modules."""


import logging
import warnings

import numpy as np
from scipy.optimize import fsolve
from numba import config, njit
from numba.extending import overload
from smew.equations import (
    water_alk_equations, biogeochem_equations, cec_calcium_equation,
    f_cec_alk_to_conc_equation, total_to_cec_equations,
)


# Integer tags passed to the native dispatcher to select the corresponding C residual
# function pointer inside the compiled cminpack wrapper.
WATER_ALK_SYSTEM = 1
BIOGEOCHEM_SYSTEM = 2
CEC_CALCIUM_SYSTEM = 3
F_CEC_ALK_TO_CONC_SYSTEM = 4
TOTAL_TO_CEC_SYSTEM = 5


@njit
def _solver_work_size(n_unknowns):
    """Number of float64 values required by MINPACK hybrd's work buffer."""
    return n_unknowns * (3 * n_unknowns + 13) // 2


# Execution contract: operates as the pure-Python fallback using SciPy's fsolve.
# When running under the compiled backend, this implementation is replaced by Numba's @overload.
def _solve_system(system, state, parameters, residual, work, xtol):
    """Return MINPACK status codes in both Python and Numba.

    Compiled callers must return failures (status <= 0) to a Python entry point
    rather than raise while they still own Numba-allocated arrays to prevent memory leaks.
    """
    if system == WATER_ALK_SYSTEM:
        equations = water_alk_equations
    elif system == BIOGEOCHEM_SYSTEM:
        equations = biogeochem_equations
    elif system == CEC_CALCIUM_SYSTEM:
        equations = cec_calcium_equation
    elif system == F_CEC_ALK_TO_CONC_SYSTEM:
        equations = f_cec_alk_to_conc_equation
    elif system == TOTAL_TO_CEC_SYSTEM:
        equations = total_to_cec_equations
    else:
        raise ValueError("Unknown equation system")
    solution, _, status, _ = fsolve(
        equations, state, args=tuple(parameters), xtol=xtol, full_output=True,
    )
    state[:] = solution
    residual[:] = equations(state, *parameters)
    if not np.isfinite(residual).all():
        return -1
    return status


def _warn_solver_status(status):
    """Report native nonconvergence outside the compiled calculation."""
    if status != 1:
        warnings.warn(
            f"MINPACK stopped with status {status}", RuntimeWarning, stacklevel=3,
        )


@overload(_solve_system)
def _compiled_solve_system(system, state, parameters, residual, work, xtol):
    return _solve_with_cminpack


# Load the native solver after defining the equation identifiers it imports.
_solve_with_cminpack = None
if not config.DISABLE_JIT:
    try:
        from smew._native._numba import _solve_with_cminpack
    except (ImportError, OSError, AttributeError, ValueError) as error:
        # Missing extensions, loader/symbol failures, and Cython ABI mismatches.
        logging.getLogger(__name__).error(
            "SMEW native solver could not be loaded (%s: %s). "
            "Falling back to Python/SciPy fsolve; chemistry calculations may be slower. "
            "Rebuild or reinstall the native extensions and restart Python to restore compiled solves.",
            type(error).__name__, error,
        )


def _solver_njit(function):
    """Compile solver callers only when the native solver can be loaded.

    Other SMEW helpers keep their normal Numba setting; do not change Numba's
    process-wide configuration when only SMEW's extensions are unavailable.
    """
    if _solve_with_cminpack is None:
        return function
    return njit(nogil=True, error_model="numpy")(function)
