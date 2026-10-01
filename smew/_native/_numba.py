# SPDX-License-Identifier: AGPL-3.0-only
"""Cminpack function pointer for Numba-compiled solver calls."""

import ctypes
from smew._native import _minpack

_library = ctypes.CDLL(_minpack.__file__)

# Binds the compiled C/cminpack solver (`smew_solve`) so Numba can issue direct C-level
# function pointer calls using array memory addresses (`.ctypes.data`) without GIL locks.
_native_solve = _library.smew_solve
_native_solve.restype = ctypes.c_int  # Solver status code
_native_solve.argtypes = (
    ctypes.c_int,     # Equation system: water, hydrogen, biogeochemistry…
    ctypes.c_int,     # Number of unknowns
    ctypes.c_void_p,  # Address of the parameter array
    ctypes.c_int,     # Number of parameters
    ctypes.c_void_p,  # Address of the state array: guess in, solution out
    ctypes.c_void_p,  # Address of the residual array
    ctypes.c_double,  # Convergence tolerance: xtol
    ctypes.c_int,     # Maximum number of equation evaluations
    ctypes.c_void_p,  # Address of the solver's working array
    ctypes.c_int,     # Length of the working array
)


def _solve_with_cminpack(system, state, parameters, residual, work, xtol):
    """Solve through the native entry point, updating caller-owned arrays."""
    status = _native_solve(
        system, state.size, parameters.ctypes.data, parameters.size,
        state.ctypes.data, residual.ctypes.data, xtol,
        200 * (state.size + 1), work.ctypes.data, work.size,
    )
    if status == -1:
        raise RuntimeError("Native residual evaluation failed")
    if status <= 0:
        raise ValueError("Native solver received invalid buffers or settings")
    return status

