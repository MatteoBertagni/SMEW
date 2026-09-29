# SPDX-License-Identifier: AGPL-3.0-only
"""Lazy Numba registration and native solver binding; no model equations here."""

import ctypes

from numba import njit
from numba.extending import overload, register_jitable

from smew import constants, vegetation, weathering
from smew._native import _minpack
from smew.biogeochem import _biogeochem_balance
from smew.ic import _conc_to_f_CEC, _total_to_f_CEC_and_conc, _Kelland
from smew._utils import _solve_system

# Registration keeps these functions ordinary Python functions in Python mode.
for helper in (
    constants.CO2_atm, constants.D_0, constants.Dw_0, constants.MM,
    constants.K_Al, constants.K_C, constants.K_GT_CEC, constants.plant_nutr_f,
    constants.carb_weath_const, constants.min_const,
    vegetation.up_act, weathering.carb_W, weathering.sil_Omega,
    weathering.sil_Wr, weathering.psd_evol, weathering.psd_number_from_mass,
    weathering.normalized_cumulative_area, weathering.wet_f_Anand,
    weathering.wetness_SA,
):
    register_jitable(helper)

# The import above connects the solver to the compiled equations. Here we access
# its C function so Numba can call it directly, keeping the library handle available.
_library = ctypes.CDLL(_minpack.__file__)
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



@overload(_solve_system)
def _native_solver_overload(system, state, parameters, residual, work, xtol):
    def solve(system, state, parameters, residual, work, xtol):
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
    return solve


# Compile the same functions used by the Python backend. The ctypes binding
# is process-local, so disk caching is deliberately disabled.
compiled_balance = njit(nogil=True, error_model="numpy", cache=False)(
    _biogeochem_balance,
)
compiled_conc_to_f_CEC = njit(nogil=True, error_model="numpy", cache=False)(
    _conc_to_f_CEC,
)
compiled_total_to_f_CEC_and_conc = njit(nogil=True, error_model="numpy", cache=False)(
    _total_to_f_CEC_and_conc,
)
compiled_Kelland = njit(nogil=True, error_model="numpy", cache=False)(
    _Kelland,
)
