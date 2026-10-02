# cython: language_level=3, boundscheck=False, wraparound=False
"""Cython callbacks connecting SMEW's compiled equations to bundled cminpack."""

import numpy as np
import warnings
from libc.float cimport DBL_EPSILON
from libc.math cimport isfinite
from smew._native._equations cimport (
    water_residual, h_residual, biogeochem_residual,
    cec_calcium_residual, total_to_cec_residual, kelland_residual,
)

from smew import _utils

# Native integers let callbacks use the shared identifiers without the GIL.
cdef int WATER_SYSTEM = _utils.WATER_SYSTEM
cdef int HYDROGEN_SYSTEM = _utils.HYDROGEN_SYSTEM
cdef int BIOGEOCHEM_SYSTEM = _utils.BIOGEOCHEM_SYSTEM
cdef int CEC_CALCIUM_SYSTEM = _utils.CEC_CALCIUM_SYSTEM
cdef int TOTAL_TO_CEC_SYSTEM = _utils.TOTAL_TO_CEC_SYSTEM
cdef int KELLAND_SYSTEM = _utils.KELLAND_SYSTEM

cdef extern from "cminpack.h":
    ctypedef int (*hybrd_callback)(void *, int, const double *, double *, int) noexcept nogil
    int hybrd(hybrd_callback callback, void *userdata, int n, double *state,
              double *residual, double xtol, int maxfev, int ml, int mu,
              double epsfcn, double *diag, int mode, double factor,
              int nprint, int *nfev, double *fjac, int ldfjac,
              double *r, int lr, double *qtf, double *wa1, double *wa2,
              double *wa3, double *wa4) nogil


cdef struct Problem:
    int system
    const double *parameters
    int callback_failed


cdef int _callback(void *userdata, int n, const double *state,
                   double *out, int iflag) noexcept nogil:
    cdef Problem *problem = <Problem *>userdata
    cdef const double *p = problem.parameters
    cdef int i
    if iflag == 0:
        return 0
    if problem.system == WATER_SYSTEM:
        out[0] = water_residual(state[0], p[0], p[1], p[2], p[3], p[4])
    elif problem.system == HYDROGEN_SYSTEM:
        out[0] = h_residual(state[0], p[0], p[1], p[2], p[3], p[4])
    elif problem.system == BIOGEOCHEM_SYSTEM:
        biogeochem_residual(
            state, p[0], p[1], p[2], p[3], p[4], p[5],
            p[6], p[7], p[8], p[9], p[10], p[11], p[12],
            p[13], p[14], p[15], p[16], p[17], p[18], p[19],
            p[20], p[21], p[22], p[23], p[24], out)
    elif problem.system == CEC_CALCIUM_SYSTEM:
        out[0] = cec_calcium_residual(
            state[0], p[0], p[1], p[2], p[3], p[4], p[5], p[6],
            p[7], p[8], p[9], p[10], p[11])
    elif problem.system == TOTAL_TO_CEC_SYSTEM:
        total_to_cec_residual(
            state, p[0], p[1], p[2], p[3], p[4], p[5],
            p[6], p[7], p[8], p[9], p[10], p[11], p[12],
            p[13], p[14], p[15], p[16], p[17], p[18], p[19], out)
    elif problem.system == KELLAND_SYSTEM:
        kelland_residual(
            state, p[0], p[1], p[2], p[3], p[4], p[5],
            p[6], p[7], p[8], p[9], p[10], p[11], p[12],
            p[13], p[14], out)
    else:
        problem.callback_failed = 1
        return -1
    for i in range(n):
        if not isfinite(out[i]):
            problem.callback_failed = 1
            return -1
    return 0


cdef bint _system_layout(int system, int *n, int *parameters) noexcept nogil:
    if system == WATER_SYSTEM or system == HYDROGEN_SYSTEM:
        n[0], parameters[0] = 1, 5
    elif system == BIOGEOCHEM_SYSTEM:
        n[0], parameters[0] = 16, 25
    elif system == CEC_CALCIUM_SYSTEM:
        n[0], parameters[0] = 1, 12
    elif system == TOTAL_TO_CEC_SYSTEM:
        n[0], parameters[0] = 13, 20
    elif system == KELLAND_SYSTEM:
        n[0], parameters[0] = 5, 15
    else:
        return False
    return True


cdef int _solve_native(int system, int n, const double *parameters,
                              int parameter_count, double *state,
                              double *residual, double xtol, int maxfev,
                              double *work, int work_len, int *nfev) noexcept nogil:
    """Shared cminpack invocation for the Python and Numba native interfaces.

    Return MINPACK's status, -1 for an invalid residual, or -2 for an invalid
    layout/settings. Parameter and work buffers belong to the caller.
    """
    cdef Problem problem
    cdef int expected_n, expected_parameters, status
    cdef int lr
    cdef double *diag
    cdef double *fjac
    cdef double *r
    cdef double *qtf
    cdef double *wa1
    cdef double *wa2
    cdef double *wa3
    cdef double *wa4
    if not _system_layout(system, &expected_n, &expected_parameters):
        return -2
    lr = n * (n + 1) // 2
    if (n != expected_n or parameter_count != expected_parameters or
            state == NULL or residual == NULL or parameters == NULL or
            work == NULL or xtol < 0 or maxfev <= 0 or
            work_len < n * n + lr + 6 * n):
        return -2
    diag = work
    fjac = diag + n
    r = fjac + n * n
    qtf = r + lr
    wa1 = qtf + n
    wa2 = wa1 + n
    wa3 = wa2 + n
    wa4 = wa3 + n
    problem.system = system
    problem.parameters = parameters
    problem.callback_failed = 0
    status = hybrd(
        _callback, &problem, n, state, residual, xtol, maxfev,
        n - 1, n - 1, DBL_EPSILON, diag, 1, 100.0, 0, nfev,
        fjac, n, r, lr, qtf, wa1, wa2, wa3, wa4,
    )
    if problem.callback_failed:
        return -1
    if _callback(&problem, n, state, residual, 1) != 0:
        return -1
    return status


cdef public int smew_solve(int system, int n, const double *parameters,
                           int parameter_count, double *state,
                           double *residual, double xtol, int maxfev,
                           double *work, int work_len) noexcept nogil:
    """C entry point called directly by Numba; buffers belong to the caller."""
    cdef int nfev = 0
    return _solve_native(system, n, parameters, parameter_count, state,
                         residual, xtol, maxfev, work, work_len, &nfev)


class Workspace:
    """Caller-owned cminpack scratch arrays for repeated solves of one size."""

    def __init__(self, dimension):
        if not isinstance(dimension, int) or dimension < 1:
            raise ValueError("workspace dimension must be a positive integer")
        self.dimension = dimension
        self.residual = np.empty(dimension, dtype=np.float64)
        self.work = np.empty(
            dimension * dimension + dimension * (dimension + 1) // 2
            + 6 * dimension,
            dtype=np.float64,
        )


def solve(int system_id, initial, parameters,
           *, xtol=1.4901161193847656e-8, maxfev=None, workspace=None):
    """Call full hybrd with the same defaults used by SciPy fsolve.

    The parameter order is the order of the matching residual function after
    its state argument, excluding its output argument. When a workspace is
    supplied, the returned residual array is reused by its next solve.
    """
    cdef int n, parameter_count
    if not _system_layout(system_id, &n, &parameter_count):
        raise ValueError("Unknown equation system")
    state = np.array(initial, dtype=np.float64, order="C", copy=True, ndmin=1)
    params = np.asarray(parameters, dtype=np.float64)
    if state.ndim != 1 or state.size != n:
        raise ValueError(f"solver needs {n} state values")
    if params.ndim != 1 or params.size != parameter_count:
        raise ValueError(f"solver needs {parameter_count} parameters")
    if not params.flags.c_contiguous:
        params = np.ascontiguousarray(params)
    if xtol < 0:
        raise ValueError("xtol must be nonnegative")
    if maxfev is None:
        maxfev = 200 * (n + 1)
    if maxfev <= 0:
        raise ValueError("maxfev must be positive")

    if workspace is None:
        workspace = Workspace(n)
    if not isinstance(workspace, Workspace) or workspace.dimension != n:
        raise ValueError(f"solver needs a workspace for {n} unknowns")
    residual = workspace.residual
    work = workspace.work
    if residual.size != n or work.size < n * n + n * (n + 1) // 2 + 6 * n:
        raise ValueError("workspace buffers have the wrong size")

    cdef double[::1] x_view = state
    cdef double[::1] p_view = params
    cdef double[::1] f_view = residual
    cdef double[::1] work_view = work
    cdef int status, nfev = 0
    cdef int evaluation_budget = maxfev
    cdef double tolerance = xtol
    cdef int work_len = work_view.shape[0]
    with nogil:
        status = _solve_native(
            system_id, n, &p_view[0], parameter_count, &x_view[0],
            &f_view[0], tolerance, evaluation_budget, &work_view[0],
            work_len, &nfev,
        )
    if status == -1:
        raise RuntimeError("native residual evaluation failed during solve or final evaluation")
    if status <= 0:
        raise ValueError("native solver received invalid buffers or settings")
    if status != 1:
        warnings.warn(
            f"MINPACK stopped with status {status}", RuntimeWarning,
            stacklevel=3,
        )
    return state, residual, status, nfev
