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
    if problem.system == 0:
        out[0] = water_residual(state[0], p[0], p[1], p[2], p[3], p[4])
    elif problem.system == 1:
        out[0] = h_residual(state[0], p[0], p[1], p[2], p[3], p[4])
    elif problem.system == 2:
        biogeochem_residual(
            state, p[0], p[1], p[2], p[3], p[4], p[5],
            p[6], p[7], p[8], p[9], p[10], p[11], p[12],
            p[13], p[14], p[15], p[16], p[17], p[18], p[19],
            p[20], p[21], p[22], p[23], p[24], out)
    elif problem.system == 3:
        out[0] = cec_calcium_residual(
            state[0], p[0], p[1], p[2], p[3], p[4], p[5], p[6],
            p[7], p[8], p[9], p[10], p[11])
    elif problem.system == 4:
        total_to_cec_residual(
            state, p[0], p[1], p[2], p[3], p[4], p[5],
            p[6], p[7], p[8], p[9], p[10], p[11], p[12],
            p[13], p[14], p[15], p[16], p[17], p[18], p[19], out)
    elif problem.system == 5:
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


def _solve(int system_id, int n, int parameter_count, initial, parameters,
           *, xtol=1.4901161193847656e-8, maxfev=None, workspace=None):
    """Call full hybrd with the same defaults used by SciPy fsolve.

    The parameter order is the order of the matching residual function after
    its state argument, excluding its output argument. When a workspace is
    supplied, the returned residual array is reused by its next solve.
    """
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
    cdef int status, nfev = 0, final_status = 0
    cdef int evaluation_budget = maxfev
    cdef double tolerance = xtol
    cdef Problem problem
    cdef double *diag = &work_view[0]
    cdef double *fjac = diag + n
    cdef double *r = fjac + n * n
    cdef double *qtf = r + n * (n + 1) // 2
    cdef double *wa1 = qtf + n
    cdef double *wa2 = wa1 + n
    cdef double *wa3 = wa2 + n
    cdef double *wa4 = wa3 + n
    problem.system = system_id
    problem.parameters = &p_view[0]
    problem.callback_failed = 0
    with nogil:
        status = hybrd(
            _callback, &problem, n, &x_view[0], &f_view[0], tolerance,
            evaluation_budget, n - 1, n - 1, DBL_EPSILON, diag, 1,
            100.0, 0, &nfev, fjac, n, r, n * (n + 1) // 2,
            qtf, wa1, wa2, wa3, wa4,
        )
    if problem.callback_failed:
        raise RuntimeError("native residual evaluation failed during solve")
    with nogil:
        final_status = _callback(&problem, n, &x_view[0], &f_view[0], 1)
    if final_status != 0:
        raise RuntimeError("final native residual evaluation failed")
    if status == 0:
        raise ValueError("native solver received invalid settings")
    if status != 1:
        warnings.warn(
            f"MINPACK stopped with status {status}", RuntimeWarning,
            stacklevel=3,
        )
    return state, residual, status, nfev


def solve_water(initial, parameters, **options):
    return _solve(0, 1, 5, initial, parameters, **options)


def solve_hydrogen(initial, parameters, **options):
    return _solve(1, 1, 5, initial, parameters, **options)


def solve_biogeochem(initial, parameters, **options):
    return _solve(2, 16, 25, initial, parameters, **options)


def solve_cec_calcium(initial, parameters, **options):
    return _solve(3, 1, 12, initial, parameters, **options)


def solve_total_to_cec(initial, parameters, **options):
    return _solve(4, 13, 20, initial, parameters, **options)


def solve_kelland(initial, parameters, **options):
    return _solve(5, 5, 15, initial, parameters, **options)
