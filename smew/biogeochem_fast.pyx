# cython: language_level=3
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport pow, sqrt, exp, log, fabs
from libc.string cimport memset, memcpy

# 1. Import Cython to use the high-performance decorators
cimport cython

# Ensure the NumPy C-API is initialized
np.import_array()

# Define the c-functions from cminpack
cdef extern from "cminpack.h":
    ctypedef int (*cminpack_func_nn)(void *p, int n, const double *x, double *fvec, int iflag) nogil
    int hybrd1(cminpack_func_nn fcn, void *p, int n, double *x, double *fvec, double tol, double *wa, int lwa) nogil


cdef struct EquationArgs:
    double Alk_tot, n, Zr, s, IC_tot, k1, k2, k_H, k_w, CEC_tot, conv_Al, Al_tot
    double K1, K2, K3, K4
    double Mg_tot, Ca_tot, Na_tot, K_tot
    double K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H
    int func_calls
    int jac_calls


# @cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int biogeochem_equations_c_fd(void *p, int n, const double *state, double *fvec, int iflag) noexcept nogil:
    cdef EquationArgs* args = <EquationArgs*>p

    cdef double Alk   = state[0]
    cdef double CO2_w = state[1]
    cdef double H     = state[2]
    # Prevents floating-point overflow when dividing by H^2 in basic soils
    if H < 1e-15:
        H = 1e-15
    cdef double R_alk = state[3]
    cdef double Al_w  = state[4]
    cdef double Al    = state[5]
    cdef double Mg    = state[6]
    cdef double Ca    = state[7]
    cdef double Na    = state[8]
    cdef double K     = state[9]
    cdef double f_Al  = state[10]
    cdef double f_Mg  = state[11]
    cdef double f_Na  = state[12]
    cdef double f_K   = state[13]
    cdef double f_H   = state[14]
    cdef double f_Ca  = state[15]

    # Precompute common terms
    cdef double nZrs1000 = args.n * args.Zr * args.s * 1000.0
    cdef double H2 = H * H
    cdef double H3 = H2 * H
    cdef double H4 = H3 * H

    # Equation 1: Alkalinity balance
    fvec[0] = (args.Alk_tot - R_alk) - Alk * nZrs1000
    # Equation 2: Inorganic Carbon balance
    fvec[1] = args.IC_tot - (CO2_w * (1.0 + args.k1/H + args.k2 * args.k1 / H2) * args.s + (CO2_w / args.k_H) * (1.0 - args.s)) * (args.n * args.Zr * 1000.0)
    # Equation 3: Carbonate system equilibrium
    fvec[2] = (args.k1 * CO2_w / H + 2.0 * args.k1 * args.k2 * CO2_w / H2 - H + args.k_w / H) - Alk
    # Equation 4: Reserve alkalinity from CEC
    fvec[3] = R_alk - (f_Mg + f_Ca + f_Na + f_K) * args.CEC_tot
    # Equation 5: Aluminum mass balance
    fvec[4] = Al_w * nZrs1000 + (f_Al / 3.0) * args.CEC_tot * args.conv_Al - args.Al_tot
    # Equation 6: Aluminum speciation
    fvec[5] = Al - (H4 / (H4 + H3 * args.K1 + H2 * args.K1 * args.K2 + H * args.K1 * args.K2 * args.K3 + args.K1 * args.K2 * args.K3 * args.K4)) * Al_w
    # Equations 7-10: Cation mass balances (Mg, Ca, Na, K)
    fvec[6] = Mg * nZrs1000 + (f_Mg / 2.0) * args.CEC_tot - args.Mg_tot
    fvec[7] = Ca * nZrs1000 + (f_Ca / 2.0) * args.CEC_tot - args.Ca_tot
    fvec[8] = Na * nZrs1000 + f_Na * args.CEC_tot - args.Na_tot
    fvec[9] = K * nZrs1000 + f_K * args.CEC_tot - args.K_tot
    # Equations 11-15: Gaines-Thomas exchange equations
    # Equation 11: Aluminum exchange
    fvec[10] = f_Al - (Al / args.conv_Al) * sqrt( (fabs(f_Ca) * fabs(f_Ca) * fabs(f_Ca)) / (args.K_Ca_Al * fabs(Ca) * fabs(Ca) * fabs(Ca)) )
    # Equation 12: Magnesium exchange (No fractional power, so no fabs needed)
    fvec[11] = f_Mg - Mg * (f_Ca / (args.K_Ca_Mg * Ca))
    # Equation 13: Sodium exchange
    fvec[12] = f_Na - Na * sqrt( fabs(f_Ca) / (args.K_Ca_Na * fabs(Ca)) )
    # Equation 14: Potassium exchange
    fvec[13] = f_K - K * sqrt( fabs(f_Ca) / (args.K_Ca_K * fabs(Ca)) )
    # Equation 15: Hydrogen exchange
    fvec[14] = f_H - H * sqrt( fabs(f_Ca) / (args.K_Ca_H * fabs(Ca)) )
    # Equation 16: Sum of exchange fractions must be 1
    fvec[15] = 1.0 - (f_Ca + f_Al + f_Mg + f_Na + f_K + f_H)

    return 0


cdef public int solve_biogeochem_eq_fd(
    double *x0,          # Input/Output pointer (Numba passes .ctypes.data)
    double *residuals,   # Output pointer for errors
    double Alk_tot_in, double n_in, double Zr_in, double s_in, double IC_tot_in,
    double k1_in, double k2_in, double k_H_in, double k_w_in, double CEC_tot_in,
    double conv_Al_in, double Al_tot_in, double K1_in, double K2_in, double K3_in,
    double K4_in, double Mg_tot_in, double Ca_tot_in, double Na_tot_in, double K_tot_in,
    double K_Ca_Al_in, double K_Ca_Mg_in, double K_Ca_Na_in, double K_Ca_K_in, double K_Ca_H_in
) nogil:
    cdef int n_vars = 16
    # MINPACK requires a very specific workspace array size for hybrd1.
    # The formula from the MINPACK documentation is: LWA >= N * (3*N + 13) / 2
    cdef int lwa = 488
    cdef double tol = 1e-12
    cdef int status

    cdef double *fvec = <double *>malloc(n_vars * sizeof(double))
    cdef double *wa = <double *>malloc(lwa * sizeof(double))

    memset(fvec, 0, n_vars*sizeof(double))
    memset(wa, 0, lwa*sizeof(double))

    if not fvec or not wa:
        if fvec: free(fvec)
        if wa: free(wa)
        return -999

    cdef EquationArgs args
    args.Alk_tot = Alk_tot_in
    args.n = n_in
    args.Zr = Zr_in
    args.s = s_in
    args.IC_tot = IC_tot_in
    args.k1 = k1_in
    args.k2 = k2_in
    args.k_H = k_H_in
    args.k_w = k_w_in
    args.CEC_tot = CEC_tot_in
    args.conv_Al = conv_Al_in
    args.Al_tot = Al_tot_in
    args.K1 = K1_in
    args.K2 = K2_in
    args.K3 = K3_in
    args.K4 = K4_in
    args.Mg_tot = Mg_tot_in
    args.Ca_tot = Ca_tot_in
    args.Na_tot = Na_tot_in
    args.K_tot = K_tot_in
    args.K_Ca_Al = K_Ca_Al_in
    args.K_Ca_Mg = K_Ca_Mg_in
    args.K_Ca_Na = K_Ca_Na_in
    args.K_Ca_K = K_Ca_K_in
    args.K_Ca_H = K_Ca_H_in

    # when cminpack runs, it will update in place x0 so we don't need to retrieve and save its values
    status = hybrd1(biogeochem_equations_c_fd, <void*>&args, n_vars, x0, fvec, tol, wa, lwa)

    # Don't copy fvec for the residuals it's unreliable.
    # Re-evaluate at the actual returned solution:
    biogeochem_equations_c_fd(<void*>&args, n_vars, x0, residuals, 1)

    free(fvec)
    free(wa)

    return status


# ==============================================================================
# 1D Solver for the Rainwater H+ Equilibrium
# ==============================================================================
cdef struct EqWaterArgs:
    double Alk_rain, k1, k2, CO2_w_rain, k_w

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int eq_water_c(void *p, int n, const double *state, double *fvec, int iflag) noexcept nogil:
    cdef EqWaterArgs* args = <EqWaterArgs*>p
    cdef double H_rain = state[0]

    fvec[0] = args.Alk_rain - (args.k1 * args.CO2_w_rain / H_rain + 2.0 * args.k1 * args.k2 * args.CO2_w_rain / (H_rain * H_rain) - H_rain + args.k_w / H_rain)
    return 0

cdef public int solve_water_eq_fd(
    double *x0,
    double Alk_rain_in, double k1_in, double k2_in, double CO2_w_rain_in, double k_w_in
) nogil:
    cdef int n_vars = 1
    cdef int lwa = 10
    cdef double tol = 1e-12
    cdef int status

    cdef double *fvec = <double *>malloc(n_vars * sizeof(double))
    cdef double *wa = <double *>malloc(lwa * sizeof(double))

    if not fvec or not wa:
        if fvec: free(fvec)
        if wa: free(wa)
        return -999

    cdef EqWaterArgs args
    args.Alk_rain = Alk_rain_in
    args.k1 = k1_in
    args.k2 = k2_in
    args.CO2_w_rain = CO2_w_rain_in
    args.k_w = k_w_in

    status = hybrd1(eq_water_c, <void*>&args, n_vars, x0, fvec, tol, wa, lwa)

    free(fvec)
    free(wa)

    return status


# ==============================================================================
# 1D Solver for the H0 Fallback Guess
# ==============================================================================
cdef struct EqHArgs:
    double k1, k2, CO2_w0, k_w, Alk0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int eq_H_c(void *p, int n, const double *state, double *fvec, int iflag) noexcept nogil:
    cdef EqHArgs* args = <EqHArgs*>p
    cdef double H0 = state[0]

    fvec[0] = (args.k1 * args.CO2_w0 / H0 + 2.0 * args.k1 * args.k2 * args.CO2_w0 / (H0 * H0) - H0 + args.k_w / H0) - args.Alk0
    return 0

cdef public int solve_H_eq_fd(
    double *x0,
    double k1_in, double k2_in, double CO2_w0_in, double k_w_in, double Alk0_in
) nogil:
    cdef int n_vars = 1
    cdef int lwa = 10
    cdef double tol = 1e-12
    cdef int status

    cdef double *fvec = <double *>malloc(n_vars * sizeof(double))
    cdef double *wa = <double *>malloc(lwa * sizeof(double))

    if not fvec or not wa:
        if fvec: free(fvec)
        if wa: free(wa)
        return -999

    cdef EqHArgs args
    args.k1 = k1_in
    args.k2 = k2_in
    args.CO2_w0 = CO2_w0_in
    args.k_w = k_w_in
    args.Alk0 = Alk0_in

    status = hybrd1(eq_H_c, <void*>&args, n_vars, x0, fvec, tol, wa, lwa)

    free(fvec)
    free(wa)

    return status
