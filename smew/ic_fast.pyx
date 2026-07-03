# cython: language_level=3
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport pow, sqrt, exp, log, fabs
from libc.string cimport memset, memcpy
cimport cython

np.import_array()

cdef extern from "cminpack.h":
    ctypedef int (*cminpack_func_nn)(void *p, int n, const double *x, double *fvec, int iflag) nogil
    int hybrd1(cminpack_func_nn fcn, void *p, int n, double *x, double *fvec, double tol, double *wa, int lwa) nogil

# ==============================================================================
# 1D Solver for eqf_Ca (conc_to_f_CEC)
# ==============================================================================
cdef struct EqfCaArgs:
    double Al, conv_Al, K_Ca_Al, Ca, Mg, K_Ca_Mg, Na, K_Ca_Na, K_K, K_Ca_K, H, K_Ca_H

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int eqf_Ca_c(void *p, int n, const double *state, double *fvec, int iflag) noexcept nogil:
    cdef EqfCaArgs* args = <EqfCaArgs*>p
    cdef double f_Ca = state[0]

    # To avoid negative base for power, use fabs where fractional powers exist
    cdef double term_Al = (args.Al / args.conv_Al) * sqrt((f_Ca * f_Ca * f_Ca) / (args.K_Ca_Al * args.Ca * args.Ca * args.Ca))
    cdef double term_Mg = args.Mg * (f_Ca / (args.K_Ca_Mg * args.Ca))
    cdef double term_Na = args.Na * sqrt(f_Ca / (args.K_Ca_Na * args.Ca))
    cdef double term_K = args.K_K * sqrt(f_Ca / (args.K_Ca_K * args.Ca))
    cdef double term_H = args.H * sqrt(f_Ca / (args.K_Ca_H * args.Ca))

    fvec[0] = 1.0 - (f_Ca + term_Al + term_Mg + term_Na + term_K + term_H)
    return 0

cdef public int solve_eqf_Ca_fd(
    double *x0,
    double Al_in, double conv_Al_in, double K_Ca_Al_in, double Ca_in,
    double Mg_in, double K_Ca_Mg_in, double Na_in, double K_Ca_Na_in,
    double K_in, double K_Ca_K_in, double H_in, double K_Ca_H_in
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

    cdef EqfCaArgs args
    args.Al = Al_in
    args.conv_Al = conv_Al_in
    args.K_Ca_Al = K_Ca_Al_in
    args.Ca = Ca_in
    args.Mg = Mg_in
    args.K_Ca_Mg = K_Ca_Mg_in
    args.Na = Na_in
    args.K_Ca_Na = K_Ca_Na_in
    args.K_K = K_in
    args.K_Ca_K = K_Ca_K_in
    args.H = H_in
    args.K_Ca_H = K_Ca_H_in

    status = hybrd1(eqf_Ca_c, <void*>&args, n_vars, x0, fvec, tol, wa, lwa)
    free(fvec)
    free(wa)
    return status

# ==============================================================================
# 13D Solver for total_to_f_CEC_and_conc
# ==============================================================================
cdef struct EqTotalArgs:
    double n_p, Zr, s0, CEC_tot, conv_Al
    double K1, K2, K3, K4
    double Mg_tot, Ca_tot, Na_tot, K_tot
    double K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H
    double H, f_acid

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int eq_total_c(void *p, int n, const double *state, double *fvec, int iflag) noexcept nogil:
    cdef EqTotalArgs* args = <EqTotalArgs*>p
    cdef double Al_w = state[0]
    cdef double Al = state[1]
    cdef double Al_tot = state[2]
    cdef double Mg = state[3]
    cdef double Ca = state[4]
    cdef double Na = state[5]
    cdef double K = state[6]
    cdef double f_Mg = state[7]
    cdef double f_Na = state[8]
    cdef double f_K = state[9]
    cdef double f_Ca = state[10]
    cdef double f_Al = state[11]
    cdef double f_H = state[12]

    cdef double nZrs1000 = args.n_p * args.Zr * args.s0 * 1000.0
    cdef double H2 = args.H * args.H
    cdef double H3 = H2 * args.H
    cdef double H4 = H3 * args.H

    fvec[0] = Al_w * nZrs1000 + (f_Al/3.0) * args.CEC_tot * args.conv_Al - Al_tot
    fvec[1] = Al - (H4 / (H4 + H3*args.K1 + H2*args.K1*args.K2 + args.H*args.K1*args.K2*args.K3 + args.K1*args.K2*args.K3*args.K4)) * Al_w
    fvec[2] = Mg * nZrs1000 + (f_Mg/2.0) * args.CEC_tot - args.Mg_tot
    fvec[3] = Ca * nZrs1000 + (f_Ca/2.0) * args.CEC_tot - args.Ca_tot
    fvec[4] = Na * nZrs1000 + f_Na * args.CEC_tot - args.Na_tot
    fvec[5] = K * nZrs1000 + f_K * args.CEC_tot - args.K_tot
    fvec[6] = f_Al - (Al/args.conv_Al) * sqrt((f_Ca*f_Ca*f_Ca) / (args.K_Ca_Al * Ca*Ca*Ca))
    fvec[7] = f_H - args.H * sqrt(f_Ca / (args.K_Ca_H * Ca))
    fvec[8] = f_H + f_Al - args.f_acid
    fvec[9] = f_Mg - Mg * (f_Ca / (args.K_Ca_Mg * Ca))
    fvec[10] = f_Na - Na * sqrt(f_Ca / (args.K_Ca_Na * Ca))
    fvec[11] = f_K - K * sqrt(f_Ca / (args.K_Ca_K * Ca))
    fvec[12] = 1.0 - (f_Ca + f_Al + f_Mg + f_Na + f_K + f_H)

    return 0

cdef public int solve_total_eq_fd(
    double *x0,
    double n_in, double Zr_in, double s0_in, double CEC_tot_in, double conv_Al_in,
    double K1_in, double K2_in, double K3_in, double K4_in,
    double Mg_tot_in, double Ca_tot_in, double Na_tot_in, double K_tot_in, double H_in,
    double K_Ca_Al_in, double K_Ca_Mg_in, double K_Ca_Na_in, double K_Ca_K_in, double K_Ca_H_in,
    double f_acid_in
) nogil:
    cdef int n_vars = 13
    cdef int lwa = 338
    cdef double tol = 1e-14
    cdef int status

    cdef double *fvec = <double *>malloc(n_vars * sizeof(double))
    cdef double *wa = <double *>malloc(lwa * sizeof(double))
    if not fvec or not wa:
        if fvec: free(fvec)
        if wa: free(wa)
        return -999

    cdef EqTotalArgs args
    args.n_p = n_in
    args.Zr = Zr_in
    args.s0 = s0_in
    args.CEC_tot = CEC_tot_in
    args.conv_Al = conv_Al_in
    args.K1 = K1_in; args.K2 = K2_in; args.K3 = K3_in; args.K4 = K4_in
    args.Mg_tot = Mg_tot_in; args.Ca_tot = Ca_tot_in; args.Na_tot = Na_tot_in; args.K_tot = K_tot_in
    args.H = H_in
    args.K_Ca_Al = K_Ca_Al_in; args.K_Ca_Mg = K_Ca_Mg_in; args.K_Ca_Na = K_Ca_Na_in
    args.K_Ca_K = K_Ca_K_in; args.K_Ca_H = K_Ca_H_in
    args.f_acid = f_acid_in

    status = hybrd1(eq_total_c, <void*>&args, n_vars, x0, fvec, tol, wa, lwa)
    free(fvec)
    free(wa)
    return status

# ==============================================================================
# 5D Solver for Kelland
# ==============================================================================
cdef struct EqKellandArgs:
    double Al_w, n_p, Zr, s0, CEC_tot, conv_Al
    double Ca, Ca_tot, Al, K_Ca_Al, H, K_Ca_H
    double f_Mg, f_Na, f_K

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int eq_kelland_c(void *p, int n, const double *state, double *fvec, int iflag) noexcept nogil:
    cdef EqKellandArgs* args = <EqKellandArgs*>p
    cdef double Al_tot = state[0]
    cdef double CaCO3 = state[1]
    cdef double f_Al = state[2]
    cdef double f_H = state[3]
    cdef double f_Ca = state[4]

    fvec[0] = args.Al_w * args.n_p * args.Zr * args.s0 * 1000.0 + (f_Al/3.0)*args.CEC_tot*args.conv_Al - Al_tot
    fvec[1] = args.Ca * args.n_p * args.Zr * args.s0 * 1000.0 + (f_Ca/2.0)*args.CEC_tot + CaCO3 - args.Ca_tot
    fvec[2] = f_Al - (args.Al/args.conv_Al) * sqrt((f_Ca*f_Ca*f_Ca)/(args.K_Ca_Al * args.Ca*args.Ca*args.Ca))
    fvec[3] = f_H - args.H * sqrt(f_Ca/(args.K_Ca_H * args.Ca))
    fvec[4] = 1.0 - (f_Ca + f_Al + args.f_Mg + args.f_Na + args.f_K + f_H)

    return 0

cdef public int solve_kelland_eq_fd(
    double *x0,
    double Al_w_in, double n_in, double Zr_in, double s0_in, double CEC_tot_in, double conv_Al_in,
    double Ca_in, double Ca_tot_in, double Al_in, double K_Ca_Al_in, double H_in, double K_Ca_H_in,
    double f_Mg_in, double f_Na_in, double f_K_in
) nogil:
    cdef int n_vars = 5
    cdef int lwa = 70
    cdef double tol = 1e-14
    cdef int status

    cdef double *fvec = <double *>malloc(n_vars * sizeof(double))
    cdef double *wa = <double *>malloc(lwa * sizeof(double))
    if not fvec or not wa:
        if fvec: free(fvec)
        if wa: free(wa)
        return -999

    cdef EqKellandArgs args
    args.Al_w = Al_w_in; args.n_p = n_in; args.Zr = Zr_in; args.s0 = s0_in
    args.CEC_tot = CEC_tot_in; args.conv_Al = conv_Al_in
    args.Ca = Ca_in; args.Ca_tot = Ca_tot_in; args.Al = Al_in
    args.K_Ca_Al = K_Ca_Al_in; args.H = H_in; args.K_Ca_H = K_Ca_H_in
    args.f_Mg = f_Mg_in; args.f_Na = f_Na_in; args.f_K = f_K_in

    status = hybrd1(eq_kelland_c, <void*>&args, n_vars, x0, fvec, tol, wa, lwa)
    free(fvec)
    free(wa)
    return status
