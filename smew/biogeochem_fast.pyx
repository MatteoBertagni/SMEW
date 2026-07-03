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

    # ctypedef int (*cminpack_func_der_nn)(void *p, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag) nogil
    # int hybrj1(cminpack_func_der_nn fcn, void *p, int n, double *x, double *fvec, double *fjac, int ldfjac, double tol, double *wa, int lwa) nogil

cdef struct EquationArgs:
    double Alk_tot, n, Zr, s, IC_tot, k1, k2, k_H, k_w, CEC_tot, conv_Al, Al_tot
    double K1, K2, K3, K4
    double Mg_tot, Ca_tot, Na_tot, K_tot
    double K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H
    int func_calls
    int jac_calls


# # equation with Jacobian and log space
# @cython.boundscheck(False)
# @cython.wraparound(False)
# cdef int biogeochem_equations_c_jac_log(void *p, int n, const double *state, double *fvec, double *fjac, int ldfjac, int iflag) nogil noexcept:
#     cdef EquationArgs* args = <EquationArgs*>p
#
#     cdef double nZr1000 = args.n * args.Zr * 1000.0
#     cdef double nZrs1000 = nZr1000 * args.s
#
#     # Alkalinity can be negative, so it stays entirely in linear space.
#     cdef double Alk   = state[0]
#     # unpack all the other variables from log space to linear space
#     cdef double CO2_w = exp(state[1])
#     cdef double H     = exp(state[2])
#     cdef double R_alk = exp(state[3])
#     cdef double Al_w  = exp(state[4])
#     cdef double Al    = exp(state[5])
#     cdef double Mg    = exp(state[6])
#     cdef double Ca    = exp(state[7])
#     cdef double Na    = exp(state[8])
#     cdef double K     = exp(state[9])
#     cdef double f_Al  = exp(state[10])
#     cdef double f_Mg  = exp(state[11])
#     cdef double f_Na  = exp(state[12])
#     cdef double f_K   = exp(state[13])
#     cdef double f_H   = exp(state[14])
#     cdef double f_Ca  = exp(state[15])
#
#     # Store linear values in an array for the Jacobian chain rule later
#     cdef double C[16]
#     C[0]=Alk; C[1]=CO2_w; C[2]=H; C[3]=R_alk; C[4]=Al_w; C[5]=Al; C[6]=Mg; C[7]=Ca;
#     C[8]=Na; C[9]=K; C[10]=f_Al; C[11]=f_Mg; C[12]=f_Na; C[13]=f_K; C[14]=f_H; C[15]=f_Ca
#
#     cdef double H2 = H * H
#     cdef double H3 = H2 * H
#     cdef double H4 = H3 * H
#
#     cdef double denom, root_Ca_ratio
#     cdef double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
#     cdef double x11, x12, x13, x14, x15, x16, x17, x18, x19, x20
#     cdef double x21, x22, x23, x24, x25, x26, x27, x28, x29, x30
#     cdef double x31, x32, x33, x34, x35, x36, x37
#
#     cdef int i_jac, i_eq, j_var
#
#     if iflag == 1:
#         args.func_calls += 1
#
#         fvec[0] = (args.Alk_tot - R_alk) - Alk * nZrs1000
#         fvec[1] = args.IC_tot - (CO2_w * (1.0 + args.k1/H + args.k2 * args.k1 / H2) * args.s + (CO2_w / args.k_H) * (1.0 - args.s)) * nZr1000
#         fvec[2] = (args.k1 * CO2_w / H + 2.0 * args.k1 * args.k2 * CO2_w / H2 - H + args.k_w / H) - Alk
#         fvec[3] = R_alk - (f_Mg + f_Ca + f_Na + f_K) * args.CEC_tot
#         fvec[4] = Al_w * nZrs1000 + (f_Al / 3.0) * args.CEC_tot * args.conv_Al - args.Al_tot
#
#         denom = H4 + H3 * args.K1 + H2 * args.K1 * args.K2 + H * args.K1 * args.K2 * args.K3 + args.K1 * args.K2 * args.K3 * args.K4
#         fvec[5] = Al - (H4 / denom) * Al_w
#
#         fvec[6] = Mg * nZrs1000 + (f_Mg / 2.0) * args.CEC_tot - args.Mg_tot
#         fvec[7] = Ca * nZrs1000 + (f_Ca / 2.0) * args.CEC_tot - args.Ca_tot
#         fvec[8] = Na * nZrs1000 + f_Na * args.CEC_tot - args.Na_tot
#         fvec[9] = K * nZrs1000 + f_K * args.CEC_tot - args.K_tot
#
#         root_Ca_ratio = sqrt(f_Ca / Ca)
#         fvec[10] = f_Al - (Al / args.conv_Al) * sqrt( (f_Ca * f_Ca * f_Ca) / (args.K_Ca_Al * Ca * Ca * Ca) )
#         fvec[11] = f_Mg - Mg * (f_Ca / (args.K_Ca_Mg * Ca))
#         fvec[12] = f_Na - Na * (root_Ca_ratio / sqrt(args.K_Ca_Na))
#         fvec[13] = f_K  - K  * (root_Ca_ratio / sqrt(args.K_Ca_K))
#         fvec[14] = f_H  - H  * (root_Ca_ratio / sqrt(args.K_Ca_H))
#
#         fvec[15] = 1.0 - (f_Ca + f_Al + f_Mg + f_Na + f_K + f_H)
#
#     elif iflag == 2:
#         args.jac_calls += 1
#
#         # initialize fjac values
#         for i_jac in range(256):
#             fjac[i_jac] = 0.0
#
#         x0 = 1000.0*args.Zr*args.n
#         x1 = args.s*x0
#         x2 = args.k1/H
#         x3 = (H*H)
#         x4 = 1.0/x3
#         x5 = args.k1*x4
#         x6 = args.k2*x5
#         x7 = (H*H*H)
#         x8 = args.k1*args.k2/x7
#         x9 = -args.CEC_tot
#         x10 = (H*H*H*H)
#         x11 = args.K1*args.K2
#         x12 = args.K3*x11
#         x13 = args.K1*x3
#         x14 = H*x12 + args.K1*x7 + args.K2*x13 + args.K4*x12 + x10
#         x15 = 1.0/x14
#         x16 = 4*x7
#         x17 = 0.5*args.CEC_tot
#         x18 = pow(Ca, -3.0/2.0)
#         x19 = 1.0/args.conv_Al
#         x20 = sqrt(1.0/args.K_Ca_Al)
#         x21 = pow(f_Ca, 3.0/2.0)*x19*x20
#         x22 = (3.0/2.0)*Al
#         x23 = sqrt(f_Ca)
#         x24 = x18*x23
#         x25 = 1.0/args.K_Ca_Mg
#         x26 = x25/Ca
#         x27 = sqrt(1.0/args.K_Ca_Na)
#         x28 = (1.0/2.0)*Na*x27
#         x29 = pow(Ca, -1.0/2.0)
#         x30 = x23*x29
#         x31 = x29/x23
#         x32 = (1.0/2.0)*x24
#         x33 = sqrt(1.0/args.K_Ca_K)
#         x34 = K*x33
#         x35 = (1.0/2.0)*x31
#         x36 = sqrt(1.0/args.K_Ca_H)
#         x37 = H*x36
#
#         # --- Jacobian Matrix Assignment ---
#         fjac[0] = -x1
#         fjac[48] = -1
#         fjac[17] = -x0*(args.s*(x2 + x6 + 1.0) + (1.0 - args.s)/args.k_H)
#         fjac[33] = -CO2_w*x1*(-x5 - 2*x8)
#         fjac[2] = -1
#         fjac[18] = x2 + 2.0*x6
#         fjac[34] = -CO2_w*x5 - 4.0*CO2_w*x8 - args.k_w*x4 - 1
#         fjac[51] = 1
#         fjac[179] = x9
#         fjac[195] = x9
#         fjac[211] = x9
#         fjac[243] = x9
#         fjac[68] = x1
#         fjac[164] = 0.33333333333333331*args.CEC_tot*args.conv_Al
#         fjac[37] = -Al_w*x10*(-2*H*x11 - x12 - 3*x13 - x16)/(x14*x14) - Al_w*x15*x16
#         fjac[69] = -x10*x15
#         fjac[85] = 1
#         fjac[102] = x1
#         fjac[182] = x17
#         fjac[119] = x1
#         fjac[247] = x17
#         fjac[136] = x1
#         fjac[200] = args.CEC_tot
#         fjac[153] = x1
#         fjac[217] = args.CEC_tot
#         fjac[90] = -x18*x21
#         fjac[122] = x21*x22/pow(Ca, 5.0/2.0)
#         fjac[170] = 1
#         fjac[250] = -x19*x20*x22*x24
#         fjac[107] = -f_Ca*x26
#         fjac[123] = Mg*f_Ca*x25/(Ca*Ca)
#         fjac[187] = 1
#         fjac[251] = -Mg*x26
#         fjac[124] = x24*x28
#         fjac[140] = -x27*x30
#         fjac[204] = 1
#         fjac[252] = -x28*x31
#         fjac[125] = x32*x34
#         fjac[157] = -x30*x33
#         fjac[221] = 1
#         fjac[253] = -x34*x35
#         fjac[46] = -x30*x36
#         fjac[126] = x32*x37
#         fjac[238] = 1
#         fjac[254] = -x35*x37
#         fjac[175] = -1
#         fjac[191] = -1
#         fjac[207] = -1
#         fjac[223] = -1
#         fjac[239] = -1
#         fjac[255] = -1
#
#         # TODO: implement in the fjac assignment directly?
#         # --- 2. THE CHAIN RULE ADJUSTMENT ---
#         # Convert dF/dC into dF/d(lnC).
#         # Since d(lnC) = dC/C, we multiply column j of the Jacobian by C[j].
#         # In MINPACK (column-major), the index is (equation_row + variable_col * 16)
#         # We start the loop at 1 instead of 0 as column 0 (Alk) is linear
#         for j_var in range(1, 16):
#             for i_eq in range(16):
#                 fjac[i_eq + j_var * 16] *= C[j_var]
#
#     return 0
#
#
# # equation with Jacobian
# @cython.boundscheck(False)
# @cython.wraparound(False)
# cdef int biogeochem_equations_c_jac(void *p, int n, const double *state, double *fvec, double *fjac, int ldfjac, int iflag) nogil noexcept:
#     cdef EquationArgs* args = <EquationArgs*>p
#
#     cdef double nZr1000 = args.n * args.Zr * 1000.0
#     cdef double nZrs1000 = args.n * args.Zr * args.s * 1000.0
#
#     cdef double Alk   = state[0]
#     cdef double CO2_w = state[1]
#     cdef double H     = state[2]
#     cdef double R_alk = state[3]
#     cdef double Al_w  = state[4]
#     cdef double Al    = state[5]
#     cdef double Mg    = state[6]
#     cdef double Ca    = state[7]
#     cdef double Na    = state[8]
#     cdef double K     = state[9]
#     cdef double f_Al  = state[10]
#     cdef double f_Mg  = state[11]
#     cdef double f_Na  = state[12]
#     cdef double f_K   = state[13]
#     cdef double f_H   = state[14]
#     cdef double f_Ca  = state[15]
#
#     cdef double H2 = H * H
#     cdef double H3 = H2 * H
#     cdef double H4 = H3 * H
#
#     cdef double denom, root_Ca_ratio
#     cdef double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
#     cdef double x11, x12, x13, x14, x15, x16, x17, x18, x19, x20
#     cdef double x21, x22, x23, x24, x25, x26, x27, x28, x29, x30
#     cdef double x31, x32, x33, x34, x35, x36, x37
#
#     cdef int i_jac, j_var
#
#     if iflag == 1:
#         args.func_calls += 1
#
#         fvec[0] = (args.Alk_tot - R_alk) - Alk * nZrs1000
#         fvec[1] = args.IC_tot - (CO2_w * (1.0 + args.k1/H + args.k2 * args.k1 / H2) * args.s + (CO2_w / args.k_H) * (1.0 - args.s)) * nZr1000
#         fvec[2] = (args.k1 * CO2_w / H + 2.0 * args.k1 * args.k2 * CO2_w / H2 - H + args.k_w / H) - Alk
#         fvec[3] = R_alk - (f_Mg + f_Ca + f_Na + f_K) * args.CEC_tot
#         fvec[4] = Al_w * nZrs1000 + (f_Al / 3.0) * args.CEC_tot * args.conv_Al - args.Al_tot
#
#         denom = H4 + H3 * args.K1 + H2 * args.K1 * args.K2 + H * args.K1 * args.K2 * args.K3 + args.K1 * args.K2 * args.K3 * args.K4
#         fvec[5] = Al - (H4 / denom) * Al_w
#
#         fvec[6] = Mg * nZrs1000 + (f_Mg / 2.0) * args.CEC_tot - args.Mg_tot
#         fvec[7] = Ca * nZrs1000 + (f_Ca / 2.0) * args.CEC_tot - args.Ca_tot
#         fvec[8] = Na * nZrs1000 + f_Na * args.CEC_tot - args.Na_tot
#         fvec[9] = K * nZrs1000 + f_K * args.CEC_tot - args.K_tot
#
#         root_Ca_ratio = sqrt(f_Ca / Ca)
#         fvec[10] = f_Al - (Al / args.conv_Al) * sqrt( (f_Ca * f_Ca * f_Ca) / (args.K_Ca_Al * Ca * Ca * Ca) )
#         fvec[11] = f_Mg - Mg * (f_Ca / (args.K_Ca_Mg * Ca))
#         fvec[12] = f_Na - Na * (root_Ca_ratio / sqrt(args.K_Ca_Na))
#         fvec[13] = f_K  - K  * (root_Ca_ratio / sqrt(args.K_Ca_K))
#         fvec[14] = f_H  - H  * (root_Ca_ratio / sqrt(args.K_Ca_H))
#
#         fvec[15] = 1.0 - (f_Ca + f_Al + f_Mg + f_Na + f_K + f_H)
#
#     elif iflag == 2:
#         args.jac_calls += 1
#
#         # initialize fjac values
#         for i_jac in range(256):
#             fjac[i_jac] = 0.0
#
#         x0 = 1000.0*args.Zr*args.n
#         x1 = args.s*x0
#         x2 = args.k1/H
#         x3 = (H*H)
#         x4 = 1.0/x3
#         x5 = args.k1*x4
#         x6 = args.k2*x5
#         x7 = (H*H*H)
#         x8 = args.k1*args.k2/x7
#         x9 = -args.CEC_tot
#         x10 = (H*H*H*H)
#         x11 = args.K1*args.K2
#         x12 = args.K3*x11
#         x13 = args.K1*x3
#         x14 = H*x12 + args.K1*x7 + args.K2*x13 + args.K4*x12 + x10
#         x15 = 1.0/x14
#         x16 = 4*x7
#         x17 = 0.5*args.CEC_tot
#         x18 = pow(Ca, -3.0/2.0)
#         x19 = 1.0/args.conv_Al
#         x20 = sqrt(1.0/args.K_Ca_Al)
#         x21 = pow(f_Ca, 3.0/2.0)*x19*x20
#         x22 = (3.0/2.0)*Al
#         x23 = sqrt(f_Ca)
#         x24 = x18*x23
#         x25 = 1.0/args.K_Ca_Mg
#         x26 = x25/Ca
#         x27 = sqrt(1.0/args.K_Ca_Na)
#         x28 = (1.0/2.0)*Na*x27
#         x29 = pow(Ca, -1.0/2.0)
#         x30 = x23*x29
#         x31 = x29/x23
#         x32 = (1.0/2.0)*x24
#         x33 = sqrt(1.0/args.K_Ca_K)
#         x34 = K*x33
#         x35 = (1.0/2.0)*x31
#         x36 = sqrt(1.0/args.K_Ca_H)
#         x37 = H*x36
#
#         # --- Jacobian Matrix Assignment ---
#         fjac[0] = -x1
#         fjac[48] = -1
#         fjac[17] = -x0*(args.s*(x2 + x6 + 1.0) + (1.0 - args.s)/args.k_H)
#         fjac[33] = -CO2_w*x1*(-x5 - 2*x8)
#         fjac[2] = -1
#         fjac[18] = x2 + 2.0*x6
#         fjac[34] = -CO2_w*x5 - 4.0*CO2_w*x8 - args.k_w*x4 - 1
#         fjac[51] = 1
#         fjac[179] = x9
#         fjac[195] = x9
#         fjac[211] = x9
#         fjac[243] = x9
#         fjac[68] = x1
#         fjac[164] = 0.33333333333333331*args.CEC_tot*args.conv_Al
#         fjac[37] = -Al_w*x10*(-2*H*x11 - x12 - 3*x13 - x16)/(x14*x14) - Al_w*x15*x16
#         fjac[69] = -x10*x15
#         fjac[85] = 1
#         fjac[102] = x1
#         fjac[182] = x17
#         fjac[119] = x1
#         fjac[247] = x17
#         fjac[136] = x1
#         fjac[200] = args.CEC_tot
#         fjac[153] = x1
#         fjac[217] = args.CEC_tot
#         fjac[90] = -x18*x21
#         fjac[122] = x21*x22/pow(Ca, 5.0/2.0)
#         fjac[170] = 1
#         fjac[250] = -x19*x20*x22*x24
#         fjac[107] = -f_Ca*x26
#         fjac[123] = Mg*f_Ca*x25/(Ca*Ca)
#         fjac[187] = 1
#         fjac[251] = -Mg*x26
#         fjac[124] = x24*x28
#         fjac[140] = -x27*x30
#         fjac[204] = 1
#         fjac[252] = -x28*x31
#         fjac[125] = x32*x34
#         fjac[157] = -x30*x33
#         fjac[221] = 1
#         fjac[253] = -x34*x35
#         fjac[46] = -x30*x36
#         fjac[126] = x32*x37
#         fjac[238] = 1
#         fjac[254] = -x35*x37
#         fjac[175] = -1
#         fjac[191] = -1
#         fjac[207] = -1
#         fjac[223] = -1
#         fjac[239] = -1
#         fjac[255] = -1
#
#     return 0

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


# # Jacobian based solver in log space
# def solve_biogeochem_eq_jac_log(
#     double[:] x0_arr,
#     double Alk_tot_in, double n_in, double Zr_in, double s_in, double IC_tot_in,
#     double k1_in, double k2_in, double k_H_in, double k_w_in, double CEC_tot_in,
#     double conv_Al_in, double Al_tot_in, double K1_in, double K2_in, double K3_in,
#     double K4_in, double Mg_tot_in, double Ca_tot_in, double Na_tot_in, double K_tot_in,
#     double K_Ca_Al_in, double K_Ca_Mg_in, double K_Ca_Na_in, double K_Ca_K_in, double K_Ca_H_in
# ):
#     cdef int n_vars = 16
#     cdef int lwa = 232
#     cdef double tol = 1e-8
#     cdef int status
#
#     cdef double x[16]
#     cdef double fvec[16]
#     cdef double fjac[256]
#     cdef double wa[232]
#
#     cdef int i
#
#     # --- TRANSFORM INPUT TO LOG SPACE (Except Alk) ---
#     x[0] = x0_arr[0]  # Alk stays linear
#     for i in range(1, n_vars):
#         x[i] = log(x0_arr[i])
#
#     cdef EquationArgs args
#     args.Alk_tot = Alk_tot_in
#     args.n = n_in
#     args.Zr = Zr_in
#     args.s = s_in
#     args.IC_tot = IC_tot_in
#     args.k1 = k1_in
#     args.k2 = k2_in
#     args.k_H = k_H_in
#     args.k_w = k_w_in
#     args.CEC_tot = CEC_tot_in
#     args.conv_Al = conv_Al_in
#     args.Al_tot = Al_tot_in
#     args.K1 = K1_in
#     args.K2 = K2_in
#     args.K3 = K3_in
#     args.K4 = K4_in
#     args.Mg_tot = Mg_tot_in
#     args.Ca_tot = Ca_tot_in
#     args.Na_tot = Na_tot_in
#     args.K_tot = K_tot_in
#     args.K_Ca_Al = K_Ca_Al_in
#     args.K_Ca_Mg = K_Ca_Mg_in
#     args.K_Ca_Na = K_Ca_Na_in
#     args.K_Ca_K = K_Ca_K_in
#     args.K_Ca_H = K_Ca_H_in
#     args.func_calls = 0
#     args.jac_calls = 0
#
#     with nogil:
#         status = hybrj1(biogeochem_equations_c_jac_log, <void*>&args, n_vars, x, fvec, fjac, n_vars, tol, wa, lwa)
#
#     if status != 1:
#         print(f"JAC_LOG: Status: {status} | F-Evals: {args.func_calls} | J-Evals: {args.jac_calls}")
#
#     cdef double[:] result = np.zeros(n_vars, dtype=np.float64)
#     cdef double[:] residuals = np.zeros(n_vars, dtype=np.float64)
#
#     # --- TRANSFORM RESULT BACK TO LINEAR SPACE ---
#     result[0] = x[0] # Alk stays linear
#     for i in range(1, n_vars):
#         result[i] = exp(x[i])
#         residuals[i] = fvec[i]  # Capture the final error of equation i
#
#     return result, residuals, status
#
#
# # Jacobian based solver
# def solve_biogeochem_eq_jac(
#     double[:] x0_arr,
#     double Alk_tot_in, double n_in, double Zr_in, double s_in, double IC_tot_in,
#     double k1_in, double k2_in, double k_H_in, double k_w_in, double CEC_tot_in,
#     double conv_Al_in, double Al_tot_in, double K1_in, double K2_in, double K3_in,
#     double K4_in, double Mg_tot_in, double Ca_tot_in, double Na_tot_in, double K_tot_in,
#     double K_Ca_Al_in, double K_Ca_Mg_in, double K_Ca_Na_in, double K_Ca_K_in, double K_Ca_H_in
# ):
#     cdef int n_vars = 16
#     cdef int lwa = 232
#     cdef double tol = 1e-8
#     cdef int status
#
#     cdef double x[16]
#     cdef double fvec[16]
#     cdef double fjac[256]
#     cdef double wa[232]
#
#     cdef int i
#     for i in range(n_vars):
#         x[i] = x0_arr[i]
#
#     cdef EquationArgs args
#     args.Alk_tot = Alk_tot_in
#     args.n = n_in
#     args.Zr = Zr_in
#     args.s = s_in
#     args.IC_tot = IC_tot_in
#     args.k1 = k1_in
#     args.k2 = k2_in
#     args.k_H = k_H_in
#     args.k_w = k_w_in
#     args.CEC_tot = CEC_tot_in
#     args.conv_Al = conv_Al_in
#     args.Al_tot = Al_tot_in
#     args.K1 = K1_in
#     args.K2 = K2_in
#     args.K3 = K3_in
#     args.K4 = K4_in
#     args.Mg_tot = Mg_tot_in
#     args.Ca_tot = Ca_tot_in
#     args.Na_tot = Na_tot_in
#     args.K_tot = K_tot_in
#     args.K_Ca_Al = K_Ca_Al_in
#     args.K_Ca_Mg = K_Ca_Mg_in
#     args.K_Ca_Na = K_Ca_Na_in
#     args.K_Ca_K = K_Ca_K_in
#     args.K_Ca_H = K_Ca_H_in
#     args.func_calls = 0
#     args.jac_calls = 0
#
#     with nogil:
#         status = hybrj1(biogeochem_equations_c_jac, <void*>&args, n_vars, x, fvec, fjac, n_vars, tol, wa, lwa)
#
#     if status != 1:
#         print(f"JAC Status: {status} | F-Evals: {args.func_calls} | J-Evals: {args.jac_calls}")
#
#     cdef double[:] result = np.zeros(n_vars, dtype=np.float64)
#     cdef double[:] residuals = np.zeros(n_vars, dtype=np.float64)
#
#     for i in range(0, n_vars):
#         result[i] = x[i]
#         residuals[i] = fvec[i]  # Capture the final error of equation i
#
#     return result, residuals, status


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
