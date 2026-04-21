# cython: language_level=3
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport pow, sqrt, fabs

# Ensure the NumPy C-API is initialized
np.import_array()

# Define the c-functions from cminpack
cdef extern from "cminpack.h":
    # The callback function signature required by cminpack
    ctypedef int (*cminpack_func_nn)(void *p, int n, const double *x, double *fvec, int iflag) nogil

    # The hybrd1 solver signature
    int hybrd1(cminpack_func_nn fcn, void *p, int n, double *x, double *fvec, double tol, double *wa, int lwa) nogil


# cminpack allows to pass a void pointer (void *p) to the equations.
# We use a struct so we don't have to use global variables.
cdef struct EquationArgs:
    double Alk_tot
    double n
    double Zr
    double s
    double IC_tot
    double k1
    double k2
    double k_H
    double k_w
    double CEC_tot
    double conv_Al
    double Al_tot
    double K1
    double K2
    double K3
    double K4
    double Mg_tot
    double Ca_tot
    double Na_tot
    double K_tot
    double K_Ca_Al
    double K_Ca_Mg
    double K_Ca_Na
    double K_Ca_K
    double K_Ca_H


cdef int biogeochem_equations_c(void *p, int n, const double *x, double *fvec, int iflag) noexcept nogil:
    # Cast the void pointer back to our struct to access constants [cite: 5]
    cdef EquationArgs* args = <EquationArgs*>p

    # Unpack the state variables from the solver's x array [cite: 5]
    cdef double Alk = x[0]
    cdef double CO2_w = x[1]
    cdef double H = x[2]
    cdef double R_alk = x[3]
    cdef double Al_w = x[4]
    cdef double Al = x[5]
    cdef double Mg = x[6]
    cdef double Ca = x[7]
    cdef double Na = x[8]
    cdef double K = x[9]
    cdef double f_Al = x[10]
    cdef double f_Mg = x[11]
    cdef double f_Na = x[12]
    cdef double f_K = x[13]
    cdef double f_H = x[14]
    cdef double f_Ca = x[15]

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
    fvec[10] = f_Al - (Al / args.conv_Al) * sqrt( (f_Ca * f_Ca * f_Ca) / (args.K_Ca_Al * Ca * Ca * Ca) )
    # Equation 12: Magnesium exchange
    fvec[11] = f_Mg - Mg * (f_Ca / (args.K_Ca_Mg * Ca))
    # Equation 13: Sodium exchange
    fvec[12] = f_Na - Na * sqrt( f_Ca / (args.K_Ca_Na * Ca) )
    # Equation 14: Potassium exchange
    fvec[13] = f_K - K * sqrt( f_Ca / (args.K_Ca_K * Ca) )
    # Equation 15: Hydrogen exchange
    fvec[14] = f_H - H * sqrt( f_Ca / (args.K_Ca_H * Ca) )
    # Equation 16: Sum of exchange fractions must be 1
    fvec[15] = 1.0 - (f_Ca + f_Al + f_Mg + f_Na + f_K + f_H)

    return 0

# The python exposed wrapper
def solve_biogeochem_eq(
    double[:] x0_arr,
    double Alk_tot_in, double n_in, double Zr_in, double s_in, double IC_tot_in,
    double k1_in, double k2_in, double k_H_in, double k_w_in, double CEC_tot_in,
    double conv_Al_in, double Al_tot_in, double K1_in, double K2_in, double K3_in,
    double K4_in, double Mg_tot_in, double Ca_tot_in, double Na_tot_in, double K_tot_in,
    double K_Ca_Al_in, double K_Ca_Mg_in, double K_Ca_Na_in, double K_Ca_K_in, double K_Ca_H_in
):
    cdef int n_vars = 16

    # MINPACK requires a very specific workspace array size for hybrd1.
    # The formula from the MINPACK documentation is: LWA >= N * (3*N + 13) / 2
    cdef int lwa = (n_vars * (3 * n_vars + 13)) // 2

    cdef double tol = 1e-12
    cdef int status

    # --- ALLOCATE C-ARRAYS ---
    # We cast the result of malloc to a double pointer (<double *>)
    # sizeof(double) ensures we reserve exactly 8 bytes per variable
    cdef double *x = <double *>malloc(n_vars * sizeof(double))
    cdef double *fvec = <double *>malloc(n_vars * sizeof(double))
    cdef double *wa = <double *>malloc(lwa * sizeof(double))

    # Safety check: If your system runs out of RAM, malloc returns a Null pointer.
    if not x or not fvec or not wa:
        if x: free(x)
        if fvec: free(fvec)
        if wa: free(wa)
        raise MemoryError("Failed to allocate C-arrays for cminpack.")

    # --- COPY PYTHON DATA TO C-MEMORY ---
    cdef int i
    for i in range(n_vars):
        x[i] = x0_arr[i]

    # Pack the Python arguments into the C-Struct
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

    # Call CMINPACK with the packed struct
    with nogil:
        status = hybrd1(biogeochem_equations_c, <void*>&args, n_vars, x, fvec, tol, wa, lwa)

    # Copy results back to a new numpy array
    cdef double[:] result = np.zeros(n_vars, dtype=np.float64)
    for i in range(n_vars):
        result[i] = x[i]

    # Free memory to prevent leaks
    free(x)
    free(fvec)
    free(wa)

    return result, status