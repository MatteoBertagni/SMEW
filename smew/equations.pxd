# Native declarations for the arithmetic maintained in equations.py.
# Scalar arguments are double precision; vector state/output use contiguous
# double memoryviews. Callers must validate vector lengths before entry.

import cython

cdef double water_residual(double H_rain, double Alk_rain, double k1,
                           double k2, double CO2_w_rain, double k_w) except? -1.0 nogil
cdef double h_residual(double H0, double k1, double k2,
                       double CO2_w0, double k_w, double Alk0) except? -1.0 nogil
cdef double cec_calcium_residual(double f_Ca, double Al, double conv_Al,
                                 double Ca, double Mg, double Na, double K,
                                 double H, double K_Ca_Al, double K_Ca_Mg,
                                 double K_Ca_Na, double K_Ca_K,
                                 double K_Ca_H) except? -1.0 nogil

@cython.locals(Alk=cython.double, CO2_w=cython.double, H=cython.double,
               R_alk=cython.double, Al_w=cython.double, Al=cython.double,
               Mg=cython.double, Ca=cython.double, Na=cython.double,
               K=cython.double, f_Al=cython.double, f_Mg=cython.double,
               f_Na=cython.double, f_K=cython.double, f_H=cython.double,
               f_Ca=cython.double, nZrs1000=cython.double)
cdef int biogeochem_residual(double[::1] p, double Alk_tot, double n,
        double Zr, double s, double IC_tot, double k1, double k2,
        double k_H, double k_w, double CEC_tot, double conv_Al,
        double Al_tot, double K1, double K2, double K3, double K4,
        double Mg_tot, double Ca_tot, double Na_tot, double K_tot,
        double K_Ca_Al, double K_Ca_Mg, double K_Ca_Na, double K_Ca_K,
        double K_Ca_H, double[::1] out) except -1 nogil

@cython.locals(Al_w=cython.double, Al=cython.double,
               Al_tot=cython.double, Mg=cython.double, Ca=cython.double,
               Na=cython.double, K=cython.double, f_Mg=cython.double,
               f_Na=cython.double, f_K=cython.double, f_Ca=cython.double,
               f_Al=cython.double, f_H=cython.double)
cdef int total_to_cec_residual(double[::1] p, double H, double n,
        double Zr, double s, double CEC_tot, double conv_Al,
        double K1, double K2, double K3, double K4, double Ca_tot,
        double Mg_tot, double K_tot, double Na_tot, double f_acid,
        double K_Ca_Al, double K_Ca_Mg, double K_Ca_Na,
        double K_Ca_K, double K_Ca_H, double[::1] out) except -1 nogil

@cython.locals(Al_tot=cython.double, CaCO3=cython.double,
               f_Al=cython.double, f_H=cython.double, f_Ca=cython.double)
cdef int kelland_residual(double[::1] p, double Al_w, double Al,
        double H, double Ca, double Ca_tot, double f_Mg, double f_K,
        double f_Na, double n, double Zr, double s, double CEC_tot,
        double conv_Al, double K_Ca_Al, double K_Ca_H,
        double[::1] out) except -1 nogil
