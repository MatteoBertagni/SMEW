# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""
import ctypes
import os
import glob

from numba import njit
import numpy as np

from smew.constants import K_GT_CEC as smew_K_GT_CEC
from smew.constants import K_Al as smew_K_Al

# ------------------------------------------------------------------------------
# Locate and setup the compiled Cython shared object
# ------------------------------------------------------------------------------
search_pattern = os.path.join(os.path.dirname(__file__), "ic_fast.*.so")
so_files = glob.glob(search_pattern)

if not so_files:
    raise ImportError(
        "Could not find the compiled Cython solver for ic.py. "
        "Ensure 'ic_fast.pyx' has been properly compiled into a .so or .pyd file."
    )

solver_lib = ctypes.CDLL(so_files[0])

# 1. 1D Solver Bridge (eqf_Ca)
solve_eqf_Ca_c = solver_lib.solve_eqf_Ca_fd
solve_eqf_Ca_c.restype = ctypes.c_int
solve_eqf_Ca_c.argtypes = [
    ctypes.c_void_p,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double
]

# 2. 13D Solver Bridge (total_to_f_CEC_and_conc)
solve_total_eq_c = solver_lib.solve_total_eq_fd
solve_total_eq_c.restype = ctypes.c_int
solve_total_eq_c.argtypes = [
    ctypes.c_void_p,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double
]

# 3. 5D Solver Bridge (Kelland)
solve_kelland_eq_c = solver_lib.solve_kelland_eq_fd
solve_kelland_eq_c.restype = ctypes.c_int
solve_kelland_eq_c.argtypes = [
    ctypes.c_void_p,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double
]


#------------------------------------------------------------------------------
 # conc to CEC fractions

@njit(nogil=True)
def conc_to_f_CEC(conc_in,pH_in,soil,conv_mol,conv_Al):

    #constants
    K_CEC = smew_K_GT_CEC(soil,conv_mol) #CEC Gaines-Thomas
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H = K_CEC
    K1, K2, K3, K4 = smew_K_Al(conv_mol) #Al speciation

    # cations (mol-conv/l)
    Ca, Mg, K, Na, Al_w = conc_in
    H = 10**(-pH_in)*conv_mol

    # aluminium speciation
    Al=(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w #mol-conv/l

    #CEC fractions
    x0 = np.array([0.2], dtype=np.float64)
    solve_eqf_Ca_c(x0.ctypes.data, Al, conv_Al, K_Ca_Al, Ca, Mg, K_Ca_Mg, Na, K_Ca_Na, K, K_Ca_K, H, K_Ca_H)
    f_Ca = x0[0]

    f_Al = (Al/conv_Al)*((f_Ca**3/(K_Ca_Al*Ca**3))**(1/2))
    f_Mg = Mg*(f_Ca/(K_Ca_Mg*Ca))
    f_Na = Na*((f_Ca/(K_Ca_Na*Ca))**(1/2))
    f_K = K*((f_Ca/(K_Ca_K*Ca))**(1/2))
    f_H = H*((f_Ca/(K_Ca_H*Ca))**(1/2))

    f_CEC_in = (f_Ca, f_Mg, f_K, f_Na, f_Al, f_H)

    return f_CEC_in, K_CEC

#------------------------------------------------------------------------------
 # CEC fractions to conc

@njit(nogil=True)
def f_CEC_to_conc(f_CEC_in, pH_in, soil, conv_mol,conv_Al):

    #pH
    H = 10**(-pH_in)*conv_mol

    #constants
    K_CEC = smew_K_GT_CEC(soil,conv_mol) #CEC Gaines-Thomas
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H = K_CEC
    K1, K2, K3, K4 = smew_K_Al(conv_mol) #Al speciation

    #f_CEC [-]
    f_Ca, f_Mg, f_K, f_Na, f_Al, f_H = f_CEC_in

    # FIX: Prevent Division by Zero ---
    eps = 1e-15
    f_H_safe = max(f_H, eps)
    f_Ca_safe = max(f_Ca, eps)

    #estimates of concentrations (division by zero error)
    # Ca = (f_Ca/K_Ca_H)*(H/f_H)**2
    # Mg = (Ca/f_Ca)*K_Ca_Mg*f_Mg
    # K = ((Ca/f_Ca)*K_Ca_K)**(1/2)*f_K
    # Na = ((Ca/f_Ca)*K_Ca_Na)**(1/2)*f_Na
    # Al = ((Ca/f_Ca)**3*K_Ca_Al)**(1/2)*f_Al*conv_Al
    # safe estimates of concentrations
    Ca = (f_Ca/K_Ca_H)*(H/f_H_safe)**2
    Mg = (Ca/f_Ca_safe)*K_Ca_Mg*f_Mg
    K = ((Ca/f_Ca_safe)*K_Ca_K)**(1/2)*f_K
    Na = ((Ca/f_Ca_safe)*K_Ca_Na)**(1/2)*f_Na
    Al = ((Ca/f_Ca_safe)**3*K_Ca_Al)**(1/2)*f_Al*conv_Al

    Al_w=Al/(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))

    conc_in = (Ca, Mg, K, Na, Al_w)
    K_CEC = (K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H)

    return conc_in, K_CEC

#------------------------------------------------------------------------------
 # Input: Total (Ca, Mg, K, Na) and CEC base saturation (or acid saturation, f_H+f_Al)

@njit(nogil=True)
def total_to_f_CEC_and_conc(total_in, pH_in, f_acid, s, soil, n,Zr,CEC_tot,conv_mol,conv_Al):

    #constants
    K_CEC = smew_K_GT_CEC(soil, conv_mol) #CEC Gaines-Thomas
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H = K_CEC
    K1, K2, K3, K4 = smew_K_Al(conv_mol) #Al speciation

    # total (mol-conv/m2)
    Ca_tot, Mg_tot, K_tot, Na_tot = total_in
    H = 10**(-pH_in)*conv_mol

    #initial guess
    f_Ca0 = 0.8*(1 - f_acid)
    f_Mg0 = 0.1*(1 - f_acid)
    f_K0 = 0.05*(1 - f_acid)
    f_Na0 = 0.05*(1 - f_acid)
    f_H0 = f_acid*0.3
    f_Al0 = f_acid*0.7
    Mg0 = (Mg_tot-f_Mg0/2*CEC_tot)/(n*Zr*s[0]*1000)
    Ca0 = (Ca_tot-f_Ca0/2*CEC_tot)/(n*Zr*s[0]*1000)
    Na0 = (Na_tot-f_Na0*CEC_tot)/(n*Zr*s[0]*1000)
    K0 = (K_tot-f_K0*CEC_tot)/(n*Zr*s[0]*1000)
    Al0 = (f_Al0*conv_Al)*((K_Ca_Al*Ca0**3)/f_Ca0**3)**(1/2)
    Al_w0 = Al0/(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))
    Al_tot0 = (f_Al0/3)*CEC_tot*conv_Al+Al_w0*n*Zr*s[0]*1000

    x0 = np.array((Al_w0, Al0, Al_tot0, Mg0, Ca0, Na0, K0, f_Mg0, f_Na0, f_K0, f_Ca0, f_Al0, f_H0), dtype=np.float64)

    solve_total_eq_c(x0.ctypes.data, n, Zr, s[0], CEC_tot, conv_Al, K1, K2, K3, K4,
                          Mg_tot, Ca_tot, Na_tot, K_tot, H, K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H, f_acid)

    Al_w, Al, Al_tot, Mg, Ca, Na, K, f_Mg, f_Na, f_K, f_Ca, f_Al, f_H = x0

    #K_Ca_Al = (Al/conv_Al/f_Al)**2*(f_Ca/Ca)**3
    #K_Ca_H = (f_Ca/Ca)*(H/f_H)**2
    K_CEC = (K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H)
    conc_in = (Ca, Mg, K, Na, Al_w)
    f_CEC_in = (f_Ca, f_Mg, f_K, f_Na, f_Al, f_H)

    return conc_in, f_CEC_in, K_CEC

#------------------------------------------------------------------------------
 # Calibration of K constants on coupled f_CEC and conc measurements

@njit(nogil=True)
def f_CEC_and_conc_to_K(f_CEC_in, conc_in, pH_in, soil, conv_mol,conv_Al):

    #pH
    H = 10**(-pH_in)*conv_mol

    #f_CEC [-]
    f_Ca, f_Mg, f_K, f_Na, f_Al, f_H = f_CEC_in

    #conc [-]
    Ca, Mg, K, Na, Al_w = conc_in

    # Al spec
    K1, K2, K3, K4 = smew_K_Al(conv_mol)
    Al=(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w

    #estimates of K constants
    K_Ca_Mg = (f_Ca/Ca)*(Mg/f_Mg)
    K_Ca_K = (f_Ca/Ca)*(K/f_K)**2
    K_Ca_Na = (f_Ca/Ca)*(Na/f_Na)**2
    K_Ca_Al = (f_Ca/Ca)**3*(f_Al/(Al/conv_Al))**2
    K_Ca_H = (f_Ca/Ca)*(H/f_H)**2

    K_CEC = (K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H)

    return K_CEC

#------------------------------------------------------------------------------
 # Input: Total (Ca, Mg, K, Na) and Al_w

@njit(nogil=True)
def Kelland(total_in, pH_in, conc_in, s, soil, n,Zr,CEC_tot,conv_mol,conv_Al):

    #constants
    K_CEC = smew_K_GT_CEC(soil, conv_mol) #CEC Gaines-Thomas
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H = K_CEC
    K1, K2, K3, K4 = smew_K_Al(conv_mol) #Al speciation

    # known (mol-conv/m2)
    Ca_tot, Mg_tot, K_tot, Na_tot = total_in
    Ca, Mg, K, Na, Al_w = conc_in
    H = 10**(-pH_in)*conv_mol
    Al=(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w

    #CEC
    f_Mg = (Mg_tot - Mg*n*Zr*s[0]*1000)*2/CEC_tot
    f_K = (K_tot - K*n*Zr*s[0]*1000)*2/CEC_tot
    f_Na = (Na_tot - Na*n*Zr*s[0]*1000)/CEC_tot

    #initial guess
    f_H0 = 1e-3
    f_Al0 = 1e-3
    f_Ca0 = 0.8*(1 - f_H0 + f_Al0)
    Al_tot0 = (f_Al0/3)*CEC_tot*conv_Al+Al_w*n*Zr*s[0]*1000
    CaCO30 = Ca_tot-Ca*n*Zr*s[0]*1000+f_Ca0/2*CEC_tot

    x0 = np.array([Al_tot0, CaCO30, f_Al0, f_H0, f_Ca0], dtype=np.float64)

    solve_kelland_eq_c(x0.ctypes.data, Al_w, n, Zr, s[0], CEC_tot, conv_Al, Ca, Ca_tot, Al, K_Ca_Al, H, K_Ca_H, f_Mg, f_Na, f_K)
    Al_tot, CaCO3,  f_Al, f_H, f_Ca = x0

    #estimating soil-dependent K_CEC
    K_Ca_Mg = (f_Ca/Ca)*(Mg/f_Mg)
    K_Ca_K = (f_Ca/Ca)*(K/f_K)**2
    K_Ca_Na = (f_Ca/Ca)*(Na/f_Na)**2
    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]

    conc_in = [Ca, Mg, K, Na, Al_w]
    f_CEC_in = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H]

    return(conc_in, f_CEC_in, K_CEC, CaCO3)

#------------------------------------------------------------------------------
 # Amann et al., fractions and Mg conc

@njit(nogil=True)
def Amann(f_CEC_in, pH_in, Mg_in, soil, conv_mol,conv_Al):

    #pH
    H = 10**(-pH_in)*conv_mol

    #constants
    K_CEC = smew_K_GT_CEC(soil,conv_mol) #CEC Gaines-Thomas
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H = K_CEC
    K1, K2, K3, K4 = smew_K_Al(conv_mol) #Al speciation

    #f_CEC [-]
    f_Ca, f_Mg, f_K, f_Na, f_Al, f_H = f_CEC_in

    #conc from measurements
    Mg = Mg_in

    #estimates of concentrations
    Ca = (Mg/f_Mg)/K_Ca_Mg*f_Ca
    K = ((Ca/f_Ca)*K_Ca_K)**(1/2)*f_K
    Na = ((Ca/f_Ca)*K_Ca_Na)**(1/2)*f_Na
    Al = ((Ca/f_Ca)**3*K_Ca_Al)**(1/2)*f_Al*conv_Al

    Al_w=Al/(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))

    #estimates of K_CEC
    K_Ca_H = (f_Ca/Ca)*(H/f_H)**2

    conc_in = (Ca, Mg, K, Na, Al_w)
    K_CEC = (K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H)

    return conc_in, K_CEC
