# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-

"""Scientific equations for SMEW's nonlinear solves.

To change an equation, edit its ``*_residual`` function below. Functions with
one unknown return a value; larger systems fill ``out``. The matching
``*_equations`` functions pass the results to SciPy's ``fsolve``. Python uses
this file directly; compiled runs require a Cython rebuild after edits.
"""

import numpy as np


def water_residual(H_rain, Alk_rain, k1, k2, CO2_w_rain, k_w):
    return Alk_rain-(k1*CO2_w_rain/H_rain+2*k1*k2*CO2_w_rain/(H_rain**2)-H_rain+k_w/H_rain)


def water_equations(p, Alk_rain, k1, k2, CO2_w_rain, k_w):
    """Rainwater alkalinity residual; ``p`` contains H"""
    return water_residual(p[0], Alk_rain, k1, k2, CO2_w_rain, k_w)


def h_residual(H0, k1, k2, CO2_w0, k_w, Alk0):
    return (k1*CO2_w0/H0+2*k1*k2*CO2_w0/(H0**2)-H0+k_w/H0)-Alk0


def h_equations(p, k1, k2, CO2_w0, k_w, Alk0):
    """Alternative hydrogen guess residual; ``p`` contains H"""
    return h_residual(p[0], k1, k2, CO2_w0, k_w, Alk0)


def biogeochem_residual(
        p, Alk_tot, n, Zr, s, IC_tot, k1, k2, k_H, k_w, CEC_tot, conv_Al, Al_tot, K1, K2, K3, K4, Mg_tot, Ca_tot,
        Na_tot, K_tot, K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H, out
):
    """Main equilibrium residuals for aqueous/adsorbed species.

    State and residual order: Alk, CO2_w, H, R_alk, Al_w, Al, Mg, Ca,
    Na, K, f_Al, f_Mg, f_Na, f_K, f_H, f_Ca.
    """
    Alk = p[0]
    CO2_w = p[1]
    H = p[2]
    R_alk = p[3]
    Al_w = p[4]
    Al = p[5]
    Mg = p[6]
    Ca = p[7]
    Na = p[8]
    K = p[9]
    f_Al = p[10]
    f_Mg = p[11]
    f_Na = p[12]
    f_K = p[13]
    f_H = p[14]
    f_Ca = p[15]

    # Precompute
    nZrs1000 = n * Zr * s * 1000

    out[0] = (Alk_tot-R_alk)-Alk*nZrs1000
    out[1] = IC_tot-(CO2_w*(1+k1/H+k2*k1/(H**2))*s+(CO2_w/k_H)*(1-s))*(n*Zr*1000)
    out[2] = (k1*CO2_w/H+2*k1*k2*CO2_w/(H**2)-H+k_w/H)-Alk
    out[3] = R_alk-(f_Mg+f_Ca+f_Na+f_K)*CEC_tot
    out[4] = Al_w*nZrs1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot
    out[5] = Al-(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w
    out[6] = Mg*nZrs1000+f_Mg/2*CEC_tot-Mg_tot
    out[7] = Ca*nZrs1000+f_Ca/2*CEC_tot-Ca_tot
    out[8] = Na*nZrs1000+f_Na*CEC_tot-Na_tot
    out[9] = K*nZrs1000+f_K*CEC_tot-K_tot
    out[10] = f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2)
    out[11] = f_Mg - Mg*(f_Ca/(K_Ca_Mg*Ca))
    out[12] = f_Na - Na*(f_Ca/(K_Ca_Na*Ca))**(1/2)
    out[13] = f_K - K*(f_Ca/(K_Ca_K*Ca))**(1/2)
    out[14] = f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2)
    out[15] = 1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H)
    return 0


def biogeochem_equations(
        p, Alk_tot, n, Zr, s, IC_tot, k1, k2, k_H, k_w, CEC_tot, conv_Al, Al_tot, K1, K2, K3, K4, Mg_tot, Ca_tot,
        Na_tot, K_tot, K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H
):
    out = np.empty(16)
    biogeochem_residual(
        p, Alk_tot, n, Zr, s, IC_tot, k1, k2, k_H, k_w, CEC_tot, conv_Al, Al_tot, K1, K2, K3, K4,
        Mg_tot, Ca_tot, Na_tot, K_tot, K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H, out
    )
    return tuple(out)


def cec_calcium_residual(f_Ca, Al, conv_Al, Ca, Mg, Na, K, H,
                         K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H):
    return 1-(f_Ca+(Al/conv_Al)*((f_Ca**3/(K_Ca_Al*Ca**3))**(1/2))
              +Mg*(f_Ca/(K_Ca_Mg*Ca))
              +Na*((f_Ca/(K_Ca_Na*Ca))**(1/2))
              +K*((f_Ca/(K_Ca_K*Ca))**(1/2))
              +H*((f_Ca/(K_Ca_H*Ca))**(1/2)))


def cec_calcium_equation(p, Al, conv_Al, Ca, Mg, Na, K, H,
                         K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H):
    """CEC saturation residual for the calcium fraction ``p[0]``."""
    return cec_calcium_residual(p[0], Al, conv_Al, Ca, Mg, Na, K, H,
                                K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H)


def total_to_cec_residual(p, H, n, Zr, s, CEC_tot, conv_Al,
                          K1, K2, K3, K4, Ca_tot, Mg_tot, K_tot,
                          Na_tot, f_acid, K_Ca_Al, K_Ca_Mg, K_Ca_Na,
                          K_Ca_K, K_Ca_H, out):
    """Residuals for the 13 initial concentrations and CEC fractions."""
    Al_w = p[0]
    Al = p[1]
    Al_tot = p[2]
    Mg = p[3]
    Ca = p[4]
    Na = p[5]
    K = p[6]
    f_Mg = p[7]
    f_Na = p[8]
    f_K = p[9]
    f_Ca = p[10]
    f_Al = p[11]
    f_H = p[12]
    out[0] = Al_w*n*Zr*s*1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot
    out[1] = Al-(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w
    out[2] = Mg*n*Zr*s*1000+f_Mg/2*CEC_tot-Mg_tot
    out[3] = Ca*n*Zr*s*1000+f_Ca/2*CEC_tot-Ca_tot
    out[4] = Na*n*Zr*s*1000+f_Na*CEC_tot-Na_tot
    out[5] = K*n*Zr*s*1000+f_K*CEC_tot-K_tot
    out[6] = f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2)
    out[7] = f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2)
    out[8] = f_H + f_Al - f_acid
    out[9] = f_Mg - Mg*(f_Ca/(K_Ca_Mg*Ca))
    out[10] = f_Na - Na*(f_Ca/(K_Ca_Na*Ca))**(1/2)
    out[11] = f_K - K*(f_Ca/(K_Ca_K*Ca))**(1/2)
    out[12] = 1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H)
    return 0


def total_to_cec_equations(p, H, n, Zr, s, CEC_tot, conv_Al,
                           K1, K2, K3, K4, Ca_tot, Mg_tot, K_tot,
                           Na_tot, f_acid, K_Ca_Al, K_Ca_Mg, K_Ca_Na,
                           K_Ca_K, K_Ca_H):
    out = np.empty(13)
    total_to_cec_residual(p, H, n, Zr, s, CEC_tot, conv_Al,
                          K1, K2, K3, K4, Ca_tot, Mg_tot, K_tot,
                          Na_tot, f_acid, K_Ca_Al, K_Ca_Mg, K_Ca_Na,
                          K_Ca_K, K_Ca_H, out)
    return tuple(out)


def kelland_residual(p, Al_w, Al, H, Ca, Ca_tot,
                      f_Mg, f_K, f_Na, n, Zr, s, CEC_tot, conv_Al,
                      K_Ca_Al, K_Ca_H, out):
    """Residuals for the five Kelland initial-condition unknowns."""
    Al_tot = p[0]
    CaCO3 = p[1]
    f_Al = p[2]
    f_H = p[3]
    f_Ca = p[4]
    out[0] = Al_w*n*Zr*s*1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot
    out[1] = Ca*n*Zr*s*1000+f_Ca/2*CEC_tot+CaCO3-Ca_tot
    out[2] = f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2)
    out[3] = f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2)
    out[4] = 1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H)
    return 0


def kelland_equations(p, Al_w, Al, H, Ca, Ca_tot,
                      f_Mg, f_K, f_Na, n, Zr, s, CEC_tot, conv_Al,
                      K_Ca_Al, K_Ca_H):
    out = np.empty(5)
    kelland_residual(p, Al_w, Al, H, Ca, Ca_tot,
                     f_Mg, f_K, f_Na, n, Zr, s, CEC_tot, conv_Al,
                     K_Ca_Al, K_Ca_H, out)
    return tuple(out)
