# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-

"""Pure Python residuals for SMEW's nonlinear systems.

Each function reads a solver state and parameters and returns residuals in solver
order. Solver settings, initial guesses, and history updates belong to callers.
"""


def water_equations(p, Alk_rain, k1, k2, CO2_w_rain, k_w):
    """Rainwater alkalinity residual; ``p`` contains H in mol/l."""
    H_rain = p[0]
    return Alk_rain-(k1*CO2_w_rain/H_rain+2*k1*k2*CO2_w_rain/(H_rain**2)-H_rain+k_w/H_rain)


def h_equations(p, k1, k2, CO2_w0, k_w, Alk0):
    """Alternative hydrogen guess residual; ``p`` contains H in mol/l."""
    H0 = p[0]
    return (k1*CO2_w0/H0+2*k1*k2*CO2_w0/(H0**2)-H0+k_w/H0)-Alk0


def biogeochem_equations(
        p, Alk_tot, n, Zr, s, IC_tot, k1, k2, k_H, k_w, CEC_tot, conv_Al, Al_tot, K1, K2, K3, K4, Mg_tot, Ca_tot,
        Na_tot, K_tot, K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H
):
    """Main equilibrium residuals for aqueous/adsorbed species.

    State and residual order: Alk, CO2_w, H, R_alk, Al_w, Al, Mg, Ca,
    Na, K, f_Al, f_Mg, f_Na, f_K, f_H, f_Ca. Concentrations are in
    mol/l; totals are in mol/m², except charged pools in mol_c/m².
    """
    Alk, CO2_w, H, R_alk, Al_w, Al, Mg, Ca, Na, K, f_Al, f_Mg, f_Na, f_K, f_H, f_Ca = p

    # Precompute
    nZrs1000 = n * Zr * s * 1000

    return (
        (Alk_tot-R_alk)-Alk*nZrs1000,
        IC_tot-(CO2_w*(1+k1/H+k2*k1/(H**2))*s+(CO2_w/k_H)*(1-s))*(n*Zr*1000),
        (k1*CO2_w/H+2*k1*k2*CO2_w/(H**2)-H+k_w/H)-Alk,
        R_alk-(f_Mg+f_Ca+f_Na+f_K)*CEC_tot,
        Al_w*nZrs1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot,
        Al-(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w,
        Mg*nZrs1000+f_Mg/2*CEC_tot-Mg_tot,
        Ca*nZrs1000+f_Ca/2*CEC_tot-Ca_tot,
        Na*nZrs1000+f_Na*CEC_tot-Na_tot,
        K*nZrs1000+f_K*CEC_tot-K_tot,
        f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2),
        f_Mg - Mg*(f_Ca/(K_Ca_Mg*Ca)),
        f_Na - Na*(f_Ca/(K_Ca_Na*Ca))**(1/2),
        f_K - K*(f_Ca/(K_Ca_K*Ca))**(1/2),
        f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2),
        1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H)
    )


def cec_calcium_equation(p, Al, conv_Al, Ca, Mg, Na, K, H,
                         K_Ca_Al, K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H):
    """CEC saturation residual for the calcium fraction ``p[0]``."""
    f_Ca = p[0]
    return 1-(f_Ca+(Al/conv_Al)*((f_Ca**3/(K_Ca_Al*Ca**3))**(1/2))
              +Mg*(f_Ca/(K_Ca_Mg*Ca))
              +Na*((f_Ca/(K_Ca_Na*Ca))**(1/2))
              +K*((f_Ca/(K_Ca_K*Ca))**(1/2))
              +H*((f_Ca/(K_Ca_H*Ca))**(1/2)))


def total_to_cec_equations(p, H, n, Zr, s, CEC_tot, conv_Al,
                           K1, K2, K3, K4, Ca_tot, Mg_tot, K_tot,
                           Na_tot, f_acid, K_Ca_Al, K_Ca_Mg, K_Ca_Na,
                           K_Ca_K, K_Ca_H):
    """Residuals for the 13 initial concentrations and CEC fractions."""
    Al_w, Al, Al_tot, Mg, Ca, Na, K, f_Mg, f_Na, f_K, f_Ca, f_Al, f_H = p
    return (Al_w*n*Zr*s*1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot,
            Al-(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w,
            Mg*n*Zr*s*1000+f_Mg/2*CEC_tot-Mg_tot,
            Ca*n*Zr*s*1000+f_Ca/2*CEC_tot-Ca_tot,
            Na*n*Zr*s*1000+f_Na*CEC_tot-Na_tot,
            K*n*Zr*s*1000+f_K*CEC_tot-K_tot,
            f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2),
            f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2),
            f_H + f_Al - f_acid,
            f_Mg - Mg*(f_Ca/(K_Ca_Mg*Ca)),
            f_Na - Na*(f_Ca/(K_Ca_Na*Ca))**(1/2),
            f_K - K*(f_Ca/(K_Ca_K*Ca))**(1/2),
            1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H))


def kelland_equations(p, Al_w, Al, H, Ca, Ca_tot,
                      f_Mg, f_K, f_Na, n, Zr, s, CEC_tot, conv_Al,
                      K_Ca_Al, K_Ca_H):
    """Residuals for the five Kelland initial-condition unknowns."""
    Al_tot, CaCO3, f_Al, f_H, f_Ca = p
    return (Al_w*n*Zr*s*1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot,
            Ca*n*Zr*s*1000+f_Ca/2*CEC_tot+CaCO3-Ca_tot,
            f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2),
            f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2),
            1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H))
