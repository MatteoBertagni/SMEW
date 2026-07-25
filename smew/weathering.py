# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import numpy as np
from scipy.integrate import simpson
    
#------------------------------------------------------------------------------
 # psd evolution (based on Beerling et al., 2020)

def psd_evol(d, delta_d, d_0, delta_d0, psd_0, n_d_cl, a, b, rho_rock):  

    lamb = np.zeros(n_d_cl)
    SSA = np.zeros(n_d_cl)
    psd = np.zeros(n_d_cl)
    
    for k in range(n_d_cl):
                        if d[k]>0:
                            lamb[k] = a*d[k]**b #[-]
                            SSA[k] = 6/(d[k]*rho_rock)*lamb[k] # [m2/g]
                            psd[k] = psd_0[k]*(d[k]/d_0[k])**3*(delta_d0[k]/delta_d[k]) # [g/m]
  
    SA = np.sum(SSA[:]*psd[:]*delta_d[:]) # [m2]
    
    return lamb, SSA, psd, SA
#------------------------------------------------------------------------------

# from psd by mass to psd by number

def psd_number_from_mass(mass_distribution_rock, diameter_rock, density_rock):
    volume_rock = np.pi * diameter_rock**3 / 6
    number_distribution = mass_distribution_rock / (density_rock * volume_rock)
    return np.nan_to_num(number_distribution, nan=0.0, posinf=0.0, neginf=0.0)

#------------------------------------------------------------------------------

 # Wetness factor for the rock surface area in the soil [0-1]

def wetness_SA(s, keyword_ssa, pore_d, pore_pdf, d, psd_rock, mixalf, lmax):

    # no scaling of the surface area with moisture (e.g., Beerling et al., 2020, Nature)
    if keyword_ssa == 'constant': 
        wet_f = 1

    # linear scaling of the surface area with moisture (e.g., Cipolla et al., 2021, WRR)
    elif keyword_ssa == 'linear': 
        wet_f = s

    # nonlinear scaling of the surface area with moisture (Anand et al., 2026, WRR)
    elif keyword_ssa == 'nonlinear': 
        wet_f = wet_f_Anand(pore_d, pore_pdf, s, d, psd_rock, mixalf, lmax)

    return wet_f

#------------------------------------------------------------------------------

    # nonlinear scaling of the surface area with moisture (Anand et al., 2026, WRR)

def wet_f_Anand(pore_d, pore_pdf, s, d, psd_rock, mixalf=1.0, lmax=None):

    if lmax is None:
        lmax = d[-1]

    if s >= 1:
        rw = pore_d[-1]
    else:
        rw = cumulative_area_index(pore_d, pore_pdf, s)[0]

    c = (1.0 - mixalf) * lmax
    d_eff = mixalf * d + c

    if rw >= d_eff[-1]:
        wet_f = 1.0
    else:
        idx = np.searchsorted(d_eff, rw)

        area_cum = cumulative_area(
            d_eff,
            d**2 * (1.0 / mixalf) * psd_rock
        )

        wet_f = area_cum[idx] / area_cum[-1]

    return wet_f

#------------------------------------------------------------------------------

def cumulative_area(xg, yg):
    cum_area = np.zeros_like(xg, dtype=float)

    for i in range(len(xg)):
        x_sub = xg[:i+1]
        y_sub = yg[:i+1]

        valid_idx = np.where(np.diff(x_sub, prepend=np.nan) != 0)[0]
        x_valid = x_sub[valid_idx]
        y_valid = y_sub[valid_idx]

        if len(x_valid) > 1:
            cum_area[i] = simpson(y=y_valid, x=x_valid)

    return cum_area


def cumulative_area_index(xg, yg, aint):
    cum_area = cumulative_area(xg, yg)
    idx = np.argmax(cum_area >= aint)
    return xg[idx], idx

#------------------------------------------------------------------------------
 
# Carbonate weathering [mol-conv/d]
#In soil, precipitates form as discontinuous coatings on the surfaces of soil pores, so the precipitation surface area and geometry are indeterminate. https://nora.nerc.ac.uk/id/eprint/511084/1/Kirk%20et%20al%202015%20Geochmica%20et%20Cosmochimica%20Acta.pdf
   
def carb_W(CaCO3, MgCO3, Omega_CaCO3, Omega_MgCO3, s, Zr, r_CaCO3, r_MgCO3, tau_CaCO3, tau_MgCO3):
        
    #CaCO3
    if Omega_CaCO3 <= 1:
        W_CaCO3 = s*CaCO3*(1-Omega_CaCO3)/tau_CaCO3 # dissolution
    else:
        W_CaCO3 = r_CaCO3*Zr*(1-Omega_CaCO3)        # precipitation 
    
    #MgCO3
    if Omega_MgCO3 <= 1:
        W_MgCO3 = s*MgCO3*(1-Omega_MgCO3)/tau_MgCO3 # dissolution
    else:
        W_MgCO3 = r_MgCO3*Zr*(1-Omega_MgCO3)        # precipitation
                                      
    return (W_CaCO3, W_MgCO3)

#------------------------------------------------------------------------------

 # Silicate weathering rate (based on Palandri et al., 2004)

def sil_Wr(mineral, Omega, H, k_H_T, k_w_T, k_OH_T, n_H, n_OH, diss_f, conv_mol):
    
    #weathering rate [mol-conv/ m2 d]
    Wr = diss_f*(k_H_T*(H/conv_mol)**n_H + k_w_T + k_OH_T*(H/conv_mol)**n_OH)*(1-Omega)
    
    return Wr
    
#------------------------------------------------------------------------------
 # Silicate saturation index (Omega)
    
def sil_Omega(mineral, Ca, Mg, K, Na, Al, AlOH4, Si, H, K_sp, conv_mol, conv_Al):
        
    if mineral == 'albite':
        #NaAlSi3O8 + 4(H+) + 4 H2O -> Na+ + Al3+ + 3H4SiO4 + (Na+)
        Omega = min(1,(Na/conv_mol)*(Al/(conv_mol*conv_Al))*(Si/conv_mol)**3/(H/conv_mol)**4/K_sp)
    
    elif mineral == 'anorthite':
        #1/2 CaAl2Si2O8 + (H+) + 1/2 H2O -> 1/2 kaolinite + 1/2 (Ca++)
        Omega = min(1,(Ca/conv_mol)**(1/2)/(H/conv_mol)/K_sp) 
    
    elif mineral == 'analcime':
        #NaAlSi2O6(H2O) + H2O -> Na+ + Al(OH)4- + SiO2
        Omega = min(1,(Na/conv_mol)*(AlOH4/(conv_mol*conv_Al))*(Si/conv_mol)**2/K_sp) 
    
    elif mineral == 'forsterite':
        Omega = min(1,(Mg/conv_mol)**(1/2)*(Si/conv_mol)**(1/4)/(H/conv_mol)/K_sp)
    
    elif mineral == 'wollastonite':
        Omega = min(1,(Ca/conv_mol)**(1/2)*(Si/conv_mol)**(1/2)/(H/conv_mol)/K_sp)
    
    elif mineral == 'diopside':
        #CaMgSi2O6 + (H+) -> 1/4 Ca++ + 1/4Mg++ 1/2 H2SiO3
        Omega = min(1,(Ca/conv_mol)**(1/4)*(Mg/conv_mol)**(1/4)*(Si/conv_mol)**(1/2)/(H/conv_mol)/K_sp)
        
    elif mineral == 'muscovite':
        #KAl3Si3O10(OH)2 + (H+) + 3/2 H2O -> 3/2 kaolinite + K+
        Omega = min(1,(K/conv_mol)/(H/conv_mol)/K_sp)
    
    elif mineral in ['labradorite', 'augite', 'alkali_feldspar', 'Fe_forsterite', 'nepheline','apatite','leucite']:
        Omega = 0
    
    else:
        raise ValueError("Unknown mineral")
        
    return Omega
