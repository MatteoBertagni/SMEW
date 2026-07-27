# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
    
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

def wet_f_Anand(pore_d, pore_pdf, s, d, psd_rock, mixalf=1.0, lmax=None):
    """
    Calculate the nonlinear scaling of the surface area with moisture
    by computing the fraction of rock-particle surface area in contact with soil water.

    Based on Anand et al. (2026), Water Resources Research,
    doi:10.1029/2025WR041479.

    The largest water-filled pore size is obtained by inverting the
    cumulative soil pore-size distribution at relative soil water 
    saturation `s` (Equation 2). Rock-powder sizes are mapped to effective
    pore locations according to the mixing parameter `mixalf` (Equation 4).
    The wet surface fraction is then calculated as the normalized
    cumulative surface area of particles located in water-filled pores
    (Equation 3).

    Continuous interpolation is used when inverting the pore-size
    distribution and evaluation of the cumulative rock surface area.

    Parameters
    ----------
    pore_d : array_like
        Soil pore-size grid.
    pore_pdf : array_like
        Probability density associated with the pore-size grid.
    s : float
        Relative soil moisture or degree of pore saturation [0-1].
    d : array_like
        Rock-powder size.
    psd_rock : array_like
        Rock-powder number distribution.
    mixalf : float, optional
        Rock-soil mixing parameter, where 1 represents perfect mixing
        and smaller values shift rock-powder particles toward larger soil pores.
    lmax : float, optional
        Maximum effective pore location. If not given, `max(d)` is used.

    Returns
    -------
    float
        Fraction of total rock-particle surface area that is wet [0-1].
    """

    if pore_d.size != pore_pdf.size:
        raise ValueError("pore_d and pore_pdf must have the same length.")

    if d.size != psd_rock.size:
        raise ValueError("d and psd_rock must have the same length.")

    if not 0 < mixalf <= 1:
        raise ValueError("mixalf must be greater than 0 and no larger than 1.")

    if lmax is None:
        lmax = np.max(d)

    s = float(np.clip(s, 0.0, 1.0))

    # Convert soil moisture into largest water-filled pore
    pore_grid, pore_cdf = normalized_cumulative_area( pore_d, pore_pdf)

    # Continuous inverse of Equation 2: F_p(rw) = s (Anand et al., 2026, WRR)
    rw = np.interp(s, pore_cdf, pore_grid)

    # Convert particle diameter into its effective location
    c = (1.0 - mixalf) * lmax
    d_eff = mixalf * d + c

    # Equation 3 integrand after transforming d to pore location (Anand et al., WRR, 2026)
    rock_area_density = ( d**2 * psd_rock / mixalf)

    rock_grid, rock_area_cdf = normalized_cumulative_area(d_eff, rock_area_density)

    # Evaluate cumulative wet surface area continuously at rw
    wet_f = np.interp( rw, rock_grid, rock_area_cdf, left=0.0, right=1.0)

    return float(np.clip(wet_f, 0.0, 1.0))
 
#------------------------------------------------------------------------------

def normalized_cumulative_area(x, y):
    """
    Construct a normalized cumulative distribution by numerically
    integrating a nonnegative density over an increasing grid.

    The returned cumulative values range from 0 to 1 and are used to
    represent the pore-size CDF in Equation 2 and the normalized
    cumulative rock surface area in Equation 3 of Anand et al. WRR (2026).
    """
    
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Ensure increasing x and remove repeated x values
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    x_unique, unique_idx = np.unique(x, return_index=True)
    y_unique = y[unique_idx]

    # The distributions should not have negative density
    y_unique = np.maximum(y_unique, 0.0)

    cumulative = cumulative_trapezoid( y_unique, x_unique, initial=0.0)

    total = cumulative[-1]

    if total <= 0:
        raise ValueError("The distribution has zero total area.")

    return x_unique, cumulative / total
    
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
