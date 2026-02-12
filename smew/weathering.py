# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import numpy as np
    
#------------------------------------------------------------------------------
 # psd evolution 
 # Beerling et al., 2020

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

 # Silicate weathering [mol-conv/d]
 # Palandri et al., 2004

def sil_Wr(mineral, Omega, s, H, k_H_T, k_w_T, k_OH_T, n_H, n_OH, diss_f, conv_mol):
    
    #weathering rate [mol-conv/ m2 d]
    Wr = s*diss_f*(k_H_T*(H/conv_mol)**n_H + k_w_T + k_OH_T*(H/conv_mol)**n_OH)*(1-Omega)
    
    return Wr
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
