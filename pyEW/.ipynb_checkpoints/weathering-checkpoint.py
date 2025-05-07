# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import numpy as np
import pyEW

#------------------------------------------------------------------------------
 # Carbonate weathering [mol-conv/d]
 #In soil,precipitates form as discontinuous coatings on the surfaces of soil pores, so the precipitation surface area and geometry are indeterminate. https://nora.nerc.ac.uk/id/eprint/511084/1/Kirk%20et%20al%202015%20Geochmica%20et%20Cosmochimica%20Acta.pdf
   
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
    
def Omega_sil(mineral, Ca, Mg, K, Na, Si, H, K_sp, conv_mol):
        
    if mineral == 'forsterite':
        Omega = min(1,(Mg/conv_mol)**(1/2)*(Si/conv_mol)**(1/4)/(H/conv_mol)/K_sp)
    
    elif mineral == 'wollastonite':
        Omega = min(1,(Ca/conv_mol)**(1/2)*(Si/conv_mol)**(1/2)/(H/conv_mol)/K_sp)
    
    elif mineral == 'diopside':
        #CaMgSi2O6 + (H+) -> 1/4 Ca++ + 1/4Mg++ 1/2 H2SiO3
        Omega = min(1,(Ca/conv_mol)**(1/4)*(Mg/conv_mol)**(1/4)*(Si/conv_mol)**(1/2)/(H/conv_mol)/K_sp)
    
    elif mineral == 'albite':
        #NaAlSiO3 + (H+) + 5/2 H2O -> 1/2 kaolinite + 2HSiO3 + (Na+)
        Omega = min(1,(Na/conv_mol)*(Si/conv_mol)**2/(H/conv_mol)/K_sp)
    
    elif mineral == 'anorthite':
        #1/2 CaAl2Si2O8+ (H+) + 1/2 H2O -> 1/2 kaolinite + 1/2 (Ca++)
        Omega = min(1,(Ca/conv_mol)**(1/2)/(H/conv_mol)/K_sp)
        
    elif mineral == 'muscovite':
        #KAl3Si3O10(OH)2 + (H+) + 3/2 H2O -> 3/2 kaolinite + K+
        Omega = min(1,(K/conv_mol)/(H/conv_mol)/K_sp)
    
    elif mineral in ['labradorite', 'augite', 'alkali_feldspar', 'Fe_forsterite', 'nepheline','apatite','leucite']:
        Omega = 0
    
    else:
        raise ValueError("Unknown mineral")
        
    return Omega
