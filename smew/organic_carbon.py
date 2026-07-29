# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 12:19:42 2019
"""
import numpy as np
import smew
from statistics import mean

def respiration(ADD, SOC_in, CO2_air_in, ratio_aut_het, soil, s, v, k_v, Zr, temp_soil,dt,conv_mol, tau_OC=None):
      
    # Preallocating the variables
    f_s = np.zeros(len(s))
    f_T = np.zeros(len(s))
    DEC = np.zeros(len(s))
    SOC = np.zeros(len(s))
    
    #constants
    [MM_Mg, MM_Ca, MM_Na, MM_K, MM_Si, MM_C, MM_Anions, MM_Al]=smew.MM(conv_mol)
    [s_h, s_w, s_i, b, K_s, n] = smew.soil_const(soil) 
    r = 0.7 # [-]: Fraction of carbon that goes into respiration
    CO2_atm = smew.CO2_atm(conv_mol) # [mol-conv/l]
    D_0 = smew.D_0() #free-air diffusion [m2/d]
    D = D_0*(1-s)**(10/3)*n**(4/3) #Mill-Quirk (1961)
      
    #moisture impact on decomposition
    for i in range(0, len(s)):
        if s[i]<=s_h:
            f_s[i] = 0
        elif s[i]>s_h and s[i]<=s_w:
            f_s[i] = (s[i]-s_h)/(s_w-s_h)
        elif s[i]>s_w and s[i]<=s_i:
            f_s[i] = 1
        elif s[i]>s_i and s[i]<=1:
            f_s[i] = (1-s[i])/(1-s_i)

    # temperature impact on decomposition; decomposition stops below 0 degC
    temp_active = np.maximum(temp_soil, 0.0)

    if np.any(temp_active > 0.0):
        f_T = temp_active / np.mean(temp_active[temp_active > 0.0])
    else:
        f_T = 0.0
       
    #CO2 gas-diffusion baricenter
    if Zr <= 0.3:
        Z_CO2 = Zr/2
    else:
        Z_CO2 = 0.15

    # input data: SOC_in and either tau_OC or CO2_air_in 
    # missing data (None) estimated via qs-state approximation
    if SOC_in is not None:
        SOC[0] = SOC_in # [gOC/m3]
        if tau_OC is not None:
            k_dec = 1/tau_OC # [1/d]
            #Fs_in = k_dec/MM_C*(r*Zr*f_s[0]*f_T[0]*SOC[0]*(1 + ratio_aut_het * v / k_v)) # [mol-conv/m2] (resp_het + resp_aut = Fs)
        elif CO2_air_in is not None:
            Fs_in = (D[0]*1000/(Z_CO2))*(CO2_air_in - CO2_atm) # [mol-conv/m2] 
            activity0 = f_s[0] * f_T[0] * (1.0 + ratio_aut_het * v[0] / k_v)
            if activity0 <= 0.0:
                raise ValueError("Cannot estimate k_dec from CO2_air_in when the initial soil is "
                                 "frozen or biologically inactive. Provide tau_OC or start from an "
                                 "unfrozen active timestep.")
            k_dec = MM_C * Fs_in / (r * Zr * activity0 * SOC[0])

    # mean decomposition activity
    f_dec = f_T * f_s
    mean_f_dec = np.mean(f_dec)

    if mean_f_dec <= 0.0:
        raise ValueError("Mean decomposition activity is zero. Cannot estimate ADD or SOC "
            "from steady-state balance.")

    # ADD estimate for qs-equilibrium (in absence of data)
    if ADD is None and SOC[0] is not None:
        ADD = r * Zr * k_dec * mean_f_dec * SOC[0]  # [gOC/(m2*d)]

    # SOC estimate for qs-equilibrium (in absence of data)
    if SOC_in is None and ADD is not None:
        SOC[0] = (ADD / Zr) / (r * k_dec * mean_f_dec)  # [gOC/m3]     
           
    # OC equation
    DEC[0] = k_dec*f_T[0]*f_s[0]*SOC[0]
    for i in range(1, len(s)):
        SOC[i] = SOC[i-1]+(ADD/Zr-r*DEC[i-1])*dt               
        DEC[i] = k_dec*f_T[i]*f_s[i]*SOC[i] # [gOC/(m3*d)]
        
    #CO2 respiration 
    r_het = r*DEC*Zr/MM_C #mol-conv/ m2 d 
    r_aut = ratio_aut_het*r_het*v/k_v #if this changes, the initial equilibrium condition above must be changed
                         
    return(SOC, r_het, r_aut, D)