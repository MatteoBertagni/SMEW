# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 12:19:42 2019
"""
import numpy as np
import pyEW
from statistics import mean

def respiration(ADD, SOC_in, CO2_air_in, ratio_aut_het, soil, s, v, k_v, Zr, temp_soil,dt,conv_mol):
      
    # Preallocating the variables
    f_d = np.zeros(len(s))
    DEC = np.zeros(len(s))
    SOC = np.zeros(len(s))
    k_dec_T = np.zeros(len(s))
    
    #constants
    [MM_Mg, MM_Ca, MM_Na, MM_K, MM_Si, MM_C, MM_Anions, MM_Al]=pyEW.MM(conv_mol)
    [s_h, s_w, s_i, b, K_s, n] = pyEW.soil_const(soil) 
    r = 0.7 # [-]: Fraction of carbon that goes into respiration
    CO2_atm = pyEW.CO2_atm(conv_mol) # [mol-conv/l] 
    
    # SOC initial condition 
    if SOC_in>0:
        SOC[0] = SOC_in # [gOC/m3] prescribed
    else: 
        SOC[0] = (ADD/Zr)/(r*mean(k_dec_T)*mean(f_d)) #qs-equilibrium 
    
    #moisture impact on decomposition
    for i in range(0, len(s)):
        if s[i]<=s_h:
            f_d[i] = 0
        elif s[i]>s_h and s[i]<=s_w:
            f_d[i] = (s[i]-s_h)/(s_w-s_h)
        elif s[i]>s_w and s[i]<=s_i:
            f_d[i] = 1
        elif s[i]>s_i and s[i]<=1:
            f_d[i] = (1-s[i])/(1-s_i)
           
    #CO2 diffusive flux
    if Zr <= 0.3:
        Z_CO2 = Zr/2
    else:
        Z_CO2 = 0.15
        
    D_0 = pyEW.D_0() #free-air diffusion [m2/d]
    D = D_0*(1-s)**(10/3)*n**(4/3) #Mill-Quirk (1961)
    Fs_in = (D[0]*1000/(Z_CO2))*(CO2_air_in - CO2_atm) #mol-conv/m2
    
    #CO2 balance (resp_het + resp_aut = Fs)
    k_dec = MM_C*Fs_in/ (r*Zr*f_d[0]*SOC[0]*(1 + ratio_aut_het * v / k_v)) # [1/d]
    
    #temperature influence
    k_dec_T = k_dec*temp_soil/temp_soil[0]   
    k_dec_T[k_dec_T<0] = 0
                 
    # OC equation
    DEC[0] = k_dec_T[0]*f_d[0]*SOC[0]
    for i in range(1, len(s)):
        SOC[i] = SOC[i-1]+(ADD/Zr-r*DEC[i-1])*dt               
        DEC[i] = k_dec_T[i]*f_d[i]*SOC[i] # [gOC/(m3*d)]
        
    #CO2 respiration 
    r_het = r*DEC*Zr/MM_C #mol-conv/ m2 d 
    r_aut = ratio_aut_het*r_het*v/k_v #if this changes, the initial equilibrium condition above must be changed
                         
    return(SOC, r_het, r_aut, D)