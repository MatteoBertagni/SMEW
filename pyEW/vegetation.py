# -*- coding: utf-8 -*-

import numpy as np
from statistics import mean

#------------------------------------------------------------------------------
 # vegetation growth
    
def veg(v_in, T_v, k_v, t0_v, temp_soil,dt):
    
    v = np.zeros(len(temp_soil))
    i_in = int(t0_v/dt)
    v[i_in] = v_in
    
    for i in range(i_in+1, len(temp_soil)):
        
        v[i] = v[i-1] + (6/(T_v*k_v))*v[i-1]*(k_v-v[i-1])*dt           
          
    return(v)

#------------------------------------------------------------------------------
 # active uptake [Ca, Mg, K, Si] inspired by Porporato et al (2003, AWR)  and Porporato (2021, ecohydrology book)
def up_act(v, delta_v, xi, dt, T, Ca, Mg, K, Si, Dw, Zr, k_v, RAI, root_d):
    
    UP_act = np.zeros(len(xi))
    
    #demand for growth
    dem = xi*delta_v/dt # [mol-conv/d] 
    
    #passive uptake
    conc = np.array([Ca, Mg, K, Si])
    UP_p = conc*T*1000 # [mol-conv/d] passive uptake
    
    #active uptake
    for j in range(0, len(dem)):
        if UP_p[j] > dem[j] or UP_p[j]==0:
            UP_act[j] = 0
        else: 
            UP_act[j] = min(v/k_v*RAI*Dw*conc[j]*1000/(root_d*Zr/(v/k_v*RAI))**(1/2), dem[j] - UP_p[j])
        
    return UP_act