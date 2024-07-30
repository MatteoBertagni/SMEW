# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:30:44 2019
"""
import numpy as np
import scipy.stats
import pyEW


def moisture_balance(rain, Zr, soil, ET0, v, k_v, keyword_wb, s_in,t_end,dt):
    
    #constants
    [s_h, s_w, s_i, b, K_s, n] = pyEW.soil_const(soil)       
    
    # Initialization
    #--------------------------------------------------------------------------      
    if keyword_wb == 1:
        s=np.zeros([len(rain)])
        s[0] = s_in # initial value
    elif keyword_wb == 0:
        s = s_in*np.ones(round(t_end/dt))
        
    L = np.zeros([len(s)])
    E = np.zeros([len(s)])
    T = np.zeros([len(s)])
    Q = np.zeros([len(s)])
    Irr = np.zeros([len(s)])
    
    # moisture dynamics
    #--------------------------------------------------------------------------      
    if keyword_wb == 1:
        
        for i in range(0, len(rain)-1):

            # Evaporation [m/d]
            E0 = 0.5*ET0
            if s[i]<=s_h:
                E[i] = 0
            elif s[i]<=s_i:
                E[i] = (s[i]-s_h)/(s_i-s_h)*E0[i]*(1-v[i]/k_v)
            elif s[i]<=1:
                E[i] = E0[i]*(1-v[i]/k_v)

            # Transpiration [m/d]
            if s[i]<=s_w:
                T[i] = 0
            elif s[i]<=s_i:
                T[i] = (s[i]-s_w)/(s_i-s_w)*ET0[i]*v[i]/k_v
            elif s[i]<=1:
                T[i] = ET0[i]*v[i]/k_v

            # Leakage [m/d]
            L[i] = K_s*s[i]**(3+2*b)
            
            # Moisture dynamics
            s[i+1] = s[i]+rain[i+1]/(n*Zr)-((E[i]+T[i]+L[i])/(n*Zr)*dt)

            #runoff [m]
            if s[i+1] >= 1:
                Q[i+1] = (s[i+1]-1)*(n*Zr)
                s[i+1] = 1  

            #to avoid numerical issues
            if s[i+1] >= 0.98:
                s[i+1] = 0.98

        #infiltration [m]
        I=rain-Q
    
    # constant moisture
    #-------------------------------------------------------------------------- 
    if keyword_wb == 0:
          
        # leakage [m/d]
        L = K_s*s**(3+2*b)
        # Evaporation [m/d]
        if s_in>=s_h:
            E = (s_in-s_h)/((s_i+1)/2-s_h)*ET0*(1-v/k_v)
        # Transpiration [m/d]
        if s_in>=s_w and s_in<=s_i:
            T = (s-s_w)/(s_i-s_w)*ET0*v/k_v
        elif s_in>=s_i:
            T = ET0*v/k_v
       
        rain = E+T+L #[m]
        I = E+T+L 
        Q = np.zeros(len(s))       
    
    return(s, s_w, s_i, I, L, T, E, Q, Irr, n)