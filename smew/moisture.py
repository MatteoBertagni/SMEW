# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:30:44 2019
"""
import numpy as np
import smew
from numba import njit


@njit
def moisture_balance(rain, Zr, soil, ET0, v, k_v, keyword_wb, s_in,t_end,dt,
                     temp_soil=None,
                     melt_rate = 0.005):
    
    """
    temp_soil : array_like, optional
    Soil temperature [degC]. If None, freezing is ignored and the
    original water balance is used. If provided, precipitation is stored
    as snow/ice when temp_soil <= 0 and melts when temp_soil > 0.
    melt rate: float, optional
    degree-day melt rate [m water equivalent / d / °C]
    """
    
    #constants
    [s_h, s_w, s_i, b, K_s, n] = smew.soil_const(soil) 
    
    # Initialization
    #--------------------------------------------------------------------------      
    if keyword_wb == 1:
        s=np.zeros((len(rain)))
        s[0] = s_in # initial value
    elif keyword_wb == 0:
        s = s_in*np.ones(round(t_end/dt))
        
    L = np.zeros((len(s)))
    E = np.zeros((len(s)))
    T = np.zeros((len(s)))
    Q = np.zeros((len(s)))
    Irr = np.zeros((len(s)))
    I = np.zeros((len(s)))
    snowpack = np.zeros((len(s))) # snow/ice storage [m water equivalent]
    
    # moisture dynamics
    #--------------------------------------------------------------------------      
    if keyword_wb == 1:

        # Potential evaporation [m/d]
        E0 = 0.5*ET0
        
        for i in range(0, len(rain)-1):

            if temp_soil is None:
                frozen = False
            else:
                frozen = temp_soil[i] <= 0

            # Frozen soil: precipitation stored as snow/ice, and no water fluxes
            if frozen:
                E[i] = 0.0
                T[i] = 0.0
                L[i] = 0.0
                Q[i + 1] = 0.0
                I[i + 1] = 0.0

                snowpack[i + 1] = snowpack[i] + rain[i + 1]
                s[i + 1] = s[i]

                continue

            # Unfrozen step: stored snow/ice can melt and contribute to liquid input.
            if temp_soil is None:
                melt = 0.0
            else:
                melt = min(snowpack[i], melt_rate * max(temp_soil[i] - 0.0, 0.0) * dt)
            
            snowpack[i + 1] = snowpack[i] - melt

            liquid_input = rain[i + 1] + melt
            
            # Evaporation [m/d]
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
            s[i+1] = s[i]+liquid_input/(n*Zr)-((E[i]+T[i]+L[i])/(n*Zr)*dt)

            # runoff [m], using 0.98 as numerical saturation cap
            if s[i+1] >= 0.98:
                Q[i+1] = (s[i+1] - 0.98) * (n * Zr)
                s[i+1] = 0.98

            # actual liquid infiltration [m]
            I[i+1] = liquid_input - Q[i+1]

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

        if temp_soil is not None:
            for i in range(0, len(s)):
                if temp_soil[i] <= 0:
                    E[i] = 0.0
                    T[i] = 0.0
                    L[i] = 0.0
       
        rain = E+T+L #[m]
        I = E+T+L 
        Q = np.zeros(len(s))       
    
    return s, s_w, s_i, I, L, T, E, Q, Irr, n