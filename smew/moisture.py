# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:30:44 2019
"""
from typing import Any

from numba import njit
import numpy as np
from numpy.typing import NDArray

import smew


@njit
def moisture_balance(
    rain: NDArray[Any],
    Zr: float,
    soil: str,
    ET0: NDArray[Any] | None,
    v: NDArray[Any],
    k_v: float,
    keyword_wb: int,
    s_in: float,
    t_end: float,
    dt: float,
    ETc: NDArray[Any] | None = None,
) -> tuple[
    NDArray[Any],  # s
    float,    # s_w
    float,    # s_i
    NDArray[Any],  # I
    NDArray[Any],  # L
    NDArray[Any],  # T
    NDArray[Any],  # E
    NDArray[Any],  # Q
    NDArray[Any],  # Irr
    float     # n
]:
    """
    Computes the moisture balance.

    If ET0 is provided:

        * When vegetation = 0, E = 0.5 * ET0 and T = 0.
        * When vegetation = 1, E = 0 and T = ET0.

    If ETc is provided:

        * When vegetation = 0, E = ET0 and T = 0.
        * When vegetation = 1, E = 0 and T = ET0.

    :param rain: Array of rainfall values over time [m/d].
    :param Zr: Root zone depth of the soil [m].
    :param soil: Identifier for the soil type used for defining soil constants (USDA soil types classification).
    :param ET0: Reference evapotranspiration values [m/d], required if ETc is None.
    :param v: Biomass mass over time.
    :param k_v: Biomass carrying capacity (aka maximum possible biomass).
    :param keyword_wb: Indicator for water balance mode (1 for dynamic moisture dynamics, 0 for constant moisture).
    :param s_in: Initial soil moisture content (dimensionless fraction).
    :param t_end: Simulation end time [days].
    :param dt: Time step for the simulations [days].
    :param ETc: Crop evapotranspiration values [m/d], required if ET0 is None.

    :returns: A tuple containing:

              * **s** Soil moisture content over time (dimensionless fractions).
              * **s_w** Soil moisture wilting point (dimensionless fractions).
              * **s_i** Soil moisture field capacity (dimensionless fractions).
              * **I** Infiltration over time [m].
              * **L** Leakage over time [m/d].
              * **T** Transpiration over time [m/d].
              * **E** Evaporation over time [m/d].
              * **Q** Runoff over time [m].
              * **Irr** Irrigation over time [m].
              * **n** Soil porosity (volume fraction).

    :raises ValueError: If `keyword_wb` is not 0 or 1, or both `ET0` and `ETc` are None or not None.
    """
    
    #constants
    [s_h, s_w, s_i, b, K_s, n] = smew.soil_const(soil)       
    
    # Initialization
    #--------------------------------------------------------------------------      
    if keyword_wb == 1:
        s = np.zeros((len(rain)), dtype=rain.dtype)
        s[0] = s_in # initial value
    elif keyword_wb == 0:
        s = s_in*np.ones(round(t_end/dt), dtype=rain.dtype)
    else:
        raise ValueError(f"Invalid keyword wb: {keyword_wb}")

    if ET0 is not None and ETc is None:
        ET_pot = ET0
    elif ET0 is None and ETc is not None:
        ET_pot = ETc
    else:
        raise ValueError(f"Exactly one of ET0 or ETc must be None")
        
    L = np.zeros((len(s)), dtype=rain.dtype)
    E = np.zeros((len(s)), dtype=rain.dtype)
    T = np.zeros((len(s)), dtype=rain.dtype)
    Q = np.zeros((len(s)), dtype=rain.dtype)
    Irr = np.zeros((len(s)), dtype=rain.dtype)
    
    # moisture dynamics
    #--------------------------------------------------------------------------      
    if keyword_wb == 1:

        # Evaporation [m/d]
        if ETc is None:
            E0 = 0.5*ET0
        else:
            E0 = ETc
        
        for i in range(0, len(rain)-1):

            if s[i]<=s_h:
                E[i] = 0
            elif s[i]<=s_i:
                E[i] = (s[i]-s_h)/(s_i-s_h)*E0[i]*(1-v[i]/k_v)
            # elif s[i]<=1:
            else:
                E[i] = E0[i]*(1-v[i]/k_v)

            # Transpiration [m/d]
            if s[i]<=s_w:
                T[i] = 0
            elif s[i]<=s_i:
                T[i] = (s[i]-s_w)/(s_i-s_w)*ET_pot[i]*v[i]/k_v
            # elif s[i]<=1:
            else:
                T[i] = ET_pot[i]*v[i]/k_v

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
            elif s[i+1] <= 0.01:
                s[i+1] = 0.01

        #infiltration [m]
        I=rain-Q
    
    # constant moisture
    #-------------------------------------------------------------------------- 
    # if keyword_wb == 0:
    else:
        # leakage [m/d]
        L = (K_s*s**(3+2*b)).astype(rain.dtype)
        # Evaporation [m/d]
        if s_in>=s_h:
            E = ((s_in-s_h)/((s_i+1)/2-s_h)*ET_pot*(1-v/k_v)).astype(rain.dtype)
        # Transpiration [m/d]
        if s_in>=s_w and s_in<=s_i:
            T = ((s-s_w)/(s_i-s_w)*ET_pot*v/k_v).astype(rain.dtype)
        elif s_in>=s_i:
            T = (ET_pot*v/k_v).astype(rain.dtype)
       
        rain = E+T+L #[m]
        I = E+T+L 
        Q = np.zeros(len(s), dtype=rain.dtype)
    
    return s, s_w, s_i, I, L, T, E, Q, Irr, n
