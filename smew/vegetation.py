# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
from numba import njit
from numba.typed import List
import numpy as np
from numpy.typing import NDArray


@njit
def get_stage_boundaries(
    stage_array: NDArray[np.int8],
    target_stage: int
) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """
    Extracts start and stop indices ``target_stage`` of the list ``stage_array``
    """
    # Using Numba's typed lists for fast appending during the loop
    starts = List()
    stops = List()

    in_block = False

    for i in range(len(stage_array)):
        if stage_array[i] == target_stage:
            if not in_block:
                starts.append(i)
                in_block = True
        else:
            if in_block:
                stops.append(i)
                in_block = False

    # If the array ends while still inside a block, cap it off
    if in_block:
        stops.append(len(stage_array))

    # Convert typed lists back to standard NumPy arrays before returning
    return np.array(starts, dtype=np.int32), np.array(stops, dtype=np.int32)


@njit
def get_season_boundaries(
    stage_array: NDArray[np.int8]
) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """
    Extracts start and stop indices for the entire growing season (any stage > 0)
    """
    starts = List()
    stops = List()

    in_block = False

    for i in range(len(stage_array)):
        if stage_array[i] > 0:
            if not in_block:
                starts.append(i)
                in_block = True
        else:
            if in_block:
                stops.append(i)
                in_block = False

    if in_block:
        stops.append(len(stage_array))

    return np.array(starts, dtype=np.int32), np.array(stops, dtype=np.int32)


#------------------------------------------------------------------------------
 # vegetation growth
    
def veg(v_in, T_v, k_v, t0_v, temp_soil,dt):
    
    v = np.zeros(len(temp_soil))
    i_in = int(t0_v/dt)
    v[i_in] = v_in
    
    for i in range(i_in+1, len(temp_soil)):
        
        v[i] = v[i-1] + (6/(T_v*k_v))*v[i-1]*(k_v-v[i-1])*dt           
          
    return(v)


@njit
def veg_seasonal(
        growing_stage: NDArray[np.int8],
        k_v: float,
        f_v_in: float
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """
    Calculates the seasonal vegetation growth and harvest with a logistic growth curve based on the carrying capacity
    ``k_v`` and the growing seasons provided by ``growing_stage``.

    :param growing_stage: A 1D NumPy array representing the sequential stages of the crop season (0 for off-season, 1 -
        4 for active stages). The function uses stage 2 to define the growth period and emergence and the end of stage 4
        for the harvesting.
    :param k_v: The vegetation carrying capacity (maximum vegetation).
    :param f_v_in: The initial vegetation fraction (v_in = f_v_in * k_v).
    :return: A tuple containing two 1D NumPy arrays, both matching the length of the input ``growing_stage`` mask:
        * **v**: The vegetation value/biomass amount calculated via the logistic growth curve.
        * **harvest_out**: The harvested vegetation values, recorded exactly on the time step following the end of each
        growing season (all other time steps are 0.0).
    """
    v = np.zeros(len(growing_stage), dtype=np.float32)
    h = np.zeros(len(growing_stage), dtype=np.float32)

    starts_season, ends_season = get_season_boundaries(growing_stage)
    starts_growing, ends_growing = get_stage_boundaries(growing_stage, 2)  # stage 2 is the growing period

    for emergence, growth, harvest in zip(starts_growing, ends_growing, ends_season):
        v[emergence] = k_v * f_v_in  # emergence from soil
        T_v = growth - emergence  # typical growth time
        for i in range(emergence + 1, harvest + 1):
            v[i] = v[i-1] + (6/(T_v*k_v))*v[i-1]*(k_v-v[i-1])
        # check if harvest is included in the time axis (for the last harvest)
        if harvest + 1 < len(growing_stage):
            h[harvest + 1] = v[harvest]

    return v, h  # type: ignore


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