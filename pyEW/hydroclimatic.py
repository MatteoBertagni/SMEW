# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:30:44 2019
"""
import numpy as np
import math
import pyeto as eto
import scipy.stats
from scipy.signal import savgol_filter


#--------------------------------------------------------------------------------------------------
# air and soil temperature

def temp(latitude,temp_av, temp_ampl_yr, temp_ampl_d, Zr,t_end,dt,day1):
         
    #air temperature
    x0 = 2*np.pi*day1/365 #initial day
    TT = (365/dt)
    x = x0 + (2*np.pi)*np.arange(0, t_end/dt)/TT #argument of sin and cos
    if latitude<0:
        x_max = 20*2*np.pi/365 #day of max temperature
    elif latitude>0:
        x_max = 200*2*np.pi/365
    phi = np.pi/2-x_max # Phase lag
    temp_air = temp_av + temp_ampl_yr*np.sin(x+phi) # [°C]
    temp_min = temp_air-temp_ampl_d/2 
    temp_max = temp_air + temp_ampl_d/2 
    
    # soil temperature
    TD = 0.203    # Thermal diffusivity [mm^2/s]
    TD = TD / 10**6 * 86400 # convert to m^2/d
    dd = (365*TD/np.pi)**(1/2) # Damping depth [m]
    temp_soil = temp_av + dd*temp_ampl_yr/(2*Zr)*(np.sin(x+phi)-np.cos(x+phi)+np.exp(-Zr/dd)*(np.cos(x-Zr/dd+phi)-np.sin(x-Zr/dd+phi)))
                                                        
    return(temp_air,temp_soil,temp_min,temp_max)


#--------------------------------------------------------------------------------------------------    
# ET0 with Penman-Monteith approach (standard grass ref - FAO - Allen et al., 1998, https://www.fao.org/3/X0490E/x0490e00.htm#Contents)

def ET0(latitude,altitude,temp_air,temp_soil,temp_min,temp_max, wind,albedo,Zr,coastal,t_end,dt,day1):
    
    # initialization
    ET0 = np.zeros(len(temp_air))
        
    for i in range(0, len(ET0)):

        #atmosphere
        atm_p = eto.atm_pressure(altitude) # atmospheric pressure [kPa]
        svp = eto.svp_from_t(temp_air[i]) # saturation vapor pressure [kPa]
        avp = eto.avp_from_tmin(temp_min[i]) # actual vapor pressure [kPa]
        delta_svp = eto.delta_svp(temp_air[i]) # Slope of saturation vapour pressure curve [kPa C-1]
        psy = eto.psy_const(atm_p) #psychrometric constant [kPa degC-1].
        
        #sun/light
        j = day1+np.floor(dt*i)
        j = j - np.floor((j-1)/365)*365
        sol_dec = eto.sol_dec(j) # solar declination [radians]
        ird = eto.inv_rel_dist_earth_sun(j) #Inverse relative distance between earth and the sun
        sha = eto.sunset_hour_angle(latitude, sol_dec) # Sunset hour angle [rad]
        
        #radiation
        et_rad = eto.et_rad(latitude, sol_dec, sha, ird) # Extraterrestrial radiation [MJ m-2 d-1]
        cs_rad = eto.cs_rad(altitude, et_rad) # Clear sky radiation [MJ m-2 day-1]
        sol_rad = eto.sol_rad_from_t(et_rad, cs_rad, temp_min[i], temp_max[i], coastal) # Gross incoming solar radiation [MJ m-2 d-1]
        ni_sw_rad = eto.net_in_sol_rad(sol_rad, albedo) # net-incoming shortwave rad [MJ m-2 day-1]
        no_lw_rad = eto.net_out_lw_rad(273+temp_min[i], 273+temp_max[i], sol_rad, cs_rad, avp) # net outgoing long wave rad [MJ m-2 d-1]
        net_rad = eto.net_rad(ni_sw_rad,no_lw_rad) # net incoming solar radiation [MJ m-2 day-1]
        shf = 0 # Soil heat flux (G) [MJ m-2 day-1]

        #potential ET
        ET0[i] = eto.fao56_penman_monteith(net_rad, temp_air[i]+273.15, wind[i], svp, avp, delta_svp, psy, shf) #[mm d-1]
    
    ET0 = ET0/1000 #[m d-1]
    
    return(ET0)


#--------------------------------------------------------------------------------------------------
# stochastic rain 

def rain_stoc(lamda, alfa, t_end, dt):
    
    # simulated event number
    nb_ev = int(2*lamda*t_end) 
    
    # interarrival time [d]
    tau = scipy.stats.expon.rvs(scale = 1/lamda, loc = 0, size = int(nb_ev))

    # intensity [m]
    h = scipy.stats.expon.rvs(scale = alfa, loc = 0, size = int(nb_ev))

    # rainfall [m]
    rain = np.zeros(int(t_end/dt))
    t_event = 1 #initialization
    for i in range(0, nb_ev):
        t_event = int(t_event+tau[i]/dt)
        if t_event<=len(rain):
            rain[t_event] = h[i]
        else:
            break
                
    return rain


#--------------------------------------------------------------------------------------------------
# stochastic rain with seasonality 

def rain_stoc_season(lamda, alfa, t_end, dt):
    
    days = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) # month days
    rain = np.array

    for j in range(0, int(t_end/365)):#year
        for i in range(0, len(days)):#month
            
            # simulated rainfall events per month
            nb_ev = int(2*lamda[i]*days[i]) 
    
            # interarrival time [d]
            tau = scipy.stats.expon.rvs(scale = 1/lamda[i], loc = 0, size = int(nb_ev))

            # intensity [m]
            h = scipy.stats.expon.rvs(scale = alfa[i], loc = 0, size = int(nb_ev))

            # rainfall array [m]
            rain_month = np.zeros(int(days[i]/dt))
            t_event = 1 #initialization
            for ii in range(0, nb_ev):
                t_event = int(t_event+tau[ii]/dt)
                if t_event<(days[i]/dt):
                    rain_month[t_event] = h[ii]
            if i == 0 and j == 0:
                rain = rain_month
            else:
                rain = np.append(rain, rain_month)
    return rain