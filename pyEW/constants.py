# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import numpy as np


#------------------------------------------------------------------------------
 # atmospheric CO2
    
def CO2_atm(conv_mol):
    
    CO2_atm = 412E-6/22.41*conv_mol # [mol-conv/l] 
                                                        
    return CO2_atm

#------------------------------------------------------------------------------
 # soil CO2 diffusivity 
    
def D_0():
    
    D_0 = 1.6E-5*3600*24 #free-air diffusion [m2/d]
                                                        
    return D_0

#------------------------------------------------------------------------------
 # solute diffusivity in soil water
    
def Dw_0():
    
    Dw_0 = 1e-9*3600*24 # [m2/d]
                                                        
    return Dw_0

#------------------------------------------------------------------------------
 # % of nutrients in plant dry matter (ideally plant dependent)
    
def plant_nutr_f():
    
    #current values are from Kelland's (2020) measurement for sorghum
    
    v_f_Ca = 0.005 #suggested 2%, range 0.1 - 5 %, Weil and Bredy 2017
    v_f_Mg = 0.001 #suggested, 0.5%, range 0.1 - 1 %, Weil and Bredy 2017
    v_f_K = 0.01 #suggested 2%, range 1 - 5 %, Weil and Bredy 2017
    v_f_Si = 0.01 #suggested 5%, range 1 - 10%, Epstein 1994, PNAS
    
    v_f = [v_f_Ca, v_f_Mg, v_f_K, v_f_Si]
    
                                                      
    return v_f

#------------------------------------------------------------------------------
 # soil constants

def soil_const(soil):
    
    # values are from Laio et al., (2001, Adv. Water Resources)
    
    if soil == 'sand':
        s_h = 0.08 #hygroscopic point
        s_w = 0.11 #wilting
        s_i = 0.33 #max transpiration 
        #s_fc = 0.35 #field capacity - not needed in the current formulation
        b = 4.05 #power law exponent of the rentention curve
        K_s = 14  # (m/d) 
        n = 0.35  # porosity
    elif soil == 'loamy sand':
        s_h = 0.08
        s_w = 0.11
        s_i = 0.31
        #s_fc = 0.52
        b = 4.4
        K_s = 13
        n = 0.42
    elif soil == 'sandy loam':
        s_h = 0.14
        s_w = 0.18
        s_i = 0.46
        #s_fc = 0.56
        b = 4.9
        K_s = 3
        n = 0.43
    elif soil == 'loam':
        s_h = 0.19
        s_w = 0.24
        s_i = 0.57
        #s_fc = 0.65
        b = 5.4
        K_s = 0.6
        n = 0.45
    elif soil == 'clay loam':
        s_h = 0.39
        s_w = 0.45
        s_i = 0.68
        b = 8.5
        K_s = 0.2 
        n = 0.47
    elif soil == 'clay':
        s_h = 0.47
        s_w = 0.52
        s_i = 0.78
        b = 11.4
        K_s = 0.11 
        n = 0.5
    else:
        print("Invalid soil type!")
                                                      
    return(s_h, s_w, s_i, b, K_s, n)

#------------------------------------------------------------------------------
 # EW mineral constants 

def min_const(mineral, T_K, conv_mol):
    
    # All values are from Palandri (2004)
    # The exceptions are the solubility product constants (K_sp), which, when available, are from Morel et al., (1993) -- https://booksite.elsevier.com/9780120885305/appendices/Web_Appendices.pdf
    
    k_diss_H_t = np.zeros(len(T_K))
    k_diss_w_t = np.zeros(len(T_K))
    k_diss_OH_t = np.zeros(len(T_K))
    
    R = 8.314/conv_mol # [J mol-1 K-1]: universal gas constant 
    T_ref = 25+273.15 # [K]: temperature standard conditions
        
    if mineral == 'forsterite': #Mg2SiO4
            #rho_min = 2.8*1e6 # [g/m3]: density
            MM_min = 140/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-6.85)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-10.64)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 0
            E_H = 67.2/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 79/conv_mol# [kJ/mol-conv]
            E_OH = 0
            n_H = 0.47 # reaction order
            n_OH = 1
            min_st = [0, 2, 0, 0, 0, 1]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = 10**(7.11) 
            
    elif mineral == 'Fe_forsterite': #FeMgSiO4
            MM_min = 172/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-6.85)*24*3600*conv_mol #(-5.37) [mol-conv m-2 d-1]
            k_diss_w = 10**(-10.64)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 0
            E_H = 67/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 79/conv_mol# [kJ/mol-conv]
            E_OH = 0
            n_H = 0.47 # reaction order
            n_OH = 1
            min_st = [0, 1, 0, 0, 0, 1]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = np.nan      
            
    elif mineral == 'wollastonite': #CaSiO3
            #rho_min = 2.9*1e6 # [g/m3]: density
            MM_min = 116/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-5.37)*24*3600*conv_mol #(-5.37) [mol-conv m-2 d-1]
            k_diss_w = 10**(-8.88)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 0
            E_H = 54.7/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 54.7/conv_mol# [kJ/mol-conv]
            E_OH = 0
            n_H = 0.4 # reaction order
            n_OH = 1
            min_st = [1, 0, 0, 0, 0, 1]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = 10**(6.82)
            
    elif mineral == 'albite': #NaAlSiO3
            #rho_min = 2.6*1e6 # [g/m3]: density
            MM_min =  263/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-10.16)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-12.56)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 10**(-15.6)*24*3600*conv_mol
            E_H = 65/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 70/conv_mol# [kJ/mol-conv]
            E_OH = 71/conv_mol
            n_H = 0.457 # reaction order
            n_OH = - 0.572
            min_st = [0, 0, 0, 1, 1, 3]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = 10**(-0.68)
            
    elif mineral == 'anorthite': #CaAl2Si2O8
            #rho_min = 2.73*1e6 # [g/m3]: density
            MM_min =  278/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-3.5)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-9.12)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 0
            E_H = 16.6/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 17.8/conv_mol# [kJ/mol-conv]
            E_OH = 0
            n_H = 1.4 # reaction order
            n_OH = 1
            min_st = [1, 0, 0, 0, 2, 2]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = 10**(9.83)
            
    elif mineral == 'labradorite': #Na0.4Ca0.6Al1.6Si2.4O8 (http://webmineral.com/data/Labradorite.shtml)
            #rho_min = 2.73*1e6 # [g/m3]: density
            MM_min =  272/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-7.87)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-10.91)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 0
            E_H = 42/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 45/conv_mol# [kJ/mol-conv]
            E_OH = 0
            n_H =  0.6 # reaction order
            n_OH = 1
            min_st = [0.6, 0, 0, 0.4, 1.6, 2.4]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = np.nan 
            
    elif mineral == 'augite': #Ca0.9Na0.1Mg0.9Fe0.2Al0.4Ti0.1Si1.9O6 (http://webmineral.com)
            #rho_min = 2.73*1e6 # [g/m3]: density
            MM_min =  236/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-6.82)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-11.97)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 0
            E_H = 78/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 78/conv_mol# [kJ/mol-conv]
            E_OH = 0
            n_H =  0.7 # reaction order
            n_OH = 1
            min_st = [0.9, 0.9, 0, 0.1, 0.4, 1.9]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = np.nan  
            
    elif mineral == 'leucite': #K(AlSi2O6)
            MM_min =  218/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-6.00)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-9.20)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 10**(-10.66)*24*3600*conv_mol
            E_H = 132/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 75/conv_mol# [kJ/mol-conv]
            E_OH = 56.6/conv_mol
            n_H =  0.7 # reaction order
            n_OH = -0.2
            min_st = [0, 0, 1, 0, 1, 2]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = np.nan  
            
    elif mineral == 'diopside': #MgCaSi2O6
            MM_min =  216/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-6.36)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-11.11)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 0
            E_H = 96/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 40/conv_mol# [kJ/mol-conv]
            E_OH = 0
            n_H =  0.71 # reaction order
            n_OH = 1
            min_st = [1, 1, 0, 0, 0, 2]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = 10**(5.30)
    
    elif mineral == 'alkali_feldspar': #K0.41Na0.56Ca0.03Al1.03Si2.97O8 (Kelland et al., 2020)
            MM_min =  156/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-10.06)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-12.41)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_OH = 10**(-21.2)*24*3600*conv_mol
            E_H = 51/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 38/conv_mol# [kJ/mol-conv]
            E_OH = 94/conv_mol
            n_H =  0.5 # reaction order
            n_OH = -0.82
            min_st = [0.03, 0, 0.41, 0.56, 1.03, 2.97]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = np.nan
            
    elif mineral == 'apatite': #Ca5(PO4)3(OH)
            MM_min =  422/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-4.29)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-6.)*24*3600*conv_mol 
            k_diss_OH = 0
            E_H = 250/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 250/conv_mol
            E_OH = 1
            n_H =  0.17 # reaction order
            n_OH = 1
            min_st = [5, 0, 0, 0, 0, 0]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = np.nan
            
    elif mineral == 'nepheline': #Na0.75K0.25(AlSiO4)
            MM_min =  146/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-2.73)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-8.56)*24*3600*conv_mol 
            k_diss_OH = 10**(-10.76)*24*3600*conv_mol 
            E_H = 63/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 65/conv_mol
            E_OH = 38/conv_mol
            n_H =  1.13 # reaction order
            n_OH = -0.2
            min_st = [0, 0, 0.25, 0.75, 1, 1]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = np.nan
            
    elif mineral == 'muscovite': #KAl3Si3O10(OH)2
            MM_min =   398.71/conv_mol # [g/mol-conv]
            k_diss_H = 10**(-11.85)*24*3600*conv_mol # [mol-conv m-2 d-1]
            k_diss_w = 10**(-13.55)*24*3600*conv_mol 
            k_diss_OH = 10**(-14.55)*24*3600*conv_mol
            E_H = 22/conv_mol# [kJ/mol-conv]: activation energies 
            E_w = 22/conv_mol
            E_OH = 22/conv_mol
            n_H =  0.37 # reaction order
            n_OH = -0.22
            min_st = [0, 0, 1, 0, 3, 3]# Stochiometric coefficients [Ca, Mg, K, Na, Al, Si]
            K_sp = 10**(4.54) #Soil and Water Chemistry An Integrative Approach (Michael E. Essington)
            
    else:
        raise ValueError("No data for this mineral")
            
    for i in range(len(T_K)): 
        k_diss_H_t[i] = k_diss_H*np.exp(-E_H*1000/R*(1/T_K[i]-1/T_ref))
        k_diss_w_t[i] = k_diss_w*np.exp(-E_w*1000/R*(1/T_K[i]-1/T_ref))
        k_diss_OH_t[i] = k_diss_OH*np.exp(-E_OH*1000/R*(1/T_K[i]-1/T_ref))
                                                      
    return(MM_min,k_diss_H_t,k_diss_w_t,k_diss_OH_t,n_H,n_OH,min_st,K_sp)

#------------------------------------------------------------------------------
 # Carbonate weathering constants  
    

def carb_weath_const(conv_mol):
    
    # Carbonate solutibility products
    #https://booksite.elsevier.com/9780120885305/appendices/Web_Appendices.pdf
    K_CaCO3 = 10**(-8.35)*conv_mol**2 # Calcite
    K_MgCO3 = 10**(-7.46)*conv_mol**2 # Magnesite
    
    #precipitation rate
    # modeling formulation follows https://nora.nerc.ac.uk/id/eprint/511084/1/Kirk%20et%20al%202015%20Geochmica%20et%20Cosmochimica%20Acta.pdf
    # parameter values are defined to reproduce experimental results
    r_CaCO3 = 3*1e7*(1e-9*1e6/(24*3600))*conv_mol # [mol-conv/d] 
    r_MgCO3 = 1e7*(1e-9*1e6/(24*3600))*conv_mol
    
    #dissolution timescale
    tau_CaCO3 = 30 # [d] 
    tau_MgCO3 = 40 
    
    return(K_CaCO3,K_MgCO3,r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3)
    
#------------------------------------------------------------------------------
 # CEC constants (Gaines-Thomas) 

def K_GT_CEC(soil, conv_mol):
            
    # Current values are from a meta-analysis of Dutch soils (0-30 cm, https://edepot.wur.nl/31605)
    # Site-specific coefficient estimates can be obtained with soil-water coupled measurements
    
    if soil in ['sand', 'loamy sand', 'sandy loam']:
        K_Ca_Mg = 10**(0.56); #[-]
        K_Ca_K = 10**(-1.16)*conv_mol #[conc]
        K_Ca_Na = 10**(0.75)*conv_mol #[conc]
        K_Ca_H = 10**(-5)*conv_mol #[conc]
        K_Ca_Al = 10**(-1.7)/conv_mol #[conc^-1] #heterovalent 
        #K_Ca_AlOH = 10**(-1.7)
        #K_Ca_AlOH2 = 10**(-1.7)*conv_mol
    
    elif soil in ['loam', 'silty loam', 'silt']:
        K_Ca_Mg = 10**(0.1); #[-]
        K_Ca_K = 10**(-2)*conv_mol
        K_Ca_Na = 10**(0.38)*conv_mol
        K_Ca_H = 10**(-5.4)*conv_mol
        K_Ca_Al = 10**(-0.86)/conv_mol
        #K_Ca_AlOH = 10**(-0.86)
        #K_Ca_AlOH2 = 10**(-0.86)*conv_mol
        
    elif soil in ['clay', 'clay loam', 'silty clay']:
        K_Ca_Mg = 10**(0.39); #[-]
        K_Ca_K = 10**(-2.42)*conv_mol
        K_Ca_Na = 10**(0.774)*conv_mol
        K_Ca_H = 10**(-6.67)*conv_mol
        K_Ca_Al = 10**(-0.2)/conv_mol
        #K_Ca_AlOH = 10**(-0.2)
        #K_Ca_AlOH2 = 10**(-0.2)*conv_mol
    
    else:
        raise ValueError("Unknown soil type")
        
    # K_Na_K = (K_Ca_K/K_Ca_Na)**(1/2)]
    # K_Na_Al = (K_Ca_Al/K_Ca_Na**3)**(1/2)
   
    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]
                                          
    return(K_CEC)

#------------------------------------------------------------------------------
 # Aluminium speciation (pag. 398 Weil and Brady)

def K_Al(conv_mol):
              
    pK1 = 5 
    pK2 = 5.1 
    pK3 = 6.7 
    pK4 = 6.2
    K1 = 10**(-pK1)*conv_mol
    K2 = 10**(-pK2)*conv_mol
    K3 = 10**(-pK3)*conv_mol
    K4 = 10**(-pK4)*conv_mol
   
    K_Al = [K1, K2, K3, K4]
                                          
    return(K_Al)

#------------------------------------------------------------------------------
 # carbonate speciation [Stumm and Morgan, 1996]

def K_C(T_K,conv_mol):
    
    T_ref = 25+273.15 # [K]: temperature standard conditions
    
    #carbonates
    pk1 = -(-356.309 - 0.0609 * T_K + 21834.37/T_K + 126.8339 * np.log10(T_K) - 1684915/(T_K)**2)
    pk2 = -(-107.887 - 0.032528 * (T_K) + 5151.79/(T_K) + 38.92561 * np.log10(T_K) - 563713.9/(T_K)**2)
    pk_w = -(-283.971 + 13323/(T_K) - 0.0507 * (T_K) + 102.24447 * np.log10(T_K) - 1119669/(T_K)**2) 
    k1 = 10**(-pk1)*conv_mol
    k2 = 10**(-pk2)*conv_mol
    
    #water
    k_w = 10**(-pk_w)*conv_mol**2
    
    #Henry
    k_H = 0.83*np.exp(2400*(1/T_K-1/T_ref))
                                          
    return(k1, k2, k_w, k_H)

#------------------------------------------------------------------------------
 # molar masses [g/mol]

def MM(conv_mol):
    
    MM_Mg = 24/conv_mol 
    MM_Ca = 40/conv_mol 
    MM_Na = 23/conv_mol
    MM_K = 39/conv_mol 
    MM_Si = 28/conv_mol 
    MM_C = 12/conv_mol
    MM_Anions = 62/conv_mol
    MM_Al = 27/conv_mol
                                          
    return(MM_Mg, MM_Ca, MM_Na, MM_K, MM_Si, MM_C, MM_Anions, MM_Al)