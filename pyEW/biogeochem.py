# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import numpy as np
import pyEW
from scipy import optimize
from scipy.optimize import fsolve
#minimize, least_squares, newton_krylov, broyden1, root, broyden2

def biogeochem_balance(n, s, L, T, I, v, k_v, RAI, root_d, Zr, r_het, r_aut, D, temp_soil, pH_in, conc_in, f_CEC_in, K_CEC, CEC_tot, Si_in, CaCO3_in, MgCO3_in, M_rock_in, t_app, mineral, rock_f_in, d_in, psd_perc_in, SSA_in, diss_f, dt, conv_Al, conv_mol, keyword_add):
            
    # Preallocating the variables
    pH = np.zeros(len(s))
    H = np.zeros(len(s))
    f_H = np.zeros(len(s))
    
    Ca_tot = np.zeros(len(s))
    Ca = np.zeros(len(s))
    f_Ca = np.zeros(len(s))
    UP_Ca = np.zeros(len(s))
    
    Omega_CaCO3 = np.zeros(len(s))
    CaCO3 = np.zeros(len(s))
    W_CaCO3 = np.zeros(len(s))
    
    Mg_tot = np.zeros(len(s))
    Mg = np.zeros(len(s))
    f_Mg = np.zeros(len(s))
    UP_Mg = np.zeros(len(s))
    
    Omega_MgCO3 = np.zeros(len(s))
    MgCO3 = np.zeros(len(s))
    W_MgCO3 = np.zeros(len(s))
    
    K_tot = np.zeros(len(s))
    K = np.zeros(len(s))
    f_K = np.zeros(len(s))
    UP_K = np.zeros(len(s))
    
    Na_tot = np.zeros(len(s))
    Na = np.zeros(len(s))
    f_Na = np.zeros(len(s))
    
    Si = np.zeros(len(s))
    Si_tot = np.zeros(len(s))
    UP_Si = np.zeros(len(s))
    
    root_ex = np.zeros([len(s)])
    
    An = np.zeros([len(s)])
    An_tot = np.zeros([len(s)])
    R_alk = np.zeros(len(s))
    Alk_tot = np.zeros(len(s))
    Alk = np.zeros(len(s))
    
    CO2_air = np.zeros(len(s))
    IC_tot = np.zeros(len(s))
    CO2_w = np.zeros(len(s))
    HCO3 = np.zeros(len(s))
    CO3 = np.zeros(len(s))
    DIC = np.zeros(len(s))
    Fs = np.zeros(len(s))
    ADV = np.zeros(len(s))
       
    H_rain = np.zeros(len(s))
    DIC_rain = np.zeros(len(s))
    
    Al = np.zeros(len(s))
    f_Al = np.zeros(len(s))
    AlOH = np.zeros(len(s))
    AlOH2 = np.zeros(len(s))
    AlOH3 = np.zeros(len(s)) 
    AlOH4 = np.zeros(len(s))
    Al_w = np.zeros(len(s))
    Al_tot = np.zeros(len(s))
    
    M_rock = np.zeros(len(s))
    SA = np.zeros(len(s))
    EW = np.zeros([1, len(s)])
    min_st = np.zeros([1, 6])
    
    d = np.zeros([1, len(s)])
    delta_d = np.zeros([1, len(s)])
    lamb = np.zeros([1, len(s)])
    SSA = np.zeros([1, len(s)])
    psd = np.zeros([1, len(s)])
     
    if M_rock_in > 0:
        number_min = len(mineral)
        rho_min = np.zeros(number_min)
        MM_min = np.zeros(number_min)
        E_H = np.zeros(number_min) 
        E_w = np.zeros(number_min)
        E_OH = np.zeros(number_min)
        n_H = np.zeros(number_min)
        n_OH = np.zeros(number_min)
        K_sp = np.zeros(number_min)
        min_st = np.zeros([number_min, 6])

        k_diss_H_t = np.zeros([number_min, len(s)])
        k_diss_w_t = np.zeros([number_min, len(s)])
        k_diss_OH_t = np.zeros([number_min, len(s)])
        Omega =  np.zeros([number_min, len(s)])
        M_min = np.zeros([number_min, len(s)])
        rock_f = np.zeros([number_min, len(s)])
        Wr = np.zeros([number_min, len(s)])
        EW = np.zeros([number_min, len(s)])
        
        n_d_cl = len(d_in)
        if n_d_cl > 1:
            d = np.zeros([n_d_cl, len(s)])
            delta_d = np.zeros([n_d_cl, len(s)]) 
            lamb = np.zeros([n_d_cl, len(s)])
            SSA = np.zeros([n_d_cl, len(s)])
            psd = np.zeros([n_d_cl, len(s)])
        
    errors = np.zeros([16, len(s)]) 
    
#------------------------------------------------------------------------------
    # Constants
    CO2_atm = pyEW.CO2_atm(conv_mol) # [mol_CO2/l_air] Atmospheric CO2 concentration
    T_K = temp_soil + 273.15
    
    # soil CO2 diffusivity 
    D_0 = pyEW.D_0() #free-air diffusion [m2/d]
    D = D_0*(1-s)**(10/3)*n**(4/3) #Mill-Quirk (1961)
    
    #solute diffusivity in soil water
    Dw_0 = pyEW.Dw_0()
    Dw = Dw_0*(n*s)**2 # Archie 1942, Grathwohl 1998 (book)
    
    # [g/mol-conv]: Molar masses    
    [MM_Mg, MM_Ca, MM_Na, MM_K, MM_Si, MM_C, MM_Anions, MM_Al]=pyEW.MM(conv_mol) 
    
    # Aluminium speciation
    [K1, K2, K3, K4] = pyEW.K_Al(conv_mol) 
    
    # carbonate spec  
    [k1, k2, k_w, k_H] = pyEW.K_C(T_K,conv_mol)  
    
    #CEC Gaines-Thomas constants
    [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]  = K_CEC
    
    #nutrient uptake by plants
    [v_f_Ca, v_f_Mg, v_f_K, v_f_Si] = pyEW.plant_nutr_f()
    dry_perc = 0.1 #percent of dry mass
    xi = dry_perc*np.array([v_f_Ca/MM_Ca, v_f_Mg/MM_Mg, v_f_K/MM_K, v_f_Si/MM_Si]) # [mol-conv/g_biomass]
    
    #carb weathering constants
    [K_CaCO3,K_MgCO3,r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3] = pyEW.carb_weath_const(conv_mol)
    
    #mineral constants
    if M_rock_in > 0: 
        for j in range(0,number_min):
            MM_min[j], k_diss_H_t[j,:], k_diss_w_t[j,:], k_diss_OH_t[j,:], n_H[j], n_OH[j], min_st[j,:], K_sp[j] = pyEW.min_const(mineral[j], T_K, conv_mol)  
    
    #rock density
    rho_rock = 3*1e6 # [g/m3]
    
    #rock surface fractality (Beerling 2020)
    b = 0.35 #[-]
    a = (1/(2*1e-10))**b #[1/m^b]
    
#------------------------------------------------------------------------------
    # RAINWATER
    
    Alk_rain = 0 #alk
    CO2_w_rain = k_H*CO2_atm # [mol/l] Henry's law
    
    for i in range(0, len(s)):
        def equations(p):
            H_rain[i] = p
            return(Alk_rain-(k1[i]*CO2_w_rain[i]/H_rain[i]+2*k1[i]*k2[i]*CO2_w_rain[i]/(H_rain[i]**2)-H_rain[i]+k_w[i]/H_rain[i]))
        
        H_rain[i] = fsolve(equations, 10**-6*conv_mol) # [mol/l]
        DIC_rain[i]=CO2_w_rain[i]+k1[i]*CO2_w_rain[i]/H_rain[i]+k2[i]*k1[i]*CO2_w_rain[i]/(H_rain[i]**2)
    
#------------------------------------------------------------------------------            
    # INITIAL CONDITIONS
    
    #pH
    pH[0] = pH_in 
    H[0] = 10**(-pH[0])*conv_mol 
           
    #pCO2
    if Zr <= 0.3:
        Z_CO2 = Zr/2
    else:
        Z_CO2 = 0.15
    CO2_air[0] = (r_het[0]+r_aut[0])/(D[0]*1000/(Z_CO2))+CO2_atm #mol-conv/l (Fs = resp_het + resp_aut)
    Fs[0] = D[0]/(Z_CO2)*(CO2_air[0]-CO2_atm)*1000 # [mol-conv/d]
    CO2_w[0] = k_H[0]*CO2_air[0] # [mol-conv/l] Henry's law
    
    #carbonate system
    HCO3[0] = k1[0]*CO2_w[0]/H[0] # [mol/l]
    CO3[0] = k2[0]*k1[0]*CO2_w[0]/(H[0]**2) # [mol/l]
    DIC[0] = HCO3[0]+CO3[0]+CO2_w[0]
    IC_tot[0] = (DIC[0]*s[0]+CO2_air[0]*(1-s[0]))*(n*Zr*1000) # [mol]
        
    #Alk
    Alk[0]=HCO3[0]+2*CO3[0]-H[0]+k_w[0]/H[0]    
    
    # cations (mol/l)
    [Ca[0], Mg[0], K[0], Na[0], Al_w[0]] = conc_in
           
    # anions (mol_c/l)
    An[0] = 2*Mg[0]+2*Ca[0]+Na[0]+K[0]-Alk[0] #[mol_c/l]
    
    if An[0]<0:
        print(An[0])
        raise ValueError("Not enough cations for this alkalinity")
        
    # aluminium speciation
    Al[0]=(H[0]**4/(H[0]**4+H[0]**3*K1+H[0]**2*K1*K2+H[0]*K1*K2*K3+K1*K2*K3*K4))*Al_w[0] #mol/l
    AlOH[0]=(H[0]**3*K1/(H[0]**4+H[0]**3*K1+H[0]**2*K1*K2+H[0]*K1*K2*K3+K1*K2*K3*K4))*Al_w[0]
    AlOH2[0]=(H[0]**2*K1*K2/(H[0]**4+H[0]**3*K1+H[0]**2*K1*K2+H[0]*K1*K2*K3+K1*K2*K3*K4))*Al_w[0]
    AlOH3[0]=(H[0]*K1*K2*K3/(H[0]**4+H[0]**3*K1+H[0]**2*K1*K2+H[0]*K1*K2*K3+K1*K2*K3*K4))*Al_w[0]
    AlOH4[0]=Al_w[0]-(Al[0]+AlOH[0]+AlOH2[0]+AlOH3[0])
    
    # Silicon
    Si[0] = Si_in
    
    #Background inputs (rain, litterfall, background weathering..)
    if keyword_add == 1:
        I_An = np.mean(T+L)*1000*An[0]*s[0]/np.mean(s) #[mol_c d-1]
        I_Ca = np.mean(T+L)*1000*Ca[0]*s[0]/np.mean(s) #[mol d-1]
        I_Mg = np.mean(T+L)*1000*Mg[0]*s[0]/np.mean(s)
        I_Na = np.mean(T+L)*1000*Na[0]*s[0]/np.mean(s)
        I_K = np.mean(T+L)*1000*K[0]*s[0]/np.mean(s)
        I_Si = np.mean(T+L)*1000*Si[0]*s[0]/np.mean(s)
    elif keyword_add == 0:
        I_An = 0
        I_Ca = 0
        I_Mg = 0 
        I_K = 0
        I_Na = 0 
        I_Si = 0
    
    #CEC adsorbed species
    [f_Ca[0], f_Mg[0], f_K[0], f_Na[0], f_Al[0], f_H[0]] = f_CEC_in
    
    #reserve of alkalinity
    R_alk[0] = (f_Mg[0]+f_Ca[0]+f_Na[0]+f_K[0])*CEC_tot # [mol_c]
    
    #total amounts (solution and adsorbed)
    Ca_tot[0] = Ca[0]*n*s[0]*Zr*1000+f_Ca[0]/2*CEC_tot # [mol] 
    Mg_tot[0] = Mg[0]*n*s[0]*Zr*1000+f_Mg[0]/2*CEC_tot # [mol] 
    K_tot[0] = K[0]*n*s[0]*Zr*1000+f_K[0]*CEC_tot # [mol]
    Na_tot[0] = Na[0]*n*s[0]*Zr*1000+f_Na[0]*CEC_tot # [mol]
    Alk_tot[0] = 2*Mg_tot[0]+2*Ca_tot[0]+Na_tot[0]+K_tot[0]-An[0]*(n*s[0]*Zr*1000) # [mol_c]
    An_tot[0] = An[0]*n*s[0]*Zr*1000 #[mol_c]
    Al_tot[0] = Al_w[0]*n*Zr*s[0]*1000+(f_Al[0]/3)*CEC_tot*conv_Al # [mol]
    Si_tot[0] = Si[0]*n*Zr*s[0]*1000
    
    #Carbonate minerals (considered as an additional pool)
    CaCO3[0] = CaCO3_in # [mol-conv]
    MgCO3[0] = MgCO3_in
    
    #Carbonate weathering
    Omega_CaCO3[0] = Ca[0]*CO3[0]/K_CaCO3 # [-]
    Omega_MgCO3[0] = Mg[0]*CO3[0]/K_MgCO3
    [W_CaCO3[0], W_MgCO3[0]] = pyEW.carb_W(CaCO3[0], MgCO3[0], Omega_CaCO3[0], Omega_MgCO3[0], s[0], Zr, r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3) # [mol-conv/ m2 d]
        
    #Silicate weathering
    if M_rock_in > 0:
        
        #application timestep
        tt_app = int(t_app/dt)
        
        #rock composition
        M_rock[tt_app] = M_rock_in #[g/m2]
        rock_f[:,tt_app] = rock_f_in
        M_min[:,tt_app] = rock_f[:,tt_app]*M_rock[tt_app] #[g/m2]
        M_iner = M_rock[tt_app]*(1-np.sum(rock_f[:,tt_app])) #[g/m2]
        
        #diameter classes
        d[:,tt_app] = d_in #[m]
        delta_d[:,tt_app] = np.insert(np.diff(d[:,tt_app]),0,d[0,tt_app])
        
        #particle size distribution
        psd[:,tt_app] = psd_perc_in*M_rock[tt_app]/delta_d[:,tt_app] #[g/m]
        
        #refinment of fractal constant based on measured SSA 
        if SSA_in > 0:
            a = (SSA_in*rho_rock*M_rock[tt_app]/6)/np.sum(d[:,tt_app]**(b-1)*psd[:,tt_app]*delta_d[:,tt_app]) #[m**-b]
        
        #surface area
        lamb[:,tt_app] = a*d[:,tt_app]**b #[-]
        SSA[:,tt_app] = 6/(d[:,tt_app]*rho_rock)*lamb[:,tt_app] # [m2/g]
        SA[tt_app] = np.sum(SSA[:,tt_app]*psd[:,tt_app]*delta_d[:,tt_app]) #[m2]
                    
        #mineral weathering
        if t_app == 0:
            for j in range(0, number_min):
                Omega[j,0] = pyEW.Omega_sil(mineral[j], Ca[0], Mg[0], K[0], Na[0], Si[0], H[0], K_sp[j], conv_mol) #[-]
                Wr[j,0] = s[0]*diss_f*(k_diss_H_t[j,0]*(H[0]/conv_mol)**n_H[j]+k_diss_w_t[j,0]+k_diss_OH_t[j,0]*(H[0]/conv_mol)**n_OH[j])*(1-Omega[j,0]) # [mol-conv/ m2 d]
                EW[j,0] = Wr[j,0]*SA[0]*rock_f[j,0] # [mol-conv/d]

#------------------------------------------------------------------------------
    #SYSTEM RESOLUTION
        
    printcounter = 0
    numb_print = 0
    
    for i in range(1, len(s)): 
        
            #CO2 advection due to moisture variation    
            if s[i]<s[i-1]:
                ADV[i] = n*Zr*1000*(s[i]-s[i-1])*CO2_atm # [mol] 
            elif s[i]>s[i-1]:
                ADV[i] = n*Zr*1000*(s[i]-s[i-1])*CO2_air[i-1]
            
            #active uptake [Ca, Mg, K, Si] 
            UP_act = pyEW.up_act(v[i], (v[i]-v[i-1]), xi, dt, T[i-1], Ca[i-1], Mg[i-1], K[i-1], Si[i-1], Dw[i-1], Zr, k_v, RAI, root_d)
            UP_Ca[i-1], UP_Mg[i-1], UP_K[i-1], UP_Si[i-1] = UP_act # [mol-conv/d] 
                                  
            #explicit mass balances # [mol]
            Ca_tot[i] = Ca_tot[i-1]+(I_Ca+np.sum(min_st[:,0]*EW[:,i-1])+W_CaCO3[i-1]-(L[i-1]+T[i-1])*1000*Ca[i-1]-UP_Ca[i-1])*dt 
            Mg_tot[i] = Mg_tot[i-1]+(I_Mg+np.sum(min_st[:,1]*EW[:,i-1])+W_MgCO3[i-1]-(L[i-1]+T[i-1])*1000*Mg[i-1]-UP_Mg[i-1])*dt
            K_tot[i] = K_tot[i-1]+(I_K+np.sum(min_st[:,2]*EW[:,i-1])-(L[i-1]+T[i-1])*1000*K[i-1]-UP_K[i-1])*dt
            Na_tot[i] = Na_tot[i-1]+(I_Na+np.sum(min_st[:,3]*EW[:,i-1])-(L[i-1]+T[i-1])*1000*Na[i-1])*dt
            Al_tot[i] = Al_tot[i-1]+(np.sum(min_st[:,4]*EW[:,i-1])*conv_Al-L[i-1]*1000*(Al[i-1]+AlOH4[i-1]))*dt
            Si_tot[i] = Si_tot[i-1]+(I_Si+np.sum(min_st[:,5]*EW[:,i-1])-(L[i-1]+T[i-1])*1000*Si[i-1]-UP_Si[i-1])*dt
            An_tot[i] = An_tot[i-1]+(I_An - (L[i-1]+T[i-1])*An[i-1]*1000)*dt # [mol_c]
            Alk_tot[i] = 2*Mg_tot[i]+2*Ca_tot[i]+Na_tot[i]+K_tot[i]-An_tot[i] # [mol_c]
            IC_tot[i] = IC_tot[i-1]+I[i]*1000*DIC_rain[i]-ADV[i]+(W_CaCO3[i-1]+W_MgCO3[i-1]+r_het[i-1]+r_aut[i-1]-Fs[i-1]-L[i-1]*1000*DIC[i-1])*dt 
                       
            #implicit system
            def equations(p):
                Alk[i], CO2_w[i], H[i], R_alk[i], Al_w[i], Al[i], Mg[i], Ca[i], Na[i], K[i], f_Al[i], f_Mg[i], f_Na[i], f_K[i], f_H[i], f_Ca[i] = p
            
                return((Alk_tot[i]-R_alk[i])-Alk[i]*(n*Zr*s[i]*1000),\
                       IC_tot[i]-(CO2_w[i]*(1+k1[i]/H[i]+k2[i]*k1[i]/(H[i]**2))*s[i]+(CO2_w[i]/k_H[i])*(1-s[i]))*(n*Zr*1000),\
                       (k1[i]*CO2_w[i]/H[i]+2*k1[i]*k2[i]*CO2_w[i]/(H[i]**2)-H[i]+k_w[i]/H[i])-Alk[i],\
                       R_alk[i]-(f_Mg[i]+f_Ca[i]+f_Na[i]+f_K[i])*CEC_tot,\
                       Al_w[i]*n*Zr*s[i]*1000+(f_Al[i]/3)*CEC_tot*conv_Al-Al_tot[i],\
                       Al[i]-(H[i]**4/(H[i]**4+H[i]**3*K1+H[i]**2*K1*K2+H[i]*K1*K2*K3+K1*K2*K3*K4))*Al_w[i],\
                       Mg[i]*n*Zr*s[i]*1000+f_Mg[i]/2*CEC_tot-Mg_tot[i],\
                       Ca[i]*n*Zr*s[i]*1000+f_Ca[i]/2*CEC_tot-Ca_tot[i],\
                       Na[i]*n*Zr*s[i]*1000+f_Na[i]*CEC_tot-Na_tot[i],\
                       K[i]*n*Zr*s[i]*1000+f_K[i]*CEC_tot-K_tot[i],\
                       f_Al[i] - (Al[i]/conv_Al)*(f_Ca[i]**3/(K_Ca_Al*Ca[i]**3))**(1/2),\
                       f_Mg[i] - Mg[i]*(f_Ca[i]/(K_Ca_Mg*Ca[i])),\
                       f_Na[i] - Na[i]*(f_Ca[i]/(K_Ca_Na*Ca[i]))**(1/2),\
                       f_K[i] - K[i]*(f_Ca[i]/(K_Ca_K*Ca[i]))**(1/2),\
                       f_H[i] - H[i]*(f_Ca[i]/(K_Ca_H*Ca[i]))**(1/2),\
                       1-(f_Ca[i]+f_Al[i]+f_Mg[i]+f_Na[i]+f_K[i]+f_H[i]))
                       
            #initial guess
            Alk0 = (Alk_tot[i]-R_alk[i-1])/(n*Zr*s[i]*1000)
            CO2_w0 = IC_tot[i]/(n*Zr*1000)*1/(s[i]*(1+k1[i]/H[i-1]+k2[i]*k1[i]/(H[i-1]**2))+(1-s[i])/k_H[i]) 
            R_alk0 = R_alk[i-1]
            Al_w0 = (Al_tot[i]-(f_Al[i-1]/3)*CEC_tot*conv_Al)/(n*Zr*s[i]*1000)#s[i-1]*Al_w[i-1]/s[i]
            Al0 = (H[i-1]**4/(H[i-1]**4+H[i-1]**3*K1+H[i-1]**2*K1*K2+H[i-1]*K1*K2*K3+K1*K2*K3*K4))*Al_w0
            Mg0 = (Mg_tot[i]-f_Mg[i-1]/2*CEC_tot)/(n*Zr*s[i]*1000) #s[i-1]*Mg[i-1]/s[i] 
            Na0 = (Na_tot[i]-f_Na[i-1]*CEC_tot)/(n*Zr*s[i]*1000) #s[i-1]*Na[i-1]/s[i] 
            Ca0 = (Ca_tot[i]-f_Ca[i-1]/2*CEC_tot)/(n*Zr*s[i]*1000) #s[i-1]*Ca[i-1]/s[i]
            K0 =  (K_tot[i]-f_K[i-1]*CEC_tot)/(n*Zr*s[i]*1000) #s[i-1]*K[i-1]/s[i]
            H0 = H[i-1]
            def eqH(p):
                    H0 = p
                    return((k1[i]*CO2_w0/H0+2*k1[i]*k2[i]*CO2_w0/(H0**2)-H0+k_w[i]/H0)-Alk0)
            H0_2 =fsolve(eqH, H[i-1])[0]
            
            #solution 1
            x0 = np.array([Alk0, CO2_w0, H0, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al[i-1],f_Mg[i-1], f_Na[i-1], f_K[i-1], f_H[i-1], f_Ca[i-1]])         
            sol = fsolve(equations,x0, xtol=1e-12)                                           
            errors[:,i] = equations(sol) #residuals
            
            #solution 2
            res_threshold = 1e-1
            if np.any(abs(errors[:,i]) > res_threshold):
                x0 = np.array([Alk0, CO2_w0, H0_2, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al[i-1],f_Mg[i-1], f_Na[i-1], f_K[i-1], f_H[i-1], f_Ca[i-1]])
                sol = fsolve(equations, x0, xtol=1e-14)
                errors[:,i] = equations(sol)
                if np.any(abs(errors[:,i]) > res_threshold):
                    print(i)
                    raise ValueError("Solution not converging")          

            #pH and C
            pH[i] = -np.log10(H[i]/conv_mol) # [-]
            CO2_air[i] = CO2_w[i]/k_H[i] #[mol/l]
            HCO3[i] = k1[i]*CO2_w[i]/H[i] 
            CO3[i] = k2[i]*k1[i]*CO2_w[i]/(H[i]**2)
            DIC[i] = CO2_w[i]+HCO3[i]+CO3[i]
              
            #Al speciation
            AlOH[i] = (H[i]**3*K1/(H[i]**4+H[i]**3*K1+H[i]**2*K1*K2+H[i]*K1*K2*K3+K1*K2*K3*K4))*Al_w[i] #[mol/l]
            AlOH2[i] = (H[i]**2*K1*K2/(H[i]**4+H[i]**3*K1+H[i]**2*K1*K2+H[i]*K1*K2*K3+K1*K2*K3*K4))*Al_w[i]
            AlOH3[i] = (H[i]*K1*K2*K3/(H[i]**4+H[i]**3*K1+H[i]**2*K1*K2+H[i]*K1*K2*K3+K1*K2*K3*K4))*Al_w[i]
            AlOH4[i] = Al_w[i]-(Al[i]+AlOH[i]+AlOH2[i]+AlOH3[i])
            
            #concentrations
            Si[i] = Si_tot[i]/(n*Zr*s[i]*1000) # [mol-conv/l]
            An[i] = An_tot[i]/(n*Zr*s[i]*1000) # [mol_c-conv/l]
                                                 
            #CO2 diff flux 
            Fs[i] = D[i]/(Z_CO2)*(CO2_air[i]-CO2_atm)*1000 # [mol/d]    
            
            #Carbonate minerals
            CaCO3[i] = CaCO3[i-1] - W_CaCO3[i-1]*dt # [mol-conv]
            MgCO3[i] = MgCO3[i-1] - W_MgCO3[i-1]*dt
            
            #Carbonate weathering
            Omega_CaCO3[i] = Ca[i]*CO3[i]/K_CaCO3 # [-]
            Omega_MgCO3[i] = Mg[i]*CO3[i]/K_MgCO3
            [W_CaCO3[i], W_MgCO3[i]] = pyEW.carb_W(CaCO3[i], MgCO3[i], Omega_CaCO3[i], Omega_MgCO3[i], s[i], Zr, r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3)
                     
            #Silicate weathering
            if M_rock_in > 0:
                
                #application time
                if i == tt_app:
                    for j in range(0, number_min):
                        Omega[j,i] = pyEW.Omega_sil(mineral[j], Ca[i], Mg[i], K[i], Na[i], Si[i], H[i], K_sp[j], conv_mol) #[-]
                        Wr[j,i] = s[i]*diss_f*(k_diss_H_t[j,i]*(H[i]/conv_mol)**n_H[j]+k_diss_w_t[j,i]+k_diss_OH_t[j,i]*(H[i]/conv_mol)**n_OH[j])*(1-Omega[j,i]) # [mol-conv/ m2 d]
                        EW[j,i] = Wr[j,i]*SA[i]*rock_f[j,i] # [mol-conv/d]
                
                #post application
                elif i > tt_app:
                
                    #diameter class variation
                    d_shrink = np.sum(rock_f[:,i-1]*Wr[:,i-1]*MM_min[:]/rho_rock)*dt # [m]
                    d[:,i] = d[:,i-1] - 2*d_shrink*lamb[:,i-1] # [m]
                    d[:,i][d[:,i] < 0] = 0
                    delta_d[:,i] = np.insert(np.diff(d[:,i]),0,d[0,i]) # [m]

                    #d-derived quantities
                    for k in range(0, n_d_cl):
                        if d[k,i]>0:
                            lamb[k,i] = a*d[k,i]**b #[-]
                            SSA[k,i] = 6/(d[k,i]*rho_rock)*lamb[k,i] # [m2/g]
                            psd[k,i] = psd[k,i-1]*(d[k,i]/d[k,i-1])**3*(delta_d[k,i-1]/delta_d[k,i]) # [g/m]

                    #mineral weathering
                    for j in range(0, number_min):
                        if M_min[j,i-1] != 0: 
                            Omega[j,i] = pyEW.Omega_sil(mineral[j], Ca[i], Mg[i], K[i], Na[i], Si[i], H[i], K_sp[j], conv_mol) #[-]
                            Wr[j,i] = s[i]*diss_f*(k_diss_H_t[j,i]*(H[i]/conv_mol)**n_H[j]+k_diss_w_t[j,i]+k_diss_OH_t[j,i]*(H[i]/conv_mol)**n_OH[j])*(1-Omega[j,i]) # [mol/ m2 d]
                            M_min[j,i] = M_min[j,i-1]-EW[j,i-1]*MM_min[j]*dt # [g]
                            if M_min[j,i] < 0:
                                M_min[j,i] = 0

                    #global quantities
                    SA[i] = np.sum(SSA[:,i]*psd[:,i]*delta_d[:,i]) # [m2]
                    M_rock[i] = np.sum(M_min[:,i]) + M_iner # [g]
                    if M_rock[i]>0:
                        rock_f[:,i] = M_min[:,i]/M_rock[i] # [-]
                    EW[:,i] = Wr[:,i]*SA[i]*rock_f[:,i] # [mol/d]
                    

    data = {k: v for k, v in locals().items()}
                               
    return data