# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import warnings

import numpy as np
import smew
from smew._utils import require_backend, _solve_system
from smew._utils import WATER_SYSTEM, HYDROGEN_SYSTEM, BIOGEOCHEM_SYSTEM


def biogeochem_balance(n, s, L, T, I, v, k_v, RAI, root_d, Zr, r_het, r_aut, D, temp_soil, pH_in, conc_in, f_CEC_in, K_CEC, CEC_tot, Si_in, CaCO3_in, MgCO3_in, M_rock_in, t_app, mineral, rock_f_in, d_in, psd_perc_in, SSA_in, diss_f, dt, conv_Al, conv_mol, keyword_add, 
                       keyword_ssa='linear', # options: 'constant', 'linear', 'nonlinear'
                       pore_d_in=None,
                       pore_pdf_in=None,
                       rho_rock_in=None,
                       mixalf_in=1.0,
                       *, backend="python"
                      ):
    """Run the shared model as Python/SciPy or Numba/Cython/cminpack."""
    require_backend(backend)
    # Forward the public arguments to the model; backend is handled here only.
    arguments = {
        "n": n,
        "s": s,
        "L": L,
        "T": T,
        "I": I,
        "v": v,
        "k_v": k_v,
        "RAI": RAI,
        "root_d": root_d,
        "Zr": Zr,
        "r_het": r_het,
        "r_aut": r_aut,
        "D": D,
        "temp_soil": temp_soil,
        "pH_in": pH_in,
        "conc_in": conc_in,
        "f_CEC_in": f_CEC_in,
        "K_CEC": K_CEC,
        "CEC_tot": CEC_tot,
        "Si_in": Si_in,
        "CaCO3_in": CaCO3_in,
        "MgCO3_in": MgCO3_in,
        "M_rock_in": M_rock_in,
        "t_app": t_app,
        "mineral": mineral,
        "rock_f_in": rock_f_in,
        "d_in": d_in,
        "psd_perc_in": psd_perc_in,
        "SSA_in": SSA_in,
        "diss_f": diss_f,
        "dt": dt,
        "conv_Al": conv_Al,
        "conv_mol": conv_mol,
        "keyword_add": keyword_add,
        "keyword_ssa": keyword_ssa,
        "pore_d_in": pore_d_in,
        "pore_pdf_in": pore_pdf_in,
        "rho_rock_in": rho_rock_in,
        "mixalf_in": mixalf_in,
    }
    if keyword_ssa == "nonlinear" and (pore_d_in is None or pore_pdf_in is None):
        raise ValueError(
            "For keyword_ssa='nonlinear', provide both pore_d_in and pore_pdf_in. "
            "They can be estimated with: pore_d_in, pore_pdf_in = smew.soil_pore_pdf(soil)."
        )
    # Standardise array inputs as contiguous float64 arrays, using empty arrays for None.
    for name in ("s", "L", "T", "I", "v", "r_het", "r_aut", "D", "temp_soil",
                 "conc_in", "f_CEC_in", "K_CEC", "rock_f_in", "d_in", "psd_perc_in",
                 "pore_d_in", "pore_pdf_in"):
        value = arguments[name]
        arguments[name] = np.ascontiguousarray(
            () if value is None else value, dtype=np.float64,
        )
    # Standardise mineral to a tuple.
    arguments["mineral"] = tuple(mineral) if mineral is not None and len(mineral) else ("",)
    model = _biogeochem_balance
    if backend == "compiled":
        from smew._simulation_compiled import compiled_balance
        model = compiled_balance
    result = dict(model(**arguments))
    for status in range(2, 6):
        if result["solver_status_counts"][status]:
            warnings.warn(f"MINPACK stopped with status {status}", RuntimeWarning, stacklevel=2)
    return result


def _biogeochem_balance(n, s, L, T, I, v, k_v, RAI, root_d, Zr, r_het, r_aut, D, temp_soil, pH_in, conc_in, f_CEC_in, K_CEC, CEC_tot, Si_in, CaCO3_in, MgCO3_in, M_rock_in, t_app, mineral, rock_f_in, d_in, psd_perc_in, SSA_in, diss_f, dt, conv_Al, conv_mol, keyword_add,
                       keyword_ssa='linear', # options: 'constant', 'linear', 'nonlinear'
                       pore_d_in=None,
                       pore_pdf_in=None,
                       rho_rock_in=None,
                       mixalf_in=1.0,
                      ):
    """Shared scientific model, executed directly or compiled by Numba."""
            
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
    
    root_ex = np.zeros(len(s))
    
    An = np.zeros(len(s))
    An_tot = np.zeros(len(s))
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
    EW = np.zeros((1, len(s)))
    min_st = np.zeros((1, 6))

    d = np.zeros((1, len(s)))
    delta_d = np.zeros((1, len(s)))
    lamb = np.zeros((1, len(s)))
    SSA = np.zeros((1, len(s)))
    psd = np.zeros((1, len(s)))
    psd_rock_num = np.zeros((1, len(s)))

    wet_f = np.zeros(len(s))
     
    number_min = len(mineral) if M_rock_in > 0 else 1
    rho_min = np.zeros(number_min)
    MM_min = np.zeros(number_min)
    E_H = np.zeros(number_min)
    E_w = np.zeros(number_min)
    E_OH = np.zeros(number_min)
    n_H = np.zeros(number_min)
    n_OH = np.zeros(number_min)
    K_sp = np.zeros(number_min)
    min_st = np.zeros((number_min, 6))
    k_diss_H = np.zeros(number_min)
    k_diss_w = np.zeros(number_min)
    k_diss_OH = np.zeros(number_min)

    k_H_T = np.zeros((number_min, len(s)))
    k_w_T = np.zeros((number_min, len(s)))
    k_OH_T = np.zeros((number_min, len(s)))
    Omega =  np.zeros((number_min, len(s)))
    M_min = np.zeros((number_min, len(s)))
    rock_f = np.zeros((number_min, len(s)))
    Wr = np.zeros((number_min, len(s)))
    EW = np.zeros((number_min, len(s)))

    n_d_cl = len(d_in) if M_rock_in > 0 else 1
    if n_d_cl > 1:
        d = np.zeros((n_d_cl, len(s)))
        delta_d = np.zeros((n_d_cl, len(s)))
        lamb = np.zeros((n_d_cl, len(s)))
        SSA = np.zeros((n_d_cl, len(s)))
        psd = np.zeros((n_d_cl, len(s)))
        psd_rock_num = np.zeros((n_d_cl, len(s)))

    errors = np.zeros((16, len(s)))
    
#------------------------------------------------------------------------------
    # Constants
    CO2_atm = smew.CO2_atm(conv_mol) # [mol_CO2/l_air] Atmospheric CO2 concentration
    T_K = temp_soil + 273.15

    #frozen soil
    frozen = temp_soil <= 0.0
    
    # soil CO2 diffusivity 
    D_0 = smew.D_0() #free-air diffusion [m2/d]
    D = D_0*(1-s)**(10/3)*n**(4/3) #Mill-Quirk (1961)
    
    #solute diffusivity in soil water
    Dw_0 = smew.Dw_0()
    Dw = Dw_0*(n*s)**2 # Archie 1942, Grathwohl 1998 (book)
    
    # [g/mol-conv]: Molar masses    
    MM_Mg, MM_Ca, MM_Na, MM_K, MM_Si, MM_C, MM_Anions, MM_Al=smew.MM(conv_mol)
    
    # Aluminium speciation
    K1, K2, K3, K4 = smew.K_Al(conv_mol)
    
    # carbonate spec  
    k1, k2, k_w, k_H = smew.K_C(T_K,conv_mol)
    
    #CEC Gaines-Thomas constants
    K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H  = K_CEC
    
    #nutrient uptake by plants
    v_f_Ca, v_f_Mg, v_f_K, v_f_Si = smew.plant_nutr_f()
    dry_perc = 0.1 #percent of dry mass
    xi = dry_perc*np.array((v_f_Ca/MM_Ca, v_f_Mg/MM_Mg, v_f_K/MM_K, v_f_Si/MM_Si)) # [mol-conv/g_biomass]
    
    #carb weathering constants
    K_CaCO3,K_MgCO3,r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3 = smew.carb_weath_const(conv_mol)
    
    #mineral constants
    if M_rock_in > 0: 
        for j in range(0,number_min):
            MM_min[j], k_diss_H[j], k_diss_w[j], k_diss_OH[j], n_H[j], n_OH[j], E_H[j], E_w[j], E_OH[j], stoichiometry, K_sp[j] = smew.min_const(mineral[j], conv_mol)
            for element in range(6):
                min_st[j,element] = stoichiometry[element]
            #temperature scaling
            k_H_T[j,:] = k_diss_H[j]*np.exp(-E_H[j]*1000/(8.314/conv_mol)*(1/T_K[:]-1/(25+273.15)))
            k_w_T[j,:] = k_diss_w[j]*np.exp(-E_w[j]*1000/(8.314/conv_mol)*(1/T_K[:]-1/(25+273.15)))
            k_OH_T[j,:] = k_diss_OH[j]*np.exp(-E_OH[j]*1000/(8.314/conv_mol)*(1/T_K[:]-1/(25+273.15)))
    
    #rock density
    if rho_rock_in is None:
        rho_rock = 3e6 # [g/m3]
    else:
        rho_rock = rho_rock_in
    
    #rock surface fractality (Beerling 2020)
    b = 0.35 #[-]
    a = (1/(2*1e-10))**b #[1/m^b]
    
#------------------------------------------------------------------------------
    # RAINWATER
    
    Alk_rain = 0 #alk
    CO2_w_rain = k_H*CO2_atm # [mol/l] Henry's law

    # Assign variable once to be reused for solver inputs and outputs
    scalar_guess = np.empty(1)
    scalar_parameters = np.empty(5)
    scalar_residual = np.empty(1)
    scalar_work = np.empty(8)
    main_residual = np.empty(16)
    main_work = np.empty(488)
    solver_status_counts = np.zeros(6, dtype=np.int64)
    for i in range(len(s)):
        scalar_guess[0] = 10**-6*conv_mol
        scalar_parameters[0] = Alk_rain
        scalar_parameters[1] = k1[i]
        scalar_parameters[2] = k2[i]
        scalar_parameters[3] = CO2_w_rain[i]
        scalar_parameters[4] = k_w[i]
        status = _solve_system(WATER_SYSTEM, scalar_guess, scalar_parameters,
                               scalar_residual, scalar_work, 1.49012e-8)  # 1.49012e-8 is the default fsolve xtol
        solver_status_counts[status] += 1
        H_rain[i] = scalar_guess[0]
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
    CO2_air[0] = (r_het[0]+r_aut[0])/(D[0]*1000/(Z_CO2))+CO2_atm #mol-conv/l_air (Fs = resp_het + resp_aut), assumption of no leaching
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
    Ca[0], Mg[0], K[0], Na[0], Al_w[0] = conc_in
           
    # anions (mol_c/l)
    An[0] = 2*Mg[0]+2*Ca[0]+Na[0]+K[0]-Alk[0] #[mol_c/l]
    
    if An[0]<0:
        print(An[0])
        raise ValueError("Not enough cations for this alkalinity")
        
    # aluminum speciation
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
    f_Ca[0], f_Mg[0], f_K[0], f_Na[0], f_Al[0], f_H[0] = f_CEC_in
    
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
    
    #Carbonate minerals (added to the soil)
    CaCO3[0] = CaCO3_in # [mol-conv]
    MgCO3[0] = MgCO3_in

    #Carbonate weathering
    Omega_CaCO3[0] = Ca[0]*CO3[0]/K_CaCO3 # [-]
    Omega_MgCO3[0] = Mg[0]*CO3[0]/K_MgCO3
    if frozen[0]:
        W_CaCO3[0] = 0.0
        W_MgCO3[0] = 0.0
    else:
        W_CaCO3[0], W_MgCO3[0] = smew.carb_W(CaCO3[0], MgCO3[0], Omega_CaCO3[0], Omega_MgCO3[0], s[0], Zr, r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3) # [mol-conv/ m2 d]
        
    #Silicate weathering
    # We define them here as they are returned
    tt_app = 0
    M_iner = 0.0
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
        delta_d[0,tt_app] = d[0,tt_app]
        delta_d[1:,tt_app] = d[1:,tt_app] - d[:-1,tt_app]
        
        #particle size distribution by mass and number
        psd[:,tt_app] = psd_perc_in*M_rock[tt_app]/delta_d[:,tt_app] #[g/m]
        psd_rock_num[:,tt_app] = smew.psd_number_from_mass(psd[:,tt_app], d[:,tt_app], rho_rock)
        
        #refinement of fractal constant based on measured SSA 
        if SSA_in > 0:
            a = (SSA_in*rho_rock*M_rock[tt_app]/6)/np.sum(d[:,tt_app]**(b-1)*psd[:,tt_app]*delta_d[:,tt_app]) #[m**-b]
        
        #surface area
        lamb[:,tt_app] = a*d[:,tt_app]**b #[-]
        SSA[:,tt_app] = 6/(d[:,tt_app]*rho_rock)*lamb[:,tt_app] # [m2/g]
        SA[tt_app] = np.sum(SSA[:,tt_app]*psd[:,tt_app]*delta_d[:,tt_app]) #[m2]

        #wet surface area fraction
        if frozen[tt_app]:
            wet_f[tt_app] = 0.0
        else:
            wet_f[tt_app] = smew.wetness_SA(s[tt_app],keyword_ssa, pore_d_in, pore_pdf_in,  d[:,tt_app], psd_rock_num[:,tt_app], mixalf_in, d[-1,tt_app])
                    
        #mineral weathering
        if t_app == 0 and not frozen[0]:
            for j in range(0, number_min):
                #saturation state [-]
                Omega[j,0] = smew.sil_Omega(mineral[j], Ca[0], Mg[0], K[0], Na[0], Al[0], AlOH4[0], Si[0], H[0], K_sp[j], conv_mol,conv_Al)
                #weathering rate [mol-conv/ m2 d]                
                Wr[j,0] = smew.sil_Wr(mineral[j], Omega[j,0], H[0], k_H_T[j,0], k_w_T[j,0],k_OH_T[j,0], n_H[j], n_OH[j], diss_f,  conv_mol) 
                #weathering flux [mol-conv/d] 
                EW[j,0] = Wr[j,0]*SA[0]*rock_f[j,0]*wet_f[0]

#------------------------------------------------------------------------------
    #frozen option
    frozen_state = (pH, H, f_H, Ca_tot, Ca, f_Ca, Mg_tot, Mg, f_Mg, K_tot, K, f_K, Na_tot, Na, f_Na, Si_tot, Si, An_tot, An,
    Alk_tot, Alk, R_alk, CO2_w, HCO3, CO3, DIC, Al_tot, Al_w, Al, AlOH, AlOH2, AlOH3, AlOH4, f_Al, CaCO3, MgCO3, Omega_CaCO3, Omega_MgCO3)

    frozen_zero = (UP_Ca, UP_Mg, UP_K, UP_Si, W_CaCO3, W_MgCO3, ADV, root_ex, wet_f)
    
    # declared even if M_rock_in == 0 as returned
    frozen_rock_state = (d, delta_d, lamb, SSA, psd, psd_rock_num, M_min, rock_f, Omega)
    
#------------------------------------------------------------------------------
    #SYSTEM RESOLUTION

    for i in range(1, len(s)): 

            # Frozen soil: aqueous chemistry and reactions pause, only gas-phase CO2 diffusion
            if frozen[i]:

                for arr in frozen_state:
                    arr[i] = arr[i - 1]


                # Gas-phase CO2 relaxation toward atmospheric CO2
                air_vol = n * Zr * (1.0 - s[i]) * 1000.0  # [L_air m-2]
                k_diff = (D[i] / Z_CO2 * 1000.0) / air_vol  # [d-1]
                CO2_air[i] = CO2_atm + (CO2_air[i - 1] - CO2_atm) * np.exp(-k_diff * dt)
                Fs[i] = air_vol * (CO2_air[i - 1] - CO2_air[i]) / dt
                IC_tot[i] = (DIC[i] * s[i] + CO2_air[i] * (1.0 - s[i])) * (n * Zr * 1000.0)

                # stop biological, hydrological, and reaction fluxes
                for arr in frozen_zero:
                    arr[i] = 0.0

                if M_rock_in > 0:
                    if i > tt_app:
                        M_rock[i] = M_rock[i - 1]
                        SA[i] = SA[i - 1]

                    for arr in frozen_rock_state:
                        arr[:, i] = arr[:, i - 1]

                    Wr[:, i] = 0.0
                    EW[:, i] = 0.0
                    wet_f[i] = 0.0

                continue
        
            #CO2 advection due to moisture variation    
            if s[i]<s[i-1]:
                ADV[i] = n*Zr*1000*(s[i]-s[i-1])*CO2_atm # [mol] 
            elif s[i]>s[i-1]:
                ADV[i] = n*Zr*1000*(s[i]-s[i-1])*CO2_air[i-1]

            #active uptake [Ca, Mg, K, Si] 
            UP_act = smew.up_act(v[i], (v[i]-v[i-1]), xi, dt, T[i-1], Ca[i-1], Mg[i-1], K[i-1], Si[i-1], Dw[i-1], Zr, k_v, RAI, root_d)
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

            scalar_guess[0] = H[i-1]
            scalar_parameters[0] = k1[i]
            scalar_parameters[1] = k2[i]
            scalar_parameters[2] = CO2_w0
            scalar_parameters[3] = k_w[i]
            scalar_parameters[4] = Alk0
            status = _solve_system(HYDROGEN_SYSTEM, scalar_guess, scalar_parameters,
                                   scalar_residual, scalar_work, 1.49012e-8)
            solver_status_counts[status] += 1
            H0_2 = scalar_guess[0]

            #solution 1
            x0 = np.array((Alk0, CO2_w0, H0, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al[i-1],f_Mg[i-1], f_Na[i-1], f_K[i-1], f_H[i-1], f_Ca[i-1]))
            parameters = np.array((
                Alk_tot[i], n, Zr, s[i], IC_tot[i], k1[i], k2[i], k_H[i],
                k_w[i], CEC_tot, conv_Al, Al_tot[i], K1, K2, K3, K4,
                Mg_tot[i], Ca_tot[i], Na_tot[i], K_tot[i], K_Ca_Al,
                K_Ca_Mg, K_Ca_Na, K_Ca_K, K_Ca_H,
            ))
            status = _solve_system(BIOGEOCHEM_SYSTEM, x0, parameters, main_residual, main_work, 1e-12)
            solver_status_counts[status] += 1
            sol = x0
            errors[:, i] = main_residual

            #solution 2
            res_threshold = 1e-1
            if np.any(np.abs(errors[:,i]) > res_threshold):
                x0 = np.array((Alk0, CO2_w0, H0_2, R_alk0, Al_w0, Al0, Mg0, Ca0, Na0, K0, f_Al[i-1],f_Mg[i-1], f_Na[i-1], f_K[i-1], f_H[i-1], f_Ca[i-1]))
                status = _solve_system(BIOGEOCHEM_SYSTEM, x0, parameters, main_residual, main_work, 1e-14)
                solver_status_counts[status] += 1
                sol = x0
                errors[:, i] = main_residual
                if np.any(np.abs(errors[:,i]) > res_threshold):
                    print(i)
                    raise ValueError("Solution not converging")          

            Alk[i], CO2_w[i], H[i], R_alk[i], Al_w[i], Al[i], Mg[i], Ca[i], Na[i], K[i], f_Al[i], f_Mg[i], f_Na[i], f_K[i], f_H[i], f_Ca[i] = sol

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
            W_CaCO3[i], W_MgCO3[i] = smew.carb_W(CaCO3[i], MgCO3[i], Omega_CaCO3[i], Omega_MgCO3[i], s[i], Zr, r_CaCO3,r_MgCO3,tau_CaCO3,tau_MgCO3)

            #Silicate weathering
            if M_rock_in > 0:

                #saturation and weathering rate
                for j in range(0, number_min):
                    Omega[j,i] = smew.sil_Omega(mineral[j], Ca[i], Mg[i], K[i], Na[i], Al[i], AlOH4[i], Si[i], H[i], K_sp[j], conv_mol,conv_Al) #[-]           
                    Wr[j,i]= smew.sil_Wr(mineral[j], Omega[j,i], H[i], k_H_T[j,i], k_w_T[j,i],k_OH_T[j,i], n_H[j], n_OH[j], diss_f,  conv_mol)

                #post application only
                if i> tt_app:

                    #mineral fractions in rock
                    M_min[:,i] = M_min[:,i-1]-EW[:,i-1]*MM_min[:]*dt # [g]
                    M_min[:, i] = np.maximum(M_min[:, i], 0)
                    M_rock[i] = np.sum(M_min[:,i]) + M_iner # [g]
                    if M_rock[i]>0:
                        rock_f[:,i] = M_min[:,i]/M_rock[i] # [-]

                    #diameter variation
                    d_shrink = np.sum(rock_f[:,i-1]*Wr[:,i-1]*MM_min[:]/rho_rock)*dt # [m]
                    d[:,i] = d[:,i-1] - 2*d_shrink*lamb[:,i-1] # [m]
                    d[:,i][d[:,i] < 0] = 0
                    delta_d[0,i] = d[0,i]
                    delta_d[1:,i] = d[1:,i] - d[:-1,i] # [m]
                    lamb[:,i], SSA[:,i], psd[:,i], SA[i] = smew.psd_evol(d[:,i], delta_d[:,i], d[:,i-1], delta_d[:,i-1], psd[:,i-1], n_d_cl, a, b, rho_rock)
                    psd_rock_num[:,i] = smew.psd_number_from_mass(psd[:,i], d[:,i], rho_rock)

                #wetness scaling of the surface area
                wet_f[i] = smew.wetness_SA(s[i],keyword_ssa, pore_d_in, pore_pdf_in, d[:,i], psd_rock_num[:, i], mixalf_in, d[-1,i])

                #weathering fluxes         
                EW[:,i] = Wr[:,i]*SA[i]*rock_f[:,i]*wet_f[i] # [mol/d]



    # Named pairs support mixed output types in Numba; the wrapper builds the dictionary.
    return (
        ("pH", pH),
        ("H", H),
        ("f_H", f_H),
        ("Ca_tot", Ca_tot),
        ("Ca", Ca),
        ("f_Ca", f_Ca),
        ("UP_Ca", UP_Ca),
        ("Omega_CaCO3", Omega_CaCO3),
        ("CaCO3", CaCO3),
        ("W_CaCO3", W_CaCO3),
        ("Mg_tot", Mg_tot),
        ("Mg", Mg),
        ("f_Mg", f_Mg),
        ("UP_Mg", UP_Mg),
        ("Omega_MgCO3", Omega_MgCO3),
        ("MgCO3", MgCO3),
        ("W_MgCO3", W_MgCO3),
        ("K_tot", K_tot),
        ("K", K),
        ("f_K", f_K),
        ("UP_K", UP_K),
        ("Na_tot", Na_tot),
        ("Na", Na),
        ("f_Na", f_Na),
        ("Si", Si),
        ("Si_tot", Si_tot),
        ("UP_Si", UP_Si),
        ("root_ex", root_ex),
        ("An", An),
        ("An_tot", An_tot),
        ("R_alk", R_alk),
        ("Alk_tot", Alk_tot),
        ("Alk", Alk),
        ("CO2_air", CO2_air),
        ("IC_tot", IC_tot),
        ("CO2_w", CO2_w),
        ("HCO3", HCO3),
        ("CO3", CO3),
        ("DIC", DIC),
        ("Fs", Fs),
        ("ADV", ADV),
        ("H_rain", H_rain),
        ("DIC_rain", DIC_rain),
        ("Al", Al),
        ("f_Al", f_Al),
        ("AlOH", AlOH),
        ("AlOH2", AlOH2),
        ("AlOH3", AlOH3),
        ("AlOH4", AlOH4),
        ("Al_w", Al_w),
        ("Al_tot", Al_tot),
        ("M_rock", M_rock),
        ("SA", SA),
        ("EW", EW),
        ("min_st", min_st),
        ("d", d),
        ("delta_d", delta_d),
        ("lamb", lamb),
        ("SSA", SSA),
        ("psd", psd),
        ("psd_rock_num", psd_rock_num),
        ("wet_f", wet_f),
        ("number_min", number_min),
        ("rho_min", rho_min),
        ("MM_min", MM_min),
        ("E_H", E_H),
        ("E_w", E_w),
        ("E_OH", E_OH),
        ("n_H", n_H),
        ("n_OH", n_OH),
        ("K_sp", K_sp),
        ("k_diss_H", k_diss_H),
        ("k_diss_w", k_diss_w),
        ("k_diss_OH", k_diss_OH),
        ("k_H_T", k_H_T),
        ("k_w_T", k_w_T),
        ("k_OH_T", k_OH_T),
        ("Omega", Omega),
        ("M_min", M_min),
        ("rock_f", rock_f),
        ("Wr", Wr),
        ("n_d_cl", n_d_cl),
        ("errors", errors),
        ("CO2_atm", CO2_atm),
        ("T_K", T_K),
        ("frozen", frozen),
        ("D_0", D_0),
        ("D", D),
        ("Dw_0", Dw_0),
        ("Dw", Dw),
        ("MM_Mg", MM_Mg),
        ("MM_Ca", MM_Ca),
        ("MM_Na", MM_Na),
        ("MM_K", MM_K),
        ("MM_Si", MM_Si),
        ("MM_C", MM_C),
        ("MM_Anions", MM_Anions),
        ("MM_Al", MM_Al),
        ("K1", K1),
        ("K2", K2),
        ("K3", K3),
        ("K4", K4),
        ("k1", k1),
        ("k2", k2),
        ("k_w", k_w),
        ("k_H", k_H),
        ("K_Ca_Mg", K_Ca_Mg),
        ("K_Ca_K", K_Ca_K),
        ("K_Ca_Na", K_Ca_Na),
        ("K_Ca_Al", K_Ca_Al),
        ("K_Ca_H", K_Ca_H),
        ("v_f_Ca", v_f_Ca),
        ("v_f_Mg", v_f_Mg),
        ("v_f_K", v_f_K),
        ("v_f_Si", v_f_Si),
        ("dry_perc", dry_perc),
        ("xi", xi),
        ("K_CaCO3", K_CaCO3),
        ("K_MgCO3", K_MgCO3),
        ("r_CaCO3", r_CaCO3),
        ("r_MgCO3", r_MgCO3),
        ("tau_CaCO3", tau_CaCO3),
        ("tau_MgCO3", tau_MgCO3),
        ("rho_rock", rho_rock),
        ("b", b),
        ("a", a),
        ("Alk_rain", Alk_rain),
        ("CO2_w_rain", CO2_w_rain),
        ("solver_status_counts", solver_status_counts),
        ("Z_CO2", Z_CO2),
        ("I_An", I_An),
        ("I_Ca", I_Ca),
        ("I_Mg", I_Mg),
        ("I_Na", I_Na),
        ("I_K", I_K),
        ("I_Si", I_Si),
        ("tt_app", tt_app),
        ("M_iner", M_iner),
        ("frozen_state", frozen_state),
        ("frozen_zero", frozen_zero),
        ("frozen_rock_state", frozen_rock_state),
    )
