# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import numpy as np
import pyEW
from scipy import optimize
from scipy.optimize import least_squares, fsolve, minimize, newton_krylov, broyden1, root, broyden2

#------------------------------------------------------------------------------
 # conc to CEC fractions
    
def conc_to_f_CEC(conc_in,pH_in,soil,conv_mol,conv_Al):
            
    #constants 
    K_CEC = pyEW.K_GT_CEC(soil,conv_mol) #CEC Gaines-Thomas
    [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H] = K_CEC
    [K1, K2, K3, K4] = pyEW.K_Al(conv_mol) #Al speciation
    
    # cations (mol-conv/l)
    [Ca, Mg, K, Na, Al_w] = conc_in
    H = 10**(-pH_in)*conv_mol 
    
    # aluminium speciation
    Al=(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w #mol-conv/l
        
    # CEC saturation, G-T convenction
    def eqf_Ca(p): #CEC Calcium (solvability eq is sum of fractions=1)
        f_Ca = p
        return(1-(f_Ca+(Al/conv_Al)*((f_Ca**3/(K_Ca_Al*Ca**3))**(1/2))+Mg*(f_Ca/(K_Ca_Mg*Ca))+Na*((f_Ca/(K_Ca_Na*Ca))**(1/2))+K*((f_Ca/(K_Ca_K*Ca))**(1/2))+H*((f_Ca/(K_Ca_H*Ca))**(1/2))))
    
    #CEC fractions
    f_Ca = fsolve(eqf_Ca, 0.2) 
    f_Al = (Al/conv_Al)*((f_Ca**3/(K_Ca_Al*Ca**3))**(1/2))
    f_Mg = Mg*(f_Ca/(K_Ca_Mg*Ca))
    f_Na = Na*((f_Ca/(K_Ca_Na*Ca))**(1/2))
    f_K = K*((f_Ca/(K_Ca_K*Ca))**(1/2))
    f_H = H*((f_Ca/(K_Ca_H*Ca))**(1/2))
    
    f_CEC_in = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H]
                                          
    return(f_CEC_in, K_CEC)

#------------------------------------------------------------------------------
 # CEC fractions to conc 
    
def f_CEC_to_conc(f_CEC_in, pH_in, soil, conv_mol,conv_Al):
            
    #pH
    H = 10**(-pH_in)*conv_mol 
    
    #constants 
    K_CEC = pyEW.K_GT_CEC(soil,conv_mol) #CEC Gaines-Thomas
    [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H] = K_CEC
    [K1, K2, K3, K4] = pyEW.K_Al(conv_mol) #Al speciation

    #f_CEC [-]
    [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H] = f_CEC_in
       
    #estimates of concentrations
    Ca = (f_Ca/K_Ca_H)*(H/f_H)**2
    Mg = (Ca/f_Ca)*K_Ca_Mg*f_Mg
    K = ((Ca/f_Ca)*K_Ca_K)**(1/2)*f_K
    Na = ((Ca/f_Ca)*K_Ca_Na)**(1/2)*f_Na
    Al = ((Ca/f_Ca)**3*K_Ca_Al)**(1/2)*f_Al*conv_Al
    
    Al_w=Al/(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))   
         
    conc_in = [Ca, Mg, K, Na, Al_w]
    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]
        
    return(conc_in, K_CEC)

#------------------------------------------------------------------------------
 # Input: Total (Ca, Mg, K, Na) and CEC base saturation (or acid saturation, f_H+f_Al)
    
def total_to_f_CEC_and_conc(total_in, pH_in, f_acid, s, soil, n,Zr,CEC_tot,conv_mol,conv_Al):

    #constants 
    K_CEC = pyEW.K_GT_CEC(soil, conv_mol) #CEC Gaines-Thomas
    [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H] = K_CEC
    [K1, K2, K3, K4] = pyEW.K_Al(conv_mol) #Al speciation

    # total (mol-conv/m2)
    [Ca_tot, Mg_tot, K_tot, Na_tot] = total_in
    H = 10**(-pH_in)*conv_mol 

    def equations(p):
        Al_w, Al, Al_tot, Mg, Ca, Na, K, f_Mg, f_Na, f_K, f_Ca, f_Al, f_H = p

        return(Al_w*n*Zr*s[0]*1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot,\
               Al-(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w,\
               Mg*n*Zr*s[0]*1000+f_Mg/2*CEC_tot-Mg_tot,\
               Ca*n*Zr*s[0]*1000+f_Ca/2*CEC_tot-Ca_tot,\
               Na*n*Zr*s[0]*1000+f_Na*CEC_tot-Na_tot,\
               K*n*Zr*s[0]*1000+f_K*CEC_tot-K_tot,\
               f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2),\
               f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2),\
               f_H + f_Al - f_acid,\
               f_Mg - Mg*(f_Ca/(K_Ca_Mg*Ca)),\
               f_Na - Na*(f_Ca/(K_Ca_Na*Ca))**(1/2),\
               f_K - K*(f_Ca/(K_Ca_K*Ca))**(1/2),\
               1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H))   
    
    #initial guess
    f_Ca0 = 0.8*(1 - f_acid)
    f_Mg0 = 0.1*(1 - f_acid)
    f_K0 = 0.05*(1 - f_acid)
    f_Na0 = 0.05*(1 - f_acid)
    f_H0 = f_acid*0.3
    f_Al0 = f_acid*0.7
    Mg0 = (Mg_tot-f_Mg0/2*CEC_tot)/(n*Zr*s[0]*1000)
    Ca0 = (Ca_tot-f_Ca0/2*CEC_tot)/(n*Zr*s[0]*1000)
    Na0 = (Na_tot-f_Na0*CEC_tot)/(n*Zr*s[0]*1000)
    K0 = (K_tot-f_K0*CEC_tot)/(n*Zr*s[0]*1000)
    Al0 = (f_Al0*conv_Al)*((K_Ca_Al*Ca0**3)/f_Ca0**3)**(1/2)
    Al_w0 = Al0/(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))
    Al_tot0 = (f_Al0/3)*CEC_tot*conv_Al+Al_w0*n*Zr*s[0]*1000

    x0 = np.array([Al_w0, Al0, Al_tot0, Mg0, Ca0, Na0, K0, f_Mg0, f_Na0, f_K0, f_Ca0, f_Al0, f_H0])

    #system solution
    sol = fsolve(equations,x0, xtol=1e-14)
    [Al_w, Al, Al_tot, Mg, Ca, Na, K, f_Mg, f_Na, f_K, f_Ca, f_Al, f_H] = sol

    #K_Ca_Al = (Al/conv_Al/f_Al)**2*(f_Ca/Ca)**3
    #K_Ca_H = (f_Ca/Ca)*(H/f_H)**2
    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]
    conc_in = [Ca, Mg, K, Na, Al_w] 
    f_CEC_in = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H]

    return(conc_in, f_CEC_in, K_CEC)

#------------------------------------------------------------------------------
 # Calibration of K constants on coupled f_CEC and conc measurements
    
def f_CEC_and_conc_to_K(f_CEC_in, conc_in, pH_in, soil, conv_mol,conv_Al):
            
    #pH
    H = 10**(-pH_in)*conv_mol 
    
    #f_CEC [-]
    [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H] = f_CEC_in
    
    #conc [-]
    [Ca, Mg, K, Na, Al_w] = conc_in
    
    # Al spec
    [K1, K2, K3, K4] = pyEW.K_Al(conv_mol)
    Al=(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w
    
    #estimates of K constants
    K_Ca_Mg = (f_Ca/Ca)*(Mg/f_Mg)
    K_Ca_K = (f_Ca/Ca)*(K/f_K)**2
    K_Ca_Na = (f_Ca/Ca)*(Na/f_Na)**2
    K_Ca_Al = (f_Ca/Ca)**3*(f_Al/(Al/conv_Al))**2
    K_Ca_H = (f_Ca/Ca)*(H/f_H)**2 
    
    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]
        
    return(K_CEC)

#------------------------------------------------------------------------------
 # Input: Total (Ca, Mg, K, Na) and Al_w
    
def Kelland(total_in, pH_in, conc_in, s, soil, n,Zr,CEC_tot,conv_mol,conv_Al):

    #constants 
    K_CEC = pyEW.K_GT_CEC(soil, conv_mol) #CEC Gaines-Thomas
    [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H] = K_CEC
    [K1, K2, K3, K4] = pyEW.K_Al(conv_mol) #Al speciation

    # known (mol-conv/m2)
    [Ca_tot, Mg_tot, K_tot, Na_tot] = total_in
    [Ca, Mg, K, Na, Al_w] = conc_in
    H = 10**(-pH_in)*conv_mol
    Al=(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))*Al_w

    #CEC
    f_Mg = (Mg_tot - Mg*n*Zr*s[0]*1000)*2/CEC_tot
    f_K = (K_tot - K*n*Zr*s[0]*1000)*2/CEC_tot
    f_Na = (Na_tot - Na*n*Zr*s[0]*1000)/CEC_tot
        
    def equations(p):
        Al_tot, CaCO3,  f_Al, f_H, f_Ca = p

        return(Al_w*n*Zr*s[0]*1000+(f_Al/3)*CEC_tot*conv_Al-Al_tot,\
               Ca*n*Zr*s[0]*1000+f_Ca/2*CEC_tot+CaCO3-Ca_tot,\
               f_Al - (Al/conv_Al)*(f_Ca**3/(K_Ca_Al*Ca**3))**(1/2),\
               f_H - H*(f_Ca/(K_Ca_H*Ca))**(1/2),\
               1-(f_Ca+f_Al+f_Mg+f_Na+f_K+f_H))   
    
    #initial guess
    f_H0 = 1e-3
    f_Al0 = 1e-3
    f_Ca0 = 0.8*(1 - f_H0 + f_Al0)
    Al_tot0 = (f_Al0/3)*CEC_tot*conv_Al+Al_w*n*Zr*s[0]*1000
    CaCO30 = Ca_tot-Ca*n*Zr*s[0]*1000+f_Ca0/2*CEC_tot   

    x0 = np.array([Al_tot0, CaCO30, f_Al0, f_H0, f_Ca0])

    #system solution
    sol = fsolve(equations,x0, xtol=1e-14)
    [Al_tot, CaCO3,  f_Al, f_H, f_Ca] = sol

    #estimating soil-dependent K_CEC
    K_Ca_Mg = (f_Ca/Ca)*(Mg/f_Mg)
    K_Ca_K = (f_Ca/Ca)*(K/f_K)**2
    K_Ca_Na = (f_Ca/Ca)*(Na/f_Na)**2
    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]
    
    conc_in = [Ca, Mg, K, Na, Al_w] 
    f_CEC_in = [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H]

    return(conc_in, f_CEC_in, K_CEC, CaCO3)

#------------------------------------------------------------------------------
 # Amann et al., fractions and Mg conc 
    
def Amann(f_CEC_in, pH_in, Mg_in, soil, conv_mol,conv_Al):
            
    #pH
    H = 10**(-pH_in)*conv_mol 
    
    #constants 
    K_CEC = pyEW.K_GT_CEC(soil,conv_mol) #CEC Gaines-Thomas
    [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H] = K_CEC
    [K1, K2, K3, K4] = pyEW.K_Al(conv_mol) #Al speciation

    #f_CEC [-]
    [f_Ca, f_Mg, f_K, f_Na, f_Al, f_H] = f_CEC_in
       
    #conc from measurements
    Mg = Mg_in
        
    #estimates of concentrations
    Ca = (Mg/f_Mg)/K_Ca_Mg*f_Ca
    K = ((Ca/f_Ca)*K_Ca_K)**(1/2)*f_K
    Na = ((Ca/f_Ca)*K_Ca_Na)**(1/2)*f_Na
    Al = ((Ca/f_Ca)**3*K_Ca_Al)**(1/2)*f_Al*conv_Al
    
    Al_w=Al/(H**4/(H**4+H**3*K1+H**2*K1*K2+H*K1*K2*K3+K1*K2*K3*K4))
    
    #estimates of K_CEC
    K_Ca_H = (f_Ca/Ca)*(H/f_H)**2      
         
    conc_in = [Ca, Mg, K, Na, Al_w]
    K_CEC = [K_Ca_Mg, K_Ca_K, K_Ca_Na, K_Ca_Al, K_Ca_H]
        
    return(conc_in, K_CEC)

