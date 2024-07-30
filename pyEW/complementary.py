# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:34:44 2019
"""

import matplotlib.pyplot as plt
import numpy as np

#-----------------------------------------------------------------------------------
#fig CEC 
#-----------------------------------------------------------------------------------

def fig_CEC(t, f_Ca_list, f_Mg_list, f_K_list, f_Na_list, f_Al_list, f_H_list, titles):
    
    plots = [
        (f_H_list, 'f$_\mathrm{H}$'),
        ([f_H + f_Na for f_H, f_Na in zip(f_H_list, f_Na_list)], 'f$_\mathrm{Na}$'),
        ([f_H + f_Na + f_K for f_H, f_Na, f_K in zip(f_H_list, f_Na_list, f_K_list)], 'f$_\mathrm{K}$'),
        ([f_H + f_Na + f_K + f_Ca for f_H, f_Na, f_K, f_Ca in zip(f_H_list, f_Na_list, f_K_list, f_Ca_list)], 'f$_\mathrm{Ca}$'),
        ([f_H + f_Na + f_K + f_Ca + f_Mg for f_H, f_Na, f_K, f_Ca, f_Mg in zip(f_H_list, f_Na_list, f_K_list, f_Ca_list, f_Mg_list)], 'f$_\mathrm{Mg}$'),
        ([f_H + f_Na + f_K + f_Ca + f_Mg + f_Al for f_H, f_Na, f_K, f_Ca, f_Mg, f_Al in zip(f_H_list, f_Na_list, f_K_list, f_Ca_list, f_Mg_list, f_Al_list)], 'f$_\mathrm{Al}$')
    ]    
    
    if isinstance(f_H_list, list):
        
        num_panels = len(f_H_list)
        fig_out, axes = plt.subplots(1, num_panels, sharey=True, figsize=(4*num_panels, 3))
        
        for i, ax in enumerate(axes):
            data_pre = 0
            
            for data, label in plots:
                ax.plot(t, data[i], label=label)
                ax.fill_between(t, data[i], data_pre, alpha=0.5)
                data_pre = data[i]
            
            ax.legend().set_visible(False)  # Hide legends for individual plots
            ax.set_xlabel('t (days)')
            ax.set_title(titles[i])
         
        axes[-1].legend(loc='center right', bbox_to_anchor=(1.3, 0.5))  # Plot the legend at the end  
    
    else:
        num_panels = 1
        fig_out, ax = plt.subplots(figsize=(6, 4))
        data_pre = 0
        
        for data, label in plots:
            ax.plot(t, data, label=label)
            ax.fill_between(t, data, data_pre, alpha=0.5)
            data_pre = data
                    
        ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))
        ax.set_xlabel('t (days)')
        ax.set_title(titles)

   
    plt.tight_layout()

    return fig_out
    
#-----------------------------------------------------------------------------------
#fig IC
#-----------------------------------------------------------------------------------
def fig_IC(t, DIC_list, CO2_air_list,s,n,Zr):
    num_panels = len(DIC_list)
    fig_out, axes = plt.subplots(1, num_panels, sharey=True, figsize=(4*num_panels, 4))

    plots = [
            (DIC_list*s*(n*Zr*1000), 'IC_w'),
            ([(DIC*s+CO2_air*(1-s))*(n*Zr*1000) for DIC, CO2_air in zip(DIC_list, CO2_air_list)], 'IC_a')
        ]
    
    for i, ax in enumerate(axes):
            data_pre = 0
            for data, label in plots:

                ax.plot(t, data[i], label=label)
                ax.fill_between(t, data[i], data_pre, alpha=0.5)
                data_pre = data[i]                         

            ax.legend().set_visible(False)  # Hide legends for individual plots
            
    axes[-1].legend(loc='center right', bbox_to_anchor=(1.2, 0.5))  # Plot the legend at the end
    
    axes[0].set_xlabel('t')
    axes[0].set_ylabel('IC [micromol]')
    
    plt.tight_layout()
        
    return fig_out

plt.show()

#-----------------------------------------------------------------------------------
# moving average
#-----------------------------------------------------------------------------------

def mov_avg(data, window_size):
    data_avg = np.zeros(len(data))
    
    #same ic
    data_avg[0] = data[0]
    
    #left part
    for i in range(1, window_size): 
        data_avg[i] = np.mean(data[0:(i+window_size)])
    
    #central part
    for i in range(window_size, len(data)-window_size): 
        data_avg[i] = np.mean(data[(i-window_size):(i+window_size)])
    
    #right part
    for i in range(len(data)-window_size,len(data)): 
        data_avg[i] = np.mean(data[(i-window_size):len(data)])

    return data_avg

