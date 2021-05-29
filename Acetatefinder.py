# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:34:12 2021

@author: Lara
"""

import copy
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from scipy.integrate import odeint, solve_ivp
from scipy import stats
from scipy.optimize import curve_fit
import scipy.io as sio
import numpy as np
from Importmatlab import load_matlab  
from scipy.signal import savgol_filter
#from SimpleIN import Cdec # importiert das eigentliche mathematische Model 



os.chdir('C:/Users/Lara/Desktop/simple model')
         
            
superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe, Rep_ohne_Fe,superdata_mit_Fe, Rep_mit_Fe = load_matlab()

print(superdata_Kuru.keys())
for replica in superdata_Kuru.keys():
    fig, axs = plt.subplots(3)
    #plt.yscale("log")
    First_Carex= int(superdata_Kuru[replica]['First_Carex_index'])
    axs[0].plot(superdata_Kuru[replica]['measured_time'][: First_Carex], superdata_Kuru[replica]['CH4'][: First_Carex], 'g')
    axs[0].plot(superdata_Kuru[replica]['measured_time'][ First_Carex:], superdata_Kuru[replica]['CH4'][ First_Carex:], 'r')
    green = axs[0].plot(superdata_Kuru[replica]['measured_time'][: First_Carex], superdata_Kuru[replica]['CH4'][: First_Carex], 'g')
    red= axs[0].plot(superdata_Kuru[replica]['measured_time'][ First_Carex:], superdata_Kuru[replica]['CH4'][ First_Carex:], 'r')
    axs[0].set_ylabel("CH4")
   
    #plt.figure()
    axs[1].plot(superdata_Kuru[replica]['measured_time'][: First_Carex], superdata_Kuru[replica]['CH4'][: First_Carex], 'g')
    axs[1].plot(superdata_Kuru[replica]['measured_time'][ First_Carex:], superdata_Kuru[replica]['CH4'][ First_Carex:], 'r')
    axs[1].set_yscale("log")
    axs[1].set_ylabel("Logplot of CH4", fontsize=7)

    
    axs[2].plot(superdata_Kuru[replica]['measured_time'][: First_Carex], superdata_Kuru[replica]['CH4'][: First_Carex], 'g')
    axs[2].set_yscale("log")
    axs[2].set_ylabel("Logplot of CH4", fontsize=7)
    
    plt.suptitle("Kurunakh" + '___' + replica, fontsize=14) 
    line_labels = ['Pre Carex', 'post Carex']
    
    fig.legend([green, red],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="upper right",   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box
           )
    save_path = os.path.join('C:/Users/Lara/Desktop/simple model/Figs/Kurunakh',replica +'.png')
    plt.savefig(save_path)
 #=========================  
for replica in superdata_Sam.keys():
     
    fig, axs = plt.subplots(3,1)
    #plt.yscale("log")
    First_Carex= int(superdata_Sam[replica]['First_Carex_index'])
    axs[0].plot(superdata_Sam[replica]['measured_time'][: First_Carex], superdata_Sam[replica]['CH4'][: First_Carex], 'm')
    axs[0].plot(superdata_Sam[replica]['measured_time'][ First_Carex:], superdata_Sam[replica]['CH4'][ First_Carex:], 'r')
    pink = axs[0].plot(superdata_Sam[replica]['measured_time'][: First_Carex], superdata_Sam[replica]['CH4'][: First_Carex], 'm')
    red= axs[0].plot(superdata_Sam[replica]['measured_time'][ First_Carex:], superdata_Sam[replica]['CH4'][ First_Carex:], 'r')
    axs[0].set_ylabel("CH4")
    plt.legend([replica])
    #plt.figure()
    #plt.yscale("log")
    axs[1].plot(superdata_Sam[replica]['measured_time'][: First_Carex], superdata_Sam[replica]['CH4'][: First_Carex], 'm')
    axs[1].plot(superdata_Sam[replica]['measured_time'][ First_Carex:], superdata_Sam[replica]['CH4'][ First_Carex:], 'r')
    axs[1].set_yscale("log")
    axs[1].set_ylabel("Logplot of CH4", fontsize=7)
    
 #=============================================================================
    #plt.figure
    #plt.yscale("log")
    axs[2].plot(superdata_Sam[replica]['measured_time'][: First_Carex], superdata_Sam[replica]['CH4'][: First_Carex], 'm')
    axs[2].set_yscale("log")
    axs[2].set_ylabel("Logplot of CH4", fontsize=7)
    plt.suptitle("Samoylov" + '___' + replica, fontsize=14)  
    fig.legend([pink, red],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="upper right",   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box
           )
    save_path = os.path.join('C:/Users/Lara/Desktop/simple model/Figs/Samoylov', replica +'.png')
    plt.savefig(save_path)
 #=============================================================================

flachsteil = [13510,13511,13512,13670,13671, 13672, 13691,13692,13700, 13722,13731, 13750, 13751, 13752]

steil =[ 13690, 13531, 13530, 13521, 13520, 13742,13741, 13740, 13732, 13731,13730,13730,13721,13720]

steilsteil= [13782, 13781, 13780, 13772, 13771, 13770]

# Proben 13520, 13691 haben keine Carexdaten !!!!
# Probe 13701 ist iwie seltsam

#%% smoothin with gauss
from scipy.ndimage import gaussian_filter1d
for key in superdata_2021_all.keys():
    smooth = gaussian_filter1d(superdata_2021_all[key]['CH4'], 2)
    superdata_2021_all[key]['CH4_smooth'] = smooth   
    plt.figure()
    plt.plot(superdata_2021_all[key]['CH4_smooth'])
    #plt.plot(superdata_2021_all[key]['CH4'], 'r--')
    plt.legend([key])
#%%


for key in superdata_Kuru.keys():
    smooth = gaussian_filter1d(superdata_Kuru[key]['CH4'], 2)
    superdata_Kuru[key]['CH4_smooth'] = smooth   
    plt.figure()
    plt.plot(superdata_Kuru[key]['CH4_smooth'])
    #plt.plot(superdata_2021_all[key]['CH4'], 'r--')
    plt.legend([key])

#%%

for key in superdata_Sam.keys():
    smooth = gaussian_filter1d(superdata_Sam[key]['CH4'], 2)
    superdata_Sam[key]['CH4_smooth'] = smooth   
    plt.figure()
    plt.plot(superdata_Sam[key]['CH4_smooth'], 'm')
    #plt.plot(superdata_2021_all[key]['CH4'], 'r--')
    plt.legend([key])
    
#%%    

# compute second derivative
smooth_d2 = np.gradient(np.gradient(smooth))

# find switching points
infls = np.where(np.diff(np.sign(smooth_d2)))[0]



