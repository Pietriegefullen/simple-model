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
from refactored import load_matlab  
from scipy.signal import savgol_filter
#from SimpleIN import Cdec # importiert das eigentliche mathematische Model 



os.chdir('C:/Users/Lara/Desktop/simple model')
         
            
superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam = load_matlab()

print(superdata_Kuru.keys())
for replica in superdata_Kuru.keys():
    print(replica)
    plt.figure()
    plt.yscale("log")
    plt.plot(superdata_Kuru[replica]['CH4'], 'g')
    plt.legend([replica])
    
for replica in superdata_Sam.keys():
    plt.figure()
    plt.yscale("log")
    plt.plot(superdata_Sam[replica]['CH4'], 'm')    
    plt.legend([replica])


flachsteil = [13510,13511,13512,13670,13671, 13672, 13691,13692,13700, 13722,13731, 13750, 13751, 13752]

steil =[ 13690, 13531, 13530, 13521, 13520, 13742,13741, 13740, 13732, 13731,13730,13730,13721,13720]

steilsteil= [13782, 13781, 13780, 13772, 13771, 13770]

#%%
plt.plot(superdata['13510']['CO2'],'b')
plt.plot(superdata['13510']['CH4'],'r')

#%%

y =superdata['13690']['CH4']
x =superdata['13690']['measured_time']

yhat = savgol_filter(y, 11, 3)
plt.plot(x,y)
plt.plot(x,yhat, color='red')
plt.show()







