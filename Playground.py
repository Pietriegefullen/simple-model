# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 08:45:18 2021

@author: Lara
"""
import pandas as pd
import os
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
#import numpy as np

#Importieren benutzter Funktionen und Daten
#from Realdataall import load_realdata
from Readalldata import load_realdata
from SimpleOUT import simplefun, simplesolve
from SimpleOUT import optifun
from Babypascal import Mol_nach_Pa

plt.close('all')

#Realdata9 = load_realdata(24)

Data1and2and3and4and5and6and7and8and9 =[0,3,6,9,12,15,18,21,24]
#





for m in Data1and2and3and4and5and6and7and8and9:
    
    
    
    Realdata = load_realdata(m)
    
    t = Realdata[:,0]
    CH4 = Realdata[:,1]
    
    fig, axs = plt.subplots(2)

    #plt.figure()
    axs[0].plot(t,CH4, 'ro')
    #plt.show()
    
   # plt.figure()
    axs[1].plot(t,CH4, "r-")
    plt.yscale('log')
    plt.ylim([1e-8, 5e1]) 
    
    INDsmallremoved = np.where(CH4>1e-4)
    VALnewlistCH4 = CH4[INDsmallremoved]
    VALnewlistT = t[INDsmallremoved]
    #plt.plot(t[INDsmallremoved],VALnewlistCH4, 'red')

    INDlargeremoved = np.where(VALnewlistCH4 < (max(VALnewlistCH4)/2))

    axs[1].plot(VALnewlistT[INDlargeremoved],VALnewlistCH4[INDlargeremoved], 'blue')

    p = np.polyfit(VALnewlistT[INDlargeremoved], np.log(VALnewlistCH4[INDlargeremoved]), 1)
    
    #plt.plot(t[INDlargeremoved],CH4[INDlargeremoved], 'red')
    
    m = p[0]
    b = p[1]
    print(m,b)
    axs[1].plot(t, np.exp(m*t + b))

    plt.show()
    
#%%



























#for m in Data1and2and3and4and5and6and7and8and9:
#    
#    
#    
#    Realdata = load_realdata(m)
#    
#    t = Realdata[:,0]
#    CH4 = Realdata[:,1]
#    
#
#    fig, axs = plt.subplots(2)
#    
#    axs[0].plot(t,CH4)
#    #plt.show()
#    
#    #plt.figure()
#    axs[1].plot(t,CH4)
#    plt.yscale('log')
#    plt.ylim([1e-8, 5e1])
#    
# 
#    
#    m, b = np.polyfit(t, CH4, 1)
#    print(m,b)
#    axs[1].plot(t, m*t + b)
#    plt.yscale('log')
#    plt.ylim([1e-8, 5e1])
#    plt.show()
#    
    
    
 


#%% 
      
for m in Data1and2and3and4and5and6and7and8and9:
    
    
    
    Realdata = load_realdata(m)
    
    t = Realdata[:,0]
    CO2 = Realdata[:,2]
    

    fig, axs = plt.subplots(2)
    
    axs[0].plot(t,CO2)
    #plt.show()
    
    #plt.figure()
    axs[1].plot(t,CO2)
    plt.yscale('log')
    plt.ylim([1e-8, 5e1])
    plt.show()
    
    
#%%

import pandas as pd
import os
df = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Gibbs.xlsx')
    
    
    
    
    
    