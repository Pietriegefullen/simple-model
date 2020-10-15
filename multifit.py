# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:45:44 2020

@author: Lara
"""
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from Mockdata import Mockdata
from SimpleOUT import simplefun



# Fitting the parameters: 

xlist = list(range(0,101,1)) # Cpool, CH4 und CO2 sind je 100 Werte lang und
xdata = xlist + xlist + xlist # aneinandergehängt. Curvefit verlangt x Werte.
ydata = list(Mockdata) # meine Mockdata an die gefittet werden soll. 

# p0 sind die initial parameter values
optimal_parameters , _ = curve_fit(simplefun, xdata, ydata, p0 = [0.001,0.001])


k_opt = optimal_parameters[0] # abbaugeschwindigkeit (Fressgeschwindigkeit)
h_opt = optimal_parameters[1] # Microbenwachstum

# calculating the model output with optimal parameters:
CCH4CO2opt = simplefun(0,k_opt, h_opt)

plt.clf()
plt.plot(CCH4CO2opt[0:303], label = "CCH4CO2")
plt.plot(ydata[0:303])

# plotting each Model prediction and mockdata pair


for x, a in zip([0,101,202],["Cpool","CH4","CO2"]):  
    plt.figure()
    plt.plot([Mockdata[i] for i in range(x,101+x)], label = "Observed")
    plt.plot([CCH4CO2opt[i] for i in range(x,101+x)], label = "Predicted")
    plt.show()
    plt.ylabel(a)
    plt.legend()
            

#%%
# Goodness of fit with R^2
    
_ ,_ , R2, _ , _ = stats.linregress(CCH4CO2opt, y=ydata)
print("R2 is", R2)



#--calculating R^2 by hand-----------------------------------------------------
# Aufteilen von CCH4CO2opt in die jeweiligen pools und Berechnung des mean

Cobsmean   = np.mean(CCH4CO2opt[0:101])
CH4obsmean = np.mean(CCH4CO2opt[101:202])
CO2obsmean = np.mean(CCH4CO2opt[202:303])
allobsmean = np.mean(CCH4CO2opt)

#---- calculating total sum of squares for each pool and global----------------
SStotC   = []
SStotCH4 = []
SStotCO2 = []
SStotall = []

for i in range(101):
    SStotC.append((Mockdata[i]- Cobsmean)**2)
    SStotCH4.append((Mockdata[i+100]- CH4obsmean)**2)
    SStotCO2.append((Mockdata[i+200]- CO2obsmean)**2)

# list comprehension für den globalen datensatz    
SStotall.append([((Mockdata[i] - allobsmean)**2) for i in range(len(Mockdata))]) 

SStotC   = np.sum(SStotC)
SStotCH4 = np.sum(SStotCH4)
SStotCO2 = np.sum(SStotCO2)
SStotall = np.sum(SStotall)


#---- calculating total sum of residuals for each pool and global-------------- 
SSresC   = []
SSresCH4 = []
SSresCO2 = []
SSresall = []


for i in range(101):
    SSresC.append((Mockdata[i]- CCH4CO2opt[i])**2)
    SSresCH4.append((Mockdata[i+100]- CCH4CO2opt[i+100])**2)
    SSresCO2.append((Mockdata[i+200]- CCH4CO2opt[i+200])**2)
    
# list comprehension für den globalen datensatz        
SSresall.append([((Mockdata[i]- CCH4CO2opt[i])**2) for i in range(len(Mockdata))])    
   

SSresC = np.sum(SSresC)
SSresCH4 = np.sum(SSresCH4)
SSresCO2 = np.sum(SSresCO2)
SSresall = np.sum(SSresall)


R2C = 1 - (SStotC / SSresC)
print("R2C is", R2C)
R2CH4 = 1 - (SStotCH4 / SSresCH4)
print("R2CH4 is", R2CH4)
R2CO2 = 1 - (SStotCO2 / SSresCO2)
print("R2CO2 is", R2CO2)
R2all = 1 - (SStotall / SSresall)
print("R2all is", R2all)







