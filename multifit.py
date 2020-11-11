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

from Realdata import load_realdata
from SimpleOUT import simplefun
from SimpleOUT import optifun

#%%
Realdata = load_realdata()
# Fitting the parameters: 

xlist = [int(Realdata[i,0]) for i in range(len(Realdata[:,0]))] # int der Tage an denen wir Messwerte haben, Länge 44
xdata = xlist + xlist  # aneinandergehängt, weil wir die werte sowohl für CH4 als auch CO2 brauchen
ydata = list(Realdata[:,1]) + list(Realdata[:,2]) # meine Realdata an die gefittet werden soll. Länge 88

# p0 sind die initial parameter values, in der Reihenfolge f_CH4,f_CO2,f_alte, w_CH4,w_CO2, w_alte, w_CH4_heil
optimal_parameters , _ = curve_fit(optifun, xdata, ydata, p0 = [0.001,0.001,0.001,   0.01,0.01,-0.5, 0.01], bounds=((0,0,0,0,0,-1, 0), (10,10,1,10, 10,0, 1)))



f_CH4_opt = optimal_parameters[0] # abbaugeschwindigkeit (Fressgeschwindigkeit)
f_CO2_opt = optimal_parameters[1] # abbaugeschwindigkeit (Fressgeschwindigkeit)
f_alte_opt = optimal_parameters[2]
w_CH4_opt = optimal_parameters[3] # Microbenwachstum
w_CO2_opt = optimal_parameters[4] # Microbenwachstum
w_alte_opt = optimal_parameters[5]
#ac_opt = optimal_parameters[6]
w_CH4_heil_opt = optimal_parameters[6]
#c2ant_opt = 1-optimal_parameters[4]
print("f_CH4_opt is",f_CH4_opt)
print("f_CO2_opt is",f_CO2_opt)
print("f_alte_opt is",f_alte_opt)
print("w_CH4_opt is",w_CH4_opt)
print("w_CO2_opt is",w_CO2_opt)
print("w_alte_opt is",w_alte_opt)
#print("ac is", ac_opt)
print("w_CH4_heil is", w_CH4_heil_opt)
#print("c2ant is", c2ant_opt)


#%%
# calculating the model output with optimal parameters:
#k_opt = .001
#h_opt = .01
#c4ant_opt = .5
CCH4CO2opt = simplefun(xlist, f_CH4_opt, f_CO2_opt, w_alte_opt, w_CH4_opt,w_CO2_opt, w_alte_opt, w_CH4_heil_opt)
CCH4CO2optList = list(CCH4CO2opt[0]) + list(CCH4CO2opt[1])


#%% plots manuel erstellen
CH4opt = CCH4CO2opt[0]
CO2opt = CCH4CO2opt[1]
plt.figure()
plt.plot(xlist,CH4opt, label='predicted')
plt.plot(xlist,Realdata[:,1], label='observed')
plt.title('CH4')

plt.figure()
plt.plot(xlist,CO2opt, label='predicted')
plt.plot(xlist,Realdata[:,2], label='observed')
plt.title('CO2')



#%% Obs and Pred für CH4 und CO2 geplotted
plt.close('all')
for x, a in zip([0,44],["CH4","CO2"]):  
    plt.figure()
    plt.plot(xlist,[ydata[i] for i in range(x,44+x)],"ro", label = "Observed")
    plt.plot(xlist,[CCH4CO2optList[i] for i in range(x,44+x)],label = "Predicted")
    plt.show()
    plt.ylabel(a)
    plt.legend()
            
    #%%
    
    
# return order CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, Microben_CH4, Microben_CO2, Microben_AltE

plt.figure()
plt.plot( CCH4CO2opt[2], label='AltEpool')
plt.title('AltEpool')

plt.figure()
plt.plot( CCH4CO2opt[3], label='AceCO2')
plt.title('AceCO2')

plt.figure()
plt.plot( CCH4CO2opt[4], label='Acetate')
plt.title('Acetate')

plt.figure()
plt.plot( CCH4CO2opt[5], label='Cpool')
plt.title('Cpool')

plt.figure()
plt.plot( CCH4CO2opt[6], label='Microben_CH4')
plt.title('Microben_CH4')


plt.figure()
plt.plot( CCH4CO2opt[7], label='Microben_CO2')
plt.title('Microben_CO2')


plt.figure()
plt.plot( CCH4CO2opt[8], label=' Microben_AltE')
plt.title(' Microben_AltE')

#%%
# Goodness of fit with R^2
    
_ ,_ , R2, _ , _ = stats.linregress(CCH4CO2optList, y=ydata)
print("R2 is", R2)

#%%

#--calculating R^2 by hand-----------------------------------------------------
# Aufteilen von CCH4CO2opt in die jeweiligen pools und Berechnung des mean


CH4obsmean = np.mean(CCH4CO2optList[0:int(len(ydata)/2)])
CO2obsmean = np.mean(CCH4CO2optList[int(len(ydata)/2):len(ydata)])
allobsmean = np.mean(CCH4CO2optList)

#---- calculating total sum of squares for each pool and global----------------

SStotCH4 = []
SStotCO2 = []
SStotall = []

for i in range(len(xlist)):
    SStotCH4.append((ydata[i]- CH4obsmean)**2)
    SStotCO2.append((ydata[i+len(xlist)]- CO2obsmean)**2)
    
# list comprehension für den globalen datensatz    
SStotall.append([((ydata[i] - allobsmean)**2) for i in range(len(ydata))]) 


SStotCH4 = np.sum(SStotCH4)
SStotCO2 = np.sum(SStotCO2)
SStotall = np.sum(SStotall)

#%%
#---- calculating total sum of residuals for each pool and global-------------- 

SSresCH4 = []
SSresCO2 = []
SSresall = []


for i in range(len(xlist)):
    SSresCH4.append((ydata[i]- CCH4CO2optList[i])**2)
    SSresCO2.append((ydata[i+len(xlist)]- CCH4CO2optList[i+len(xlist)])**2)
    
# list comprehension für den globalen datensatz        
SSresall.append([((ydata[i]- CCH4CO2optList[i])**2) for i in range(len(ydata))])    
   

SSresCH4 = np.sum(SSresCH4)
SSresCO2 = np.sum(SSresCO2)
SSresall = np.sum(SSresall)



R2CH4 = 1 - (SStotCH4 / SSresCH4)
print("R2CH4 is", R2CH4)
R2CO2 = 1 - (SStotCO2 / SSresCO2)
print("R2CO2 is", R2CO2)
R2all = 1 - (SStotall / SSresall)
print("R2all is", R2all)







