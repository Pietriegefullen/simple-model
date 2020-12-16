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

from Realdataall import load_realdata
from SimpleOUT import simplefun
from SimpleOUT import optifun

plt.close('all')
#%% die drei verschiedenen Datensätze 
Data1 = [0]
Data1and2 = [0,3]
Data1and3and4 =[0,3,6]

for m in Data1:
    Realdata = load_realdata(m)
    # Fitting the parameters: 
    #%%    
    
    xlist = [int(Realdata[i,0]) for i in range(len(Realdata[:,0]))] # int der Tage an denen wir Messwerte haben, Länge 44
    xdata = xlist + xlist  # aneinandergehängt, weil wir die werte sowohl für CH4 als auch CO2 brauchen
    ydata = list(Realdata[:,1]) + list(Realdata[:,2]) # meine Realdata an die gefittet werden soll. Länge 88
    
    
    #  Vmax
    optimal_parameters , _ = curve_fit(optifun, xdata, ydata, p0 = [ 50.0], bounds=((0), (100)))
    
    
#    f_A_CH4_opt = optimal_parameters[0] # Abbaugeschwindigkeit (Fressgeschwindigkeit)
#    f_CO2_opt = optimal_parameters[1] # Abbaugeschwindigkeit (Fressgeschwindigkeit)
#    f_H_CH4_opt = optimal_parameters[2]
#    w_A_CH4_opt = optimal_parameters[3] # Mikrobenwachstum
#    w_CO2_opt = optimal_parameters[4] # Mikrobenwachstum
#    w_alte_opt = optimal_parameters[5]# Mikrobenwachstum
#    w_A_CH4_heil_opt = optimal_parameters[6] # Mikrobenheilrate
#    w_H_CH4_opt = optimal_parameters[7]
#    w_Homo_opt = optimal_parameters[8]
#    print("f_A_CH4_opt is",f_A_CH4_opt)
#    print("f_CO2_opt is",f_CO2_opt)
#    print("f_H_CH4_opt is",f_H_CH4_opt)
#    print("w_A_CH4_opt is",w_A_CH4_opt)
#    print("w_CO2_opt is",w_CO2_opt)
#    print("w_alte_opt is",w_alte_opt)
#    print("w_A_CH4_heil is", w_A_CH4_heil_opt)
#    print("w_H_CH4_opt is",w_H_CH4_opt)
#    print("w_Homo_opt is",w_Homo_opt)
  
#    w_A_CH4_opt = optimal_parameters[0] # Mikrobenwachstum
#    #w_CO2_opt = optimal_parameters[3] # Mikrobenwachstum
#    w_alte_opt = optimal_parameters[1]# Mikrobenwachstum
#    w_A_CH4_heil_opt = optimal_parameters[2] # Mikrobenheilrate
#    #w_H_CH4_opt = optimal_parameters[5]
#    w_Homo_opt = optimal_parameters[3]
#   
#    print("w_A_CH4_opt is",w_A_CH4_opt)
#   # print("w_CO2_opt is",w_CO2_opt)
#    print("w_alte_opt is",w_alte_opt)
#    print("w_A_CH4_heil is", w_A_CH4_heil_opt)
#   # print("w_H_CH4_opt is",w_H_CH4_opt)
#    print("w_Homo_opt is",w_Homo_opt)
    
    Vmax_opt = optimal_parameters[0] # Mikrobenwachstum
    print("Vmax_opt is", Vmax_opt)
    
    #%%
    # calculating the model output with optimal parameters:

    #CCH4CO2opt = simplefun(xlist, f_A_CH4_opt, f_CO2_opt,f_H_CH4_opt,  w_A_CH4_opt,w_CO2_opt, w_alte_opt, w_A_CH4_heil_opt, w_Homo_opt)
    CCH4CO2opt = simplefun(xlist,  Vmax_opt)
    CCH4CO2optList = list(CCH4CO2opt[0]) + list(CCH4CO2opt[1])
       
    
    #%% Obs and Pred für CH4 und CO2 geplotted
    #plt.close('all')
    for x, a in zip([0,44],["CH4","CO2"]):  
        plt.figure()
        plt.plot(xlist,[ydata[i] for i in range(x,44+x)],"ro", label = "Observed")
        plt.plot(xlist,[CCH4CO2optList[i] for i in range(x,44+x)],label = "Predicted")
        plt.ylabel(a)
        plt.legend()
                
        #%%
        
    # return order CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_CH4, M_CO2, M_AltE, M_Hyd
    
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
    plt.plot( CCH4CO2opt[6], label='Aceto_Microben_CH4')
    plt.title('Aceto_Microben_CH4')
    #axes = plt.gca()
    #axes.set_ylim([0,200000])
    #np.min(np.where(np.array(CCH4CO2opt[6][0:100]) >0))
    
    plt.figure()
    plt.plot( CCH4CO2opt[7], label='Fermentierer_CO2')
    plt.title('Fermentierer_CO2')
    
    plt.figure()
    plt.plot( CCH4CO2opt[8], label=' Microben_AltE')
    plt.title(' Microben_AltE')
    
    plt.figure()
    plt.plot( CCH4CO2opt[9], label=' Hydro_Microben_CH4')
    plt.title(' Hydro_Microben_CH4')
    
    plt.figure()
    plt.plot( CCH4CO2opt[10], label=' H2')
    plt.title('H2')
    
    #%%
    # Goodness of fit with R^2 automatic
        
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
    
plt.show()   
    
    
    
    
    
