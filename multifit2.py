# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:45:44 2020

@author: Lara
"""

#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
#import numpy as np

#Importieren benutzter Funktionen und Daten
#from Realdataall import load_realdata
from Readalldata import load_realdata
from SimpleOUT import simplefun
from SimpleOUT import optifun
from Babypascal import Pascal_Pa

plt.close('all')
#die drei verschiedenen Datensätze 
Data1 = [0]
Data1and2 = [0,3]
Data1and2and3 =[0,3,6]
Data1and2and3and4and5and6 =[0,3,6,9,12,15]
Data1and2and3and4and5and6and7and8and9 =[0,3,6,9,12,15,18,21,24]


for m in Data1:#and2and3and4and5and6and7and8and9:
    Realdata = load_realdata(m)
    # Fitting the parameters:    
    
    xlist = [int(Realdata[i,0]) for i in range(len(Realdata[:,0]))] # int der Tage an denen wir Messwerte haben
    xdata = xlist + xlist  # aneinandergehängt, weil wir die werte sowohl für CH4 als auch CO2 brauchen
    ydata = list(Realdata[:,1]) + list(Realdata[:,2]) # meine Realdata an die gefittet werden soll.
    
    
    # Fitted Parameters are: 
    # Vmax_Ferm, Stoch_ALtE,Vprod_max_AltE, Vprod_max_Homo, Vprod_max_Hydro,
    # Vprod_max_Ace, w_Ferm, w_AltE, w_Hydro, w_Homo, w_Ace, Sensenmann = Fitters
    # Boundaries für w aus bekannten YATP und ATP prod. Wachstumsfaktoren.
    
    
    p0 = [ 0.1,    # Vmax Ferm 0.01  roden2003competition wert ist 17 !
           0.8,     # Vmax AltE 0.9
           0.133,    # Vmax Homo 0.05
           0.086,    # Vmax Hydro 0.05
           0.207,    # Vmax Ace 0.15   roden2003competition wert ist 15 !
           0.05,    # w Ferm 0.04
           0.013,    # w AltE 0.02
           0.024,    # w Hydro 0.02
           0.049,    # w Homo 0.05
           0.04,    # w Ace 0.03
           0.000,   # Sensenmann 0.0001
           4,        # Stoch AltE 7  #  philben2020anaerobic nehmen einen Wert von 4 für Fe3 an
           5.75,      # AltE pool init 7
           10,      # Kmb ferm, for inverse M-M Biomass 10
           10,      # Kmh ferm, for hemmung of fermenters by acetate 10
           10,      # Kmb AltE 10
           10,      # Kmb auto 10
           10]      # Kmb hydro 10
    
    bounds = [[0.01, 0.10],     # Vmax Ferm
              [0.029, 1.9],     # Vmax AltE
              [0.005, 1],       # Vmax Homo
              [0.03,  0.2],       # Vmax Hydro
              [0.05, 0.39],     # Vmax Ace
              [0.03, 0.05],     # w Ferm
              [0.01, 0.05],     # w AltE
              [0.01, 0.05],     # w Hydro
              [0.01, 0.05],     # w Homo
              [0.01, 0.05],     # w Ace
              [-0.000000001, 0.0000000001],    # Sensenmann
              [1,    8],        # Stoch AltE
              [5,   10],        # AltE pool init
              [ 1, 10 ],        # Kmb ferm
              [ 1, 10 ],        # Kmh ferm
              [ 1, 10 ],        # Kmb AltE
              [ 1, 10 ],        # Kmb Auto
              [ 1, 10 ]]        # Kmh Hydro
           
    optimal_parameters , _ = curve_fit(optifun, xdata, ydata, #method="dogbox",
                                       p0 = p0, 
                                       bounds=tuple(zip(*bounds)))
    
   
    names = ["Vmax_Ferm","Vprod_max_AltE","Vprod_max_Homo", "Vprod_max_Hydro", "Vprod_max_Ace", "w_Ferm","w_AltE","w_Hydro","w_Homo","w_Ace", "Sensenmann", "Stoch_ALtE", "AltEpool", "KmB ferm", "Kmh ferm", "Kmb alte", "Kmb Auto", "Kmb Hydro"]
    units = ["μmol/mg",           "μmol/mg",       "μmol/mg",       "μmol/mg",          "μmol/mg",      "mg/μmol","mg/μmol","mg/μmol","mg/μmol","mg/μmol", "-",  "-","μmol", "mg" , "μmol", "mg", "mg", "mg"]
    song =  ["0.5",                "",            "0.15",           "0.15",            "0.5",         "",         "",     "",     "",         "",    "",  "",     "" ,"", "", "" , "", "", ""]
   
    #Printing the Parameter and its value
    for n, p, u in zip(names, optimal_parameters, units):
        print("{:<18} {:6.3f} {:<10}".format(n,p,u))
    
    # for comparison with Song
    SOIL_DENSITY = 1.3
    print('')
    for n, p, u, s in zip(names, optimal_parameters, units, song):
        if n.startswith("V"):
            print("{:<18} {:6.3f} {:<10} ({:})".format(n,p*SOIL_DENSITY,"mol/m^3", s))
            
            
    #Calculating the model output with optimal parameters:
    
  #  Fitters_opt = optimal_parameters
    CCH4CO2opt = simplefun(xlist,  *optimal_parameters)
    CCH4CO2optList = list(CCH4CO2opt[0]) + list(CCH4CO2opt[1])
      
    Pascal_Pa (n= max(CCH4CO2opt[9]))
    
#### Observed and Predicted für CH4 und CO2 geplotted

    data_len = len(Realdata)
    for x, a in zip([0,data_len],["CH4","CO2"]):  
        plt.figure()
        plt.plot(xlist,[ydata[i] for i in range(x,data_len+x)],"ro", label = "Observed")
        plt.plot(xlist,[CCH4CO2optList[i] for i in range(x,data_len+x)],label = "Predicted")
        plt.ylabel(a)
        plt.legend()
   
    _ ,_ , R2, _ , _ = stats.linregress(CCH4CO2optList, y=ydata)
    print("R2 is", R2)
    
             
#%%      
    # return order CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_CH4,    M_CO2, M_AltE, M_Hyd
    #              CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_Ferm, M_AltE, H2, M_H_CH4, M_Homo, AltE_init, deltaCO2_Hydro
    
    plt.figure()
    plt.plot( CCH4CO2opt[2], label='AltEpool')
    plt.title('AltEpool')
    
    plt.figure()
    plt.plot( CCH4CO2opt[3], label='CO2 aus Acetate')
    plt.title('CO2 aus Acetate')
    
    plt.figure()
    plt.plot( CCH4CO2opt[4], label='Acetate')
    plt.title('Acetate')
    
    plt.figure()
    plt.plot( CCH4CO2opt[5], label='Cpool')
    plt.title('Cpool')
    
    plt.figure()
    plt.plot( CCH4CO2opt[6], label='Acetoclastische_Mikroben')
    plt.title('Acetoclastische_Mikroben')
    
    plt.figure()
    plt.plot( CCH4CO2opt[7], label='Fermentierer')
    plt.title('Fermentierer')
    
    plt.figure()
    plt.plot( CCH4CO2opt[8], label=' AltE_Microben')
    plt.title(' AltE_Microben')
    
    plt.figure()
    plt.plot( CCH4CO2opt[9], label=' H2')
    plt.title('H2')
    
    plt.figure()
    plt.plot( CCH4CO2opt[10], label=' Hydro_Microben')
    plt.title(' Hydro_Microben')
    
    plt.figure()
    plt.plot( CCH4CO2opt[11], label=' Homo_Microben')
    plt.title(' Homo_Microben')
    
    #plt.figure()
    #plt.plot( CCH4CO2opt[12], label=' AltE_init')
    #plt.title(' AltE_init')
    
    plt.figure()
    plt.plot( CCH4CO2opt[13], label='deltaH2_Hydro')
    plt.title('deltaH2_Hydro')
    
    plt.figure()
    plt.plot( CCH4CO2opt[14], label='deltaH2_Homo')
    plt.title('deltaH2_Homo')
    
    plt.figure()
    plt.plot( CCH4CO2opt[15], label='CO2 aus Hydro')
    plt.title('CO2 aus Hydro')

    plt.figure()
    plt.plot( CCH4CO2opt[16], label='CH4 aus Hydro')
    plt.title('CH4 aus Hydro')
    
    plt.figure()
    plt.plot( CCH4CO2opt[17], label='H2 aus Ferm2')
    plt.title('H2 aus Ferm2')


    
    # Goodness of fit with R^2 automatic
        

################ ab hier ist es uninteressant #################################   
#    #--calculating R^2 by hand-----------------------------------------------------
#    # Aufteilen von CCH4CO2opt in die jeweiligen pools und Berechnung des mean
#    
#    
#    CH4obsmean = np.mean(CCH4CO2optList[0:int(len(ydata)/2)])
#    CO2obsmean = np.mean(CCH4CO2optList[int(len(ydata)/2):len(ydata)])
#    allobsmean = np.mean(CCH4CO2optList)
#    
#    #---- calculating total sum of squares for each pool and global----------------
#    
#    SStotCH4 = []
#    SStotCO2 = []
#    SStotall = []
#    
#    for i in range(len(xlist)):
#        SStotCH4.append((ydata[i]- CH4obsmean)**2)
#        SStotCO2.append((ydata[i+len(xlist)]- CO2obsmean)**2)
#        
#    # list comprehension für den globalen datensatz    
#    SStotall.append([((ydata[i] - allobsmean)**2) for i in range(len(ydata))]) 
#    
#    
#    SStotCH4 = np.sum(SStotCH4)
#    SStotCO2 = np.sum(SStotCO2)
#    SStotall = np.sum(SStotall)
#    
#    #%%
#    #---- calculating total sum of residuals for each pool and global-------------- 
#    
#    SSresCH4 = []
#    SSresCO2 = []
#    SSresall = []
#    
#    
#    for i in range(len(xlist)):
#        SSresCH4.append((ydata[i]- CCH4CO2optList[i])**2)
#        SSresCO2.append((ydata[i+len(xlist)]- CCH4CO2optList[i+len(xlist)])**2)
#        
#    # list comprehension für den globalen datensatz        
#    SSresall.append([((ydata[i]- CCH4CO2optList[i])**2) for i in range(len(ydata))])    
#       
#    
#    SSresCH4 = np.sum(SSresCH4)
#    SSresCO2 = np.sum(SSresCO2)
#    SSresall = np.sum(SSresall)
#    
#    
#    R2CH4 = 1 - (SStotCH4 / SSresCH4)
#    print("R2CH4 is", R2CH4)
#    R2CO2 = 1 - (SStotCO2 / SSresCO2)
#    print("R2CO2 is", R2CO2)
#    R2all = 1 - (SStotall / SSresall)
#    print("R2all is", R2all)
#    
plt.show()   
    
    
    
    
    
