# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:45:44 2020

@author: Lara
"""
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import numpy as np
from SimpleIN import Cdec # importiert das eigentliche mathematische Model 
from scipy.integrate import odeint

#Importieren benutzter Funktionen und Daten
#from Realdataall import load_realdata
from Readalldata import load_realdata
from Babypascal import Mol_nach_Pa



def optifun(xdata, *Fitters):    
    xtage = xdata[0:int(len(xdata)/2)]
 
    #CH4, CO2,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_= simplefun(xtage,*Fitters)
    CH4, CO2,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_= simplesolve(xtage,*Fitters)
    return CH4 + CO2 # setzt die für uns interessanten Ausgaben in curvefit format zusammen
# CH4 und CO2 werden aneinandergehängt um einen abhängigen Least Sqaure zu berechnen und nicht jeweils einen unabhängigen pro Kurve  
    
# SIMPELFUN, der Euler Forward Mechanismus

def simplesolve(xtage,*Fitters ):
    
    # Festgelegte Initialwerte
    # Die Berechnung hier für Cpool ist noch für zwei Pools gedacht. Momentan nicht relevant
    m_gluc = 180   # molar mass of glucose, g dw pro mol
    m_cell = 162   # molar mass of cellulose, g dw pro mol
    TOC = 0.04     # Knoblauchs Daten , g dw
    labile = 0.005 # Knoblauchs Daten, anteil von TOC einheitslos
    Cpool_init = (TOC*labile/m_gluc + TOC*(1-labile)/m_cell) * (10**6) # 0.00024679012345679013 * (10**6) mikromol pro g
    
    # die Werte sind alle in mikroMol pro gram Trockengewicht
    AltE_init = Fitters[-1] # 7 #146       # 146 cf. Yao, Conrad 1999 (Hyun2017 sagt werte um 100 sind extrem hoch), Philben: 6 jeweils Mikromol pro g dw
    Fitters = Fitters[:-1]
    M_A_CH4_init = 0.001     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Ferm_init = 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    M_AltE_init= 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    deltaH2_Hydro = 0        # für plot in multifit
    M_H_CH4_init = 0.2     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Homo_init = 0.2        # Monteux 2020, mg Mikrobielles C pro g dw
    CH4_init = 0
    CO2_init = 0
    Acetate_init = 0        # Philben wert ca 3, in mikromol pro g 
    AceCO2_init = 0
    H2_init = 0
    deltaH2_Homo_init = 0
    CO2_Hydro_init= 0
    CH4_Hydro_init = 0
    H2_Ferm2_init = 0
    M_Ferm2_init = 0.2
    
    xs = np.linspace(0,int(max(xtage)), int(max(xtage)+1))
    y0 = Cpool_init,  AltE_init,  M_A_CH4_init,  M_Ferm_init,  M_AltE_init, M_H_CH4_init,  M_Homo_init,  CH4_init,  CO2_init,  AceCO2_init, Acetate_init,  H2_init,  CO2_Hydro_init,  CH4_Hydro_init, H2_Ferm2_init,  M_Ferm2_init
    
    ys = odeint(Cdec, y0, xs, args = (Fitters,))
    
    Cpool, AltEpool, M_A_CH4, M_Ferm, M_AltE, M_H_CH4, M_Homo, CH4, CO2, AceCO2, Acetate, H2,  CO2_Hydro, CH4_Hydro, H2_Ferm2, M_Ferm2 = zip(*ys)
    
   # AltEpool[0::int(1/step)]

   
    CH4 = [CH4[int(i)] for i in xtage]
    CO2 = [CO2[int(i)] for i in xtage]


    return CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_Ferm, M_AltE, H2, M_H_CH4, M_Homo, AltE_init, deltaH2_Hydro,[0 for _ in range(len(CO2))], CO2_Hydro, CH4_Hydro,H2_Ferm2
    

def simplefun(xtage, *Fitters):    
    
    step = 1/24.0
    
    n = int(max(xtage)/step)
    number_steps = range(1, n+1)
    
    # Festgelegte Initialwerte
    # Die Berechnung hier für Cpool ist noch für zwei Pools gedacht. Momentan nicht relevant
    m_gluc = 180   # molar mass of glucose, g dw pro mol
    m_cell = 162   # molar mass of cellulose, g dw pro mol
    TOC = 0.04     # Knoblauchs Daten , g dw
    labile = 0.005 # Knoblauchs Daten, anteil von TOC einheitslos
    Cpool_init = (TOC*labile/m_gluc + TOC*(1-labile)/m_cell) * (10**6) # 0.00024679012345679013 * (10**6) mikromol pro g
    
    # die Werte sind alle in mikroMol pro gram Trockengewicht
    AltE_init = Fitters[-1] # 7 #146       # 146 cf. Yao, Conrad 1999 (Hyun2017 sagt werte um 100 sind extrem hoch), Philben: 6 jeweils Mikromol pro g dw
    Fitters = Fitters[:-1]
    M_A_CH4_init = 0.001     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Ferm_init = 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    M_AltE_init= 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    deltaH2_Hydro = 0        # für plot in multifit
    M_H_CH4_init = 0.2     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Homo_init = 0.2        # Monteux 2020, mg Mikrobielles C pro g dw
    CH4_init = 0
    CO2_init = 0
    Acetate_init = 0        # Philben wert ca 3, in mikromol pro g 
    AceCO2_init = 0
    H2_init = 0
    deltaH2_Homo_init = 0
    CO2_Hydro_init= 0
    CH4_Hydro_init = 0
    H2_Ferm2_init = 0
    M_Ferm2_init = 0.2
    
    Cpool    = [Cpool_init]
    AltEpool = [AltE_init]
    M_A_CH4  = [M_A_CH4_init]
    M_Ferm    = [M_Ferm_init]
    M_AltE   = [M_AltE_init]
    deltaH2_Hydro = [deltaH2_Hydro]
    M_H_CH4  = [M_H_CH4_init]
    M_Homo   = [M_Homo_init]
    CH4      = [CH4_init]
    CO2      = [CO2_init]
    Acetate  = [Acetate_init]
    AceCO2   = [AceCO2_init]
    H2       = [H2_init]
    deltaH2_Homo = [deltaH2_Homo_init]
    CO2_Hydro = [CO2_Hydro_init]
    CH4_Hydro = [CH4_Hydro_init]
    H2_Ferm2 = [H2_Ferm2_init]
    M_Ferm2 = [M_Ferm2_init]

    
    for t in number_steps: # iteriert über n+1 Zeitschritte (n = max(xtage))

        delta = Cdec(Cpool[-1], AltEpool[-1], M_A_CH4[-1],M_Ferm[-1], M_AltE[-1],M_H_CH4[-1],M_Homo[-1], CH4[-1], CO2[-1], AceCO2[-1], Acetate[-1], H2[-1],CO2_Hydro[-1], CH4_Hydro[-1], H2_Ferm2[-1], M_Ferm2[-1], *Fitters)

        
        Cpool.append(Cpool[-1] +       step*delta[0])
        AltEpool.append(AltEpool[-1] + step*delta[1])
        M_A_CH4.append(M_A_CH4[-1] +   step*delta[2])
        M_Ferm.append(M_Ferm[-1] +     step*delta[3])
        M_AltE.append(M_AltE[-1] +     step*delta[4])
        M_H_CH4.append(M_H_CH4[-1] +   step*delta[5])
        M_Homo.append(M_Homo[-1] +     step*delta[6])
        CH4.append(CH4[-1] +           step*delta[7])
        CO2.append(CO2[-1] +           step*delta[8])
        AceCO2.append(AceCO2[-1] +     step*delta[9])
        Acetate.append(Acetate[-1] +   step*delta[10])
        H2.append(H2[-1] +             step*delta[11])
        deltaH2_Hydro.append(          step*delta[12])
        deltaH2_Homo.append(           step*delta[13])
        CO2_Hydro.append(              step*delta[14])
        CH4_Hydro.append(              step*delta[15])
        H2_Ferm2.append(               step*delta[16])
        M_Ferm2.append( M_Ferm2[-1] +  step*delta[17])
    
    
   # AltEpool[0::int(1/step)]
    
    CH4 = [CH4[int(i/step)] for i in xtage]
    CO2 = [CO2[int(i/step)] for i in xtage]
    
    AltEpool = AltEpool[0::int(1/step)] # von 0 bis Ende in Schritten von , damit der Plot nicht in Stunden sondern in Tagen ist
    AceCO2 = AceCO2[0::int(1/step)]
    Acetate= Acetate[0::int(1/step)]
    Cpool= Cpool[0::int(1/step)]
    M_A_CH4= M_A_CH4[0::int(1/step)]
    M_Ferm= M_Ferm[0::int(1/step)]
    M_AltE= M_AltE[0::int(1/step)]
    H2= H2[0::int(1/step)]
    M_H_CH4= M_H_CH4[0::int(1/step)]
    M_Homo= M_Homo[0::int(1/step)]
    AltE_init= AltE_init
    deltaH2_Hydro= deltaH2_Hydro[0::int(1/step)]
    deltaH2_Homo= deltaH2_Homo[0::int(1/step)]
    CO2_Hydro= CO2_Hydro[0::int(1/step)]
    CH4_Hydro= CH4_Hydro[0::int(1/step)]
    H2_Ferm2= H2_Ferm2[0::int(1/step)]

    return CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_Ferm, M_AltE, H2, M_H_CH4, M_Homo, AltE_init, deltaH2_Hydro,deltaH2_Homo, CO2_Hydro, CH4_Hydro,H2_Ferm2


if __name__ == '__main__':
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
               0.3,     # Vmax AltE 0.9
               0.133,    # Vmax Homo 0.05
               0.086,    # Vmax Hydro 0.05
               0.207,    # Vmax Ace 0.15   roden2003competition wert ist 15 !
               0.05,    # w Ferm 0.04
               0.013,    # w AltE 0.02
               0.024,    # w Hydro 0.02
               0.049,    # w Homo 0.05
               0.04,    # w Ace 0.03
               0.0000833,   # Sensenmann 0.0001   delattre2020thermodynamic nimmt 8.33*10^-4 h^-1für alle mikroben
               4,        # Stoch AltE 7  #  philben2020anaerobic nehmen einen Wert von 4 für Fe3 an
               5.75,      # AltE pool init 7
               10,      # Kmb ferm, for inverse M-M Biomass 10
               10,      # Kmh ferm, for hemmung of fermenters by acetate 10
               10,      # Kmb AltE 10
               10,      # Kmb auto 10
               10]      # Kmb hydro 10
        
        bounds = [[0.01, 0.11],     # Vmax Ferm
                  [0.029, 1.9],     # Vmax AltE
                  [0.005, 1],       # Vmax Homo
                  [0.03,  0.2],       # Vmax Hydro
                  [0.05, 0.39],     # Vmax Ace
                  [0.03, 0.05],     # w Ferm
                  [0.01, 0.05],     # w AltE
                  [0.01, 0.05],     # w Hydro
                  [0.01, 0.05],     # w Homo
                  [0.01, 0.05],     # w Ace
                  [-0.000000001, 0.0000844],    # Sensenmann
                  [1,    8],        # Stoch AltE
                  [2,   10],        # AltE pool init
                  [ 1, 10 ],        # Kmb ferm
                  [ 1, 10 ],        # Kmh ferm
                  [ 1, 10 ],        # Kmb AltE
                  [ 1, 10 ],        # Kmb Auto
                  [ 1, 10 ]]        # Kmh Hydro
               
        optimal_parameters , _ = curve_fit(optifun, xdata, ydata, #method="dogbox",
                                           p0 = p0, 
                                           bounds=tuple(zip(*bounds)))
        
        optimal_parameters , _ = curve_fit(optifun, xdata, ydata, #method="dogbox",
                                           p0 = optimal_parameters, 
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
        #CCH4CO2opt = simplefun(xlist,  *optimal_parameters)
        CCH4CO2opt = simplesolve(xlist,  *optimal_parameters)
        CCH4CO2optList = list(CCH4CO2opt[0]) + list(CCH4CO2opt[1])
          
        Mol_nach_Pa(n= max(CCH4CO2opt[9]))
        
    #### Observed and Predicted für CH4 und CO2 geplotted
    #%%
        data_len = len(Realdata)
        for x, a, col in zip([0,data_len],["CH4","CO2"], ["r", "b"]):
            plt.figure()
            plt.plot(xlist,[ydata[i] for i in range(x,data_len+x)],col+"o", label = "Observed")
            plt.plot(xlist,[CCH4CO2optList[i] for i in range(x,data_len+x)],"k-",label = "Predicted")
            plt.ylabel(a)
            plt.legend()
            save_path = os.path.join('C:/Users/Lara/Desktop/simple model/Figs', a +'_'+str(m)+'.png')
            plt.savefig(save_path)
        _ ,_ , R2, _ , _ = stats.linregress(CCH4CO2optList, y=ydata)
        print("R2 is", R2)
        
         
    #%%      
        # return order CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_CH4,    M_CO2, M_AltE, M_Hyd
        #              CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_Ferm, M_AltE, H2, M_H_CH4, M_Homo, AltE_init, deltaCO2_Hydro
        
        plt.figure()
        plt.plot( CCH4CO2opt[2], label='AltEpool')
        plt.title('AltEpool')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/AltEpool.png')    
        plt.ylabel('mg Mikrobielles C pro g dw')
        
        plt.figure()
        plt.plot( CCH4CO2opt[3], label='CO2 aus Acetate')
        plt.title('CO2 aus Acetate')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/CO2ausAcetate.png')
        plt.ylabel('mg Mikrobielles C pro g dw')
        
        plt.figure()
        plt.plot( CCH4CO2opt[4], label='Acetate')
        plt.title('Acetate')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Acetate.png')
        plt.ylabel('μmol')
        
        plt.figure()
        plt.plot( CCH4CO2opt[5], label='Cpool')
        plt.title('Cpool')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Cpool.png')
        plt.ylabel('mikromol pro g dw')
        
        plt.figure()
        plt.plot( CCH4CO2opt[6], label='Acetoclastische_Mikroben')
        plt.title('Acetoclastische_Mikroben')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Acetoclastic_Mikroben.png')  
        plt.ylabel('mg Mikrobielles C pro g dw')
        
        plt.figure()
        plt.plot( CCH4CO2opt[7], label='Fermentierer')
        plt.title('Fermentierer')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Fermentation_Mikroben.png')  
        plt.ylabel('mg Mikrobielles C pro g dw')
        
        plt.figure()
        plt.plot( CCH4CO2opt[8], label=' AltE_Microben')
        plt.title(' AltE_Microben')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/AltE_Mikroben.png') 
        plt.ylabel('mg Mikrobielles C pro g dw')
        
        plt.figure()
        plt.plot( CCH4CO2opt[9], label=' H2') #typically in the nanomolar range (10 to 180 nM, correponding to about 8 to 140 ppmv in the gas phase [Table conrad
        plt.title('H2')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/H2.png') 
        plt.ylabel('μmol')
        
        plt.figure()
        plt.plot( CCH4CO2opt[10], label=' Hydro_Microben')
        plt.title(' Hydro_Microben')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Hydro_Mikroben.png')    
        plt.ylabel('mg Mikrobielles C pro g dw')
        
        plt.figure()
        plt.plot( CCH4CO2opt[11], label=' Homo_Microben')
        plt.title(' Homo_Microben')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Homo_Mikroben.png') 
        plt.ylabel('mg Mikrobielles C pro g dw')
        
        #plt.figure()
        #plt.plot( CCH4CO2opt[12], label=' AltE_init')
        #plt.title(' AltE_init')
        
        plt.figure()
        plt.plot( CCH4CO2opt[13], label='deltaH2_Hydro')
        plt.title('deltaH2_Hydro')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/deltaH2_Hydro.png')
        plt.ylabel('μmol')        
        
        plt.figure()
        plt.plot( CCH4CO2opt[14], label='deltaH2_Homo')
        plt.title('deltaH2_Homo')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/deltaH2_Homo.png')
        plt.ylabel('μmol')
        
        plt.figure()
        plt.plot( CCH4CO2opt[15], label='CO2 aus Hydro')
        plt.title('CO2 aus Hydro')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/CO2_aus_Hydro.png')
        plt.ylabel('μmol')    
    
        plt.figure()
        plt.plot( CCH4CO2opt[16], label='CH4 aus Hydro')
        plt.title('CH4 aus Hydro')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/CH4_aus_Hydro.png') 
        plt.ylabel('μmol')
        
        plt.figure()
        plt.plot( CCH4CO2opt[17], label='H2 aus Ferm2')
        plt.title('H2 aus Ferm2')
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/H2_aus_Ferm.png')  
        plt.ylabel('μmol')
    
    
        
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
        
        
        
    
    
