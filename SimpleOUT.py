import matplotlib.pyplot as plt
from SimpleIN import Cdec # importiert das eigentliche mathematische Model 
import numpy as np
#


# ruft Cdec für 100 Zeitschritte auf 

# Parameter:
# xdata: unrelevnt für Cdec, wird aber von Curvefit als Vorraussetzung gesehen. 

# Variablen:
# k und h für die optimierung

# Returns:
# Eine Liste mit den Werten von Cpool, CH4 und CO2. Brauchen wir um den least square
# nicht nur für jede Kurve einzeln sondern global zu berechnen. Ein k und ein h, statt 3.

def optifun(xdata,k_CH4, k_CO2,h_CH4, h_CO2, c4ant):
    #print(xdata)
    xtage = xdata[0:int(len(xdata)/2)]
    #print(type(xtage))
    CH4, CO2 = simplefun(xtage,k_CH4, k_CO2, h_CH4 ,h_CO2, c4ant)
    return CH4 + CO2 # setzt die für uns interessanten Ausgaben in curvefit format zusammen
    
    


def simplefun(xtage, k_CH4, k_CO2, h_CH4, h_CO2, c4ant): # xdata muss eingegeben werden für curve_fit in multifit
    n = max(xtage)
    
    # Festgelegte Initialwerte
    C_init = 100
    Microben_CH4_init = 0.001
    Microben_CO2_init = 0.001
    CH4_init = 0
    CO2_init = 0
    Cused_init = 0


    
    Cpool = [C_init]
    Microben_CH4 = [Microben_CH4_init]
    Microben_CO2= [Microben_CO2_init]
    CH4 = [CH4_init]
    CO2 = [CO2_init]
    Cused = [Cused_init]
    

    
    for t in range(1,n+1): # iteriert über 100 Zeitschritte
        
    
        delta = Cdec(Cpool[-1],Microben_CH4[-1],Microben_CO2[-1], CH4[-1], CO2[-1],k_CH4, k_CO2, h_CH4, h_CO2, c4ant)
        
        Cpool.append(Cpool[-1] + delta[0])# hängt den Wert aus jedem Zeitschrit aus Cdec return [0] an.
        Microben_CH4.append(Microben_CH4[-1] + delta[1])
        Microben_CO2.append(Microben_CO2[-1] + delta[2])
        CH4.append(CH4[-1] + delta[3])
        CO2.append(CO2[-1] + delta[4])
        Cused.append(Cused[-1] + delta[5])

    if False:
        plt.close('all') #closes all previous plots
        #plt.plot(Cmet)
        plt.plot(Cpool, label = "Cpool")
        #plt.ylabel('C pool')
        #plt.show()
        
        #plt.figure()
        
        plt.plot(CH4, label = "CH4")
        #plt.ylabel('Methane')
        #plt.show()
        
        #plt.figure()
        
        plt.plot(CO2, label = "CO2",linestyle='--')
        #plt.ylabel('CO2')
        #plt.show()
        
        #plt.figure()
        
        plt.plot(Microben_CH4, label = "Microben CH4")
        plt.ylabel('Menge')
        plt.show
        plt.legend()
        
        plt.plot(Microben_CO2, label = "Microben CO2")
        plt.ylabel('Menge')
        plt.show
        plt.legend()
        
        
        plt.figure()
        plt.plot(Cused, label = "Cused")
    
    CH4 = [CH4[i] for i in xtage]
    CO2 = [CO2[i] for i in xtage]
    return CH4, CO2
