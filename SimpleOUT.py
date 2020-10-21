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

def simplefun(xdata, k, h, n = 101): # xdata muss eingegeben werden für curve_fit in multifit
    
    # Festgelegte Initialwerte
    C_init = 140
    Microben_init = 0.01
    CH4_init = 0
    CO2_init = 0
    Cv_init = 0
    Temperature = np.linspace(0,20, n)

    
    Cpool = [C_init]
    Microben = [Microben_init]
    CH4 = [CH4_init]
    CO2 = [CO2_init]
    Cv = [Cv_init]
    CCH4CO2 =[]

    
    for t in range(1,n): # iteriert über 100 Zeitschritte
        
    
        delta = Cdec(Cpool[-1],Microben[-1], CH4[-1], CO2[-1],k, h)
        
        Cpool.append(Cpool[-1] + delta[0])# hängt den Wert aus jedem Zeitschrit aus Cdec return [0] an.
        Microben.append(Microben[-1] + delta[1])
        CH4.append(CH4[-1] + delta[2])
        CO2.append(CO2[-1] + delta[3])
        Cv.append(Cv[-1] + delta[4])
       
    CCH4CO2 = Cpool + CH4 + CO2 # setzt die drei für uns interessanten Ausgaben zu einer 
    #Ausgabe zusammen um in Curvefit ein least square nutzen zu können. 

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
        
        plt.plot(Microben, label = "Microben")
        plt.ylabel('Menge')
        plt.show
        plt.legend()
        
        plt.figure()
        plt.plot(Cv, label = "Cv")
    
    
    return CCH4CO2
