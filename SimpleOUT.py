import matplotlib.pyplot as plt
from SimpleIN import Cdec # importiert das eigentliche mathematische Model 

#


# ruft Cdec für 100 Zeitschritte auf 

# Parameter:
# xdata: unrelevnt für Cdec, wird aber von Curvefit als Vorraussetzung gesehen. 

# Variablen:
# k und h für die optimierung

# Returns:
# Eine Liste mit den Werten von Cpool, CH4 und CO2. Brauchen wir um den least square
# nicht nur für jede Kurve einzeln sondern global zu berechnen. Ein k und ein h, statt 3.

def optifun(xdata,f_CH4, f_CO2,  w_CH4, w_CO2, w_alte, w_CH4_heil):
    #print(xdata)
    xtage = xdata[0:int(len(xdata)/2)]
    #print(type(xtage))
    CH4, CO2, _,_,_,_,_,_,_ = simplefun(xtage,f_CH4, f_CO2,  w_CH4 , w_CO2, w_alte, w_CH4_heil)
    return CH4 + CO2 # setzt die für uns interessanten Ausgaben in curvefit format zusammen
    
    


def simplefun(xtage, f_CH4, f_CO2,  w_CH4, w_CO2, w_alte, w_CH4_heil): # xdata muss eingegeben werden für curve_fit in multifit
    n = max(xtage)
    
    # Festgelegte Initialwerte
    C_init = 1000
    AltE_init = 4
    Microben_CH4_init = 0
    Microben_CO2_init = 0.001
    Microben_AltE_init= 0.001
    Microben_CH4_krank_init = 0.001
    CH4_init = 0
    CO2_init = 0
    Cused_init = 0
    Acetate_init = 0.01
    AceCO2_init = 0
    
    Cpool = [C_init]
    AltEpool = [AltE_init]
    Microben_CH4 = [Microben_CH4_init]
    Microben_CO2= [Microben_CO2_init]
    Microben_AltE = [Microben_AltE_init]
    Microben_CH4_krank = [Microben_CH4_krank_init]
    CH4 = [CH4_init]
    CO2 = [CO2_init]
    Cused = [Cused_init]
    Acetate = [Acetate_init]
    AceCO2 = [AceCO2_init]

    
    for t in range(1,n+1): # iteriert über 100 Zeitschritte

                #Cpool,          AltEpool,  Microben_CH4,     Microben_CO2,    Microben_AltE,     FermCO2,  Acetate,  Microben_CH4_krank, CH4, CO2,
        delta = Cdec(Cpool[-1], AltEpool[-1],Microben_CH4[-1],Microben_CO2[-1], Microben_AltE[-1], CH4[-1], CO2[-1], AceCO2[-1], Acetate[-1], Microben_CH4_krank[-1], f_CH4, f_CO2, w_CH4, w_CO2, w_alte, w_CH4_heil)
        
        Cpool.append(Cpool[-1] + delta[0])# hängt den Wert aus jedem Zeitschrit aus Cdec return [0] an.
        AltEpool.append(AltEpool[-1] + delta[1])
        Microben_CH4.append(Microben_CH4[-1] + delta[2])
        Microben_CO2.append(Microben_CO2[-1] + delta[3])
        Microben_AltE.append(Microben_AltE[-1] + delta[4])
        CH4.append(CH4[-1] + delta[5])
        CO2.append(CO2[-1] + delta[6])
        #Cused.append(Cused[-1] + delta[7])
        AceCO2.append(AceCO2[-1] + delta[7])
        Acetate.append(Acetate[-1] + delta[8])
        Microben_CH4_krank.append(Microben_CH4_krank[-1] + delta[9])
        
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
    return CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, Microben_CH4, Microben_CO2, Microben_AltE
