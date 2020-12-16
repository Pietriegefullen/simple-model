import matplotlib.pyplot as plt
from SimpleIN import Cdec # importiert das eigentliche mathematische Model 

# OPTIFUN 
# optifun bringt simplefun in ein format das von scipy für Curvefit benötigt wird

def optifun(xdata, Vmax):    
    xtage = xdata[0:int(len(xdata)/2)]
 
    CH4, CO2, _,_,_,_,_,_,_,_,_ = simplefun(xtage, Vmax)
    return CH4 + CO2 # setzt die für uns interessanten Ausgaben in curvefit format zusammen
    
# SIMPELFUN 
# Simpelfun ruft Cdec für die Länge von xdata auf 
# Parameter:
# xdata: unrelevnt für Cdec, wird aber von Curvefit als Vorraussetzung gesehen. 
# Variablen:
# f_CH4, f_CO2,  w_CH4, w_CO2, w_alte, w_CH4_heil
# Returns: CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, Microben_CH4, Microben_CO2, Microben_AltE
# CH4 und CO2 brauchen wir um den least square zu berechnen.
# CH4 und CO2 werden aneinandergehängt um einen abhängigen Least Sqaure zu berechnen und nicht einen pro Kurve



def simplefun(xtage, Vmax):    
    n = max(xtage)
    
    # Festgelegte Initialwerte
    m_gluc = 180 # molar mass of glucose, g pro mol
    m_cell = 162 # molar mass of cellulose, g pro mol
    TOC = 0.04 # Knoblauchs Daten , g dw
    labile = 0.005 # Knoblauchs Daten, g dw
    #Cpool_init = 50
    Cpool_init = TOC*labile/m_gluc + TOC*(1-labile)/m_cell

    AltE_init = 146 # cf. Yao, Conrad 1999
    M_A_CH4_init = 0.0002 # Monteux 2020
    M_CO2_init = 0.0002 # Monteux 2020
    M_AltE_init= 0.0002 # Monteux 2020
    M_A_CH4_krank_init = 0.0002 # Monteux 2020
    M_H_CH4_init = 0.0002 # Monteux 2020
    M_Homo_init = 0.0002 # Monteux 2020
    CH4_init = 0
    CO2_init = 0
    Acetate_init = 0.01
    AceCO2_init = 0
    H2_init = 0
    
    Cpool = [Cpool_init]
    AltEpool = [AltE_init]
    M_A_CH4 = [M_A_CH4_init]
    M_CO2= [M_CO2_init]
    M_AltE = [M_AltE_init]
    M_A_CH4_krank = [M_A_CH4_krank_init]
    M_H_CH4 =[M_H_CH4_init]
    M_Homo = [M_Homo_init]
    CH4 = [CH4_init]
    CO2 = [CO2_init]
    Acetate = [Acetate_init]
    AceCO2 = [AceCO2_init]
    H2 = [H2_init]

    
    for t in range(1,n+1): # iteriert über n+1 Zeitschritte (n = max(xtage))

            
      
        delta = Cdec(Cpool[-1], AltEpool[-1],M_A_CH4[-1],M_CO2[-1], M_AltE[-1],M_H_CH4[-1],M_Homo[-1], CH4[-1], CO2[-1], AceCO2[-1], Acetate[-1], M_A_CH4_krank[-1], H2[-1], Vmax)
        
        
        Cpool.append(Cpool[-1] + delta[0])
        AltEpool.append(AltEpool[-1] + delta[1])
        M_A_CH4.append(M_A_CH4[-1] + delta[2])
        M_CO2.append(M_CO2[-1] + delta[3])
        M_AltE.append(M_AltE[-1] + delta[4])
        M_H_CH4.append(M_H_CH4[-1] + delta[5])
        M_Homo.append(M_Homo[-1] + delta[6])
        CH4.append(CH4[-1] + delta[7])
        CO2.append(CO2[-1] + delta[8])
        #Cused.append(Cused[-1] + delta[7])
        AceCO2.append(AceCO2[-1] + delta[9])
        Acetate.append(Acetate[-1] + delta[10])
        M_A_CH4_krank.append(M_A_CH4_krank[-1] + delta[11])
        H2.append(H2[-1] + delta[12])

    
    CH4 = [CH4[i] for i in xtage]
    CO2 = [CO2[i] for i in xtage]
    return CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_CO2, M_AltE, M_H_CH4, H2
