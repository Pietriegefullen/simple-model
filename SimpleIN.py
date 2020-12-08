
# Das eigentliche mathematische Model, durch das der Verlauf der Kurven berechnet wird. 

# Variabeln:
# Cpool: der frei Verfügbare C Gehalt zu jedem Zeitschritt im Boden
# Microbe: die Größe des Microbenpools zu jedem Zeitschritt
# CH4: cummulatives CH4 zu jedem Zeitschritt
# CO2: cummulatives CO2 zu jedem Zeitschritt
# Cv: die Abbaugeschwindigkeit abhängig von der Cpool größe und der Microbenmenge 
# und dem Abbauparameter k
# Cused: das im jeweiligen Zeitschritt abgebaute C

# Parameter:
# k (Abbaugeschwindigkeit) und h (Microbenwachstumsgeschwindigkeit), c4ant (anteil an C der zu CH4 umgewandelt wird)
# sind die zu optimierenden Parameter.
# Das im Zeitschritt abgebaute C (Cused) wird zu c4ant nach CH4 umgewandelt und zu 1- c4ant nach CO2. Zusammen 1

# Returns:
# Cpool, Microben, CH4, CO2, Cv
from Microbe import Fermenters, Hydrotrophes
#import numpy as np
#Cpool, AltEpool, Microben_CH4, Microben_CO2, Microben_AltE, CH4, CO2, AceCO2, Acetate, Microben_CH4_krank, f_CH4, f_CO2, f_alte,  w_CH4, w_CO2, w_alte, w_CH4_heil

#def Cdec(Cpool, AltEpool, M_A_CH4, M_CO2, M_AltE, M_H_CH4, M_Homo, CH4, CO2, AceCO2, Acetate, M_A_CH4_krank, H2, f_A_CH4, f_CO2, f_H_CH4, w_A_CH4, w_CO2, w_alte, w_A_CH4_heil, w_H_CH4, w_Homo):
def Cdec(C_lab, C_stab , AltEpool, M_A_CH4, M_CO2, M_AltE, M_H_CH4, M_Homo, CH4, CO2, AceCO2, Acetate, M_A_CH4_krank, H2,  w_A_CH4,  w_alte, w_A_CH4_heil, w_Homo):    
     
    #Ferm atmen einen Teil Cpool und machen einen Teil zu Acetate, CO2 und H2
    deltaM_CO2, deltaC_lab, deltaC_stab =   Fermenters(M_CO2, C_lab, C_stab)
    Ace_Ferm_prod = deltaC_lab + deltaC_stab
    CO2_Ferm_prod = Ace_Ferm_prod * 0.5 
    H2_Ferm_prod = Ace_Ferm_prod * (1/6)

    
    

    Acetate_tot  = 0 if  Acetate <= 0 else Acetate
  
#    # HOMO HOMO HOMO HOMO HOMO HOMO HOMO 4H2 + 2CO2→CH3COOH + 2H2O
#
    deltaCO2_Homo = 0
    deltaH2_Homo  = 0
    Ace_Homo_prod = 0
    
#    
    #HYDRO HYDRO HYDRO HYDRO 4H2 + CO2 ==> CH4 + 2H2O Fenchel -131kj/mol
    # RhCO2-C + 0:67RhH2--->CH4-C + 3RhH2O (Grant 1998)
   
    deltaM_H_CH4, deltaH2_Hydro, deltaCO2_Hydro =   Hydrotrophes(M_H_CH4, H2, CO2)
    CH4_M_H_CH4_prod = deltaCO2_Hydro # weil pro CO2 verbraucht, ein CH4 entsteht
    
    
    # ALT E ALT E ALT E ALT E 
    # IN: Acetat und Alt e 1:1
    # Out: CO2, Mikrobenwachstum
    # nur solange Alt e UND Acetat vorhanden
    

    f_alte = max(Acetate_tot * AltEpool, 0)
    Ace_used_AltE_resp = f_alte * M_AltE
    #Ace_used_AltE_resp = Acetate_tot if temp >= Acetate_tot  else  temp #1 Acetate wird zu zwei CO2
    #Ace_used_AltE_resp = f_alte * Microben_AltE * Acetate_tot  * 2 *AltEpool
    deltaM_AltE = min(M_AltE * w_alte * f_alte, 0 )
    deltaAltEpool = - min(Ace_used_AltE_resp, AltEpool) # pro Acetate 1 AltE- aufgebraucht?
    deltaAcetate_AltE = -min(Ace_used_AltE_resp, Acetate)

    Ace_used_Aceto = 0
    Ace_used_Aceto_resp_CH4 = 0
    Ace_used_Aceto_resp_CO2 = 0
    deltaM_A_CH4 = 0
    deltaM_A_CH4_krank = 0
    M_A_CH4_geheilt = 0
    deltaAcetate_Aceto = 0
    deltaM_Homo = 0  
    
    # ACETO ACETO ACETO ACETO 
    if AltEpool <= 0.01 :
        #Ace_used_AltE_resp = 0
        deltaM_AltE = - min(M_AltE * w_alte, M_AltE) # damit keine negativen Microben entstehen
        
           
        M_A_CH4_geheilt = 0 if M_A_CH4_krank <= 0 else (w_A_CH4_heil * M_A_CH4_krank)
        deltaM_A_CH4_krank = - M_A_CH4_geheilt
        deltaM_A_CH4 = w_A_CH4 *  M_A_CH4  +  M_A_CH4_geheilt
        
        
     # Aceto
    # IN: Acetat 
    # Out: CO2, CH4 1:1
    # nur solange Acetat und Microben_CH4 vorhanden 
        f_A_CH4 = 0.5 # 0.5 von Song/Grant1998
        Ace_used_Aceto = min(Acetate_tot  * f_A_CH4 * M_A_CH4, Acetate) 
        deltaAcetate_Aceto = -Ace_used_Aceto
        Ace_used_Aceto_resp_CH4 = 0.5 *  Ace_used_Aceto
        Ace_used_Aceto_resp_CO2 = 0.5 *  Ace_used_Aceto
        
        
      
        # HOMO HOMO HOMO HOMO HOMO HOMO HOMO 4H2 + 2CO2→CH3COOH + 2H2O
        # Homoacetogenese 4H2 + 2CO2→CH3COOH + 2H2O, 
        # 2CO2 + 8H2  = CH3COOH +H2O. Laut Fenchel
        #im If weil solange Alt E- hat Hydro Thermodynamische Vorteile, 
        #weil S und Fe H2 verbrauchen und der partialdruck zu gering ist (Ye 2013). Sonst 0
        f_Homo = 0.15 # Song/Grant 
        Ace_Homo_prod = min(H2/4, CO2/2) * f_Homo * M_Homo
        deltaH2_Homo = Ace_Homo_prod*4
        deltaCO2_Homo = Ace_Homo_prod*2
        deltaM_Homo = M_Homo * w_Homo

    
        
    #nur Ferm verbrauchen Cpool C
 
    
    
    # CH4 nur aus Acetate, CO2 aus Acetat und Ferm
    deltaCH4 = Ace_used_Aceto_resp_CH4 + CH4_M_H_CH4_prod
    deltaCO2 = Ace_used_Aceto_resp_CO2  + CO2_Ferm_prod + Ace_used_AltE_resp - deltaCO2_Homo  - deltaCO2_Hydro
   # deltaCused = Cused
    deltaAceCO2 = Ace_used_Aceto_resp_CO2
    deltaAcetate =  Ace_Ferm_prod + deltaAcetate_AltE + deltaAcetate_Aceto + Ace_Homo_prod
    deltaH2 = H2_Ferm_prod - deltaH2_Homo - deltaH2_Hydro
    deltaM_H_CH4 =  M_H_CH4 
   # deltaM_Homo = deltaM_Homo
    
 
    return deltaC_lab, deltaC_stab, deltaAltEpool, deltaM_A_CH4, deltaM_CO2, deltaM_H_CH4, deltaM_AltE, deltaM_Homo, deltaCH4, deltaCO2, deltaAceCO2, deltaAcetate, deltaM_A_CH4_krank, deltaH2



  



