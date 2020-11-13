
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

#import numpy as np
#Cpool, AltEpool, Microben_CH4, Microben_CO2, Microben_AltE, CH4, CO2, AceCO2, Acetate, Microben_CH4_krank, f_CH4, f_CO2, f_alte,  w_CH4, w_CO2, w_alte, w_CH4_heil

def Cdec(Cpool, AltEpool, Microben_CH4, Microben_CO2, Microben_AltE, CH4, CO2, AceCO2, Acetate, Microben_CH4_krank, H2, f_CH4, f_CO2,  w_CH4, w_CO2, w_alte, w_CH4_heil):
    
     
    #Ferm atmen einen Teil Cpool und machen einen Teil zu Acetate
  
    Cpool_Ferm_Abbau =  (f_CO2 * Microben_CO2) 
    Cused_Ferm_resp = (1-0.67) * Cpool_Ferm_Abbau
    Cused_Ferm_Ace = 0.67 * Cpool_Ferm_Abbau # ac anteil vom gefressenen der in Acetate umgewandelt wird                 
    deltaMicroben_CO2 = w_CO2 * Microben_CO2 
    
    Ace_Ferm_prod = Cused_Ferm_Ace # Acetate das in diesem Schritt entsteht  
    H2_Ferm_prod = Cused_Ferm_resp/3 # /3 weil auf 33 CO2, 11 H2 bei Ferm entstehen
    
    #Acetate_tot = 0 if Acetate + Ace_Ferm_prod <= 0 else  Acetate + Ace_Ferm_prod
    Acetate_tot  = 0 if  Acetate <= 0 else Acetate
  
    # Homoacetogenese 4H2 + 2CO2→CH3COOH + 2H2O
     
    Ace_Homo_prod = min( H2/4, CO2/2)
    deltaH2_Homo = Ace_Homo_prod*4
    deltaCO2_Homo = Ace_Homo_prod*2

    
    # Alt e
    # IN: Acetat und Alt e 1:1
    # Out: CO2, Mukrobenwachstum
    # nur solange Alt e UND Acetat vorhanden
    

    f_alte = max(Acetate_tot * AltEpool, 0)
    Ace_used_AltE_resp = f_alte * Microben_AltE
    #Ace_used_AltE_resp = Acetate_tot if temp >= Acetate_tot  else  temp #1 Acetate wird zu zwei CO2
    #Ace_used_AltE_resp = f_alte * Microben_AltE * Acetate_tot  * 2 *AltEpool
    deltaMicroben_AltE = Microben_AltE * w_alte * f_alte
    deltaAltEpool = - min(Ace_used_AltE_resp, AltEpool) # pro Acetate 1 AltE- aufgebraucht?
    deltaAcetate_AltE = -min(Ace_used_AltE_resp, Acetate)

    Ace_used_Aceto = 0
    Ace_used_Aceto_resp_CH4 = 0
    Ace_used_Aceto_resp_CO2 = 0
    deltaMicroben_CH4 = 0
    deltaMicroben_CH4_krank = 0
    Microben_CH4_geheilt = 0
    deltaAcetate_Aceto = 0
        
    if AltEpool <= 0.01 :
        #Ace_used_AltE_resp = 0
        deltaMicroben_AltE = - min(Microben_AltE * w_alte, Microben_AltE)
        
           
        Microben_CH4_geheilt = 0 if Microben_CH4_krank <= 0 else (w_CH4_heil * Microben_CH4_krank)
        deltaMicroben_CH4_krank = - Microben_CH4_geheilt
        deltaMicroben_CH4 = w_CH4 *  Microben_CH4  +  Microben_CH4_geheilt
        
        
     # Aceto
    # IN: Acetat 
    # Out: CO2, CH4 1:1
    # nur solange Acetat und Microben_CH4 vorhanden    
        Ace_used_Aceto = min(Acetate_tot  * f_CH4 * Microben_CH4, Acetate)
        deltaAcetate_Aceto = -Ace_used_Aceto
        Ace_used_Aceto_resp_CH4 = 0.5 *  Ace_used_Aceto
        Ace_used_Aceto_resp_CO2 = 0.5 *  Ace_used_Aceto
        
        
      

        

    
        
    #nur Ferm verbrauchen Cpool C
    Cused_tot = Cpool_Ferm_Abbau 
        
    # Microbenwachstum
    deltaCpool = - Cused_tot
    #deltaMicroben_CH4 = h_CH4 * Microben_CH4
    
    
    
    # CH4 nur aus Acetate, CO2 aus Acetat und Ferm
    deltaCH4 = Ace_used_Aceto_resp_CH4 
    deltaCO2 = Ace_used_Aceto_resp_CO2  + Cused_Ferm_resp + Ace_used_AltE_resp - deltaCO2_Homo
   # deltaCused = Cused
    deltaAceCO2 = Ace_used_Aceto_resp_CO2
    deltaAcetate =  Ace_Ferm_prod + deltaAcetate_AltE + deltaAcetate_Aceto + Ace_Homo_prod
    deltaH2 = H2_Ferm_prod - deltaH2_Homo 
 
    return deltaCpool, deltaAltEpool, deltaMicroben_CH4, deltaMicroben_CO2, deltaMicroben_AltE, deltaCH4, deltaCO2, deltaAceCO2, deltaAcetate, deltaMicroben_CH4_krank, deltaH2



  



