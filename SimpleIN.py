
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



def Cdec(Cpool, AltEpool, Microben_CH4, Microben_CO2, Microben_AltE, FermCO2, Acetate, CH4, CO2, k_CH4, k_CO2, k_alte,  h_CH4, h_CO2, h_alte, ac):
    
     
    #Ferm atmen einen Teil Cpool und machen einen Teil zu Acetate
    Cpool_abbau_Ferm = Cpool* (k_CO2 * Microben_CO2)
    Cused_CO2_resp = (1-ac) * Cpool_abbau_Ferm
    Cused_CO2_ac = ac * Cpool_abbau_Ferm # ac anteil vom gefressenen der in Acetate umgewandelt wird                 
   
    Acetate = Cused_CO2_ac
    
    # AltE haben vorrecht auf Acetate consumption
    if AltEpool > 1:
        Cused_AltE_resp = Acetate * k_alte * Microben_AltE
        
        deltaAltEpool = -AltEpool * k_alte * Microben_AltE
        Cused_CH4ac_resp = 0
        Cused_CO2ac_resp = 0
        deltaMicroben_CH4 = 0
    else:   
        # Acetate geht zu gleichen Teilen in CH4 und CO2 (stochiometrie)
        Acetate_Dinner = Acetate * (k_CH4 * Microben_CH4)
        Cused_CH4ac_resp = 0.5 *  Acetate_Dinner
        Cused_CO2ac_resp = 0.5 *  Acetate_Dinner
        Cused_AltE_resp = 0
        deltaAltEpool = 0
        deltaMicroben_CH4 = h_CH4 * Microben_CH4
        
    #nur Ferm verbrauchen Cpool C
    Cused = Cused_CO2_resp + Cused_CO2_ac 
        
    # Microbenwachstum
    deltaCpool = - Cused
    #deltaMicroben_CH4 = h_CH4 * Microben_CH4
    deltaMicroben_CO2 = h_CO2 * Microben_CO2
    deltaMicroben_AltE = h_alte * Microben_AltE
    
    # CH4 nur aus Acetate, CO2 aus Acetat und Ferm
    deltaCH4 = Cused_CH4ac_resp 
    deltaCO2 = Cused_CO2ac_resp  + Cused_CO2_resp + Cused_AltE_resp
    deltaCused = Cused
    deltaFermCO2 = Cused_CO2ac_resp
    deltaAcetate = Acetate
    
 
    return deltaCpool, deltaAltEpool, deltaMicroben_CH4, deltaMicroben_CO2, deltaMicroben_AltE, deltaCH4, deltaCO2, deltaCused, deltaFermCO2, deltaAcetate
















