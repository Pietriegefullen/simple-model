
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



def Cdec(Cpool,Microben_CH4, Microben_CO2, CH4, CO2, k_CH4, k_CO2,  h_CH4, h_CO2, c4ant):
    Cused = (Cpool * (k_CH4 * Microben_CH4)) + (Cpool * (k_CO2 * Microben_CO2))  
    Cused_CH4 = Cpool * (k_CH4 * Microben_CH4)
    Cused_CO2 = Cpool * (k_CO2 * Microben_CO2)
   
    deltaCpool = - Cused
    deltaMicroben_CH4 = h_CH4 * Microben_CH4
    deltaMicroben_CO2 = h_CO2 * Microben_CO2
    
    deltaCH4 = (Cused_CH4 * c4ant) 
    deltaCO2 = (Cused_CH4 * (1-c4ant)) + Cused_CO2
    deltaCused = Cused
    
 
    return deltaCpool, deltaMicroben_CH4, deltaMicroben_CO2, deltaCH4, deltaCO2, deltaCused
















