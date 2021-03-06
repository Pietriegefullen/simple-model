
# Das eigentliche mathematische Model, durch das der Verlauf der Kurven berechnet wird. 

# Parameter:
# Cpool: der frei Verfügbare C Gehalt zu jedem Zeitschritt im Boden
# Microbe: die Größe des Microbenpools zu jedem Zeitschritt
# CH4: cummulatives CH4 zu jedem Zeitschritt
# CO2: cummulatives CO2 zu jedem Zeitschritt
# Cv: die Abbaugeschwindigkeit abhängig von der Cpool größe und der Microbenmenge 
# und dem Abbauparameter k
# Cused: das im jeweiligen Zeitschritt abgebaute C

# Variabeln:
# k (Abbaugeschwindigkeit) und h (Microbenwachstumsgeschwindigkeit) sind die beiden zu optimierenden
# Parameter. Das im Zeitschritt abgebaute C (Cused) wird zu gleichen, festgesetzten Teilen in 
# CH4 und CO2 aufgeteilt. 

# Returns:
# Cpool, Microben, CH4, CO2, Cv



def Cdec(Cpool,Microben,CH4, CO2, k, h):
    Cv = Cpool * (k * Microben)
    Cused =  Cpool * Cv
    
    
   
    deltaCpool = - Cused
    deltaMicroben = Microben * (1+( h * Cpool))-Microben
    deltaCH4 = (Cused * 0.5) 
    deltaCO2 = (Cused * 0.5)
    deltaCv = Cv
    
 
    return deltaCpool, deltaMicroben, deltaCH4,deltaCO2, deltaCv
















