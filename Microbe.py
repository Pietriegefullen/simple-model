# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""
#mol/m3 = micromol/cm^3

SOIL_DENSITY = 1.3 # g/cm3 # 1.3 dry density for clay from Knoblauch data

################## Heteromikrobe ##############################################
#Yatp: wie viel ATP kann der Mikrobentyp aus dem Substrattyp gewinnen
#ATP:  in mikromol
# w : fasst Yatp und ATP zusammen, die Boundaries wurden entsprechend angepasst
#Substrat : in Mikromol
#Sensenmann: fixe Sterberate für alle Mikroben, 0.06 bei Song, ohne Begründung

def HeteroMicrobe(Biomass, Sub1, Sub2, w, Km1, Vmax, Sensenmann):  #ATPprod, Yatp, 
    MM1 = Sub1 / (Km1 + Sub1) if Sub1 > 0 else 0 # vermeidet negative Werte
    MM2 = Sub2 / (10 + Sub2) if Sub2 > 0 else 0 # die 1 ist ausgedacht
    deltaSub1 =  MM1 * Vmax  * Biomass 
    deltaBiomass1 = deltaSub1 * w * (1-MM2) #w setzt sich zusammen aus Yatp und ATPprod. die bekannten werte wurden verrechnet und als boundary conditions gewählt
    
#    #Stabil pool 
#    MM2 = Sub2 / (5*Km1 + Sub2) if Sub2 > 0 else 0 # die 5 ist erfunden!!!! 
#    deltaSub2 = MM2 * Vmax * Biomass #*100
#    ATP = deltaSub2 * ATPprod # - KostenEnzym
#    deltaBiomass2 = ATP * Yatp  #
    deltaSub2 = 0
 
    ToteMicroben = Biomass * Sensenmann
    deltaBiomass = deltaBiomass1 - ToteMicroben #+ deltaBiomass2
    
   # if Sub2 > 2:
    #    deltaBiomass1 = 0
    
    return deltaBiomass, -deltaSub1, -deltaSub2, ToteMicroben

################### Automikrobe ###############################################

def AutoMicrobe(Biomass, Sub1, Sub2, w, Km1, Km2, Vprod_max, Stoch, Sensenmann): # ATPprod, Yatp, 
    MM_1 = Sub1 / (Km1 + Sub1) if Sub1 > 0 else 0
    MM_2 = Sub2 / (Km2 + Sub2) if Sub2 > 0 else 0
    deltaSub1 = MM_1 * MM_2 * Vprod_max  * Biomass 
   
    if deltaSub1  > Sub1: # Verhindert dass aus leeren Pools geschöpft wird
        deltaSub1 = Sub1
        
    deltaSub2 = Stoch * deltaSub1
        
    if deltaSub2  > Sub2:
        deltaSub2 = Sub2
        deltaSub1 = deltaSub2/Stoch
    
    ToteMicroben = Biomass * Sensenmann
    deltaBiomass = deltaSub1 * w - ToteMicroben
        
    return deltaBiomass, -deltaSub1, -deltaSub2 , ToteMicroben

###############################################################################

def Fermenters(Biomass, Sub1, Sub2 , Vmax, w_Ferm, Sensenmann):
    Km1 = 10 / SOIL_DENSITY    # 10 from Song  mikromol pro gram 
    #Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
        
    deltaBiomass, deltaSub1, deltaSub2, ToteMicroben = HeteroMicrobe(Biomass, Sub1, Sub2, w_Ferm, Km1, Vmax, Sensenmann)
    
    return deltaBiomass, deltaSub1, deltaSub2, ToteMicroben



def AltE(Biomass, Acetate, AltEpool, Stoch_ALtE, Vprod_max,w_AltE,Sensenmann):
    Km1 = 0.005/SOIL_DENSITY   #0.01 / SOIL_DENSITY # wert nach Roden 2003 10 mal kleiner als bei Ace. Aber passt vlt nicht mehr mit den anderen werten zusammen
    Km2 = 0                     # damit AltE Pool keine michaelis menten gleichung hat 
    #Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    Stoch = Stoch_ALtE #jason 2001: ca 7 , weil acetate + 8Fe(III)  + 4H2O ---> 2HCO3 + 8Fe(II) + 9H+

    deltaBiomass, deltaAcetate, deltaAltE, ToteMicroben = AutoMicrobe(Biomass,  Acetate, AltEpool, w_AltE, Km1, Km2, Vprod_max, Stoch, Sensenmann)
        
    return deltaBiomass, deltaAcetate, deltaAltE, ToteMicroben



def Hydrotrophes(Biomass, CO2, H2, w_Hydro, Vprod_max,Sensenmann ):
    #Sensenmann = 0
    Km1 = 0.05 / SOIL_DENSITY     # 0.05 mikromol pro cm^3 from Song
    Km2 = 0.01 / SOIL_DENSITY     # 0.01 mikromol pro cm^3 from Song
    #Vprod_max = 0.15 / SOIL_DENSITY # 0.15 mikromol pro cm^3 from Song
    Stoch = 4                     # Stochiomitry: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O
    
    deltaBiomass, deltaCO2, deltaH2, ToteMicroben = AutoMicrobe(Biomass, CO2, H2, w_Hydro, Km1, Km2, Vprod_max, Stoch, Sensenmann)
    
    return deltaBiomass, deltaCO2, deltaH2, ToteMicroben


def Homo(Biomass, CO2, H2, w_Homo, Vprod_max, Sensenmann):
    #Sensenmann = 0
    Km1 = 0.05 / SOIL_DENSITY    # 0.05 from Song
    Km2 = 0.01 / SOIL_DENSITY    # 0.01 from Song
    #Vprod_max = 0.15 / SOIL_DENSITY # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    Stoch = 2                    # song/fenchel : 4H2 + 2CO2 → CH3COOH+ 2H2O, 2CO2 + 8H2  = CH3COOH +H2O.
    
    deltaBiomass, deltaCO2, deltaH2, ToteMicroben = AutoMicrobe(Biomass, CO2, H2, w_Homo, Km1, Km2, Vprod_max, Stoch,Sensenmann)
    
    return deltaBiomass, deltaCO2, deltaH2, ToteMicroben



def Acetoclast(Biomass, Acetate, w_Ace, Vprod_max, Sensenmann):
    Km1 = 0.5 / SOIL_DENSITY #0.05 / SOIL_DENSITY   # 0.05 from song, Wert sollte 10 mal höher sein als bei AltE laut roden2003
    Km2 = 0 / SOIL_DENSITY      # 
    #Vprod_max = 0.5/ SOIL_DENSITY # 0.5 from song
    Stoch = 0
    

    deltaBiomass, deltaAcetate, _, ToteMicroben = AutoMicrobe(Biomass, Acetate, 1, w_Ace, Km1, Km2, Vprod_max, Stoch,Sensenmann)
    
    return deltaBiomass, deltaAcetate, ToteMicroben






    
    