# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""

# laut van1999effect findet kein Wachstum statt, laut philben2020anaerobic schon
#mol/m3 = micromol/cm^3

# YATP ist 10g Biomasse pro 1 MOl ATP. 
# ATPprod sind ca 4 für Ferm.... 

SOIL_DENSITY = 1.3 # g/cm3 # 1.3 dry density for clay from Knoblauch data
m_C = 12.01*1e-3 # mg/micromol molar mass of carbon

################## Heteromikrobe ##############################################
#Yatp: wie viel ATP kann der Mikrobentyp aus dem Substrattyp gewinnen
#ATP:  in mikromol
# w : fasst Yatp und ATP zusammen, die Boundaries wurden entsprechend angepasst
#Substrat : in Mikromol
#Sensenmann: fixe Sterberate für alle Mikroben, 0.06 bei Song, ohne Begründung

def HeteroMicrobe(Biomass, Sub1, Sub2, w, Km1, Vmax, Sensenmann, Kmb, Kmh):  #ATPprod, Yatp, 
    MM1 = Sub1 / (Km1 + Sub1) if Sub1 > 0 else 0 # vermeidet negative Werte
    MM2 = Sub2 / (Kmh + Sub2) if Sub2 > 0 else 0 # die 10 ist ausgedacht
    MMB = Biomass / (Kmb + Biomass) if Biomass >0 else 0 

    deltaSub1Resp =  MM1 * Vmax * MMB
    deltaSub1Grow = deltaSub1Resp * w/m_C 
    deltaSub1 = deltaSub1Resp + deltaSub1Grow

 
    ToteMicroben = Biomass * Sensenmann
    deltaBiomass = deltaSub1Grow * m_C * (1-MM2) - ToteMicroben
    # schink1997energetics, Acetate hemmt ab 10 mikromol (für ein bestimmtes Bakterium)
    deltaSub2 = 0
 
    return deltaBiomass, -deltaSub1, -deltaSub2, ToteMicroben

################### Automikrobe ###############################################

def AutoMicrobe(Biomass, Sub1, Sub2, w, Km1, Km2, Vprod_max, Stoch, Sensenmann, Kmb): # ATPprod, Yatp, 
    MM_1 = Sub1 / (Km1 + Sub1) if Sub1 > 0 else 0
    MM_2 = Sub2 / (Km2 + Sub2) if Sub2 > 0 else 0
    MMB = Biomass / (Kmb + Biomass) if Biomass >0 else 0

    deltaSub1Resp = MM_1 * MM_2  * Vprod_max * Biomass #* MMB # micromol
    deltaSub1Grow = deltaSub1Resp * w/m_C     # micromol 
    deltaSub1 = deltaSub1Resp + deltaSub1Grow 
    deltaSub2 = Stoch * deltaSub1 
    
    limited = False
    if deltaSub2 > Sub2:
        deltaSub1 = deltaSub2/ Stoch
        limited = True
    
    if deltaSub1 > Sub1:
        deltaSub1 = Sub1
        limited = True
    
    if limited == True:  
        deltaSub2 = Stoch * deltaSub1
        deltaSub1Resp = deltaSub1/(1+w/m_C)
        deltaSub1Grow = deltaSub1 - deltaSub1Resp
        
    ToteMicroben = Biomass * Sensenmann
    deltaBiomass = deltaSub1Grow * m_C - ToteMicroben
        
    return deltaBiomass, -deltaSub1, -deltaSub2 , ToteMicroben

###############################################################################

def Fermenters(Biomass, Sub1, Sub2 , Vmax, w_Ferm, Sensenmann, Kmb, Kmh):
    Km1 = 10 / SOIL_DENSITY    # 10 from Song  mikromol pro gram 
    #Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
        
    deltaBiomass, deltaSub1, deltaSub2, ToteMicroben = HeteroMicrobe(Biomass, Sub1, Sub2, w_Ferm, Km1, Vmax, Sensenmann, Kmb, Kmh)
    
    return deltaBiomass, deltaSub1, deltaSub2, ToteMicroben



def AltE(Biomass, Acetate, AltEpool, Stoch_ALtE, Vprod_max,w_AltE,Sensenmann, Kmb):
    Km1 = 0.005/SOIL_DENSITY   #0.01 / SOIL_DENSITY # wert nach Roden 2003 10 mal kleiner als bei Ace (0.8). Aber passt vlt nicht mehr mit den anderen werten zusammen
    Km2 = 0                     # damit AltE Pool keine michaelis menten gleichung hat 
    #Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    Stoch = Stoch_ALtE #jason 2001: ca 7 , weil acetate + 8Fe(III)  + 4H2O ---> 2HCO3 + 8Fe(II) + 9H+

    deltaBiomass, deltaAcetate, deltaAltE, ToteMicroben = AutoMicrobe(Biomass,  Acetate, AltEpool, w_AltE, Km1, Km2, Vprod_max, Stoch, Sensenmann, Kmb)

    return deltaBiomass, deltaAcetate, deltaAltE, ToteMicroben



def Hydrotrophes(Biomass, CO2, H2, w_Hydro, Vprod_max,Sensenmann, Kmb ):
    # thauer1993reactions unter feldbedingungen 1 mol ATP pro Mol CH4, 
    #Sensenmann = 0
    Km1 = 0.05 / SOIL_DENSITY     # 0.05 mikromol pro cm^3 from Song
    Km2 = 0.01 / SOIL_DENSITY     # 0.01 mikromol pro cm^3 from Song
    #Vprod_max = 0.15 / SOIL_DENSITY # 0.15 mikromol pro cm^3 from Song
    Stoch = 4                     # Stochiomitry: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O
    
    deltaBiomass, deltaCO2, deltaH2, ToteMicroben = AutoMicrobe(Biomass, CO2, H2, w_Hydro, Km1, Km2, Vprod_max, Stoch, Sensenmann, Kmb)
    
    return deltaBiomass, deltaCO2, deltaH2, ToteMicroben


def Homo(Biomass, CO2, H2, w_Homo, Vprod_max, Sensenmann):
    #Sensenmann = 0
    Km1 = 0.05 / SOIL_DENSITY    # 0.05 from Song, laut (van1999effects) größer als Hydro, laut schink1997energetics sollte der Wert mit sinkender Temp, mit zunehmendem Acetate und sinkendem PH sinken (Im vlg zu Hydro)
    Km2 = 0.01 / SOIL_DENSITY    # 0.01 from Song
    #Vprod_max = 0.15 / SOIL_DENSITY # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    Stoch = 2                    # song/fenchel : 4H2 + 2CO2 → CH3COOH+ 2H2O, 2CO2 + 8H2  = CH3COOH +H2O., laut Thauer sind es 4 H2 + Co2 -> CH4 +2H2O
    
    deltaBiomass, deltaCO2, deltaH2, ToteMicroben = AutoMicrobe(Biomass, CO2, H2, w_Homo, Km1, Km2, Vprod_max, Stoch,Sensenmann, 0)
    # Die Reaktionsrichtung ändert sich abhängig vom H2 Druck cabrol2017microbial
    
    return deltaBiomass, deltaCO2, deltaH2, ToteMicroben



def Acetoclast(Biomass, Acetate, w_Ace, Vprod_max, Sensenmann, Kmb):
    Km1 = 0.05 / SOIL_DENSITY #0.05 / SOIL_DENSITY   # 0.05 from song, Wert sollte 10 mal höher sein als bei AltE laut roden2003 bei 12Mikromol !!!!
    Km2 = 0 / SOIL_DENSITY      # 
    #Vprod_max = 0.5/ SOIL_DENSITY # 0.5 from song
    Stoch = 0
    

    deltaBiomass, deltaAcetate, _, ToteMicroben = AutoMicrobe(Biomass, Acetate, 1, w_Ace, Km1, Km2, Vprod_max, Stoch,Sensenmann, Kmb)
    
    return deltaBiomass, deltaAcetate, ToteMicroben






    
    