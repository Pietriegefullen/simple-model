# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""

SOIL_DENSITY = 1.3 # g/cm3 # 1.3 dry density for clay from Knoblauch data

################## Heteromikrobe ##############################################
#Yatp: wie viel ATP kann der Mikrobentyp aus dem Substrattyp gewinnen
#ATP:  in mikromol
#Substrat : in Mikromol

def HeteroMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Vmax): 
    MM1 = Sub1 / (Km1 + Sub1) if Sub1 > 0 else 0 # vermeidet negative Werte
    deltaSub1 = MM1 * Vmax  * Biomass 
    ATP = deltaSub1 * ATPprod 
    deltaBiomass1 = ATP * Yatp 
    
#    #Stabil pool 
#    MM2 = Sub2 / (5*Km1 + Sub2) if Sub2 > 0 else 0 # die 5 ist erfunden!!!! 
#    deltaSub2 = MM2 * Vmax * Biomass #*100
#    ATP = deltaSub2 * ATPprod # - KostenEnzym
#    deltaBiomass2 = ATP * Yatp  #
    deltaSub2 = 0
 
    deltaBiomass = deltaBiomass1  #+ deltaBiomass2
    
    return deltaBiomass, -deltaSub1, -deltaSub2

################### Automikrobe ###############################################

def AutoMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Km2, Vprod_max, Stoch): 
    MM_1 = Sub1 / (Km1 + Sub1) if Sub1 > 0 else 0
    MM_2 = Sub2 / (Km2 + Sub2) if Sub2 > 0 else 0
    deltaSub1 = MM_1 * MM_2 * Vprod_max  * Biomass 
    ATP = deltaSub1 * ATPprod 
    deltaBiomass = ATP * Yatp 
   
    if deltaSub1  > Sub1: # Verhindert das aus leeren Pools geschöpft wird
        deltaSub1 = Sub1
        
    deltaSub2 = Stoch*deltaSub1
        
    if deltaSub2  > Sub2:
        deltaSub2 = Sub2
        deltaSub1 = deltaSub2/Stoch
       
    return deltaBiomass, -deltaSub1, -deltaSub2

###############################################################################

def Fermenters(Biomass, Sub1, Sub2 , Vmax_Ferm, ATPprod_Ferm,Yatp_Ferm):
    
    ATPprod = ATPprod_Ferm     # 3 kommt von Fenchel, 4 stochiometrie Mikromol
    Yatp = Yatp_Ferm           # 10e-3 # Fenchel mg C Biomass pro mikromol ATP 
    Km1 = 10 / SOIL_DENSITY    # 10 from Song  mikromol pro gram 
    #Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
    Vmax = Vmax_Ferm/SOIL_DENSITY
        
    deltaBiomass, deltaSub1, deltaSub2 = HeteroMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Vmax)
    
    return deltaBiomass, deltaSub1, deltaSub2



def AltE(Biomass, Acetate, AltEpool, Stoch_ALtE, Vprod_max_AltE,ATPprod_AltE,Yatp_AltE):
    
    ATPprod_AltE = ATPprod_AltE # 2 laut Fenchel maximalwert, theoretisch kleiner
    Yatp = Yatp_AltE            # 10e-3 # Fenchel mg pro mikromol ATP 
    Km1 = 0.01 / SOIL_DENSITY # 
    Km2 = 0                     # damit AltE Pool keine michaelis menten gleichung hat 
    #Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    Vprod_max = Vprod_max_AltE/ SOIL_DENSITY 
    Stoch = Stoch_ALtE

    deltaBiomass, deltaAcetate, deltaAltE = AutoMicrobe(Biomass,  Acetate, AltEpool, ATPprod_AltE, Yatp, Km1, Km2, Vprod_max, Stoch)
        
    return deltaBiomass, deltaAcetate, deltaAltE



def Hydrotrophes(Biomass, CO2, H2, ATPprod_Hydro, Yatp_Hydro):
    
    ATPprod_Hydro = ATPprod_Hydro # 2 laut Fenchel maximalwert, theoretisch kleiner
    Yatp = Yatp_Hydro             #10e-3 #Fenchel mg pro mikromol ATP  
    Km1 = 0.05 / SOIL_DENSITY     # 0.05 mikromol pro cm^3 from Song
    Km2 = 0.01 / SOIL_DENSITY     # 0.01 mikromol pro cm^3 from Song
    Vprod_max = 0.15 / SOIL_DENSITY # 0.15 mikromol pro cm^3 from Song
    Stoch = 4                     # Stochiomitry: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O
    
    deltaBiomass, deltaCO2, deltaH2 = AutoMicrobe(Biomass, CO2, H2, ATPprod_Hydro, Yatp, Km1, Km2, Vprod_max, Stoch)
    
    return deltaBiomass, deltaCO2, deltaH2



def Homo(Biomass, CO2, H2, ATPprod_Homo, Yatp_Homo):
    
    ATPprod_Homo = ATPprod_Homo  #5 Hugenholtz, eigentlich sollte es mehr sein 
    Yatp = Yatp_Homo             # 10e-3 #Fenchel mg pro mikromol ATP 
    Km1 = 0.05 / SOIL_DENSITY    # 0.05 from Song
    Km2 = 0.01 / SOIL_DENSITY    # 0.01 from Song
    Vprod_max = 0.15 / SOIL_DENSITY # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    Stoch = 2                    # song/fenchel : 4H2 + 2CO2 → CH3COOH+ 2H2O, 2CO2 + 8H2  = CH3COOH +H2O.
    
    deltaBiomass, deltaCO2, deltaH2 = AutoMicrobe(Biomass, CO2, H2, ATPprod_Homo, Yatp, Km1, Km2, Vprod_max, Stoch)
    
    return deltaBiomass, deltaCO2, deltaH2



def Acetoclast(Biomass, Acetate, ATPprod_Ace, Yatp_Ace):
    
    ATPprod_Ace = ATPprod_Ace   # weniger als Hydro, genauer Wert nicht gefunden
    Yatp =  Yatp_Ace            #10e-3 #Fenchel mg pro mikromol ATP 
    Km1 = 0.05 / SOIL_DENSITY   # 0.05 from song
    Km2 = 0 / SOIL_DENSITY      # 
    Vprod_max = 0.5/ SOIL_DENSITY # 0.5 from song
    Stoch = 1
    

    deltaBiomass, deltaAcetate, _ = AutoMicrobe(Biomass, Acetate, 1, ATPprod_Ace, Yatp, Km1, Km2, Vprod_max, Stoch)
    
    return deltaBiomass, deltaAcetate 






    
    