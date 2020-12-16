# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""

SOIL_DENSITY = 1.3 * 10**6 # g/m3 # 1.3 dry density for clay from Knoblauch data


def HeteroMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Vmax): 
    
    deltaSub1 = Sub1 / (Km1 + Sub1) * Vmax  * Biomass
    ATP = deltaSub1 * ATPprod
    deltaBiomass1 = ATP * Yatp # growth rate = deltaBiomass1 = ATPprod*Yatp*Substrate
    
    deltaSub2 = Sub2 / (5*Km1 + Sub2) * Vmax * Biomass
    ATP = deltaSub2 * ATPprod# - KostenEnzym
    deltaBiomass2 = ATP * Yatp
    
    deltaBiomass = deltaBiomass1 + deltaBiomass2
   
    return deltaBiomass, -deltaSub1, -deltaSub2

def AutoMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Km2, Vprod_max, Stoch): 
    
    deltaSub1 = Sub1 / (Km1 + Sub1) * Sub2 / (Km2 + Sub2) * Vprod_max  * Biomass
    ATP = deltaSub1 * ATPprod
    deltaBiomass = ATP * Yatp # growth rate = deltaBiomass1 = ATPprod*Yatp*Substrate
    
    deltaSub2 = Stoch*deltaSub1
           
    return deltaBiomass, -deltaSub1, -deltaSub2


def Fermenters(Biomass, Sub1, Sub2 , Vmax):
    
    ATPprod = 4 # stochiometrie
    Yatp = 10 # Fenchel 
    Km1 = 10 / SOIL_DENSITY # 10 from Song
    #Vmax = 0.5/ SOIL_DENSITY # 0.5 from Song
    
    Vmax = Vmax/SOIL_DENSITY
    
    deltaBiomass, deltaSub1, deltaSub2 = HeteroMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Vmax)
    
    
    return deltaBiomass, deltaSub1, deltaSub2







def Hydrotrophes(Biomass, CO2, H2):
    
    ATPprod = 2 # laut Fenchel maximalwert, theoretisch kleiner
    Yatp = 10 # Fenchel 
    Km1 = 0.05 / SOIL_DENSITY # 0.05 from Song
    Km2 = 0.01 / SOIL_DENSITY # 0.01 from Song
    Vprod_max = 0.15/ SOIL_DENSITY # 0.15 from Song
    Stoch = 4 # Stochiomitry: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O
    

    deltaBiomass, deltaCO2, deltaH2 = AutoMicrobe(Biomass, CO2, H2, ATPprod, Yatp, Km1, Km2, Vprod_max, Stoch)
    
    return deltaBiomass, deltaCO2, deltaH2




def Homo(Biomass, CO2, H2):
    
    ATPprod = 5 # Hugenholtz, eigentlich sollte es mehr sein 
    Yatp = 10 #Fenchel
    Km1 = 0.05 / SOIL_DENSITY # 0.05 from Song
    Km2 = 0.01 / SOIL_DENSITY # 0.01 from Song
    Vprod_max = 0.15/ SOIL_DENSITY # 0.15 from Song
    Stoch = 2 # weil nur ein substrat (Acetate),song/fenchel : 4H2 + 2CO2 → CH3COOH+ 2H2O, 2CO2 + 8H2  = CH3COOH +H2O. Laut Fenchel
    

    deltaBiomass, deltaCO2, deltaH2 = AutoMicrobe(Biomass, CO2, H2, ATPprod, Yatp, Km1, Km2, Vprod_max, Stoch)
    
    return deltaBiomass, deltaCO2, deltaH2



def Acetoclast(Biomass,Acetate):
    
    ATPprod = 1 #  weniger als Hydro
    Yatp = 10 #Fenchel
    Km1 = 0.05 / SOIL_DENSITY # song
    Km2 = 0 / SOIL_DENSITY # 
    Vprod_max = 0.5/ SOIL_DENSITY # song
    Stoch = 1
    

    deltaBiomass, deltaAcetate, _ = AutoMicrobe(Biomass, Acetate, 1, ATPprod, Yatp, Km1, Km2, Vprod_max, Stoch)
    
    return deltaBiomass, deltaAcetate 





def AltE(Biomass, Acetate, AltEpool):
    
    ATPprod = 2 # laut Fenchel maximalwert, theoretisch kleiner
    Yatp = 10 # Fenchel 
    Km1 = 0.01 / SOIL_DENSITY # 
    Km2 = 0 # damit alt e Pool keine michaelis menten gleichung hat 
    Vprod_max = 0.3/ SOIL_DENSITY # geschätzt
    Stoch = 1
    

    deltaBiomass, deltaAcetate, deltaAltE = AutoMicrobe(Biomass,  Acetate, AltEpool, ATPprod, Yatp, Km1, Km2, Vprod_max, Stoch)
    
    return deltaBiomass, deltaAcetate, deltaAltE




    
    