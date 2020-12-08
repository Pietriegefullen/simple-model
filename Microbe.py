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
   
    return deltaBiomass, deltaSub1, deltaSub2


def Fermenters(Biomass, Sub1, Sub2):
    
    ATPprod = 4 # stochiometrie
    Yatp = 10 # Fenchel 
    Km1 = 10 / SOIL_DENSITY # 10 from Song
    Vmax = 0.5/ SOIL_DENSITY # 0.5 from Song

    
    deltaBiomass, deltaSub1, deltaSub2 = HeteroMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Vmax)
    
    return deltaBiomass, deltaSub1, deltaSub2




def AutoMicrobe(Biomass, Sub1, Sub2, ATPprod, Yatp, Km1, Km2, Vprod_max, Stoch): 
    
    deltaSub1 = Sub1 / (Km1 + Sub1) * Sub2 / (Km2 + Sub2) * Vprod_max  * Biomass
    ATP = deltaSub1 * ATPprod
    deltaBiomass = ATP * Yatp # growth rate = deltaBiomass1 = ATPprod*Yatp*Substrate
    
    deltaSub2 = Stoch*deltaSub1
       
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



    
    