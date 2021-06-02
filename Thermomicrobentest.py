# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""
import math as math

# laut van1999effect findet kein Wachstum statt, laut philben2020anaerobic schon
#mol/m3 = micromol/cm^3

# YATP ist 10g Biomasse pro 1 MOl ATP. 
# ATPprod sind ca 4 für Ferm.... 

SOIL_DENSITY = 1.3 # g/cm3 # 1.3 dry density for clay from Knoblauch data
m_C = 12.01*1e-3 # mg/micromol molar mass of carbon

def thermodynamics(educt_dict, product_dict):
    Edu_Q = list()
    Edu_DGf = list()
    
    Prod_Q= list()
    Prod_DGf =list()

    for educt, edu_dict in educt_dict.items():  
       # to get Q values of Educts
        
        Concentration = edu_dict['concentration']
        Stoch = edu_dict['Stoch']
        
        Edu_Q.append(Concentration**Stoch)
        
        # to get Gibbs Energy of formation of the Educts
        Edu_DGf.append(edu_dict['DGf'])
        
    Edu_Q_total = numpy.prod(Edu_Q)
    Edu_DGf_total = sum(Edu_DGf)
      
    for product, produ_dict in product_dict.items():  
     # to get Q values of the Products
        
      
        Concentration = produ_dict['concentration']
        Stoch = produ_dict['Stoch']
        
        Prod_Q.append = Concentration**Stoch
        
        # to get Gibbs Energy of formation of the Products
        
        Prod_DGf.append(produ_dict['DGf'])
        
    Prod_Q_total = numpy.prod(Prod_Q)
    Prod_DGf_total = sum(Prod_DGf)
        
           

    Q = Prod_Q_total/ Edu_Q_total
    DGs = Prod_DGf_total - Edu_DGf_total
    
    R = 8.31446261815324 	# in J⋅K−1⋅mol−1
    T = 277.15 # 4 °C in Kelvin
    
    DGr = DGs + R * T * math.log( Q )
    
    return 1 - np.exp(min(0, (DGr - (-26)/R*T)))

def GeneralPathway(microbe_dict, educt_dict, product_dict):
    """
    Parameters
    ----------
    microbe_dict : dict
        'concentration': Menge der vorhandenen Biomasse
        'Vmax': Maximale Reaktionsgeschwindigkeit, bezogen auf den 
                Elektronenspender
        'growth_rate': Wachstumsrate der Biomasse, bezogen auf die Menge
                       der umgesetzten Elektronenspender
        'death_rate': Sterberate, in Abhängigkeit der vorhandenen Biomasse
        'Kmb': half-saturation coefficient für die inverse MM-Kinetik
        
    educt_dict : dict(dict)
        für jedes Edukt dieser Reaktion ein dict:
        educt: dict
            'concentration'
            'Stoch'
            'DGf'
            'Km'
            
    product_dict : dict(dict)
        für jedes Produkt dieser Reaktion ein dict:
        product : dict
            'concentration'
            'Stoch'
            'DGf'

    Returns
    -------
    pool_change_dict : dict
        für jeden beteiligten pool einen eintrag, z.B.
        'CO2': -45.7
        außerdem einen Eintrag
        'biomass': +34.8

    """
    # name as string like "Ferm"
    
    MM_factors = list()
    for educt, edu_dict in educt_dict.items():  
        Concentration = edu_dict['concentration']
        substance_MM = Concentration/(edu_dict['Km'] + Concentration) if Concentration > 0 else 0
        MM_factors.append(substance_MM)
        
    MM_factors_total = numpy.prod(MM_factors)
    
    # erstellen der inversen MM der Biomasse
    Biomass = microbe_dict['concentration']
    MMB = Biomass / (microbe_dict['Kmb'] + Biomass)  if Biomass > 0 else 0 
#-------------------------------------------------------------------------------------------------------------
    # die tatsächliche Stoffwechselrate, gegeben die termodynamischen und kinetischen hindernisse
    Vmax = microbe_dict['Vmax']
    
    #TODO: die fermentation hat keine thermodyn., stattdessen einen Hemmfaktor
    
    thermodynamic_factor = thermodynamics(educt_dict, product_dict)
    V = Biomass * Vmax * MM_factors_total  * MMB * thermodynamic_factor# micromol ???
#-------------------------------------------------------------------------------------------------------------    
    
    # Mikrobenzuwachs und verbrauch von C Biomasse , ACHTUNG!! DAS IST EIGENTLICH NICHT UNBEDINGT KORRECKT, weil nur in den seltensten
    # Fällen der E_Donor auch die C Quelle für Mikrobenwachstum ( nur bei Acetate und Ferm evtl. )
    #deltaSub1Grow = deltaSub1Resp * w/m_C     # micromol 
    
    # erstellen der Veränderungen der Konzentrationen  des Bezugsstoffs
    #Edu_change_dict = {Edu_Bezug_name : deltaSub1Resp + deltaSub1Grow}
        
    
    Edu_Bezug_name = list(educt_dict.keys())[0]
    Edu_Bezug_stoich = educt_dict[Edu_Bezug_name]['Stoch']
    
    # der Bezugstoff ist der erste im jeweiligen Dict und ist immer mit 1 angenommen. 
    # z.b wird 1 Acetate umgewandelt, auch wenn die stoich eigentlich 2 ist
    limiting_dict = dict()
    for educt, edu_dict in educt_dict.items():
        normalized_stoich = edu_dict['Stoch']/Edu_Bezug_stoich
        number_of_reactions = edu_dict['concentration']/normalized_stoich
        limiting_dict[educt] = number_of_reactions
    
    limiting_educt = min(limiting_dict, key = limiting_dict.get)
    limiting_reaction_rate = max(0,limiting_dict[limiting_educt])

    actual_reaction_rate = min(V,limiting_reaction_rate)
    
    pool_change_dict = dict()
    for educt, edu_dict in educt_dict.items():
        normalized_stoich = edu_dict['Stoch']/Edu_Bezug_stoich
        pool_change_dict[educt] = - normalized_stoich*actual_reaction_rate
        
    # wie viel der jeweiligen Produkte werden produziert    
    for product, produ_dict in product_dict.items():
        normalized_stoich = produ_dict['Stoch']/Edu_Bezug_stoich
        pool_change_dict[product] = normalized_stoich*actual_reaction_rate
         
    w = microbe_dict['growth_rate']
    dead_microbes = Biomass*microbe_dict['death_rate']
    biomass_change = w*actual_reaction_rate - dead_microbes

    pool_change_dict['biomass'] = biomass_change
        
    return pool_change_dict

def AltEPathway():
    biomass_dict = 
    educt_dict = 
    product_dict = 
    pool_change_dict = GeneralPathway(biomass_dict, educt_dict, product_dict)
    pool_change_dict['M_Fe'] = pool_change_dict.pop('biomass')
    return pool_change_dict








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








