# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""
import math as math
import numpy as np

np.seterr(divide='raise', over = 'warn', under = 'warn')

from order import Henrys_dict, enthalpy, Gibbs_formation

####################### Definition globaler Konstanten ########################

SOIL_DENSITY = 1.3 # g/cm3 # 1.3 dry density for clay from Knoblauch data
m_C = 12.01*1e-3 # mg/micromol, molar mass of carbon
T_C = 4 # hier die Temperatur in Celius eingeben
T = T_C + 273.15 # von Celcius nach Kelvin  
T0 = 25 + 273.15

######################## Berechnung der Reaktionen ############################

def henrys_law (H_cp_Standard, H_cp_temp):
    # H2, CO2 und CH4 steigen als Gase in den Headspace. Aus der Konzentration, der Temperatur und den H_cp - Werten 
    # berechnen wir den gelösten Anteil in der auatischen Phase, der zur Reaktion zur Verfügung steht.
    # werte aus Sander
    
    T_standard = 298.15 # muss in Kelvin sein
    
    H_cp_temp_adjusted = H_cp_temp - T # Standard Temp Wert - actual Temp
      
    H_cc = H_cp_Standard *  np.exp(H_cp_temp_adjusted * ((1/T) - (1/ T_standard)))

    return H_cc

def H_plus_funk (pool_dict):
    pH = pool_dict['pH']                                                        # in mol/L 
    total_water_volume_in_flask = pool_dict['water']                            # ml
    total_dry_soil_in_flask = pool_dict['weight'] #pool_dict['weigth']          # g
    water_volume = total_water_volume_in_flask/total_dry_soil_in_flask          # ml/gdw
    H_plus_per_liter = (10**-pH) * 1e6                                          # mikromol/ l
    H_plus = (H_plus_per_liter/1000) * water_volume                             # Mikromol / gdw Konzentration der vorhandenen H+ Ionen in Mikromol
    H_plus = H_plus*1e6 #TODO das * 1e6 ist ausgedacht damit alles funktioniert.
    return H_plus
     
def thermodynamics_Ferm(product_dict, microbe_dict):
    # 1- Acetate/(KmA + Acetate) je mehr Acetate desto weniger schnell läuft die Fermentation 
    
    # inverse MM erlaubt ? weil aceate durch enzyme abgebaut werden was einer MM folgt und Acetate unsere Hemmung auslöst. die auflösung der Hemmung 
    # ist also gleichgesetzt mit dem Abbau des acetate und die Hemmung ist deshalb der inverse abbau von Acetate

    Acetate = product_dict['Acetate']['concentration']
    #TODO je niedriger der pH (je weiter unter 7.4 (dem zellinneren))desto stärker die Hemmung durch Acetate pinhal2019acetate

    # TODO hier haben wir keine Temperaturabhängigkeit
    Inhibition = microbe_dict['Inhibition_Ferm']
   # print( '1- Acetate/(Inhibition + Acetate)', 1- Acetate/(KmA + Acetate))

    
    return 1- Acetate/(Inhibition + Acetate) # als gegengleichung zu pinhal2019acetate so-called inhibition index (i)


def GeneralPathway(microbe_dict, educt_dict, product_dict, pathway_name = ''):
    """
    Parameters
    ----------
    microbe_dict : dict
        'concentration': Menge der vorhandenen Biomasse
        'Vmax': Maximale Reaktionsgeschwindigkeit, bezogen auf den 
                Elektronenspender
        'growth_rate': Wachstumsrate der Biomasse, bezogen auf die Menge
                       der umgesetzten Elektronenspender (veraltet)
        'death_rate': Sterberate, in Abhängigkeit der vorhandenen Biomasse
        'Kmb': half-saturation coefficient für die inverse MM-Kinetik
         'Ferm_help': nur für die if condition
         'microbe' : Name der Mikrobe
         'CUE': Carbon use efficiency der Mikrobe
        
           
    educt_dict : dict(dict)
        für jedes Edukt dieser Reaktion ein dict:
        educt: dict
            'concentration'
            'Stoch'
            'DGf'
            'Km'
            'C_atoms'  : Anzahl der C Atome in der C source
            
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
        'DGr' : die DGr der jeweiligen Mikrobe

    """
    
    
#--------------- Alle Edukte die INNERHALB einer Mikrobe enzymatisch aufgespalten werden folgen MM----------------
    
    MM_factors_total = 1.0 # Fall es keine MM gibt, gibt es auch keine Hemmung. Also 1
    for educt, edu_dict in educt_dict.items():  
        # Alle Educte werden innerhalb der Mikroben durch Enyme aufgespalten
        # Enzymatische Reaktionen folgen einer MM-Gleichung, die Substratlimitiert ist. 
        # Innerhalb der Mikrobe ist genug Enzym vorhanden, nicht genug Substrat
        Concentration = edu_dict['concentration']  #mikromol/gdw
        substance_MM = Concentration/(edu_dict['Km'] + Concentration) if Concentration > 0 else 0
        MM_factors_total *= substance_MM # alle MM_Gleichungen der Edukte werden in einen MM_factor verrechnent
    # ------------ 
    # Alle Edukte (Corg), die außerhalb der Mikrobe enzymatisch aufgespalten werden sind nicht Substrat limitiert,
    # sondern Enzymlimitiert. das kommt später über die Invers MM zu tragen
    # For Ferm_help, Km must be 0 so that MM_factors evaluates to 1.0, Process is Enzyme Limited not Substrate limited
    # ------------

#------------------------------Berechnung thermodynamischer Faktor aller Pathways----------------------------------------------
#---------------------------und Berechnung des Biomassefaktors für alle Pathways   ------------------------------------------
    DGr_Ausgabe = 0
    Biomass = microbe_dict['concentration'] # Mikrobielle Biomasse in mg Mikrobielles C /gdw
    MMB = max(0,Biomass)
    ferm_factor = 1
    
    if 'Ferm_help' in microbe_dict and microbe_dict['Ferm_help'] and Biomass > 0:
        MMB = Biomass / (microbe_dict['Kmb'] + Biomass)
        # Die Biomasse dient als Proxy für die Limitierung an Exoenzymen
    
    if 'Ferm' in microbe_dict and microbe_dict['Ferm']:
        # Berechnung für den endoenzymatischen Prozess der Fermentation 
        ferm_factor = thermodynamics_Ferm(product_dict, microbe_dict)
        # Die Ferm folgt einer simplen Acetate-Hemmung ohne genaue Angaben

        
#-----------------------------Berechnung der Reaktionsgeschwindigkeit--------------------------------------------
             
    Vmax = microbe_dict['Vmax'] # die Maximale Rate, wenn alle Umweltumstände ideal sind    
#-------------------------------------------------------------------------------------------------------------            
    V = Vmax * MM_factors_total  * MMB * ferm_factor # die tatsächliche Stoffwechselrate, gegeben die termodynamischen und kinetischen Hindernisse
    
#-------------------------------------------------------------------------------------------------------------    
   
#-----------------------------Berechnung der Stoffumsatzmengen abhänging von V und der Stoichiometrie------------------
   
    Edu_Bezug_name = list(educt_dict.keys())[0]
    Edu_Bezug_stoich = educt_dict[Edu_Bezug_name]['Stoch']
    
#-------------------------------------------------------------------------------------------------------    

#-----------------------------Berechnung mit der Carbon use efficiency------------------------------------------------      
    #print('microbe----------------------------------------', microbe_dict['microbe'])
    C_for_growth = 0
    pool_change_dict = dict()
    for educt, edu_dict in educt_dict.items():
        normalized_stoich = edu_dict['Stoch']/Edu_Bezug_stoich
        CUE = 0
        
        # für die C Quelle der jeweiligen Mikrobe wird die CUE aus dem dict geholt
        if 'C_source' in microbe_dict and educt == microbe_dict['C_source']:    
            CUE = microbe_dict['CUE']
            
            respiration = normalized_stoich*V # falls Edu_Bezug nicht C_source ist
            total_metabolized = respiration / (1-CUE) # siehe nächste Zeilen 
            # CUE = growth/total = growth/ (resp+ growth) = growth = total* CUE
            # => respiration = total - growth = total-total*CUE = total* (1-CUE) 
            # => total = respiration/ (1-CUE) 
            growth_metabolized = CUE * total_metabolized # Anteil für growth von der Gesamtmenge
            
            C_for_growth = growth_metabolized * edu_dict['C_atoms'] # mikromol C für Wachstum aus der C_source
        
        pool_change_dict[educt] = - normalized_stoich*V / (1-CUE)       

    # wie viel der jeweiligen Produkte werden produziert, korrigiert nach CUE 
    for product, produ_dict in product_dict.items():
        normalized_stoich = produ_dict['Stoch']/Edu_Bezug_stoich # Produktion pro Reaktion 
        if normalized_stoich*V < 1e-30: # damit der solver nicht crashed
            pool_change_dict[product] = 0
        else:
            pool_change_dict[product] = normalized_stoich*V # Tatsächliche Produktion 
 
#-------------------------------------------------------------------------------------------------------

#-----------------------------Berechnung der Biomass change-------------------------------------------         
        
    dead_microbes = 0#Biomass*microbe_dict['death_rate'] # Tote Mikroben
    
    biomass_change = - dead_microbes  +  C_for_growth * m_C #  m_C Gewicht 1 Mol C =  12.011 g /mol 
    # das Gewicht der toten Mikroben geht ab, währen alles C_for_growth als Gewicht hinzukommt. Keine w_growth mehr
    # w_growth geht sozusagen über CUE ein 
    pool_change_dict['biomass'] = biomass_change
    
    
    pool_change_dict['MM'] = MM_factors_total
#-------------------------------------------------------------------------------------------------------    
 
    return pool_change_dict


################################# MIKROBEN ####################################


def Ferm_help_Pathway(pool_dict,model_parameter_dict):
    # Ferm_help ist keine wirklich eigenständige Mikrobe, sondern nur eine Hilfsmikrobe die die Exoenzym menge
    # der Fermentation repräsentiert
    
    microbe_dict = {'concentration' : pool_dict['M_Ferm'], 
                    'Vmax'          : model_parameter_dict['Vmax_help_Ferm'], 
                    #'growth_rate'   : 0,
                    'death_rate'    : 0,
                    'Kmb'           : model_parameter_dict['Kmb_help_Ferm'], # MM Faktor für die Exoenzyme
                    'Ferm_help'     : True       , # ist nur wichtig das es den key 'Ferm_help' gibt 
                    'microbe'       : 'Ferm_help',
                    'CUE'           :       0    , # weil Ferm_help nicht wächst (nur Ferm wächst)
                    'pH'            : pool_dict['pH']}
    
    educt_dict =  {'C'              : {'concentration'  : pool_dict['C'],
                                       'Stoch'          : 1,
                                       'Km'             : 10000/SOIL_DENSITY} # For Ferm_help, Km must be 0 so that MM_factors evaluates to 1.0, Process is Enzyme Limited not Substrate limited 
                                                                 }
    
    product_dict = { 'DOC' : {'concentration': pool_dict['DOC'],
                                  'Stoch'        : 1       }
                                                                 }
                                                                  
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Ferm')

    if 'biomass' in pool_change_dict:
         pool_change_dict['M_Ferm'] = pool_change_dict.pop('biomass') # Ferm_help benutzt die Biomasse von Ferm. 
         
    return pool_change_dict


def Ferm_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY: 1 DOC ---->  6 Acetate + 3 CO2 + 1 H2                           # Ratio aus Grant 1998

    microbe_dict = {'concentration' : pool_dict['M_Ferm'], 
                    'Vmax'          : model_parameter_dict['Vmax_Ferm'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'Inhibition_Ferm': model_parameter_dict['Inhibition_Ferm'],        # MM Faktor für die Acetatehemmung
                    'Ferm'          : True,
                    'microbe'       : "Ferm" ,
                    'CUE'           :   0.5,
                    'C_source'       : 'DOC'}
    
    
    educt_dict =  {'DOC'           : {'concentration':pool_dict['DOC'],
                                    'Stoch'          : 6                 ,                                    
                                    'Km'             : 100 / SOIL_DENSITY,      # 10 from Song  mikromol pro gram 
                                    'C_atoms'        : 6                 }}    # weil glucose (und andere Monomere) 6 C atome hat und ein momomer aus der spaltung von Coellulose ist 
                        
#------------------ für alle Gase wird die auqatische Phase berechnet   ------------------------------------ 
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_H2 = henrys_law( Henrys_dict['H2']['H_cp_Standard'],  Henrys_dict['H2']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_H2 =  pool_dict['H2']* (H_cc_H2/  (H_cc_H2+1))
    
    dissolved_CO2_total = dissolved_CO2 
    
#------------------------------------------------------------------------------------------------------------                                       
    #laut Grant A :6 , CO2:3, H2 : 1
    #laut Knoblauch A :6 , CO2:3, H2 : 6
    #laut mir DOC = A: 7 ,CO2 :6 H2:12 (teile durch 2 weil dann näher am Grant verhältnis)
    product_dict = { 'Acetate' : {'concentration': pool_dict['Acetate'],
                                  'Stoch'        : 3.5             }  ,
                    
                     'CO2'      : {'concentration': dissolved_CO2_total,
                                   'Stoch'        :  3            }  , 
                         
                      'H2'     : {'concentration': dissolved_H2   ,
                                  'Stoch'        : 6          }}   #6          # die 6 kommt aus Gesprächen von Christian und Christian
                                                                  
   
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Ferm')

    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Ferm'] = pool_change_dict.pop('biomass')  
        pool_change_dict['dissolved_CO2_total'] = dissolved_CO2_total
        pool_change_dict['dissolved_H2'] = dissolved_H2
        
        


    return pool_change_dict


def Fe3_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY:   C2H3O2 − + 4H2O + 8Fe3(III)   --->       9H+ + 2 HCO3− + 8Fe2+   Delattre 2019 , jason 2001
    # H2 + 2Fe3+ -> 2H+ +2Fe2+ aus der Vorlesung von Christian K.
    microbe_dict = {'concentration' : pool_dict['M_Fe3'], 
                    'Vmax'          : model_parameter_dict['Vmax_Fe3'],    
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Fe3',
                    'CUE'           :  0.5  ,
                    'C_source'      :  'Acetate'}
                        
    educt_dict =  {'Acetate'       : {'concentration':pool_dict['Acetate'],
                                    'Stoch'          : 1, #1                  ,
                                    'Km'             : 0.01/ SOIL_DENSITY,    # 0.01 / SOIL_DENSITY # wert nach Roden 2003 10 mal kleiner als bei Ace (0.8). Aber passt vlt nicht mehr mit den anderen werten zusammen
                                    'C_atoms'        : 2}, # kJ/mol
    
                      'Fe3'        :{'concentration':pool_dict['Fe3'] ,
                                   'Stoch'          : 8                   ,
                                   'Km'             : 2 / SOIL_DENSITY }}#, # TODO WISO 0 ? 0 wäre keine Hemmung
                      

    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))

    product_dict = { 'Fe2' : {'concentration': pool_dict['Fe2'],
                              'Stoch'        : 8               } ,
                     'CO2' : {'concentration': dissolved_CO2,
                              'Stoch'        : 2               }}
                                                                  
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Fe3')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Fe3'] = pool_change_dict.pop('biomass')
    
    return pool_change_dict


def Hydro_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O conrad2000selective, Fe3nchel -131kj/mol
    # Pathway beer2007transport: CO2 (aq) + 4H2 (aq)   - > CH4 (aq) + 2H2O (l)

    microbe_dict = {'concentration' : pool_dict['M_Hydro']              , 
                    'Vmax'          : model_parameter_dict['Vmax_Hydro'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Hydro',
                    'CUE'           :  0.5,
                    'C_source'      :  'CO2'}
                      
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_H2 = henrys_law( Henrys_dict['H2']['H_cp_Standard'],  Henrys_dict['H2']['H_cp_temp'])
    H_cc_CH4 = henrys_law(Henrys_dict['CH4']['H_cp_Standard'], Henrys_dict['CH4']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_H2 =  pool_dict['H2']* (H_cc_H2/  (H_cc_H2+1))#*1e104
    dissolved_CH4 = pool_dict['CH4']*(H_cc_CH4/ (H_cc_CH4+1))
                                       
    educt_dict = { 'H2'  : {'concentration':dissolved_H2 ,
                            'Stoch'     : 4                  , 
                            'Km'        : 0.0001 / SOIL_DENSITY},                # 0.01 mikromol pro cm^3 from Song
    
                  'CO2'  :{'concentration': dissolved_CO2,
                           'Stoch'        : 1 ,
                           'Km'           : 0.005/SOIL_DENSITY  ,               # 0.05 mikromol pro cm^3 from Song
                           'C_atoms'      : 1                 }}
    
    product_dict = {'CH4' :{'concentration': dissolved_CH4 ,                 
                            'Stoch'        : 1                }}
                                                            
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Hydro')
        
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Hydro'] = pool_change_dict.pop('biomass')
        
    pool_change_dict['hydro_diss_co2'] = dissolved_CO2
    pool_change_dict['hydro_diss_h2'] = dissolved_H2
    
        
    return pool_change_dict


def Homo_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY: song,Fe3nchel,conrad2000selective : 4H2 + 2CO2 → CH3COOH+ 2H2O. Alternativ: 2CO2 + 8H2  = CH3COOH +H2O., laut Thauer sind es 4 H2 + Co2 -> CH4 +2H2O
    # Pathway beer2007transport :CO2 (aq) + 2H2 (aq) -> 1/2CH3COO?(aq) + 1/2H+(aq) + H2O(l)
    microbe_dict = {'concentration' : pool_dict['M_Homo'], 
                    'Vmax'          : model_parameter_dict['Vmax_Homo'] , 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Homo',
                    'CUE'           :       0.5 ,
                    'C_source'      :  'CO2'}
                    
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_H2 = henrys_law( Henrys_dict['H2']['H_cp_Standard'],  Henrys_dict['H2']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_H2 =  pool_dict['H2']* (H_cc_H2/  (H_cc_H2+1))     #*1e104         
    
    educt_dict = { 'H2'  : {'concentration': dissolved_H2       ,
                             'Stoch'       : 4                     ,
                               'Km'        : 0.001 / SOIL_DENSITY   },          # 0.01 from Song
    
                  'CO2'  :{'concentration' : dissolved_CO2         ,
                           'Stoch'         : 2                     ,
                           'Km'            : 0.005 / SOIL_DENSITY   ,           # 0.05 from Song, laut (van1999efFe3cts) größer als Hydro, laut schink1997energetics sollte der Wert mit sinkender Temp, mit zunehmendem Acetate und sinkendem PH sinken (Im vlg zu Hydro)
                           'C_atoms'       : 1                 }}
    
    product_dict = {'Acetate' : {'concentration': pool_dict['Acetate'] ,
                                 'Stoch'        : 1                   }}
    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Homo')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Homo'] = pool_change_dict.pop('biomass')
   
    return pool_change_dict


def Ac_Pathway(pool_dict,model_parameter_dict):
    #PATHWAY:  CH3COO + H+ ->  CH4 + CO2 Fe3y_Conrad2000
    #CH3COO (aq) + H+ (aq) ->CO2 (aq) + CH4 (aq), beer2007transport

    microbe_dict = {'concentration' : pool_dict['M_Ac'], 
                    'Vmax'          : model_parameter_dict['Vmax_Ac'],  
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Ac',
                    'CUE'           :  0.5,
                    'C_source'      :  'Acetate'}
    
    educt_dict = { 'Acetate' : {'concentration': pool_dict['Acetate']   ,
                                'Stoch'        :  1                     ,
                                 'Km'          :  0.2 / SOIL_DENSITY   ,      #0.05 / SOIL_DENSITY   # 0.05 from song, Wert sollte 10 mal höher sein als bei AltE laut roden2003 bei 12Mikromol !!!!
                                 'C_atoms'      : 2                 }}
    
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_CH4 = henrys_law(Henrys_dict['CH4']['H_cp_Standard'], Henrys_dict['CH4']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 + 1))
    dissolved_CH4 = pool_dict['CH4']*(H_cc_CH4/ (H_cc_CH4 + 1))

    product_dict = { 'CH4' : {'concentration': dissolved_CH4, 
                              'Stoch'        : 1               },
                    
                     'CO2' : {'concentration': dissolved_CO2,
                              'Stoch'        : 1               }}
    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Ac')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Ac'] = pool_change_dict.pop('biomass')

    return pool_change_dict

