# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""
import math as math
import numpy as np

# laut van1999efFe3ct findet kein Wachstum statt, laut philben2020anaerobic schon
#mol/m3 = micromol/cm^3

# YATP ist 10g Biomasse pro 1 MOl ATP. 
# ATPprod sind ca 4 für Ferm.... 

SOIL_DENSITY = 1.3 # g/cm3 # 1.3 dry density for clay from Knoblauch data
m_C = 12.01*1e-3 # mg/micromol molar mass of carbon

def thermodynamics_Ferm(product_dict, microbe_dict):
    # inverse mm erlaubt ? weil aceate durch enzyme abgebaut werden was einer MM folgt und Acetate unsere Hemmung auslöst. die auflösung der Hemmung 
    # ist also gleichgesetzt mit dem Abbau des acetate und die Hemmung ist deshalb der inverse abbau von Acetate
    # print('thermo fermo')
    # print(product_dict)
    Acetate = product_dict['Acetate']['concentration']
    # print(microbe_dict)
    # input()
    KmA = microbe_dict['KmA_Ferm']
    # print('Fermentations ',1- Acetate/(KmA + Acetate))
    # input('...')
    return 1- Acetate/(KmA + Acetate)



def thermodynamics(educt_dict, product_dict):
    
    # print('thermononfermo')
    # print(educt_dict)
    # print(product_dict)
    # input()
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
        
    Edu_Q_total = np.prod(Edu_Q)
    Edu_DGf_total = sum(Edu_DGf)
      
    for product, produ_dict in product_dict.items():  
     # to get Q values of the Products
        
      
        Concentration = produ_dict['concentration']
        Stoch = produ_dict['Stoch']
        
        Prod_Q.append(Concentration**Stoch)
        
        # to get Gibbs Energy of formation of the Products
        
        Prod_DGf.append(produ_dict['DGf'])
        
    Prod_Q_total = np.prod(Prod_Q)
    Prod_DGf_total = sum(Prod_DGf)
        
           
    #print('Q', Prod_Q_total, Edu_Q_total)
    # input()
    if Prod_Q_total == 0:
        return 1.0
    
    Q = Prod_Q_total/ Edu_Q_total
    DGs = Prod_DGf_total - Edu_DGf_total
    
    R = 8.31446261815324 	# in J⋅K−1⋅mol−1
    T = 277.15 # 4 °C in Kelvin
    
    DGr = DGs + R * T * math.log( Q )
    # print('DGr', DGr)
    DGmin = -26
    return 1 - np.exp(min(0, DGr - DGmin)/R*T)

def GeneralPathway(microbe_dict, educt_dict, product_dict, pathway_name = ''):
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
    #print('life really sucks right now')
    for educt, edu_dict in educt_dict.items():  
        if edu_dict['concentration'] <= 0:
            return dict()
    
    #print('entering '+pathway_name)
    # input()
    MM_factors = list()
    for educt, edu_dict in educt_dict.items():  
        Concentration = edu_dict['concentration']
        substance_MM = Concentration/(edu_dict['Km'] + Concentration) if Concentration > 0 else 0
        MM_factors.append(substance_MM)
    # Note that for Ferm_help, Km must be 0 so that MM_factors evaluates to 1.0
    MM_factors_total = np.prod(MM_factors)
#-------------------------------------------------------------------------------------------------------------
    # die tatsächliche Stoffwechselrate, gegeben die termodynamischen und kinetischen hindernisse
    Vmax = microbe_dict['Vmax']
    Biomass = microbe_dict['concentration']
    if 'Ferm' in microbe_dict and microbe_dict['Ferm']==True:
        thermodynamic_factor = thermodynamics_Ferm(product_dict, microbe_dict) #TODO Hier gehört die acetatehemmung 
        MMB = Biomass  if Biomass > 0 else 0 

    elif 'Ferm_help' in microbe_dict and microbe_dict['Ferm_help']==True:
        thermodynamic_factor = 1.0
        MMB = Biomass / (microbe_dict['Kmb'] + Biomass)  if Biomass > 0 else 0 

    else:
        thermodynamic_factor = thermodynamics(educt_dict, product_dict)
        MMB = Biomass  if Biomass > 0 else 0 

    
#-------------------------------------------------------------------------------------------------------------            
    V = Vmax * MM_factors_total  * MMB * thermodynamic_factor# micromol ???
#-------------------------------------------------------------------------------------------------------------    
    print(microbe_dict['microbe']   ) 
    print(V)
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
    # print(actual_reaction_rate)
    # input('actual reactual')

    pool_change_dict = dict()
    for educt, edu_dict in educt_dict.items():
        normalized_stoich = edu_dict['Stoch']/Edu_Bezug_stoich
        pool_change_dict[educt] = - normalized_stoich*actual_reaction_rate
        #print(educt)
        #print(Edu_Bezug_stoich)
        #print(pool_change_dict[educt])
        
    # wie viel der jeweiligen Produkte werden produziert    
    for product, produ_dict in product_dict.items():
        normalized_stoich = produ_dict['Stoch']/Edu_Bezug_stoich
        pool_change_dict[product] = normalized_stoich*actual_reaction_rate
       # print(product) 
       # print(Edu_Bezug_stoich)
       # print(pool_change_dict[product] )
        
    w = microbe_dict['growth_rate']
    dead_microbes = Biomass*microbe_dict['death_rate']
    biomass_change = Biomass*w*actual_reaction_rate - dead_microbes # !TODO massenbilanz ausgleichen

    pool_change_dict['biomass'] = biomass_change
    # print(biomass_change)
    # input('biomass '+pathway_name)
    return pool_change_dict





###############################################################################
###############################################################################

def Ferm_help_Pathway(pool_dict,model_parameter_dict):
    # print('Ferm')
    # input('..')

    microbe_dict = {'concentration' : pool_dict['M_Ferm'], 
                    'Vmax'          : model_parameter_dict['Vmax_help_Ferm'],          #Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
                    'growth_rate'   : 0,
                    'death_rate'    : 0,
                    'Kmb'           : model_parameter_dict['Kmb_help_Ferm'],
                    'Ferm_help'     : True, # ist nur wichtig das es den key 'Ferm_help' gibt 
                    'microbe'       : 'Ferm_help'}
    
    educt_dict =  {'C'              : {'concentration'  : pool_dict['C'],
                                       'Stoch'          : 1,
                                       'Km'             : 0} }# 10 from Song  mikromol pro gram 
                     
    
    product_dict = { 'DOC' : {'concentration': pool_dict['DOC'],
                                  'Stoch'        : 1               }}
                                                                  
   
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Ferm')

    if 'biomass' in pool_change_dict:
         pool_change_dict['M_Ferm'] = pool_change_dict.pop('biomass')

    return pool_change_dict










# TODO Ferm pathway mit acetatehemmung schreibe


def Ferm_Pathway(pool_dict,model_parameter_dict):
    # print('Ferm')
    # input('..')

    microbe_dict = {'concentration' : pool_dict['M_Ferm'], 
                    'Vmax'          : model_parameter_dict['Vmax_Ferm'],          #Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
                    'growth_rate'   : model_parameter_dict['w_Ferm'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'KmA_Ferm'      : model_parameter_dict['KmA_Ferm'],
                    'Ferm'          : True,
                    'microbe'       : "Ferm" }
    
    
    educt_dict =  {'DOC'           : {'concentration':pool_dict['DOC'],
                                    'Stoch'          : 6                 ,                                    
                                    'Km'             : 10 / SOIL_DENSITY} }# 10 from Song  mikromol pro gram 
                     
    
    product_dict = { 'Acetate' : {'concentration': pool_dict['Fe2'],
                                  'Stoch'        : 6               }  ,
                    
                     'CO2'      : {'concentration': pool_dict['CO2'],
                                   'Stoch'        : 3              }  , 
                         
                      'H2'     : {'concentration': pool_dict['H2'],
                                  'Stoch'        : 1             }}
                                                                  
   
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Ferm')

    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Ferm'] = pool_change_dict.pop('biomass')


    return pool_change_dict










def Fe3_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY:   C2H3O2 − + 4H2O + 8Fe3(III)   --->       9H+ + 2 HCO3− + 8Fe3+2   Delattre 2019 , jason 2001
    #print('Fe3')
    microbe_dict = {'concentration' : pool_dict['M_Fe3'], 
                    'Vmax'          : model_parameter_dict['Vmax_Fe3'],          #Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
                    'growth_rate'   : model_parameter_dict['w_Fe3'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Fe3'}
    
    
    educt_dict =  {'Acetate'       : {'concentration':pool_dict['Acetate'],
                                    'Stoch'          : 1                  ,
                                    'DGf'            : -369.31            ,
                                    'Km'             : 0.01 / SOIL_DENSITY}, #0.01 / SOIL_DENSITY # wert nach Roden 2003 10 mal kleiner als bei Ace (0.8). Aber passt vlt nicht mehr mit den anderen werten zusammen
    
                      'Fe3'        :{'concentration':pool_dict['Fe3'] ,
                                   'Stoch'          : 8                   ,
                                   'DGf'            : -4.7                ,
                                   'Km'             : 0                   }} # TODO WISO 0 ? 
                     
    
    product_dict = { 'Fe2' : {'concentration': pool_dict['Fe2'],
                              'Stoch'        : 8               ,
                              'DGf'          : -78.9           }  ,
                    
                     'CO2' : {'concentration': pool_dict['CO2'],
                              'Stoch'        : 2               ,
                              'DGf'          :   0              }  , #TODO DIESE ZAHL IST FALSCH
                         
                      'H2' : {'concentration': pool_dict['H2'],
                              'Stoch'        : 0              ,
                              'DGf'          :  0              }   }#TODO DIESE ZAHL IST FALSCH
                                                                  
   
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Fe3')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Fe3'] = pool_change_dict.pop('biomass')
    
    #print(pool_change_dict)
    
    return pool_change_dict





def Hydro_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O conrad2000selective, Fe3nchel -131kj/mol
    # print('Hydro')
    microbe_dict = {'concentration' : pool_dict['M_Hydro']              , 
                    'Vmax'          : model_parameter_dict['Vmax_Hydro'], ## 0.15 mikromol pro cm^3 from Song
                    'growth_rate'   : model_parameter_dict['w_Hydro']   , 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Hydro'}
    
    educt_dict = { 'H2'  : {'concentration':pool_dict['H2']  ,
                            'Stoch'     : 4                  , 
                            'DGf'       : 0                  ,
                            'Km'        : 0.01 / SOIL_DENSITY},# 0.01 mikromol pro cm^3 from Song
    
                  'CO2'  :{'concentration': pool_dict['CO2'] ,
                           'Stoch'        : 1                                ,
                           'DGf'          : -394.36                          ,
                           'Km'           : 0.05/SOIL_DENSITY }              } # 0.05 mikromol pro cm^3 from Song
                    
    
    product_dict = {'CH4' :{'concentration': pool_dict['CH4'] ,
                            'Stoch'        : 1                ,
                            'DGf'          : -50.72}          }
   
    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Hydro')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Hydro'] = pool_change_dict.pop('biomass')
    # print(pool_change_dict['CH4'])
    
    return pool_change_dict



def Homo_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY: song,Fe3nchel,conrad2000selective : 4H2 + 2CO2 → CH3COOH+ 2H2O. Alternativ: 2CO2 + 8H2  = CH3COOH +H2O., laut Thauer sind es 4 H2 + Co2 -> CH4 +2H2O
    # print('Homo')
    microbe_dict = {'concentration' : pool_dict['M_Homo'], 
                    'Vmax'          : model_parameter_dict['Vmax_Homo'] , # # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
                    'growth_rate'   : model_parameter_dict['w_Homo']    , 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Homo'}
    
    educt_dict = { 'H2'  : {'concentration':pool_dict['H2']        ,
                             'Stoch'       : 4                     ,
                             'DGf'         : 0                     ,
                               'Km'        : 0.01 / SOIL_DENSITY   },    # 0.01 from Song
    
                  'CO2'  :{'concentration' :pool_dict['CO2']       ,
                           'Stoch'         : 2                     ,
                           'DGf'           : -394.36               ,
                           'Km'            : 0.05 / SOIL_DENSITY } }  # 0.05 from Song, laut (van1999efFe3cts) größer als Hydro, laut schink1997energetics sollte der Wert mit sinkender Temp, mit zunehmendem Acetate und sinkendem PH sinken (Im vlg zu Hydro)
                     
    
    product_dict = {'Acetate' : {'concentration': pool_dict['Acetate'] ,
                                 'Stoch'        : 1                    ,
                                 'DGf'          : -369.31}             , }
   
    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Homo')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Homo'] = pool_change_dict.pop('biomass')
    
    return pool_change_dict



def Ac_Pathway(pool_dict,model_parameter_dict):
    #PATHWAY:  CH3COO + H+ ->  CH4 + CO2 Fe3y_Conrad2000
    # print('Ac')
    microbe_dict = {'concentration' : pool_dict['M_Ac'], 
                    'Vmax'          : model_parameter_dict['Vmax_Ac'],  #Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
                    'growth_rate'   : model_parameter_dict['w_Ac'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Ac'}
    
    educt_dict = { 'Acetate' : {'concentration':pool_dict['Acetate']    ,
                                'Stoch'        :  1                     ,
                                 'DGf'         : -369.31                ,
                                 'Km'          :  0.05 / SOIL_DENSITY   }} #0.05 / SOIL_DENSITY   # 0.05 from song, Wert sollte 10 mal höher sein als bei AltE laut roden2003 bei 12Mikromol !!!!
    
                      
    product_dict = { 'CH4' : {'concentration': pool_dict['CH4'],
                              'Stoch'        : 1               ,
                              'DGf'          : 186.26          },
                    
                     'CO2' : {'concentration': pool_dict['CO2'],
                              'Stoch'        : 1               ,
                              'DGf'          :-394.36          } }

    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Ac')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Ac'] = pool_change_dict.pop('biomass')
    
    
    return pool_change_dict



# =============================================================================
# #%%
# from scipy.io import savemat
# a = [1,2,3,'dfjalkf']
# 
# test_dict = { "Moli": {"Buy": 'Apples',"Sell": a ,"Quantity": 300},
#              "Anna":  {"Buy": 55,"Sell": 83,"Quantity": 154}}
# 
# 
# 
# test_dict1 =  {"Buy": 75,"Sell": 53,"Quantity": 300}
# 
# 
# savemat("test_dict.mat", test_dict)
# 
# 
# 
# 
# test_dict2 =  {"Buy": 55,"Sell": 83,"Quantity": 154}
# test_dict3 =  {"Chicken": 55,"Sell": 83,"Quantity": 154}
# 
# test_list =[test_dict1, test_dict2, test_dict3]
# 
# 
# 
# test_summms = sum(list([dictionary['Buy'] for dictionary in test_list if 'Buy' in dictionary]))
# print(test_summms)
# #%%
# =============================================================================


# =============================================================================
# for key, item in test_dict.items():
#     if ('Buy' in item ) :
#         print('pear')
# 
# 
# if ('Moli' in test_dict):
#     print(test_dict[(list(test_dict.keys())[0])]['Buy'] )
#     print('banana')
# 
# if test_dict[(list(test_dict.keys())[0])]['Buy'] in test_dict:
#     print(test_dict[(list(test_dict.keys())[0])]['Buy'] )
#     print('banana')
# 
# 
# =============================================================================
