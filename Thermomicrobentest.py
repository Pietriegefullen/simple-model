# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:38:25 2020

@author: Lara
"""
import math as math
import numpy as np

from order import Henrys_dict

####################### Definition globaler Konstanten 

SOIL_DENSITY = 1.3 # g/cm3 # 1.3 dry density for clay from Knoblauch data
m_C = 12.01*1e-3 # mg/micromol molar mass of carbon
T = 4 # hier die Temperatur in Celius eingeben
T = T + 273.15 # von Celcius nach Kelvin  

########################

def henrys_law (H_cp_Standard, H_cp_temp):
    # H2, CO2 und CH4 steigen als Gase in den Headspace. Aus der Konzentration, der Temperatur und den H_cp - Werten 
    # berechnen wir den gelösten Anteil in der auatischen Phase, der zur Reaktion zur Verfügung steht.
    
    T_standard = 298.15 # muss in Kelvin sein
    
    H_cp_temp_adjusted = H_cp_temp - T # Standard Temp Wert - actual Temp
      
    H_cc = H_cp_Standard *  np.exp(H_cp_temp_adjusted * ((1/T) - (1/ T_standard)))

    return H_cc
    
   
def thermodynamics_Ferm(product_dict, microbe_dict):
    # 1- Acetate/(KmA + Acetate) je mehr Acetate desto weniger schnell läuft die Fermentation 
    
    # inverse MM erlaubt ? weil aceate durch enzyme abgebaut werden was einer MM folgt und Acetate unsere Hemmung auslöst. die auflösung der Hemmung 
    # ist also gleichgesetzt mit dem Abbau des acetate und die Hemmung ist deshalb der inverse abbau von Acetate

    Acetate = product_dict['Acetate']['concentration']

    # TODO hier haben wir keine Temperaturabhängigkeit
    KmA = microbe_dict['KmA_Ferm']

    
    return 1- Acetate/(KmA + Acetate)

def thermodynamics(educt_dict, product_dict, microbe_dict):
    # die Berechnung des Thermodynamischen Faktors für alle Mikroben, die klare Edukte und Produkte haben
    
    #--------------------------Reaction Quotient Q---------------------
    DGf_educt = 0
    Q_educts = 1.0
    for name, educt in educt_dict.items():
        # Berechnung von Q_educts  für den Reaction quotient Q 
        DGf_educt += educt['DGf'] # zur Berechnung von DGs
        Q_educts *= (1e-6*educt['concentration'])**educt['Stoch'] # concentration must be MOL, Zur Berechnung von DGr
        
    DGf_product = 0
    Q_products = 1.0
    for name, product in product_dict.items():
        # Berechnung von Q_products  für den Reaction quotient Q 
        DGf_product += product['DGf']# zur Berechnung von DGs
        Q_products *= (1e-6*product['concentration'])**product['Stoch'] # concentration must be MOL, Zur Berechnung von DGr
        
     #--------------------------   
    
    if Q_educts <= 0:
      #falls ein edukt fehlt, kann die Reaktion nicht stattfinden und der thermofactor muss 0 sein
       return 0.0, 110000000.0 # 110000000.0  random wert der im plot auffallen soll
   
    if Q_products <= 0:
        #TODO falls ein Produkt fehlt, liegt das Gewicht sehr weit rechts und die Reaktion ist ungehemmt. 
        return 1.0, 100000000.0 # 110000000.0  random wert der im plot auffallen soll
  
  
    #--------------------------   
    

    #-------------------------Berechnung von Delta Gibbs Energien---------------------

    DGs = DGf_product - DGf_educt # DGs (Delta Gibbs Standart der Reaktion)
    
    R = 8.31446261815324 	   # in J⋅K−1⋅mol−1, Gaskonstante ,Einheit passt zur DGmin
    #T  ist global definiert
    
    DGr = DGs + R * T * np.log( Q_products/Q_educts )  # DGr: Delta Gibbsenergie der Reaktion in J⋅mol-1
   
    DGmin = -26.*1e3          # J⋅mol-1, Einheit passt zur Gaskonstante - 26 in kJ/mol aus z.b. blodau2011thermodynamic

    if DGr >= DGmin: # wenn keine Energie gewonnen wird findet die Reaktion nicht statt
        #print('minimum not reached, return 0', DGr - DGmin)
        return 0.0, DGr
    
    
    return 1 - np.exp((DGr - DGmin)/(R*T)) , DGr # DGr_Ausgabe nur fürs plotting

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
    
    for educt, edu_dict in educt_dict.items():  
        # Falls ein Edukt < 0 ist wird ein leeres Dict zurückgegeben,
        # weil die Reaktion nicht stattfindet, alle changes sind 0
        if edu_dict['concentration'] <= 0:
            return dict()
    

    MM_factors_total = 1.0 # Fall es keine MM gibt, gibt es auch keine Hemmung. Also 1
    for educt, edu_dict in educt_dict.items():  
        # Alle Educte werden innerhalb der Mikroben durch Enyme aufgespalten
        # Enzymatische Reaktionen folgen einer MM-Gleichung, die Substratlimitiert ist. 
        Concentration = edu_dict['concentration']
        substance_MM = Concentration/(edu_dict['Km'] + Concentration) if Concentration > 0 else 0
    
        MM_factors_total *= substance_MM # alle MM_Gleichungen der Edukte werden in einen MM_factor verrechnent
        
    # For Ferm_help, Km must be 0 so that MM_factors evaluates to 1.0, Process is Enzyme Limited not Substrate limited
    
#-------------------------------------------------------------------------------------------------------------
    # die tatsächliche Stoffwechselrate, gegeben die termodynamischen und kinetischen hindernisse
    Vmax = microbe_dict['Vmax']
    DGr_Ausgabe = 0
    Biomass = microbe_dict['concentration']
    if 'Ferm' in microbe_dict and microbe_dict['Ferm']==True:
        thermodynamic_factor = thermodynamics_Ferm(product_dict, microbe_dict) #TODO Hier gehört die acetatehemmung 
        MMB = Biomass  if Biomass > 0 else 0 

    elif 'Ferm_help' in microbe_dict and microbe_dict['Ferm_help']==True:
        thermodynamic_factor = 1.0
        MMB = Biomass / (microbe_dict['Kmb'] + Biomass)  if Biomass > 0 else 0 

    else:
        #thermodynamic_factor = 1
        thermodynamic_factor, DGr_Ausgabe = thermodynamics(educt_dict, product_dict, microbe_dict)
        MMB = Biomass  if Biomass > 0 else 0 
           
#-------------------------------------------------------------------------------------------------------------            
    V = Vmax * MM_factors_total  * MMB * thermodynamic_factor# micromol ???
#-------------------------------------------------------------------------------------------------------------    
    #print(microbe_dict['microbe']   ) 
    #print(V)
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

    C_for_growth = 0
    pool_change_dict = dict()
    for educt, edu_dict in educt_dict.items():
        normalized_stoich = edu_dict['Stoch']/Edu_Bezug_stoich
        CUE = 0

        if 'C_source' in microbe_dict and educt == microbe_dict['C_source']:    
            CUE = microbe_dict['CUE']
            
            total_metabolized = normalized_stoich*actual_reaction_rate / (1-CUE)
            growth_metabolized = CUE * total_metabolized 
            
            C_for_growth = growth_metabolized * edu_dict['C_atoms'] # mikromol C für Wachstum 
        pool_change_dict[educt] = - normalized_stoich*actual_reaction_rate / (1-CUE)
       
        
     
    # wie viel der jeweiligen Produkte werden produziert    
    for product, produ_dict in product_dict.items():
        normalized_stoich = produ_dict['Stoch']/Edu_Bezug_stoich
        pool_change_dict[product] = normalized_stoich*actual_reaction_rate
       # print(product) 
       # print(Edu_Bezug_stoich)
       # print(pool_change_dict[product] )
        
    #W_max = microbe_dict['growth_rate']
    dead_microbes = Biomass*microbe_dict['death_rate']
    #biomass_change = - dead_microbes  + Biomass* W_max *actual_reaction_rate *C_for_growth* 12.011/1000
    
    biomass_change = - dead_microbes  +  C_for_growth * m_C #  m_C Gewicht 1 Mol C =  12.011 g /mol 
    pool_change_dict['biomass'] = biomass_change
    # print(biomass_change)
    # input('biomass '+pathway_name)
    #print(DGr_Ausgabe)
    pool_change_dict['DGr'] = DGr_Ausgabe*1e-3 # umrechnung von J/mol nach kJ/mol für den plot
    #print(pool_change_dict['DGr'])
    return pool_change_dict





###############################################################################
###############################################################################

def Ferm_help_Pathway(pool_dict,model_parameter_dict):
    # print('Ferm')
    # input('..')

    microbe_dict = {'concentration' : pool_dict['M_Ferm'], 
                    'Vmax'          : model_parameter_dict['Vmax_help_Ferm'],          #Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
                    #'growth_rate'   : 0,
                    'death_rate'    : 0,
                    'Kmb'           : model_parameter_dict['Kmb_help_Ferm'],
                    'Ferm_help'     : True, # ist nur wichtig das es den key 'Ferm_help' gibt 
                    'microbe'       : 'Ferm_help',
                    'CUE'           :       1 } # völlig unwichtig
            
    
    educt_dict =  {'C'              : {'concentration'  : pool_dict['C'],
                                       'Stoch'          : 1,
                                       'Km'             : 0}}# 10 from Song  mikromol pro gram 
                     
    
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
                    #'growth_rate'   : model_parameter_dict['w_Ferm'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'KmA_Ferm'      : model_parameter_dict['KmA_Ferm'],
                    'Ferm'          : True,
                    'microbe'       : "Ferm" ,
                    'CUE'           :   0.5,
                    'C_source'       : 'DOC'}
    
    
    educt_dict =  {'DOC'           : {'concentration':pool_dict['DOC'],
                                    'Stoch'          : 6                 ,                                    
                                    'Km'             : 10 / SOIL_DENSITY,# 10 from Song  mikromol pro gram 
                                    'C_atoms'      : 6                 }}    # weil glucose 6 C atome hat und ein momomer aus der spaltung von Coellulose ist 
                     
    
    
        
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_H2 = henrys_law( Henrys_dict['H2']['H_cp_Standard'],  Henrys_dict['H2']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_H2 =  pool_dict['H2']* (H_cc_H2/  (H_cc_H2+1))
                                       
    
    
    
    product_dict = { 'Acetate' : {'concentration': pool_dict['Acetate'],
                                  'Stoch'        : 6               }  ,
                    
                     'CO2'      : {'concentration': dissolved_CO2,
                                   'Stoch'        : 3              }  , 
                         
                      'H2'     : {'concentration': dissolved_H2   ,
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
                    #'growth_rate'   : model_parameter_dict['w_Fe3'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Fe3',
                    'CUE'           :  0.5  ,
                    'C_source'      :  'Acetate'}
                    #'DGs'           : -109.45 }
    
    educt_dict =  {'Acetate'       : {'concentration':pool_dict['Acetate'],
                                    'Stoch'          : 1                  ,
                                    'DGf'            : -369.31*1e3        ,
                                    'Km'             : 0.01 / SOIL_DENSITY, #0.01 / SOIL_DENSITY # wert nach Roden 2003 10 mal kleiner als bei Ace (0.8). Aber passt vlt nicht mehr mit den anderen werten zusammen
                                    'C_atoms'      : 2                    }, 
    
                      'Fe3'        :{'concentration':pool_dict['Fe3'] ,
                                   'Stoch'          : 8                   ,
                                   'DGf'            : -4.7*1e3            ,
                                   'Km'             : 0                   }} # TODO WISO 0 ? 
                     
    
                          
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_H2 = henrys_law( Henrys_dict['H2']['H_cp_Standard'],  Henrys_dict['H2']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_H2 =  pool_dict['H2']* (H_cc_H2/  (H_cc_H2+1))
                                       
    
    
    product_dict = { 'Fe2' : {'concentration': pool_dict['Fe2'],
                              'Stoch'        : 8               ,
                              'DGf'          : -89.1*1e3           }  ,#https://www.engineeringtoolbox.com/standard-state-enthalpy-formation-definition-value-Gibbs-free-energy-entropy-molar-heat-capacity-d_1978.html
                    
                     'CO2' : {'concentration': dissolved_CO2,
                              'Stoch'        : 2               ,
                              'DGf'          : -394.36*1e3                 }  , #Vaxa
                         
                      'H2' : {'concentration': dissolved_H2,
                              'Stoch'        : 0              ,
                              'DGf'          :  0              }   }
                                                                  
   
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Fe3')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Fe3'] = pool_change_dict.pop('biomass')
    # if 'DGr' in pool_change_dict:
    #     pool_change_dict['DGr_Fe3'] = pool_change_dict.pop('DGr')
    
    #print(pool_change_dict)
    
    return pool_change_dict





def Hydro_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O conrad2000selective, Fe3nchel -131kj/mol
    # print('Hydro')
    microbe_dict = {'concentration' : pool_dict['M_Hydro']              , 
                    'Vmax'          : model_parameter_dict['Vmax_Hydro'], ## 0.15 mikromol pro cm^3 from Song
                    #'growth_rate'   : model_parameter_dict['w_Hydro']   , 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Hydro',
                    'CUE'           :       0.5,
                    'C_source'      :  'CO2'}
                      #'DGs '       : 343.56 } # ja positiv
                      
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_H2 = henrys_law( Henrys_dict['H2']['H_cp_Standard'],  Henrys_dict['H2']['H_cp_temp'])
    H_cc_CH4 = henrys_law(Henrys_dict['CH4']['H_cp_Standard'], Henrys_dict['CH4']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_H2 =  pool_dict['H2']* (H_cc_H2/  (H_cc_H2+1))
    dissolved_CH4 = pool_dict['CH4']*(H_cc_CH4/ (H_cc_CH4+1))
                                       
    
    educt_dict = { 'H2'  : {'concentration':dissolved_H2 ,
                            'Stoch'     : 4                  , 
                            'DGf'       : 0                  ,
                            'Km'        : 0.01 / SOIL_DENSITY},# 0.01 mikromol pro cm^3 from Song
    
                  'CO2'  :{'concentration': dissolved_CO2,
                           'Stoch'        : 1                                ,
                           'DGf'          : -394.36*1e3                          ,
                           'Km'           : 0.05/SOIL_DENSITY  , # 0.05 mikromol pro cm^3 from Song
                           'C_atoms'      : 1                 }}
    
    product_dict = {'CH4' :{'concentration': dissolved_CH4 , # auf 0 weil es in den Headspace diffundiert
                            'Stoch'        : 1                ,
                            'DGf'          : -50.8*1e3}          }
   
    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Hydro')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Hydro'] = pool_change_dict.pop('biomass')
    # if 'DGr' in pool_change_dict:
    #     pool_change_dict['DGr_Hydro'] = pool_change_dict.pop('DGr')

    # print(pool_change_dict['CH4'])
    
    return pool_change_dict



def Homo_Pathway(pool_dict,model_parameter_dict):
    # PATHWAY: song,Fe3nchel,conrad2000selective : 4H2 + 2CO2 → CH3COOH+ 2H2O. Alternativ: 2CO2 + 8H2  = CH3COOH +H2O., laut Thauer sind es 4 H2 + Co2 -> CH4 +2H2O
    # print('Homo')
    microbe_dict = {'concentration' : pool_dict['M_Homo'], 
                    'Vmax'          : model_parameter_dict['Vmax_Homo'] , # # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
                    #'growth_rate'   : model_parameter_dict['w_Homo']    , 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Homo',
                    'CUE'           :       0.5 ,
                    'C_source'      :  'CO2'}
                    #'DGs'          : 25 }
                    
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_H2 = henrys_law( Henrys_dict['H2']['H_cp_Standard'],  Henrys_dict['H2']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_H2 =  pool_dict['H2']* (H_cc_H2/  (H_cc_H2+1))              
    
    educt_dict = { 'H2'  : {'concentration': dissolved_H2       ,
                             'Stoch'       : 4                     ,
                             'DGf'         : 0                     ,
                               'Km'        : 0.01 / SOIL_DENSITY   },    # 0.01 from Song
    
                  'CO2'  :{'concentration' : dissolved_CO2         ,
                           'Stoch'         : 2                     ,
                           'DGf'           : -394.36*1e3               ,
                           'Km'            : 0.05 / SOIL_DENSITY   ,  # 0.05 from Song, laut (van1999efFe3cts) größer als Hydro, laut schink1997energetics sollte der Wert mit sinkender Temp, mit zunehmendem Acetate und sinkendem PH sinken (Im vlg zu Hydro)
                           'C_atoms'       : 1                 }}
    
    product_dict = {'Acetate' : {'concentration': pool_dict['Acetate'] ,
                                 'Stoch'        : 1                    ,
                                 'DGf'          : -369.31*1e3}             , }
   
    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Homo')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Homo'] = pool_change_dict.pop('biomass')
    # if 'DGr' in pool_change_dict:
    #     pool_change_dict['DGr_Homo'] = pool_change_dict.pop('DGr')

    return pool_change_dict



def Ac_Pathway(pool_dict,model_parameter_dict):
    #PATHWAY:  CH3COO + H+ ->  CH4 + CO2 Fe3y_Conrad2000
    # print('Ac')
    microbe_dict = {'concentration' : pool_dict['M_Ac'], 
                    'Vmax'          : model_parameter_dict['Vmax_Ac'],  #Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
                    #'growth_rate'   : model_parameter_dict['w_Ac'], 
                    'death_rate'    : model_parameter_dict['Sensenmann'],
                    'microbe'       : 'M_Ac',
                    'CUE'           :       0.5,
                    'C_source'      :  'Acetate'}
                      #'DGs'        : -75.89 }
    
    educt_dict = { 'Acetate' : {'concentration': pool_dict['Acetate']    ,
                                'Stoch'        :  1                     ,
                                 'DGf'         : -369.31*1e3                ,
                                 'Km'          :  0.05 / SOIL_DENSITY   , #0.05 / SOIL_DENSITY   # 0.05 from song, Wert sollte 10 mal höher sein als bei AltE laut roden2003 bei 12Mikromol !!!!
                                 'C_atoms'      : 2                 }}
    
    H_cc_CO2 = henrys_law(Henrys_dict['CO2']['H_cp_Standard'], Henrys_dict['CO2']['H_cp_temp'])
    H_cc_CH4 = henrys_law(Henrys_dict['CH4']['H_cp_Standard'], Henrys_dict['CH4']['H_cp_temp'])
    
    dissolved_CO2 = pool_dict['CO2']*(H_cc_CO2/ (H_cc_CO2 +1))
    dissolved_CH4 = pool_dict['CH4']*(H_cc_CH4/ (H_cc_CH4+1))
    
                      
    product_dict = { 'CH4' : {'concentration': dissolved_CH4, #auf 0, weil es in den Headspace diffundiert
                              'Stoch'        : 1               ,
                              'DGf'          : -50.8*1e3          }, # wert aus Vaxa
                    
                     'CO2' : {'concentration': dissolved_CO2,
                              'Stoch'        : 1               ,
                              'DGf'          :-394.36*1e3         } }

    
    pool_change_dict = GeneralPathway(microbe_dict, educt_dict, product_dict, 'Ac')
    
    if 'biomass' in pool_change_dict:
        pool_change_dict['M_Ac'] = pool_change_dict.pop('biomass')
    # if 'DGr' in pool_change_dict:
    #     pool_change_dict['DGr_Ac'] = pool_change_dict.pop('DGr')

    
    return pool_change_dict


