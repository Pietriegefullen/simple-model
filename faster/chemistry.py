
import numpy as np

import CONSTANTS


HENRYS_DICT = {'CO2' : {'H_cp_Standard': 3.4*1e-4,
                         'H_cp_temp':    2400                  },

               'CH4' : {'H_cp_Standard': 1.4*1e-5,
                            'H_cp_temp':    1700                },

               'H2' : {'H_cp_Standard': 7.7*1e-6,
                         'H_cp_temp':    500              },
               
               #'H2O': {'H_cp_Standard': 0, 
               #          'H_cp_temp':   0              }
               }

GIBBS_FORMATION = {'Acetate' :  -369.46*1e3,    # Tabellenwerte üblicherweise in kJ /mol, hier J/mol
                   'Fe3'     :  -4.6*1e3 , 
                   'Fe2'     :  -78.*1e3,
                   'CO2'     :  -386.36*1e3, # Wert für aq!
                   'H2'      :   0.0,
                   'CH4'     :  -34.4*1e3  ,# Wert für aq!
                   'HCO3'    :  -586.85 * 1e3, # # Wert für aq!
                   'H2O'    :   -237.178408 *1e3,# Wert für aq!
                   'H_plus' :          0}  


GIBBS_MINIMUM = -26.*1e3      # J/mol, Einheit passt zur Gaskonstante - 26 in kJ/mol aus z.b. blodau2011thermodynamic
                              # wert DGmin z.b aus Schink 1997, ist 1/3 der Energie die für ein ATP Herstellung benötigt wird

def C_atoms(substance):
    sources = {'CO2': 1,
               'Acetate': 2,
               'DOC': 10/6}
    return sources[substance]

def henrys_law(substance):
    if not substance in HENRYS_DICT:
        return 1.0

    R = CONSTANTS.GAS_CONSTANT
    T = CONSTANTS.SPECIMEN_TEMPERATURE
    T_std = CONSTANTS.STANDARD_TEMPERATURE

    C_cp = HENRYS_DICT[substance]['H_cp_temp']
    H_cp_Standard = HENRYS_DICT[substance]['H_cp_Standard']

    H_cc = H_cp_Standard * R * T *  np.exp(C_cp * (1/T - 1/T_std))
    return H_cc
