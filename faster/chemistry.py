
import numpy as np

import CONSTANTS


HENRYS_DICT = {'CO2' : {'H_cp_Standard': 3.4*10e-4,
                         'H_cp_temp':    2400                  },

               'CH4' : {'H_cp_Standard': 1.4*10e-5,
                            'H_cp_temp':    1700                },

               'H2' : {'H_cp_Standard': 7.7*10e-6,
                         'H_cp_temp':    500              }
               }

GIBBS_FORMATION = {'Acetate' :  -396.46*1e3,    # Tabellenwerte üblicherweise in kJ /mol, hier J/mol
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

def henrys_law(substance):
    if not substance in HENRYS_DICT:
        return 1.0

    T = CONSTANTS.SPECIMEN_TEMPERATURE + 273.15

    H_cp_temp = HENRYS_DICT[substance]['H_cp_temp']
    H_cp_Standard = HENRYS_DICT[substance]['H_cp_Standard']

    H_cp_temp_adjusted = H_cp_temp - T # Standard Temp Wert - actual Temp
    H_cc = H_cp_Standard *  np.exp(H_cp_temp_adjusted * (1/T - 1/CONSTANTS.STANDARD_TEMPERATURE))

    return H_cc/(H_cc+1)
