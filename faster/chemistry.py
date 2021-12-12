
import numpy as np

import CONSTANTS


HENRYS_DICT = { 'CO2' : {'H_cp_Standard': 3.4*10e-4,
                         'H_cp_temp':    2400                  },

               'CH4' : {'H_cp_Standard': 1.4*10e-5,
                            'H_cp_temp':    1700                },

               'H2' : {'H_cp_Standard': 7.7*10e-6,
                         'H_cp_temp':    500              }
               }


def henrys_law(substance):
    if not substance in HENRYS_DICT:
        return 1.0

    T = CONSTANTS.SPECIMEN_TEMPERATURE + 273.15

    H_cp_temp = HENRYS_DICT[substance]['H_cp_temp']
    H_cp_Standard = HENRYS_DICT[substance]['H_cp_Standard']

    H_cp_temp_adjusted = H_cp_temp - T # Standard Temp Wert - actual Temp
    H_cc = H_cp_Standard *  np.exp(H_cp_temp_adjusted * (1/T - 1/CONSTANTS.STANDARD_TEMPERATURE))

    return H_cc/(H_cc+1)
