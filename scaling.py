# -*- coding: utf-8 -*-
"""
Created on Sat May 29 19:13:35 2021

@author: Lara
"""

scaling_factors = {
           # 'C':            0.001,
           # 'M_Ferm':       10,
           # 'M_Fe3':        10,
           # 'M_Ferm2':      10,
           'Vmax_help_Ferm': 100,
            'Vmax_Ferm':    100,
            'Vmax_Fe':      10,
            'Vmax_Homo':    10,
            'Vmax_Hydro':   10,
            'Vmax_Ac':      10,
           # 'w_Ferm':       100,
           # 'w_Fe3':         100,
           # 'w_Hydro':      10,
           # 'w_Homo':       10,
           # 'w_Ac':         100,
            'Sensenmann':   1e5,
           # 'Stoch_Fe3':     10,
            'Kmb_help_Ferm':     0.1,
           # 'Kmb_Ferm':     0.1,
           # 'Kmh_Ferm':     0.1,
           # 'Kmb_Fe3':       0.1,
           # 'Kmb_Ac':       0.1,
            #'Kmb_Hydro':    0.1,
            'Fe3':           0.1,
            'M_Ac':         10,
            #'Kmb_Homo' :    0.1,
            'KmA_Ferm':     0.1,
        }
scaling_factors = {k:v*100000. for k,v in scaling_factors.items()}

#scaling_factors = {k:1.0 for k,v in scaling_factors.items()} # switches off scaling



def scale(parameter_dict):
    scaled_dict = dict()
    for name, value in parameter_dict.items():
        factor = 1.
        if name in scaling_factors.keys():
            factor = scaling_factors[name]
       
        if isinstance(value, tuple):
            scaled_value = tuple([v*factor for v in value])
        else:
            scaled_value = value*factor
        scaled_dict[name] = scaled_value
    return scaled_dict
    
def rescale(parameter_dict):
    scaled_dict = dict()
    for name, value in parameter_dict.items():
        factor = 1.
        if name in scaling_factors.keys():
            factor = scaling_factors[name]
               
        if isinstance(value, tuple):
            scaled_value = tuple([v/factor for v in value])
        else:
            scaled_value = value/factor

        scaled_dict[name] = scaled_value
    return scaled_dict