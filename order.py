# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 10:16:35 2021

@author: Lara
"""


pool_order = ['C',
            'DOC',
            'CH4',
            'CO2',
            'Fe3',
            'M_Ac',
            'M_Ferm',
            'M_Fe3',
            'M_Hydro',
            'M_Homo',
            'Acetate',
            'H2',
            'Fe2',
            'pH',
            'H2O',
            'HCO3',
            'weight',
            'water'
            # 'CO2_Ac',
            # 'CO2_Hydro',
            # 'CH4_Hydro',
            # 'CO2_Ferm',
            # 'CO2_Fe3',
            # 'CO2_Homo',
            # 'H2_Homo',
            # 'H2_Hydro',
            # 'DGr_Fe3',
            # 'DGr_Hydro',
            # 'DGr_Homo',
            # 'DGr_Ac'
            ]



# durch auskommentieren bestimmen wir was optimiert wird und was nicht
changeables_order = [
                     'Vmax_help_Ferm',
                     'Vmax_Ferm',
                     'Vmax_Fe3',
                     'Vmax_Homo',
                     'Vmax_Hydro',
                     'Vmax_Ac',
                    #'w_Ferm',
                    # 'w_Fe3',
                     #'w_Hydro', 
                     #'w_Homo',
                     #'w_Ac',
                     'Sensenmann',
                     'Kmb_help_Ferm',
                     'KmA_Ferm',
                     'Fe3',
                     'M_Ac'
                     ]

parameter_units = {'Vmax_Ferm':'μmol/mg',  
                   'Vmax_help_Ferm':'μmol/mg', 
                   'Vmax_Fe3':'μmol/mg',       
                   'Vmax_Homo':'μmol/mg',      
                   'Vmax_Hydro':'μmol/mg',          
                   'Vmax_Ac':'μmol/mg',      
                  # 'w_Ferm':'mg/μmol',
                 #  'w_Fe3':'mg/μmol',
                  # 'w_Hydro':'mg/μmol',
                   #'w_Homo':'mg/μmol',
                   #'w_Ac':'mg/μmol',
                   'Sensenmann':'-', 
                   #'Stoch_Fe3':'-',
                   'Kmb_help_Ferm':'mg' ,
                   #'Kmh_Ferm':'μmol', 
                   'KmA_Ferm':'μmol',
                   'Fe3':'μmol',
                   'C':'μmol',
                   'DOC':'μmol',
                   'M_Ac':'mg',
                   'M_Ferm':'mg',
                   'M_Fe3':'mg',
                   'M_Hydro':'mg',
                   'M_Homo':'mg',
                   'CH4':'μmol',
                   'CO2':'μmol',
                   'CO2_Ac':'μmol',
                   'Acetate':'μmol',
                   'H2':'μmol',
                   'CO2_Hydro':'μmol',
                   'CH4_Hydro':'μmol'}
                  # 'H2_Ferm2':'μmol',
                   #'M_Ferm2':'mg'}



Henrys_dict = { 'CO2' : {'H_cp_Standard': 3.4*10e-4,
                         'H_cp_temp':    2400                  },
                   
               'CH4' : {'H_cp_Standard': 1.4*10e-5,
                            'H_cp_temp':    1700                },
               
               'H2' : {'H_cp_Standard': 7.7*10e-6,
                         'H_cp_temp':    500              },}



enthalpy = {'Acetate' :-484.13 *1e3,
            'Fe3'     : -48.5344  *1e3,
            'Fe2'     : -89.1192  *1e3,
            'CO2'     : -413.8 *1e3, # Wert für aq!
            'H2'      :  0.0,
            'HCO3'    : 1,
            'CH4'     : -89 *1e3}
            #'H2O'     : -285.82996 * 1e3}

Gibbs_formation = {'Acetate' :  -396.46*1e3,    # Tabellenwerte in kJ /mol, hier J/mol
                   'Fe3'     :  -4.6*1e3 , 
                   'Fe2'     :  -78.*1e3,
                   'CO2'     :  -386.36*1e3, # Wert für aq!
                   'H2'      :   0.0,
                   'CH4'     :  -34.4*1e3  ,# Wert für aq!
                   'HCO3'    :  -586.85 * 1e3, # # Wert für aq!
                   #'H2O'    :   -237.178408 *1e3,
                   'H_plus' :          0}  


def get_fixed_quantities(): 
    fixed_quantities_dict = dict()        
    #fixed_quantities_dict['M_Ac'] = 0.2 # superdata_Kuru = 0.001
    fixed_quantities_dict['M_Ferm'] =  0.5#1.3e-07 #
    fixed_quantities_dict['M_Fe3'] =0.2  # 0.2
    fixed_quantities_dict['M_Hydro'] = 0.2 # 0.2
    fixed_quantities_dict['M_Homo'] = 0.2 # 0.2
    
    return(fixed_quantities_dict)
    
def get_initial_guesses():
    # specify initial guesses and bounds for the parameters to be optimized
    initial_guess_dict = dict()         #   init    lower upper  ok guesses to start with
    initial_guess_dict['Vmax_help_Ferm'] =  (0.002,   0.009,1.71)  # 0.05
    initial_guess_dict['Vmax_Ferm'] =       (0.9,   0.01, 0.2)  # 0.011       # Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
    initial_guess_dict['Vmax_Fe3'] =        (0.9,   0.029, 3)  # 0.8         # Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    initial_guess_dict['Vmax_Homo'] =       (0.0869, 0.005, 1.)   # 0.869       # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    initial_guess_dict['Vmax_Hydro'] =      (0.0182, 0.03, 0.3)   # 0.182 1.8   # 0.15 mikromol pro cm^3 from Song
    initial_guess_dict['Vmax_Ac'] =         (0.2,   0.05, 3.0)  # 0.99           # Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
   #initial_guess_dict['w_Ferm'] =          (0.00032,  0.03, 0.05)  # 0.05
   #initial_guess_dict['w_Fe3'] =           (0.0093, 0.01, 0.05)  # 0.013
   #initial_guess_dict['w_Hydro'] =         (0.00024, 0.01, 0.05)  # 0.024
   #initial_guess_dict['w_Homo'] =          (0.00049, 0.01, 0.05)  # 0.049
   #initial_guess_dict['w_Ac'] =            (0.0036,  0.01, 0.05)  # 0.04
    initial_guess_dict['Sensenmann'] =      (8.33e-5, 0, 8.44e-5)# 0
    initial_guess_dict['Kmb_help_Ferm'] =   (1.009,    0.5,  10)      # 10
    # initial_guess_dict['Kmh_Ferm'] =        (10,    1,  10)    # 10
    initial_guess_dict['Fe3'] =             (0.15,  0,  100)    # 15.587,
    initial_guess_dict['M_Ac'] =            (0.02,  1.3e-08,  5e05) # 0.002
    initial_guess_dict['KmA_Ferm']=         (0.9315, 0.001, 20)     # 17.315 # Diese Boundaries müssen anhander Acetatekurven angepasst werden

    return(initial_guess_dict)



