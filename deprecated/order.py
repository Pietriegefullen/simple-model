# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 10:16:35 2021

@author: Lara
"""

#global_switches =


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
            #'H2O',
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
                     #'Sensenmann',
                     'Kmb_help_Ferm',
                     'Inhibition_Ferm',
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
                   'Inhibition_Ferm':'μmol',
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
    fixed_quantities_dict['M_Ferm'] =  0.7#1.3e-07 #
    fixed_quantities_dict['M_Fe3'] =0.2  # 0.2
    fixed_quantities_dict['M_Hydro'] = 0.0025 # 0.2
    fixed_quantities_dict['M_Homo'] = 0.0001 # 0.2

    return(fixed_quantities_dict)

def get_initial_guesses():
    # specify initial guesses and bounds for the parameters to be optimized
    initial_guess_dict = dict()         #   init    lower upper  ok guesses to start with
    initial_guess_dict['Vmax_help_Ferm'] =  (0.177,   0.006,0.5)  # 0.05
    initial_guess_dict['Vmax_Ferm'] =       (0.606,   0.01, 1.)  # 0.011       # Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
    initial_guess_dict['Vmax_Fe3'] =        (1.709,   0.02, 3)  # 0.8         # Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    initial_guess_dict['Vmax_Homo'] =       (0.85, 0.05, 1.)   # 0.869       # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    initial_guess_dict['Vmax_Hydro'] =      (0.423, 0.03, 1.)   # 0.182 1.8   # 0.15 mikromol pro cm^3 from Song
    initial_guess_dict['Vmax_Ac'] =         (0.159,   0.05, 1.0)  # 0.99           # Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
   #initial_guess_dict['w_Ferm'] =          (0.00032,  0.03, 0.05)  # 0.05
   #initial_guess_dict['w_Fe3'] =           (0.0093, 0.01, 0.05)  # 0.013
   #initial_guess_dict['w_Hydro'] =         (0.00024, 0.01, 0.05)  # 0.024
   #initial_guess_dict['w_Homo'] =          (0.00049, 0.01, 0.05)  # 0.049
   #initial_guess_dict['w_Ac'] =            (0.0036,  0.01, 0.05)  # 0.04
    initial_guess_dict['Sensenmann'] =      (8.33e-5, 0, 1.5)# (8.33e-5, 0, 8.44e-5)# 0
    initial_guess_dict['Kmb_help_Ferm'] =   (9.829,    0.005,  20)      # 10
    # initial_guess_dict['Kmh_Ferm'] =        (10,    1,  10)    # 10
    initial_guess_dict['Fe3'] =             (3.645,  0,  100)    # 15.587,
    initial_guess_dict['M_Ac'] =            (0.001,  1.3e-08,  0.5) # 0.002
    initial_guess_dict['Inhibition_Ferm']=  (0.412, 0.001, 20) # Je niedriger desto hemmung # Diese Boundaries müssen anhander Acetatekurven angepasst werden

    return(initial_guess_dict)


"""

25750: SSR     10.502   best     10.501
Vmax_help_Ferm      0.177 μmol/mg
Vmax_Ferm           0.606 μmol/mg
Vmax_Fe3            1.709 μmol/mg
Vmax_Homo           0.850 μmol/mg
Vmax_Hydro          0.423 μmol/mg
Vmax_Ac             0.159 μmol/mg
Kmb_help_Ferm       9.829 mg
Inhibition_Ferm     0.412 μmol
Fe3                 3.645 μmol
M_Ac                0.001 mg
calling extra info
r2 for CH4 is 0.9145413005069074
r2adj  CH4 is 0.8992808184545694
r2 for CO2 is 0.9038053555573695
r2adj  CO2 is 0.8866277404783284
D:\OneDrive - Universität Hamburg\simple model\Thermorefactored.py:599: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).
  plt.figure()




SSR (mean squared) 41.706

   Fe3                   19.3639 *
   M_Ac                   1.3211 *
model parameters:
   Vmax_help_Ferm         0.0433 *
   Vmax_Ferm              0.0331 *
   Vmax_Fe3               1.1593 *
   Vmax_Homo              0.4160 *
   Vmax_Hydro             0.0597 *
   Vmax_Ac                0.6788 *
   Kmb_help_Ferm          3.3514 *
   Inhibition_Ferm       19.4313 *


SSR (mean squared) 39.48
   initial pool values:
   Fe3                   10.1351 *
   M_Ac                   0.5125 *
model parameters:
   Vmax_help_Ferm         0.0100 *
   Vmax_Ferm              0.0348 *
   Vmax_Fe3               2.2136 *
   Vmax_Homo              0.5960 *
   Vmax_Hydro             0.0438 *
   Vmax_Ac                1.4306 *
   Kmb_help_Ferm          1.0010 *
   Inhibition_Ferm       11.0730 *


SSR(mean squared) 32.26
    initial pool values:
   Fe3                   12.0230 *
   M_Ac                   0.7038 *
model parameters:
   Vmax_help_Ferm         0.0719 *
   Vmax_Ferm              0.0377 *
   Vmax_Fe3               2.6594 *
   Vmax_Homo              0.9000 *
   Vmax_Hydro             0.1557 *
   Vmax_Ac                2.9555 *
   Kmb_help_Ferm          2.3993 *
   Inhibition_Ferm       17.5532 *

"""