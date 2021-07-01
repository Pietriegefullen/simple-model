# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 10:16:35 2021

@author: Lara
"""


pool_order = ['C',
            'DOC',
            'Fe3',
            'M_Ac',
            'M_Ferm',
            'M_Fe3',
            'M_Hydro',
            'M_Homo',
            'CH4',
            'CO2',
            'CO2_Ac',
            'Acetate',
            'H2',
            'CO2_Hydro',
            'CH4_Hydro',
            'H2_Ferm2',
            'M_Ferm2',
            'CO2_Ferm',
            'CO2_Fe3',
            'CO2_Homo',
            'H2_Homo',
            'H2_Hydro',
            'Fe2']



# durch auskommentieren bestimmen wir was optimiert wird und was nicht
changeables_order = [
                     'Vmax_Ferm',
                     'Vmax_help_Ferm',
                     'Vmax_Fe3',
                     'Vmax_Homo',
                     'Vmax_Hydro',
                     'Vmax_Ac',
                     'w_Ferm',
                     'w_Fe3',
                     'w_Hydro', 
                     'w_Homo',
                     'w_Ac',
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
                   'w_Ferm':'mg/μmol',
                   'w_Fe3':'mg/μmol',
                   'w_Hydro':'mg/μmol',
                   'w_Homo':'mg/μmol',
                   'w_Ac':'mg/μmol',
                   'Sensenmann':'-', 
                   'Stoch_Fe3':'-',
                   'Kmb_Ferm':'mg' ,
                   'Kmh_Ferm':'μmol', 
                   'Kmb_Fe3':'mg',
                   'Kmb_Ac':'mg',
                   'Kmb_Hydro':'mg',
                   'Kmb_Homo':'mg',
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
                   'CH4_Hydro':'μmol',
                   'H2_Ferm2':'μmol',
                   'M_Ferm2':'mg'}