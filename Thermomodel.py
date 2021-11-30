
import numpy as np

from order import pool_order
from Thermomicrobentest import Ferm_help_Pathway,Ferm_Pathway, Hydro_Pathway, Fe3_Pathway, Ac_Pathway, Homo_Pathway # Importiert Funktionen aus dem Skript "Microbe"

################# Aufruf der Mikroben und Aufsummieren der Pooländerungen #####


def Cdec_wrapper(model_parameter_dict, return_thermodynamics = False):
    # wrapper lädt das model_parameter_dict das für jeden durchlauf gleich ist
        
    def Cdec(t, system_state):
        pool_dict = dict(zip(pool_order, system_state))
        
        #safeguard to prevent negative and zero pools
        pool_dict = {k: v if v > 0 else 1e-30  for k,v in pool_dict.items()}
        
#---------------- Aufrur der einzelnen Mikroben -------------------------------        
     
        #FERM FERM FERM 
        pool_change_dict_Ferm_help =  Ferm_help_Pathway(pool_dict, model_parameter_dict)
        
        pool_change_dict_Ferm =  Ferm_Pathway(pool_dict, model_parameter_dict)
        
        # Fe3 Fe3 Fe3 Fe3,                                                      # Pathway:    C2H3O2 − + 4H2O + 8Fe3(III)   --->  9H+ + 2 HCO3− + 8Fe3+2   Delattre 2019 , jason 2001
        pool_change_dict_Fe3 =  Fe3_Pathway(pool_dict, model_parameter_dict)

        
        # ACETO ACETO ACETO ACETO,                                           Pathway:  CH3COO + H+ --->  CH4 + CO2 Fe3y_Conrad2000
        pool_change_dict_Ac =  Ac_Pathway(pool_dict, model_parameter_dict)

        #HYDRO HYDRO HYDRO HYDRO,                                              Pathway: 4H2 + CO2 ---> CH4 + 2H2O, conrad2000selective, Fe3nchel -131kj/mol
        pool_change_dict_Hydro =  Hydro_Pathway(pool_dict, model_parameter_dict)

        # HOMO HOMO HOMO HOMO HOMO:,                                           Pathway: 4H2 + 2CO2 ---> CH3COOH + 2H2O  conrad2000selective
        pool_change_dict_Homo = Homo_Pathway(pool_dict, model_parameter_dict)        
        
        pool_change_dict_list = [pool_change_dict_Ferm_help, 
                                 pool_change_dict_Ferm, 
                                 pool_change_dict_Fe3, 
                                 pool_change_dict_Ac,
                                 pool_change_dict_Hydro,
                                 pool_change_dict_Homo
                                 ]
        
#------------------------------------------------------------------------------        
        
    
#---------- Beitrag aller Mikroben zu einer Substanz wird aufsummiert ---------

        changes_dict = dict()
       # Aufsummieren der Gesamtpooländerungen
        for pathway_change_dict in pool_change_dict_list:
            for pool, change in pathway_change_dict.items():
                if not pool in changes_dict:
                    changes_dict[pool] = 0.0 # wenn der Pool nicht im changes_dict ist wird 0 addiert
                changes_dict[pool] += change # ansonsten wird die Änderung addiert          
              
        # Check, dass nicht mehr aus dem Pool geholt wird als drin ist
        pool_changes_array = np.array([-min(-changes_dict[pool_name], pool_dict[pool_name]) # die Änderung ist das was jeweils kleiner ist, die Änderung, oder der Pool 
                                       if pool_name in changes_dict else 0.0   # wenn der pool_name nicht in changes_dict, wird nichts rausgenommen. Der solver will aber einen Wert, also 0 
                                       for pool_name in pool_order])           # für alle pool_names 
    
#------------------------------------------------------------------------------        
        
#--------------------------Extra Info fürs Plotting ---------------------------        
                  
        if return_thermodynamics is True:# True in 'Thermorefactored' in  'compute_extra_info'
            # Der Solver kommt nicht klar mit diesen Sonder-ausgaben
           
            extra_dict =  dict()
            extra_dict['CO2_Ferm'] = pool_change_dict_Ferm['CO2']  if 'CO2' in pool_change_dict_Ferm else 0.
            extra_dict['CH4_Hydro'] = pool_change_dict_Hydro['CH4']  if 'CH4' in pool_change_dict_Hydro else 0.
            extra_dict['CO2_Ac'] = pool_change_dict_Ac['CO2'] if 'CO2' in pool_change_dict_Ac else 0.
            extra_dict['CO2_Hydro'] = pool_change_dict_Hydro['CO2'] if 'CO2' in pool_change_dict_Hydro else 0.
            extra_dict['CO2_Fe3'] = pool_change_dict_Fe3['CO2'] if 'CO2' in pool_change_dict_Fe3 else 0.
            extra_dict['CO2_Homo'] =pool_change_dict_Homo['CO2'] if 'CO2' in pool_change_dict_Homo else 0.
            extra_dict['H2_Homo'] = pool_change_dict_Homo['H2'] if 'H2' in pool_change_dict_Homo else 0.
            extra_dict['H2_Hydro'] = pool_change_dict_Hydro['H2'] if 'H2' in pool_change_dict_Hydro else 0.
            
            extra_dict['Acetate_used'] = pool_change_dict_Ac['Acetate'] if 'Acetate' in pool_change_dict_Ac else 0.
            extra_dict['Acetate_used'] =+ pool_change_dict_Fe3['Acetate'] if 'Acetate' in pool_change_dict_Fe3 else 0.
            #extra_dict['CO2dissolved_CO2_total'] = changes_dict['dissolved_CO2_total'] if 'dissolved_CO2_total' in changes_dict else 0.
            #extra_dict['dissolved_H2'] = changes_dict['dissolved_H2'] if 'dissolved_H2' in changes_dict else 0.
            extra_dict['MM'] = pool_change_dict_Hydro['MM'] if 'MM' in pool_change_dict_Hydro else 0.
            #extra_dict['hydro_diss_co2'] = pool_change_dict_Hydro['hydro_diss_co2'] if 'hydro_diss_co2' in pool_change_dict_Hydro else 0.
            #extra_dict['hydro_diss_h2'] = pool_change_dict_Hydro['hydro_diss_h2'] if 'hydro_diss_h2' in pool_change_dict_Hydro else 0.

            
            return pool_changes_array, extra_dict # Ausgabe fürs plotting 
 #------------------------------------------------------------------------------    
   
        else:
            return pool_changes_array # Ausgabe der relevanten Pools für den solver
                          
    return Cdec


  



