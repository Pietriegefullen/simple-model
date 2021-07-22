# Hier ruFe3n wir die Mikroben auf und summieren die Änderungen der Pools 


import numpy as np
from order import pool_order

from Thermomicrobentest import Ferm_help_Pathway,Ferm_Pathway, Hydro_Pathway, Fe3_Pathway, Ac_Pathway, Homo_Pathway # Importiert Funktionen aus dem Skript "Microbe"

def Cdec_wrapper(model_parameter_dict, return_thermodynamics = False):
  
        
    def Cdec(t, system_state):
        pool_dict = dict(zip(pool_order, system_state))

     
        # print('C_pool' , C_pool)
        
        
        # Ferm Ferm Ferm Ferm Ferm  Ratio aus Grant 1998
       
        #deltaM_Ferm2, deltaC_pool2, _, Tot_Ferm2 =   Fermenters(M_Ferm2_pool, C_pool, Acetate_pool, Vmax_Ferm, w_Ferm, Sensenmann, Kmb_Ferm, Kmh_Ferm)
        # Produziert:
        #deltaH2_Ferm2 =  -deltaC_pool2 * 1.2 # 1.2 aus cabrol2017microbial
         
        #FERM FERM FERM 
        pool_change_dict_Ferm_help =  Ferm_help_Pathway(pool_dict, model_parameter_dict)
        
        pool_change_dict_Ferm =  Ferm_Pathway(pool_dict, model_parameter_dict)
        
        # Fe3 Fe3 Fe3 Fe3 : C2H3O2 − + 4H2O + 8Fe3+3 → 9H+ + 2 HCO3− + 8Fe3+2   Delattre 2019
        pool_change_dict_Fe3 =  Fe3_Pathway(pool_dict, model_parameter_dict)
        #print(pool_change_dict_Fe3)
     
        # ACETO ACETO ACETO ACETO  CH3COO + H+ ->  CH4 + CO2 Fe3y_Conrad2000
        pool_change_dict_Ac =  Ac_Pathway(pool_dict, model_parameter_dict)
    
        #HYDRO HYDRO HYDRO HYDRO 4H2 + CO2 → CH4 + 2H2O, conrad2000selective, Fe3nchel -131kj/mol
        pool_change_dict_Hydro =  Hydro_Pathway(pool_dict, model_parameter_dict)
      
        # HOMO HOMO HOMO HOMO HOMO: 4H2 + 2CO2 → CH3COOH + 2H2O  conrad2000selective
        pool_change_dict_Homo = Homo_Pathway(pool_dict, model_parameter_dict)        
        
        pool_change_dict_list = [pool_change_dict_Ferm_help, 
                                 pool_change_dict_Ferm, 
                                 pool_change_dict_Fe3, 
                                 pool_change_dict_Ac,
                                 pool_change_dict_Hydro,
                                 pool_change_dict_Homo]
        
        
        #print(pool_change_dict_list)
        # DELTA DELTA DELTA                           
        
       # m_C = 12.01*1e-3 # molar mass of carbon
       # das übernimmt jetzt die spezielle mikrobe? = changes_dict['C'] = -min(-deltaC_pool, C_pool)  +  (Tot_Hyd + Tot_Ac + Tot_ALtE + Tot_Ferm)/m_C
        
        # hier werden die Beitr�ge aller Mikroben zur Poolgrößenänderung zusammengefasst
        changes_dict = dict()
        # for pool in pool_dict:
        #     changes_dict[pool] = sum(list([dictionary[pool] for dictionary in pool_change_dict_list if pool in dictionary]))
        for pathway_change_dict in pool_change_dict_list:
            for pool, change in pathway_change_dict.items():
                if not pool in changes_dict:
                    changes_dict[pool] = 0.0
                changes_dict[pool] += change

        if 'Fe2' in changes_dict and changes_dict['Fe2'] < 0:
           print('!!!', changes_dict['Fe2'])        
            #input('!!!')
        #changes_dict['H2_Ferm2'] = 0
        # print('changes CO2_Ferm', changes_dict['CO2_Ferm'])
        # input('')
      
        # for k,v in changes_dict.items():
        #     print(k,v)
        # print('DGr_Fe3 is',changes_dict['DGr_Fe3'])
        # print('CO2_Homo is',type(changes_dict['CO2_Homo']))
       # changes_dict['CO2'] = changes_dict['CO2'] /2
        #changes_dict['CH4'] =0
        
        pool_changes_array = np.array([changes_dict[pool_name] if pool_name in changes_dict else 0.0 for pool_name in pool_order])
        # input('dcdec done..')
        #print(pool_changes_array[:][0])       
        #print(dict(zip(pool_order, pool_changes_array)))
        #input('pth..')

        # print('')
        # for p,v in zip(pool_order, pool_changes_array):
        #     print(p,v)
        #print(pool_changes_array)
        #input('input reached')
        if not return_thermodynamics:
            return pool_changes_array
        
        else:
            
            # speziellere aufteilung der Pools for plotting
            extra_dict =  dict()
            extra_dict['CO2_Ferm'] = pool_change_dict_Ferm['CO2']  if 'CO2' in pool_change_dict_Ferm else 0.
            extra_dict['CH4_Hydro'] = pool_change_dict_Hydro['CH4']  if 'CH4' in pool_change_dict_Hydro else 0.
            extra_dict['CO2_Ac'] = pool_change_dict_Ac['CO2'] if 'CO2' in pool_change_dict_Ac else 0.
            extra_dict['CO2_Hydro'] = pool_change_dict_Hydro['CO2'] if 'CO2' in pool_change_dict_Hydro else 0.
            extra_dict['CO2_Fe3'] = pool_change_dict_Fe3['CO2'] if 'CO2' in pool_change_dict_Fe3 else 0.
            extra_dict['CO2_Homo'] =pool_change_dict_Homo['CO2'] if 'CO2' in pool_change_dict_Homo else 0.
            extra_dict['H2_Homo'] = pool_change_dict_Homo['H2'] if 'H2' in pool_change_dict_Homo else 0.
            extra_dict['H2_Hydro'] = pool_change_dict_Hydro['H2'] if 'H2' in pool_change_dict_Hydro else 0.
            #extra_dict['Fe2_Fe3'] = pool_change_dict_Fe3['Fe2'] if 'Fe2' in pool_change_dict_Fe3 else 0.
            
            #
            extra_dict['DGr_Fe3'] = pool_change_dict_Fe3['DGr'] if 'DGr' in pool_change_dict_Fe3 else 0
            extra_dict['DGr_Hydro'] = pool_change_dict_Hydro['DGr'] if 'DGr' in pool_change_dict_Hydro else 0.
            extra_dict['DGr_Homo'] = pool_change_dict_Homo['DGr'] if 'DGr' in pool_change_dict_Homo else 0.
            extra_dict['DGr_Ac'] = pool_change_dict_Ac['DGr'] if 'DGr' in pool_change_dict_Ac else 0.
                    
            return pool_changes_array, extra_dict
        
    return Cdec


  



