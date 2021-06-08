# Hier ruFe3n wir die Mikroben auf und summieren die Änderungen der Pools 


import numpy as np
from order import pool_order

from Thermomicrobentest import Ferm_Pathway, Hydro_Pathway, Fe3_Pathway, Ac_Pathway, Homo_Pathway # Importiert Funktionen aus dem Skript "Microbe"

def Cdec_wrapper(model_parameter_dict):
    Vmax_Ferm = model_parameter_dict['Vmax_Ferm']
    Vmax_Fe3 = model_parameter_dict['Vmax_Fe3']
    Vmax_Ac = model_parameter_dict['Vmax_Ac']
    Vmax_Hydro = model_parameter_dict['Vmax_Hydro']
    Vmax_Homo = model_parameter_dict['Vmax_Homo']
    w_Ferm = model_parameter_dict['w_Ferm']
    w_Fe3 = model_parameter_dict['w_Fe3']
    w_Ac = model_parameter_dict['w_Ac']
    w_Hydro = model_parameter_dict['w_Hydro']
    w_Homo = model_parameter_dict['w_Homo']
    Kmb_Ferm = model_parameter_dict['Kmb_Ferm']
    Kmh_Ferm = model_parameter_dict['Kmh_Ferm']
    Kmb_Fe3 = model_parameter_dict['Kmb_Fe3']
    Kmb_Ac = model_parameter_dict['Kmb_Ac']
    Kmb_Hydro = model_parameter_dict['Kmb_Hydro']
    Stoch_Fe3 = model_parameter_dict['Stoch_Fe3'] # ist 8 und sollte rausgenommen werden. 
    Sensenmann = model_parameter_dict['Sensenmann']
    KmA_Herm = model_parameter_dict['KmA_Ferm']
    Kmb_Homo = model_parameter_dict['Kmb_Homo']
        
    def Cdec(system_state, t):
        pool_dict = dict(zip(pool_order, system_state))

        C_pool = pool_dict['C']
        Fe3_pool = pool_dict['Fe3']
        Acetate_pool = pool_dict['Acetate']
        M_Ac_pool = pool_dict['M_Ac']
        M_Ferm_pool = pool_dict['M_Ferm']
        M_Fe3_pool = pool_dict['M_Fe3']
        M_Hydro_pool = pool_dict['M_Hydro']
        M_Homo_pool = pool_dict['M_Homo']
        M_Ferm2_pool = pool_dict['M_Ferm2']
        CO2_pool = pool_dict['CO2']
        H2_pool = pool_dict['H2']
        Fe2_pool = pool_dict['Fe2']
        # print('C_pool' , C_pool)
        
        
        # Ferm Ferm Ferm Ferm Ferm  Ratio aus Grant 1998
       
        #deltaM_Ferm2, deltaC_pool2, _, Tot_Ferm2 =   Fermenters(M_Ferm2_pool, C_pool, Acetate_pool, Vmax_Ferm, w_Ferm, Sensenmann, Kmb_Ferm, Kmh_Ferm)
        # Produziert:
        #deltaH2_Ferm2 =  -deltaC_pool2 * 1.2 # 1.2 aus cabrol2017microbial
        deltaH2_Ferm2 = 0
        deltaM_Ferm2 = 0
        
        #FERM FERM FERM 
        pool_change_dict_Ferm =  Ferm_Pathway(pool_dict, model_parameter_dict)

        
        # Fe3 Fe3 Fe3 Fe3 : C2H3O2 − + 4H2O + 8Fe3+3 → 9H+ + 2 HCO3− + 8Fe3+2   Delattre 2019
        pool_change_dict_Fe =  Fe3_Pathway(pool_dict, model_parameter_dict)
     
        # ACETO ACETO ACETO ACETO  CH3COO + H+ ->  CH4 + CO2 Fe3y_Conrad2000
        pool_change_dict_Ac =  Ac_Pathway(pool_dict, model_parameter_dict)
    
        #HYDRO HYDRO HYDRO HYDRO 4H2 + CO2 → CH4 + 2H2O, conrad2000selective, Fe3nchel -131kj/mol
        pool_change_dict_Hydro =  Hydro_Pathway(pool_dict, model_parameter_dict)
      
        # HOMO HOMO HOMO HOMO HOMO: 4H2 + 2CO2 → CH3COOH + 2H2O  conrad2000selective
        pool_change_dict_Homo = Homo_Pathway(pool_dict, model_parameter_dict)

        
        
        pool_change_dict_list = [pool_change_dict_Ferm, pool_change_dict_Fe, pool_change_dict_Ac,pool_change_dict_Hydro,pool_change_dict_Homo, pool_change_dict_Ferm]
        
        
        
        # DELTA DELTA DELTA                           
        
        m_C = 12.01*1e-3 # molar mass of carbon
        changes_dict = dict()
       # das übernimmt jetzt die spezielle mikrobe? = changes_dict['C'] = -min(-deltaC_pool, C_pool)  +  (Tot_Hyd + Tot_Ac + Tot_ALtE + Tot_Ferm)/m_C
        
        # hier werden die Beitr�ge aller Mikroben zur Poolgrößenänderung zusammengefasst
        for pool in pool_dict:
            changes_dict[pool] = sum(list([dictionary[pool] for dictionary in pool_change_dict_list if pool in dictionary]))
        
        # print(changes_dict)
        # input('pth..')

       # speziellere aufteilung der Pools fürs plotting
        changes_dict['CO2_Ferm'] = pool_change_dict_Ferm['CO2']  if 'CO2' in pool_change_dict_Ferm else 0.
        changes_dict['CH4_Hydro'] = pool_change_dict_Hydro['CH4']  if 'CH4' in pool_change_dict_Hydro else 0.
        changes_dict['CO2_Ac'] = pool_change_dict_Ac['CO2'] if 'CO2' in pool_change_dict_Ac else 0.
        changes_dict['CO2_Hydro'] = pool_change_dict_Hydro['CO2'] if 'CO2' in pool_change_dict_Hydro else 0.
        changes_dict['CO2_Fe3'] = pool_change_dict_Fe['Fe3'] if 'Fe3' in pool_change_dict_Fe else 0.
        changes_dict['CO2_Homo'] =pool_change_dict_Homo['CO2'] if 'CO2' in pool_change_dict_Homo else 0.
        changes_dict['H2_Homo'] = pool_change_dict_Homo['H2'] if 'H2' in pool_change_dict_Homo else 0.
        changes_dict['H2_Hydro'] = pool_change_dict_Hydro['H2'] if 'H2' in pool_change_dict_Hydro else 0.
        changes_dict['Fe2_Fe'] = pool_change_dict_Fe['Fe2'] if 'Fe2' in pool_change_dict_Fe else 0.
        
        changes_dict['H2_Ferm2'] = 0
        # print('changes CO2_Ferm', changes_dict['CO2_Ferm'])
        # input('')
        
        pool_changes_array = np.array([changes_dict[pool_name] if pool_name in changes_dict else 0.0 for pool_name in pool_order])
        # input('dcdec done..')

        return pool_changes_array
    return Cdec


  



