
# Das eigentliche  Model, durch das der Verlauf der Kurven berechnet wird. 

################## Variabeln:##################################################
# C_pool:       C_pool der von den Fermentierern genutzt werden kann (Mikromol pro g dw)
# Fe_pool :     Noch zu reduzierende Alternative E Akzeptoren (Mikromol pro g dw)
# M_Ac_pool:    Mikrobielle Biomasse der Acetoclastischen Methanogenese
# M_Ferm_pool:  Mikrobielle Biomasse der Fermentierer
# M_Fe_pool:    Mikrobielle Biomasse der Mikroben die Alternative E nutzen
# M_Hydro_pool: Mikrobielle Biomasse der Hydrogenotrophen Methanogenen
# M_Homo_pool:  Mikrobielle Biomasse der Homoacetoclasten
# CH4:          CH4 Pool
# CO2:          CO2 Pool
# CO2_A:        CO2 dass aus der Acetoclastischen Methanogene ensteht (Nur zum Plotten)
# Acetate_pool: Acetate_pool Pool
# H2:           H2 Pool
# CO2_Hydro':   CO2 verbraucht durch Hydro
# CH4_Hydro':   CH4 erzeugt durch Hydro
# H2_Ferm2':    H2 aus Ferm2
# M_Ferm2':     Mikroben Ferm2


#################### Parameter:################################################
#Vmax_Ferm',:       Maximale Reaktionsraten
#Vmax_Fe',:
#Vmax_Homo',:
#Vmax_Hydro',:
#Vmax_Ac',:
#w_Ferm',:          Maximale Wachstumsraten
#w_Fe',:
#w_Hydro',: 
#w_Homo',:
#w_Ac',:
#Sensenmann',:
#Stoch_Fe',:
#Kmb_Ferm',:        MM werte
#Kmh_Ferm',:
#Kmb_Fe',:
#Kmb_Ac',:
#Kmb_Hydro',:
#Fe':               Initialer Fe Pool
###################### Returns:################################################
#delta"Pool" :  Die Änderungen in den entsprechenden Pools

import numpy as np
from order import pool_order

from Microbe import Fermenters, Hydrotrophes, AltE, Acetoclast, Homo # Importiert Funktionen aus dem Skript "Microbe"

def Cdec_wrapper(model_parameter_dict):
    Vmax_Ferm = model_parameter_dict['Vmax_Ferm']
    Vmax_Fe = model_parameter_dict['Vmax_Fe']
    Vmax_Ac = model_parameter_dict['Vmax_Ac']
    Vmax_Hydro = model_parameter_dict['Vmax_Hydro']
    Vmax_Homo = model_parameter_dict['Vmax_Homo']
    w_Ferm = model_parameter_dict['w_Ferm']
    w_Fe = model_parameter_dict['w_Fe']
    w_Ac = model_parameter_dict['w_Ac']
    w_Hydro = model_parameter_dict['w_Hydro']
    w_Homo = model_parameter_dict['w_Homo']
    Kmb_Ferm = model_parameter_dict['Kmb_Ferm']
    Kmh_Ferm = model_parameter_dict['Kmh_Ferm']
    Kmb_Fe = model_parameter_dict['Kmb_Fe']
    Kmb_Ac = model_parameter_dict['Kmb_Ac']
    Kmb_Hydro = model_parameter_dict['Kmb_Hydro']
    Stoch_Fe = model_parameter_dict['Stoch_Fe']
    Sensenmann = model_parameter_dict['Sensenmann']
        
    def Cdec(system_state, t):
        pool_dict = dict(zip(pool_order, system_state))
        
        C_pool = pool_dict['C']
        Fe_pool = pool_dict['Fe']
        Acetate_pool = pool_dict['Acetate']
        M_Ac_pool = pool_dict['M_Ac']
        M_Ferm_pool = pool_dict['M_Ferm']
        M_Fe_pool = pool_dict['M_Fe']
        M_Hydro_pool = pool_dict['M_Hydro']
        M_Homo_pool = pool_dict['M_Homo']
        M_Ferm2_pool = pool_dict['M_Ferm2']
        CO2_pool = pool_dict['CO2']
        H2_pool = pool_dict['H2']
        
        
        
        # FERM FERM FERM FERM FERM  Ratio aus Grant 1998
        #Ferm atmen einen Teil C_pool und machen einen Teil zu Acetate_pool, CO2 und H2
        # Die Aufteilung folgt Conrad und wird so auch von Song genutzt
        #deltaM_Ferm, deltaC_pool, _ =   Fermenters(M_Ferm_pool, C_pool, Acetate_pool, Vmax_Ferm, ATPprod_Ferm, Yatp_Ferm)
        deltaM_Ferm, deltaC_pool, _, Tot_Ferm =   Fermenters(M_Ferm_pool, C_pool, Acetate_pool, Vmax_Ferm, w_Ferm, Sensenmann, Kmb_Ferm, Kmh_Ferm)
        deltaAce_Ferm = - (deltaC_pool)
        # Produziert:
        deltaCO2_Ferm = deltaAce_Ferm * 0.5 
        deltaH2_Ferm = deltaAce_Ferm *(1/6) 
      
        
        
        #deltaM_Ferm2, deltaC_pool2, _, Tot_Ferm2 =   Fermenters(M_Ferm2_pool, C_pool, Acetate_pool, Vmax_Ferm, w_Ferm, Sensenmann, Kmb_Ferm, Kmh_Ferm)
        # Produziert:
        #deltaH2_Ferm2 =  -deltaC_pool2 * 1.2 # 1.2 aus cabrol2017microbial
        deltaH2_Ferm2 = 0
        deltaM_Ferm2 = 0
        
        # FE FE FE FE : C2H3O2 − + 4H2O + 8Fe+3 → 9H+ + 2 HCO3− + 8Fe+2   Delattre 2019
        # nur solange AltE UND Acetat vorhanden
        deltaM_Fe, deltaAcetate_Fe, deltaAltEpool, Tot_ALtE =   AltE(M_Fe_pool, Acetate_pool, Fe_pool, Stoch_Fe, Vmax_Fe, w_Fe, Sensenmann, Kmb_Fe)
        deltaCO2_Alte = - deltaAcetate_Fe * 2 # pro 1 Acetate_pool entstehen zwei CO2
        deltaH2_Alte = - deltaAcetate_Fe * 0 # was ist der wirkliche Wert? 
     
      

        # ACETO ACETO ACETO ACETO  CH3COO + H+ ->  CH4 + CO2 Fey_Conrad2000
    
        deltaM_A, deltaAcetate_A,Tot_Ac =   Acetoclast(M_Ac_pool, Acetate_pool,w_Ac, Vmax_Ac, Sensenmann, Kmb_Ac)
        deltaCH4_A = - deltaAcetate_A * 0.5 # pro mol Acetate_pool entsteht 0.5 Mol CH4
        deltaCO2_A = - deltaAcetate_A * 0.5 # pro mol Acetate_pool entsteht 0.5 Mol CO2
    
    
        #HYDRO HYDRO HYDRO HYDRO 4H2 + CO2 → CH4 + 2H2O, conrad2000selective, Fenchel -131kj/mol

       
        deltaM_Hydro, deltaCO2_Hydro, deltaH2_Hydro, Tot_Hyd =   Hydrotrophes(M_Hydro_pool, CO2_pool, H2_pool, w_Hydro, Vmax_Hydro, Sensenmann, Kmb_Hydro)
        deltaCH4_Hydro = - deltaCO2_Hydro # pro mol CO2 entsteht 1 mol CH4 
        
        
      
        # HOMO HOMO HOMO HOMO HOMO: 4H2 + 2CO2 → CH3COOH + 2H2O  conrad2000selective

        
        deltaM_Homo, deltaCO2_Homo ,deltaH2_Homo, Tot_Homo =   Homo(M_Homo_pool, CO2_pool, H2_pool, w_Homo, Vmax_Homo, Sensenmann)
       # deltaM_Homo_CH4, deltaCO2_Homo ,deltaH2_Homo =   0,0,0
        deltaAcetate_Homo  = - deltaH2_Homo * 4 # aus 4 mol H2 wird ein mol Acetate_pool
        #deltaM_Homo = 0   
        #deltaH2_Homo = 0
        
        
        
        
        
        # DELTA DELTA DELTA
        deltaCO2 =      deltaCO2_Ferm   + deltaCO2_Alte         + deltaCO2_A       + deltaCO2_Hydro                   + deltaCO2_Homo  
        deltaAcetate =  deltaAce_Ferm   + deltaAcetate_Fe     + deltaAcetate_A                                      + deltaAcetate_Homo 
        deltaH2 =       deltaH2_Ferm    + deltaH2_Alte                             + deltaH2_Hydro   + deltaH2_Ferm2  + deltaH2_Homo          
        deltaCH4 =                                                deltaCH4_A       + deltaCH4_Hydro 
        
        m_C = 12.01*1e-3 # molar mass of carbon
        changes_dict = dict()
        changes_dict['C'] = -min(-deltaC_pool, C_pool)  +  (Tot_Hyd + Tot_Ac + Tot_ALtE + Tot_Ferm)/m_C #+ Tot_Homo
        changes_dict['Acetate'] = -min(-deltaAcetate, Acetate_pool)
        changes_dict['H2'] = -min(-deltaH2, H2_pool)
        changes_dict['CO2'] = -min(-deltaCO2, CO2_pool)
        changes_dict['CH4'] = deltaCH4
        
        changes_dict['Fe'] = deltaAltEpool
        changes_dict['M_Ac'] = deltaM_A
        changes_dict['M_Ferm'] = deltaM_Ferm
        changes_dict['M_Fe'] = deltaM_Fe
        changes_dict['M_Hydro'] = deltaM_Hydro
        changes_dict['M_Homo'] = deltaM_Homo
        changes_dict['M_Ferm2'] = deltaM_Ferm2
        changes_dict['CH4_Hydro'] = deltaCH4_Hydro
        changes_dict['CO2_Ac'] = deltaCO2_A
        changes_dict['CO2_Hydro'] = deltaCO2_Hydro
        changes_dict['CO2_Ferm'] = deltaCO2_Ferm
        changes_dict['CO2_Fe'] = deltaCO2_Alte
        changes_dict['CO2_Homo'] = deltaCO2_Homo  
        changes_dict['H2_Homo'] = deltaH2_Homo
        changes_dict['H2_Ferm2'] = deltaH2_Ferm2
        changes_dict['H2_Hydro'] = deltaH2_Hydro
        
        pool_dict_changes = np.array([changes_dict[pool_name] for pool_name in pool_order])
        return pool_dict_changes
    return Cdec


  



