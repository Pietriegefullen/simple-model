
# Das eigentliche  Model, durch das der Verlauf der Kurven berechnet wird. 

################## Variabeln:##################################################
# C_pool:     C_pool der von den Fermentierern genutzt werden kann (Mikromol pro g dw)
# Fe_pool : Noch zu reduzierende Alternative E Akzeptoren (Mikromol pro g dw)
# M_Ac_pool:       Mikrobielle Biomasse der Acetoclastischen Methanogenese
# M_Ferm_pool:    Mikrobielle Biomasse der Fermentierer
# M_Fe_pool:    Mikrobielle Biomasse der Mikroben die Alternative E nutzen
# M_Hydro_pool:   Mikrobielle Biomasse der Hydrogenotrophen Methanogenen
# M_Homo_pool:    Mikrobielle Biomasse der Homoacetoclasten
# CH4:       CH4 Pool
# CO2:       CO2 Pool
# CO2_A:     CO2 dass aus der Acetoclastischen Methanogene ensteht (Nur zum Plotten)
# Acetate_pool:   Acetate_pool Pool
# M_A_CH4_krank: Menge an "kranken" Methanogenen (Relikt aus Vorgängermodel)
# H2:         H2 Pool
# *Fitters:   Zu optimierende Parameter

#################### Parameter:################################################
# Vmax_Ferm:      Maximale Abbaurate der Fermentierer
# Stoch_Fe:     Stochiometrischer Wert für AltE und Acetate_pool
# Vmax_Fe: Maximale Abbaurate der Mikroben, die AltE nutzen
# ATPprod_Ferm:   Mol ATP produziert pro mol Glukose von Fermentierern 
# ATPprod_Fe:   Mol ATP produziert pro mol Acetate_pool von Mikroben, die AltE nutzen
# ATPprod_Hydro:  Mol ATP produziert pro mol H2 & CO2 von Hydrogenotropen
# ATPprod_Homo:   Mol ATP produziert pro mol H2 & CO2 von Homoacetogenen
# ATPprod_Ac:    Mol ATP produziert pro mol Acetate_pool von Acetoclasten
# Yatp_Ferm:      g C Biomasse die Fermentierer aus 1 Mol ATP aufbauen
# Yatp_Fe:      g C Biomasse die Mikroben, die Alt E nutzen,  aus 1 Mol ATP aufbauen
# Yatp_Hydro:     g C Biomasse die Hydrogenotrophe aus 1 Mol ATP aufbauen
# Yatp_Homo:      g C Biomasse die Homoacetogene aus 1 Mol ATP aufbauen
# Yatp_Ac        g C Biomasse die Acetoclasten aus 1 Mol ATP aufbauen

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
        #CH4_pool = pool_dict['CH4']
        #CH4_Hydro_pool = pool_dict['CH4_Hydro']
        CO2_pool = pool_dict['CO2']
        #CO2_Ac_pool = pool_dict['CO2_Ac']
        #CO2_Hydro_pool = pool_dict['CO2_Hydro']
        #CO2_Ferm_pool = pool_dict['CO2_Ferm']
        #CO2_Fe_pool = pool_dict['CO2_Fe']
        #CO2_Homo_pool = pool_dict['CO2_Homo']   
        H2_pool = pool_dict['H2']
        #H2_Homo_pool = pool_dict['H2_Homo']
        #H2_Ferm2_pool = pool_dict['H2_Ferm2']
        #H2_Hydro_pool = pool_dict['H2_Hydro']
        
        
         # FERM FERM FERM FERM FERM 
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
        
        # ALT E ALT E ALT E ALT E 
        # nur solange AltE UND Acetat vorhanden
        deltaM_Fe, deltaAcetate_Fe, deltaAltEpool, Tot_ALtE =   AltE(M_Fe_pool, Acetate_pool, Fe_pool, Stoch_Fe, Vmax_Fe, w_Fe, Sensenmann, Kmb_Fe)
        deltaCO2_Alte = - deltaAcetate_Fe * 2 # pro 1 Acetate_pool entstehen zwei CO2
        deltaH2_Alte = - deltaAcetate_Fe #* 0 # was ist der wirkliche Wert? 
     
    #    manueller ALtE Abbau (um die anderen Prozesse früher in Gang zu bringen, eigentlich kein echter Teil des Models
    #    deltaAltEpool = - min(Fe_pool, 2)  
    #    deltaAcetate_Fe = - min(-deltaAltEpool *0.01, Acetate_pool)
    #    deltaCO2_Alte = - deltaAcetate_Fe * 2
    #    deltaM_Fe = 0
    #    
    #    # HOMO HOMO HOMO HOMO  
    #    deltaCO2_Homo = 0
    #    deltaH2_Homo  = 0
    #    deltaAcetate_Homo  = 0   
    #    deltaM_Homo = 0  
    #
    #    # ACETO ACETO ACETO ACETO 
    #    deltaAcetate_A = 0
    #    deltaCH4_A = 0
    #    deltaCO2_A = 0
    #    deltaM_A = 0
    #    
    #    #HYDRO HYDRO HYDRO HYDRO
    #    deltaCH4_Hydro = 0
    #    deltaCO2_Hydro = 0
    #    deltaH2_Hydro = 0
    #    deltaM_Hydro = 0
    #
    #    Tot_Homo = 0
    #    Tot_Hyd = 0 
    #    Tot_Ac = 0
      
        
    #    if Fe_pool <= 0: # termodyn, AltE besser als Methanogen (Gao 2019), AltE müssen leer sein bevor Ace anfängt
            # Sterben muss noch ins MiKrobenskript 
            
        # Relikt aus Vorgängermodel, im Moment nicht funktionsfähig
           # M_A_CH4_geheilt = 0 if M_A_CH4_krank <= 0 else (w_A_CH4_heil * M_A_CH4_krank)
           # deltaM_A_CH4_krank = - M_A_CH4_geheilt
           # deltaM_A_CH4 = w_A_CH4 *  M_A_CH4  +  M_A_CH4_geheilt
            
            
    # ACETO ACETO ACETO ACETO  
    
        deltaM_A, deltaAcetate_A,Tot_Ac =   Acetoclast(M_Ac_pool, Acetate_pool,w_Ac, Vmax_Ac, Sensenmann, Kmb_Ac)
        deltaCH4_A = - deltaAcetate_A * 0.5 # pro mol Acetate_pool entsteht 0.5 Mol CH4
        deltaCO2_A = - deltaAcetate_A * 0.5 # pro mol Acetate_pool entsteht 0.5 Mol CO2
    
    
        #HYDRO HYDRO HYDRO HYDRO 4H2 + CO2 → CH4 + 2H2O, Fenchel -131kj/mol
        #Evtl aus dem If Statement raus
        #Hydrogenotrophy hat Thermodynamische Vorteile vor Homo, bei hohen Temp
       
        deltaM_Hydro, deltaCO2_Hydro, deltaH2_Hydro, Tot_Hyd =   Hydrotrophes(M_Hydro_pool, CO2_pool, H2_pool, w_Hydro, Vmax_Hydro, Sensenmann, Kmb_Hydro)
        deltaCH4_Hydro = - deltaCO2_Hydro # pro mol CO2 entsteht 1 mol CH4 
        
        
      
        # HOMO HOMO HOMO HOMO HOMO: 4H2 + 2CO2 → CH3COOH + 2H2O
        #im If weil, Alt E- verbrauchen H2 und senken den H2 Partialdruck, den Homo braucht(Ye 2013) 
        
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


  



