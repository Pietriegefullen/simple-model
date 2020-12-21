
# Das eigentliche mathematische Model, durch das der Verlauf der Kurven berechnet wird. 

# Variabeln:

# Parameter:

# Returns:

from Microbe import Fermenters, Hydrotrophes, AltE, Acetoclast, Homo


def Cdec(Cpool, AltEpool, M_A, M_Ferm, M_AltE, M_Hydro, M_Homo, CH4, CO2, CO2_A, Acetate, M_A_CH4_krank, H2, *Fitters):    
     
    
    Vmax_Ferm, Stoch_ALtE,Vprod_max_AltE, ATPprod_Ferm, ATPprod_AltE, ATPprod_Hydro,ATPprod_Homo,ATPprod_Ace = Fitters
    
    
    # FERM FERM FERM FERM FERM 
    #Ferm atmen einen Teil Cpool und machen einen Teil zu Acetate, CO2 und H2
    deltaM_Ferm, deltaCpool, _ =   Fermenters(M_Ferm, Cpool, 0, Vmax_Ferm, ATPprod_Ferm)
    deltaAce_Ferm =-(deltaCpool)
    # Produziert:
    deltaCO2_Ferm = deltaAce_Ferm * 0.5 
    deltaH2_Ferm = deltaAce_Ferm * (1/6)

    Acetate_tot  = 0 if  Acetate <= 0 else Acetate
    
    # ALT E ALT E ALT E ALT E 
    # nur solange Alt e UND Acetat vorhanden
    deltaM_AltE, deltaAcetate_AltE, deltaAltEpool =   AltE(M_AltE, Acetate_tot, AltEpool, Stoch_ALtE, Vprod_max_AltE,ATPprod_AltE)
    deltaCO2_Alte = - deltaAcetate_AltE * 2 #1 Acetate wird zu zwei CO2
    deltaAltEpool  = - min(-deltaAltEpool, AltEpool)
##    
#    
#    deltaAltEpool = - min(AltEpool, 2)  
#    deltaAcetate_AltE = - min(-deltaAltEpool *0.01, Acetate)
#    deltaCO2_Alte = - deltaAcetate_AltE * 2
#    deltaM_AltE = 0
    
    # HOMO HOMO HOMO HOMO  

    deltaCO2_Homo = 0
    deltaH2_Homo  = 0
    deltaAcetate_Homo  = 0   
    deltaM_Homo = 0  

    # ACETO ACETO ACETO ACETO 
    deltaAcetate_A = 0
    deltaCH4_A = 0
    deltaCO2_A = 0
    deltaM_A = 0
    
    #HYDRO HYDRO HYDRO HYDRO
    deltaCH4_Hydro = 0
    deltaCO2_Hydro = 0
    deltaH2_Hydro = 0
    deltaM_Hydro = 0

  
    
    if AltEpool <= 0:
        #Ace_used_AltE_resp = 0
        # STerben muss noch ins microbenskript deltaM_AltE = - min(M_AltE * w_alte, M_AltE) # damit keine negativen Microben entstehen
        
           
       # M_A_CH4_geheilt = 0 if M_A_CH4_krank <= 0 else (w_A_CH4_heil * M_A_CH4_krank)
       # deltaM_A_CH4_krank = - M_A_CH4_geheilt
       # deltaM_A_CH4 = w_A_CH4 *  M_A_CH4  +  M_A_CH4_geheilt
        
        
    # ACETO ACETO ACETO ACETO 
    # IN: Acetat 
    # Out: CO2, CH4 1:1
    # nur solange Acetat und Microben_CH4 vorhanden 
    
        deltaM_A, deltaAcetate_A =   Acetoclast(M_A, Acetate,ATPprod_Ace)
        deltaCH4_A = - deltaAcetate_A * 0.5
        deltaCO2_A = - deltaAcetate_A * 0.5
    

        #HYDRO HYDRO HYDRO HYDRO 4H2 + CO2 → CH4 + 2H2O, Fenchel -131kj/mol
   
        deltaM_Hydro, deltaCO2_Hydro, deltaH2_Hydro =   Hydrotrophes(M_Hydro, CO2, H2,ATPprod_Hydro)
        deltaCH4_Hydro = - deltaCO2_Hydro # weil pro CO2 verbraucht, ein CH4 entsteht
        
        
      
        # HOMO HOMO HOMO HOMO HOMO: 4H2 + 2CO2 → CH3COOH + 2H2O

        #im If weil solange Alt E- hat Hydro Thermodynamische Vorteile, 
        #weil S und Fe H2 verbrauchen und der partialdruck zu gering ist (Ye 2013). Sonst 0
        
        deltaM_Homo_CH4, deltaCO2_Homo ,deltaH2_Homo =   Homo(M_Homo, CO2, H2,ATPprod_Homo)
        #deltaM_H_CH4, deltaAcetate_A , =   Homo(M_Homo, CO2, H2)
        deltaAcetate_Homo  = - deltaH2_Homo * 4 # aus 4 H2 wird ein Acetate
 
        
    
 
    
    # DELTA DELTA DELTA
    deltaCO2 =      deltaCO2_Ferm   + deltaCO2_Alte         + deltaCO2_A       + deltaCO2_Hydro    + deltaCO2_Homo  
    deltaAcetate =  deltaAce_Ferm   + deltaAcetate_AltE     + deltaAcetate_A                       + deltaAcetate_Homo 
    deltaH2 =       deltaH2_Ferm                                               + deltaH2_Hydro     - deltaH2_Homo          
    deltaCH4 =                                                deltaCH4_A       + deltaCH4_Hydro 
    
    
    
    deltaM_A_CH4_krank = 0
    
    deltaCpool  = 0 if  Cpool + deltaCpool < 0 else deltaCpool
   
 
    return deltaCpool, deltaAltEpool, deltaM_A, deltaM_Ferm, deltaM_Hydro, deltaM_AltE, deltaM_Homo, deltaCH4, deltaCO2, deltaCO2_A, deltaAcetate, deltaM_A_CH4_krank, deltaH2



  



