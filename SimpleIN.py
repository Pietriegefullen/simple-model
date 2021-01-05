
# Das eigentliche  Model, durch das der Verlauf der Kurven berechnet wird. 

################## Variabeln:##################################################
# Cpool:     Cpool der von den Fermentierern genutzt werden kann (Mikromol pro g dw)
# AltEpool : Noch zu reduzierende Alternative E Akzeptoren (Mikromol pro g dw)
# M_A:       Mikrobielle Biomasse der Acetoclastischen Methanogenese
# M_Ferm:    Mikrobielle Biomasse der Fermentierer
# M_AltE:    Mikrobielle Biomasse der Mikroben die Alternative E nutzen
# M_Hydro:   Mikrobielle Biomasse der Hydrogenotrophen Methanogenen
# M_Homo:    Mikrobielle Biomasse der Homoacetoclasten
# CH4:       CH4 Pool
# CO2:       CO2 Pool
# CO2_A:     CO2 dass aus der Acetoclastischen Methanogene ensteht (Nur zum Plotten)
# Acetate:   Acetate Pool
# M_A_CH4_krank: Menge an "kranken" Methanogenen (Relikt aus Vorgängermodel)
# H2:         H2 Pool
# *Fitters:   Zu optimierende Parameter

#################### Parameter:################################################
# Vmax_Ferm:      Maximale Abbaurate der Fermentierer
# Stoch_ALtE:     Stochiometrischer Wert für AltE und Acetate
# Vprod_max_AltE: Maximale Abbaurate der Mikroben, die AltE nutzen
# ATPprod_Ferm:   Mol ATP produziert pro mol Glukose von Fermentierern 
# ATPprod_AltE:   Mol ATP produziert pro mol Acetate von Mikroben, die AltE nutzen
# ATPprod_Hydro:  Mol ATP produziert pro mol H2 & CO2 von Hydrogenotropen
# ATPprod_Homo:   Mol ATP produziert pro mol H2 & CO2 von Homoacetogenen
# ATPprod_Ace:    Mol ATP produziert pro mol Acetate von Acetoclasten
# Yatp_Ferm:      g C Biomasse die Fermentierer aus 1 Mol ATP aufbauen
# Yatp_AltE:      g C Biomasse die Mikroben, die Alt E nutzen,  aus 1 Mol ATP aufbauen
# Yatp_Hydro:     g C Biomasse die Hydrogenotrophe aus 1 Mol ATP aufbauen
# Yatp_Homo:      g C Biomasse die Homoacetogene aus 1 Mol ATP aufbauen
# Yatp_Ace        g C Biomasse die Acetoclasten aus 1 Mol ATP aufbauen

###################### Returns:################################################
#delta"Pool" :  Die Änderungen in den entsprechenden Pools


from Microbe import Fermenters, Hydrotrophes, AltE, Acetoclast, Homo # Importiert Funktionen aus dem Skript "Microbe"


def Cdec(Cpool, AltEpool, M_A, M_Ferm, M_AltE, M_Hydro, M_Homo, CH4, CO2, CO2_A, Acetate, H2, M_A_CH4_krank, *Fitters):    
     
    
    Vmax_Ferm, Stoch_ALtE,Vprod_max_AltE, ATPprod_Ferm, ATPprod_AltE, ATPprod_Hydro,ATPprod_Homo,ATPprod_Ace, Yatp_Ferm, Yatp_AltE, Yatp_Hydro,Yatp_Homo,Yatp_Ace = Fitters
    
    
    # FERM FERM FERM FERM FERM 
    #Ferm atmen einen Teil Cpool und machen einen Teil zu Acetate, CO2 und H2
    # Die Aufteilung folgt Conrad und wird so auch von Song genutzt
    deltaM_Ferm, deltaCpool, _ =   Fermenters(M_Ferm, Cpool, 0, Vmax_Ferm, ATPprod_Ferm, Yatp_Ferm)
    deltaAce_Ferm = - (deltaCpool)
    # Produziert:
    deltaCO2_Ferm = deltaAce_Ferm * 0.5 
    deltaH2_Ferm = deltaAce_Ferm * (1/6)
    
    # ALT E ALT E ALT E ALT E 
    # nur solange AltE UND Acetat vorhanden
    deltaM_AltE, deltaAcetate_AltE, deltaAltEpool =   AltE(M_AltE, Acetate, AltEpool, Stoch_ALtE, Vprod_max_AltE, ATPprod_AltE, Yatp_AltE)
    deltaCO2_Alte = - deltaAcetate_AltE * 2 # pro 1 Acetate entstehen zwei CO2
##    
#    manueller ALtE Abbau (um die anderen Prozesse früher in Gang zu bringen, eigentlich kein echter Teil des Models
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

  
    
    if AltEpool <= 0: # termodyn, AltE besser als Methanogen (Gao 2019), AltE müssen leer sein bevor Ace anfängt
        # Sterben muss noch ins MiKrobenskript 
        
    # Relikt aus Vorgängermodel, im Moment nicht funktionsfähig
       # M_A_CH4_geheilt = 0 if M_A_CH4_krank <= 0 else (w_A_CH4_heil * M_A_CH4_krank)
       # deltaM_A_CH4_krank = - M_A_CH4_geheilt
       # deltaM_A_CH4 = w_A_CH4 *  M_A_CH4  +  M_A_CH4_geheilt
        
        
    # ACETO ACETO ACETO ACETO  
    
        deltaM_A, deltaAcetate_A =   Acetoclast(M_A, Acetate,ATPprod_Ace,Yatp_Ace)
        deltaCH4_A = - deltaAcetate_A * 0.5 # pro mol Acetate entsteht 0.5 Mol CH4
        deltaCO2_A = - deltaAcetate_A * 0.5 # pro mol Acetate entsteht 0.5 Mol CO2
    

        #HYDRO HYDRO HYDRO HYDRO 4H2 + CO2 → CH4 + 2H2O, Fenchel -131kj/mol
        #Evtl aus dem If Statement raus
        #Hydrogenotrophy hat Thermodynamische Vorteile vor Homo, bei hohen Temp
   
        deltaM_Hydro, deltaCO2_Hydro, deltaH2_Hydro =   Hydrotrophes(M_Hydro, CO2, H2,ATPprod_Hydro,Yatp_Hydro)
        deltaCH4_Hydro = - deltaCO2_Hydro # pro mol CO2 entsteht 1 mol CH4 
        
        
      
        # HOMO HOMO HOMO HOMO HOMO: 4H2 + 2CO2 → CH3COOH + 2H2O
        #im If weil, Alt E- verbrauchen H2 und senken den H2 Partialdruck, den Homo braucht(Ye 2013) 
        
        deltaM_Homo_CH4, deltaCO2_Homo ,deltaH2_Homo =   Homo(M_Homo, CO2, H2,ATPprod_Homo,Yatp_Homo)
        deltaAcetate_Homo  = - deltaH2_Homo * 4 # aus 4 mol H2 wird ein mol Acetate
        
       
    
 
    
    # DELTA DELTA DELTA
    deltaCO2 =      deltaCO2_Ferm   + deltaCO2_Alte         + deltaCO2_A       + deltaCO2_Hydro    + deltaCO2_Homo  
    deltaAcetate =  deltaAce_Ferm   + deltaAcetate_AltE     + deltaAcetate_A                       + deltaAcetate_Homo 
    deltaH2 =       deltaH2_Ferm                                               + deltaH2_Hydro     + deltaH2_Homo          
    deltaCH4 =                                                deltaCH4_A       + deltaCH4_Hydro 
    
    deltaCpool =   -min(-deltaCpool, Cpool)
    deltaAcetate = -min(-deltaAcetate, Acetate)
    deltaH2 =      -min(-deltaH2, H2)
    deltaCO2 =     -min(-deltaCO2, CO2)
    
    deltaM_A_CH4_krank = 0 # Relikt
 
    return deltaCpool, deltaAltEpool, deltaM_A, deltaM_Ferm, deltaM_Hydro, deltaM_AltE, deltaM_Homo, deltaCH4, deltaCO2, deltaCO2_A, deltaAcetate, deltaH2 ,deltaM_A_CH4_krank



  



