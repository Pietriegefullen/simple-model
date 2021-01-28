from SimpleIN import Cdec # importiert das eigentliche mathematische Model 

# OPTIFUN 
# optifun bringt simplefun in ein format (added xdata), das von scipy für Curvefit benötigt wird

def optifun(xdata, *Fitters):    
    xtage = xdata[0:int(len(xdata)/2)]
 
    CH4, CO2,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_= simplefun(xtage,*Fitters)
    return CH4 + CO2 # setzt die für uns interessanten Ausgaben in curvefit format zusammen
# CH4 und CO2 werden aneinandergehängt um einen abhängigen Least Sqaure zu berechnen und nicht jeweils einen unabhängigen pro Kurve  
    
# SIMPELFUN, der Euler Forward Mechanismus

def simplefun(xtage, *Fitters):    
    
    step = 1/24.0
    
    n = int(max(xtage)/step)
    number_steps = range(1, n+1)
    
    # Festgelegte Initialwerte
    # Die Berechnung hier für Cpool ist noch für zwei Pools gedacht. Momentan nicht relevant
    m_gluc = 180   # molar mass of glucose, g dw pro mol
    m_cell = 162   # molar mass of cellulose, g dw pro mol
    TOC = 0.04     # Knoblauchs Daten , g dw
    labile = 0.005 # Knoblauchs Daten, anteil von TOC einheitslos
    Cpool_init = (TOC*labile/m_gluc + TOC*(1-labile)/m_cell) * (10**6) # 0.00024679012345679013 * (10**6) mikromol pro g
    
    # die Werte sind alle in mikroMol pro gram Trockengewicht
    AltE_init = Fitters[-1] # 7 #146       # 146 cf. Yao, Conrad 1999 (Hyun2017 sagt werte um 100 sind extrem hoch), Philben: 6 jeweils Mikromol pro g dw
    Fitters = Fitters[:-1]
    M_A_CH4_init = 0.001     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Ferm_init = 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    M_AltE_init= 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    deltaH2_Hydro = 0        # für plot in multifit
    M_H_CH4_init = 0.2     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Homo_init = 0.2        # Monteux 2020, mg Mikrobielles C pro g dw
    CH4_init = 0
    CO2_init = 0
    Acetate_init = 0        # Philben wert ca 3, in mikromol pro g 
    AceCO2_init = 0
    H2_init = 0
    deltaH2_Homo_init = 0
    CO2_Hydro_init= 0
    CH4_Hydro_init = 0
    H2_Ferm2_init = 0
    M_Ferm2_init = 0.2
    
    Cpool    = [Cpool_init]
    AltEpool = [AltE_init]
    M_A_CH4  = [M_A_CH4_init]
    M_Ferm    = [M_Ferm_init]
    M_AltE   = [M_AltE_init]
    deltaH2_Hydro = [deltaH2_Hydro]
    M_H_CH4  = [M_H_CH4_init]
    M_Homo   = [M_Homo_init]
    CH4      = [CH4_init]
    CO2      = [CO2_init]
    Acetate  = [Acetate_init]
    AceCO2   = [AceCO2_init]
    H2       = [H2_init]
    deltaH2_Homo = [deltaH2_Homo_init]
    CO2_Hydro = [CO2_Hydro_init]
    CH4_Hydro = [CH4_Hydro_init]
    H2_Ferm2 = [H2_Ferm2_init]
    M_Ferm2 = [M_Ferm2_init]

    
    for t in number_steps: # iteriert über n+1 Zeitschritte (n = max(xtage))

        delta = Cdec(Cpool[-1], AltEpool[-1], M_A_CH4[-1],M_Ferm[-1], M_AltE[-1],M_H_CH4[-1],M_Homo[-1], CH4[-1], CO2[-1], AceCO2[-1], Acetate[-1], H2[-1],CO2_Hydro[-1], CH4_Hydro[-1], H2_Ferm2[-1], M_Ferm2[-1], *Fitters)

        
        Cpool.append(Cpool[-1] +       step*delta[0])
        AltEpool.append(AltEpool[-1] + step*delta[1])
        M_A_CH4.append(M_A_CH4[-1] +   step*delta[2])
        M_Ferm.append(M_Ferm[-1] +     step*delta[3])
        M_AltE.append(M_AltE[-1] +     step*delta[4])
        M_H_CH4.append(M_H_CH4[-1] +   step*delta[5])
        M_Homo.append(M_Homo[-1] +     step*delta[6])
        CH4.append(CH4[-1] +           step*delta[7])
        CO2.append(CO2[-1] +           step*delta[8])
        AceCO2.append(AceCO2[-1] +     step*delta[9])
        Acetate.append(Acetate[-1] +   step*delta[10])
        H2.append(H2[-1] +             step*delta[11])
        deltaH2_Hydro.append(          step*delta[12])
        deltaH2_Homo.append(           step*delta[13])
        CO2_Hydro.append(              step*delta[14])
        CH4_Hydro.append(              step*delta[15])
        H2_Ferm2.append(               step*delta[16])
        M_Ferm2.append( M_Ferm2[-1] +  step*delta[17])
    
    
   # AltEpool[0::int(1/step)]
    
    CH4 = [CH4[int(i/step)] for i in xtage]
    CO2 = [CO2[int(i/step)] for i in xtage]
    AltEpool = AltEpool[0::int(1/step)] # von 0 bis Ende in Schritten von 24, damit der Plot nicht in Stunden sondern in Tagen ist
    AceCO2 = AceCO2[0::int(1/step)]
    Acetate= Acetate[0::int(1/step)]
    Cpool= Cpool[0::int(1/step)]
    M_A_CH4= M_A_CH4[0::int(1/step)]
    M_Ferm= M_Ferm[0::int(1/step)]
    M_AltE= M_AltE[0::int(1/step)]
    H2= H2[0::int(1/step)]
    M_H_CH4= M_H_CH4[0::int(1/step)]
    M_Homo= M_Homo[0::int(1/step)]
    AltE_init= AltE_init
    deltaH2_Hydro= deltaH2_Hydro[0::int(1/step)]
    deltaH2_Homo= deltaH2_Homo[0::int(1/step)]
    CO2_Hydro= CO2_Hydro[0::int(1/step)]
    CH4_Hydro= CH4_Hydro[0::int(1/step)]
    H2_Ferm2= H2_Ferm2[0::int(1/step)]

    return CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_Ferm, M_AltE, H2, M_H_CH4, M_Homo, AltE_init, deltaH2_Hydro,deltaH2_Homo, CO2_Hydro, CH4_Hydro,H2_Ferm2