

import data
import CONSTANTS
import OPTIMIZATION_PARAMETERS

def default_model_parameters(specimen_index = None, site = 'all'):

        #Startparameter für das slidertool
    model_parameters =     {'M_Fe3':            0.2,
                            'M_Ferm':           0.02,
                            'M_Hydro':          0.0025,
                            'M_Homo':           0.0001,
                            'M_Ac':             0.01010,

                            'Sensenmann':       8.33e-5,

                            'Vmax_Fe3':        0.0053274970925869525,
                            'Vmax_help_Ferm':   0.347557284254151,
                            'Vmax_Ferm':       2.1826192982834263,
                            'Vmax_Homo':        0.3842183674485621,
                            'Vmax_Hydro':       0.6787923256204944,
                            'Vmax_Ac':         0.6185701320273594,

                            'Kmb_help_Ferm':    1322.0602068699052,
                            'Inhibition_Ferm':  10.118875996129882,
                            
                            "Km_help_Ferm":    5188.51170010838/CONSTANTS.SOIL_DENSITY,
                            "Km_Ac_Acetate" :  7.030348185396872/ CONSTANTS.SOIL_DENSITY, 
                            "Km_Homo_CO2"  :   848.9503500185294 / CONSTANTS.SOIL_DENSITY,
                            "Km_Homo_H2"    :  584.5061484463681 / CONSTANTS.SOIL_DENSITY ,
                            "Km_Hydro_CO2"  :  869.0991529683715/CONSTANTS.SOIL_DENSITY ,
                            "Km_Hydro_H2"    : 489.1731054166966 / CONSTANTS.SOIL_DENSITY,
                            "Km_Fe3_Fe3"   :   533.3685302693583/CONSTANTS.SOIL_DENSITY,
                            "Km_Fe3_Acetate"  :818.2462137022104/CONSTANTS.SOIL_DENSITY,
                            "Km_Ferm"       :  630.4916492204161 / CONSTANTS.SOIL_DENSITY,
                                                
                            

                            'Fe3':              117.37574379080921,
                            
                            'Acetate':          1,
                            
                            'temperature':      4.0,
                            
                           "CUE_Ferm": 0.3698723525408035,
                            "CUE_Fe3": 0.5572385007984572,
                            "CUE_Ac": 0.08453242432990166,
                            "CUE_Homo": 0.3602422303950234,
                            "CUE_Hydro": 0.7433100141164364
                            }
    

    model_parameters_2 = {}        #   init    lower upper  ok guesses to start with
    model_parameters_2['Vmax_help_Ferm'] =  0.9999306161130834  # 0.05
    model_parameters_2['Vmax_Ferm'] =      4.1380785394837005  # 0.011       # Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
    model_parameters_2['Vmax_Fe3'] =       0.89817626250035  # 0.8         # Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    model_parameters_2['Vmax_Homo'] =       0.5457394394384627  # 0.869       # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    model_parameters_2['Vmax_Hydro'] =      0.24706244797942814   # 0.182 1.8   # 0.15 mikromol pro cm^3 from Song
    model_parameters_2['Vmax_Ac'] =         0.5665115584644631  # 0.99           # Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
    model_parameters_2['Sensenmann'] =      8.33e-5# (8.33e-5, 0, 8.44e-5)# 0
    model_parameters_2['Kmb_help_Ferm'] =   891.7732649811559      # 10
    model_parameters_2["Km_help_Ferm"] =    3618.5403636705505
    model_parameters_2["Km_Ac_Acetate"] =   77.45512661594466
    model_parameters_2["Km_Homo_CO2"] =     267.17305425499063
    model_parameters_2["Km_Homo_H2"] =      752.0140435559757
    model_parameters_2["Km_Hydro_CO2"] =    357.87351679143444 
    model_parameters_2["Km_Hydro_H2"] =     424.58131001540846 
    model_parameters_2["Km_Fe3_Fe3"] =      738.8197995059913
    model_parameters_2["Km_Fe3_Acetate"] =  458.2426314432057
    model_parameters_2["Km_Ferm"] =         780.2347021088556 
    model_parameters_2['Fe3'] =             119.09123130773733    # 15.587,
    model_parameters_2['Acetate'] =         1
    model_parameters_2['M_Ac'] =            0.015401704858753227 # 0.002
    model_parameters_2['Inhibition_Ferm']=  7.350034345938164 # Je niedriger desto hemmung # Diese Boundaries müssen anhander Acetatekurven angepasst werden
    model_parameters_2['CUE_Ferm']=         0.22251799068790443
    model_parameters_2['CUE_Fe3']=           0.308964760641128
    model_parameters_2['CUE_Ac']=           0.26964564053185136
    model_parameters_2['CUE_Homo']=         0.7182020298979199
    model_parameters_2['CUE_Hydro']=        0.6157316796638491
    model_parameters.update(model_parameters_2)

    specimen_model_parameters = {}
    if not specimen_index is None:
        specimen_model_parameters = data.model_parameters_from_data(specimen_index,
                                                                    site = site)
        model_parameters.update(specimen_model_parameters)

    print('model parameters:')
    print('================')
    for k, v in model_parameters.items():
        source = 'specimen data' if k in specimen_model_parameters else ''
        print(f'   {k[:20]:20} {v:10g} {source}')
    print('')

    return model_parameters


def Ferm_help(model_parameters):
    microbe = {'name' :         'M_Ferm',
               'vmax':          model_parameters['Vmax_help_Ferm'],
               'death_rate':    0,
               'Kmb':           model_parameters['Kmb_help_Ferm'], # MM Faktor für die Exoenzyme
               'CUE':           0, # weil Ferm_help nicht wächst (nur Ferm wächst)
               'thermodynamics': False
               }

    educts =  [{'name':         'C',
                'stoich':       1,
                'Km':           model_parameters['Km_help_Ferm']} #10000/CONSTANTS.SOIL_DENSITY} # For Ferm_help, Km must be 0 so that MM_factors evaluates to 1.0, Process is Enzyme Limited not Substrate limited
               ]

    products = [{'name':        'DOC',
                 'stoich':      1}
                ]

    return microbe, educts, products

def Ferm(model_parameters):
    # laut wikipedia: C6H12O6 + 2 H2O → 2 CH3CO2H + 2 CO2 + 4 H2

    microbe = {'name':          'M_Ferm',
               'vmax':          model_parameters['Vmax_Ferm'],
               'death_rate':    model_parameters['Sensenmann'],
               'microbe'       : "Ferm" ,
               'CUE'           :   model_parameters['CUE_Ferm'],
               'C_source'       : 'DOC',
               'thermodynamics': False}


    educts =  [{'name':         'DOC',
                'stoich':       6, # warum 6 ? warum nicht 1? 
                'Km':           model_parameters['Km_Ferm'], #100 / CONSTANTS.SOIL_DENSITY,      # 10 from Song  mikromol pro gram
                'C_atoms':      6}    # weil glucose (und andere Monomere) 6 C atome hat und ein momomer aus der spaltung von Coellulose ist
               ]

    #------------------------------------------------------------------------------------------------------------
    #laut Grant A :6 , CO2:3, H2 : 1
    #laut Knoblauch A :6 , CO2:3, H2 : 6
    #laut mir DOC = A: 7 ,CO2 :6 H2:12 (teile durch 2 weil dann näher am Grant verhältnis)
    #laut neumann2016modeling  Glucose = A:2 CO2:2 H2:4  C6H12O6 +2H2O -> 2CH3COOH+2CO2 + 4H2
    products = [{'name':        'Acetate',
                 'stoich':       3.5, # 3.5
                 'inhibition':   model_parameters['Inhibition_Ferm']}  , # MM faktor für Acetathemmung

                {'name':        'CO2',
                 'stoich':       3}  ,

                {'name':        'H2',
                 'stoich':      2}   # laut berechnung eine 6
                ]
    return microbe, educts, products

def Fe3(model_parameters):
    # PATHWAY:   C2H3O2 − + 4H2O + 8Fe3(III)   --->       9H+ + 2 HCO3− + 8Fe2+   Delattre 2019 , jason 2001
    # H2 + 2Fe3+ -> 2H+ +2Fe2+ aus der Vorlesung von Christian K.
    microbe = {'name':          'M_Fe3',
               'C_source':      'Acetate',
               'CUE':           model_parameters['CUE_Fe3'],
               'death_rate':    model_parameters['Sensenmann'],
               'vmax':          model_parameters['Vmax_Fe3']}

    educts = [{'name':          'Acetate',
               'stoich':        1,
               'Km':            model_parameters['Km_Fe3_Acetate'] ,#0.01/CONSTANTS.SOIL_DENSITY,
               'C_atoms':       2},
              
              {'name':          'H2O',
               'stoich':        1,
               'Km':            0 }, # kein Einfluss von H2O
               
              {'name':          'Fe3',
               'stoich':        8,
               'Km':            model_parameters['Km_Fe3_Fe3'],}] #2/CONSTANTS.SOIL_DENSITY}]

    products = [{'name':        'Fe2',
                 'stoich':      8},
                

                {'name':        'CO2', # proxy für HCO3
                 'stoich':      2}]

    return microbe, educts, products

def Hydro(model_parameters):
    # PATHWAY: 4 H2 + 1 CO2 -> 1 CH4 + 2 H2O conrad2000selective, Fe3nchel -131kj/mol
    # Pathway beer2007transport: CO2 (aq) + 4H2 (aq)   - > CH4 (aq) + 2H2O (l)
    microbe = {'name':          'M_Hydro',
               'vmax':          model_parameters['Vmax_Hydro'],
               'death_rate':    model_parameters['Sensenmann'],
               'CUE':           model_parameters['CUE_Hydro'],
               'C_source':      'CO2'}

    educts = [{'name':          'H2',
               'stoich':        4,
               'Km':            model_parameters['Km_Hydro_H2']}, #0.001 / CONSTANTS.SOIL_DENSITY},                # 0.01 mikromol pro cm^3 from Song

              {'name':          'CO2',
               'stoich':        1,
               'Km':            model_parameters['Km_Hydro_CO2'], #0.005/CONSTANTS.SOIL_DENSITY  ,               # 0.05 mikromol pro cm^3 from Song
               'C_atoms':       1}
              ]

    products = [{'name':        'CH4',
                 'stoich':      1},
                
                {'name':        'H2O',
                 'stoich':      2}
                ]

    return microbe, educts, products

def Homo(model_parameters):
    # PATHWAY: song,Fe3nchel,conrad2000selective : 4H2 + 2CO2 → CH3COOH+ 2H2O. Alternativ: 2CO2 + 8H2  = CH3COOH +H2O., laut Thauer sind es 4 H2 + Co2 -> CH4 +2H2O
    # Pathway beer2007transport :CO2 (aq) + 2H2 (aq) -> 1/2CH3COO?(aq) + 1/2H+(aq) + H2O(l)
    microbe = {'name':          'M_Homo',
               'vmax':          model_parameters['Vmax_Homo'],
               'death_rate':    model_parameters['Sensenmann'],
               'CUE':           model_parameters['CUE_Homo'] ,
               'C_source':      'CO2'}

    educts = [{'name':          'H2',
               'stoich':        4,
               'Km':            model_parameters['Km_Homo_H2']}, #0.0001 / CONSTANTS.SOIL_DENSITY   },          # 0.01 from Song

               {'name':         'CO2',
                'stoich':       2,
                'Km':           model_parameters['Km_Homo_CO2'], #0.0005 / CONSTANTS.SOIL_DENSITY   ,           # 0.05 from Song, laut (van1999efFe3cts) größer als Hydro, laut schink1997energetics sollte der Wert mit sinkender Temp, mit zunehmendem Acetate und sinkendem PH sinken (Im vlg zu Hydro)
                'C_atoms':      1}
               ]

    products = [{'name':        'Acetate',
                 'stoich':       1},
                
                 {'name':        'H2O',
                 'stoich':      2}
                ]

    return microbe, educts, products

def Ac(model_parameters):
     #PATHWAY:  CH3COO + H+ ->  CH4 + CO2 Fe3y_Conrad2000
    #CH3COO (aq) + H+ (aq) ->CO2 (aq) + CH4 (aq), beer2007transport
    microbe = {'name':          'M_Ac',
               'vmax':          model_parameters['Vmax_Ac'],
               'death_rate':    model_parameters['Sensenmann'],
               'CUE':           model_parameters['CUE_Ac'],
               'C_source':      'Acetate'}

    educts = [{'name':          'Acetate',
               'stoich':        1,
               'Km':            model_parameters['Km_Ac_Acetate'] ,#0.2 / CONSTANTS.SOIL_DENSITY   ,      #0.05 / SOIL_DENSITY   # 0.05 from song, Wert sollte 10 mal höher sein als bei AltE laut roden2003 bei 12Mikromol !!!!
               'C_atoms':       2}
              ]

    products = [{'name':        'CH4',
                 'stoich':      1},

                {'name':        'CO2',
                 'stoich':      1}
                ]

    return microbe, educts, products
