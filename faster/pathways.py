

import data
import CONSTANTS
import OPTIMIZATION_PARAMETERS

def default_model_parameters(specimen_index = None, site = 'all'):

    model_parameters =     {'M_Fe3':            0.2,
                            'M_Ferm':           0.02,
                            'M_Hydro':          0.0025,
                            'M_Homo':           0.0001,
                            'M_Ac':             0.0002,

                            'Sensenmann':       0,#8.33e-5,

                            'Vmax_Fe3':         20.27,
                            'Vmax_help_Ferm':   0.030966,
                            'Vmax_Ferm':        0.0943461,
                            'Vmax_Homo':        0.066571,
                            'Vmax_Hydro':       0.070321,
                            'Vmax_Ac':          0.040890,

                            'Kmb_help_Ferm':    11.748221,
                            'Inhibition_Ferm':  9.424312,

                            'Fe3':              15.090783}

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
               }

    educts =  [{'name':         'C',
                'stoich':       1,
                'Km':           10000/CONSTANTS.SOIL_DENSITY} # For Ferm_help, Km must be 0 so that MM_factors evaluates to 1.0, Process is Enzyme Limited not Substrate limited
               ]

    products = [{'name':        'DOC',
                 'stoich':      1}
                ]

    return microbe, educts, products

def Ferm(model_parameters):

    microbe = {'name':          'M_Ferm',
               'vmax':          model_parameters['Vmax_Ferm'],
               'death_rate':    model_parameters['Sensenmann'],
               'microbe'       : "Ferm" ,
               'CUE'           :   0.5,
               'C_source'       : 'DOC'}


    educts =  [{'name':         'DOC',
                'stoich':       6,
                'Km':           100 / CONSTANTS.SOIL_DENSITY,      # 10 from Song  mikromol pro gram
                'C_atoms':      6}    # weil glucose (und andere Monomere) 6 C atome hat und ein momomer aus der spaltung von Coellulose ist
               ]

    #------------------------------------------------------------------------------------------------------------
    #laut Grant A :6 , CO2:3, H2 : 1
    #laut Knoblauch A :6 , CO2:3, H2 : 6
    #laut mir DOC = A: 7 ,CO2 :6 H2:12 (teile durch 2 weil dann näher am Grant verhältnis)
    #laut neumann2016modeling  Glucose = A:2 CO2:2 H2:4  C6H12O6 +2H2O -> 2CH3COOH+2CO2 + 4H2
    products = [{'name':        'Acetate',
                 'stoich':       3.5,
                 'inhibition':   model_parameters['Inhibition_Ferm']}  , # MM faktor für Acetathemmung

                {'name':        'CO2',
                 'stoich':       3}  ,

                {'name':        'H2',
                 'stoich':      6}
                ]
    return microbe, educts, products

def Fe3(model_parameters):
    microbe = {'name':          'M_Fe3',
               'C_source':      'Acetate',
               'CUE':           0.5,
               'death_rate':    model_parameters['Sensenmann'],
               'vmax':          model_parameters['Vmax_Fe3']}

    educts = [{'name':          'Acetate',
               'stoich':        1,
               'Km':            0.01/CONSTANTS.SOIL_DENSITY,
               'C_atoms':       2},

              {'name':          'Fe3',
               'stoich':        8,
               'Km':            2/CONSTANTS.SOIL_DENSITY}]

    products = [{'name':        'Fe2',
                 'stoich':      8},

                {'name':        'CO2',
                 'stoich':      2}]

    return microbe, educts, products

def Hydro(model_parameters):
    microbe = {'name':          'M_Hydro',
               'vmax':          model_parameters['Vmax_Hydro'],
               'death_rate':    model_parameters['Sensenmann'],
               'CUE':           0.5,
               'C_source':      'CO2'}

    educts = [{'name':          'H2',
               'stoich':        4,
               'Km':            0.001 / CONSTANTS.SOIL_DENSITY},                # 0.01 mikromol pro cm^3 from Song

              {'name':          'CO2',
               'stoich':        1,
               'Km':            0.005/CONSTANTS.SOIL_DENSITY  ,               # 0.05 mikromol pro cm^3 from Song
               'C_atoms':       1}
              ]

    products = [{'name':        'CH4',
                 'stoich':      1}
                ]

    return microbe, educts, products

def Homo(model_parameters):
    microbe = {'name':          'M_Homo',
               'vmax':          model_parameters['Vmax_Homo'],
               'death_rate':    model_parameters['Sensenmann'],
               'CUE':           0.5 ,
               'C_source':      'CO2'}

    educts = [{'name':          'H2',
               'stoich':        4,
               'Km':            0.0001 / CONSTANTS.SOIL_DENSITY   },          # 0.01 from Song

               {'name':         'CO2',
                'stoich':       2,
                'Km':           0.0005 / CONSTANTS.SOIL_DENSITY   ,           # 0.05 from Song, laut (van1999efFe3cts) größer als Hydro, laut schink1997energetics sollte der Wert mit sinkender Temp, mit zunehmendem Acetate und sinkendem PH sinken (Im vlg zu Hydro)
                'C_atoms':      1}
               ]

    products = [{'name':        'Acetate',
                 'stoich':       1}
                ]

    return microbe, educts, products

def Ac(model_parameters):
    microbe = {'name':          'M_Ac',
               'vmax':          model_parameters['Vmax_Ac'],
               'death_rate':    model_parameters['Sensenmann'],
               'CUE':           0.5,
               'C_source':      'Acetate'}

    educts = [{'name':          'Acetate',
               'stoich':        1,
               'Km':            0.2 / CONSTANTS.SOIL_DENSITY   ,      #0.05 / SOIL_DENSITY   # 0.05 from song, Wert sollte 10 mal höher sein als bei AltE laut roden2003 bei 12Mikromol !!!!
               'C_atoms':       2}
              ]

    products = [{'name':        'CH4',
                 'stoich':      1},

                {'name':        'CO2',
                 'stoich':      1}
                ]

    return microbe, educts, products
