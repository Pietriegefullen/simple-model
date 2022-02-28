
MEASURE_DAYS_WEIGHTING = True

ALGORITHM = 'PSO' #'differential_evolution'#'gradient' #
WORKERS = -1 # -1 uses all that are available on this machine

PLOT_LIVE_FIT = False

# das ist was ich optimieren will
CHANGEABLES = ['Vmax_help_Ferm',
                'Vmax_Ferm',
                'Vmax_Fe3',
                'Vmax_Homo',
                'Vmax_Hydro',
                'Vmax_Ac',
                'Kmb_help_Ferm',
                'Inhibition_Ferm',
                'Fe3',
                'M_Ac',
                'CUE_Ferm',
                'CUE_Fe3',
                'CUE_Ac',
                'CUE_Homo',
                'CUE_Hydro'
                ]

DIFF_EVOL_PARAMETERS = {'strategy':         'best1bin',
                        'updating':         'immediate'}

GRADIENT_PARAMETERS = {'method':            'L-BFGS-B',
                       'options':           {'maxiter':     2000,
                                             'disp':        True}
                       }

PSO_PARAMETERS = {'options':{'c1': 0.5, 'c2': 0.3, 'w':0.9},
                  'particles':100,
                  'iterations':5000
                  }

def get_initial_guesses():
    # specify initial guesses and bounds for the parameters to be optimized
    initial_guess_dict = {}        #   init    lower upper  ok guesses to start with
    initial_guess_dict['Vmax_help_Ferm'] =  (0.0210966,   0.006,     0.1)  # 0.05
    initial_guess_dict['Vmax_Ferm'] =       (0.43461,   0.1,      0.5)  # 0.011       # Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
    initial_guess_dict['Vmax_Fe3'] =        (0.797,   0.02,      3)  # 0.8         # Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    initial_guess_dict['Vmax_Homo'] =       (0.0466571,    0.05,      1.)   # 0.869       # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    initial_guess_dict['Vmax_Hydro'] =      (0.370321,   0.03,      1.)   # 0.182 1.8   # 0.15 mikromol pro cm^3 from Song
    initial_guess_dict['Vmax_Ac'] =         (0.240890,   0.05,      1.0)  # 0.99           # Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
    initial_guess_dict['Sensenmann'] =      (8.33e-4,    0.0000001,      0.001)# (8.33e-5, 0, 8.44e-5)# 0
    initial_guess_dict['Kmb_help_Ferm'] =   (2.6174,   0.005,    20)      # 10
    initial_guess_dict['Fe3'] =             (25.09,      10,    30)    # 15.587,
    initial_guess_dict['Acetate'] =         (0.1, 0, 10 )
    initial_guess_dict['M_Ac'] =            (0.002,   1.3e-08,   0.02) # 0.002
    initial_guess_dict['Inhibition_Ferm']=  (9.4,   0.001,    20) # Je niedriger desto hemmung # Diese Boundaries müssen anhander Acetatekurven angepasst werden
    initial_guess_dict['CUE_Ferm']=  (.5,   0.,    0.1)
    initial_guess_dict['CUE_Fe3']=  (.5,   0.,    0.1)
    initial_guess_dict['CUE_Ac']=  (.5,   0.,    0.1)
    initial_guess_dict['CUE_Homo']=  (.5,   0.,    0.1)
    initial_guess_dict['CUE_Hydro']=  (.5,   0.,    0.1)
    return(initial_guess_dict)


  