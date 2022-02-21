
MEASURE_DAYS_WEIGHTING = True

ALGORITHM = 'PSO' #'differential_evolution'#'gradient' #
WORKERS = 64 # -1 uses all that are available on this machine

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
                'M_Ac'
                ]

DIFF_EVOL_PARAMETERS = {'strategy':         'best1bin',
                        'updating':         'immediate'}

GRADIENT_PARAMETERS = {'method':            'L-BFGS-B',
                       'options':           {'maxiter':     2000,
                                             'disp':        True}
                       }

PSO_PARAMETERS = {'options':{'c1': 0.5, 'c2': 0.3, 'w':0.9},
                  'particles':10,#0,
                  'iterations':100#500#0
                  }

def get_initial_guesses():
    # specify initial guesses and bounds for the parameters to be optimized
    initial_guess_dict = {}        #   init    lower upper  ok guesses to start with
    initial_guess_dict['Vmax_help_Ferm'] =  (0.0210966,   0.006,     0.1)  # 0.05
    initial_guess_dict['Vmax_Ferm'] =       (0.143461,   0.1,      0.5)  # 0.011       # Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
    initial_guess_dict['Vmax_Fe3'] =        (1.527,   0.02,      3)  # 0.8         # Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    initial_guess_dict['Vmax_Homo'] =       (0.1466571,    0.05,      1.)   # 0.869       # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    initial_guess_dict['Vmax_Hydro'] =      (0.70321,   0.03,      1.)   # 0.182 1.8   # 0.15 mikromol pro cm^3 from Song
    initial_guess_dict['Vmax_Ac'] =         (0.240890,   0.05,      1.0)  # 0.99           # Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
    initial_guess_dict['Sensenmann'] =      (8.33e-4,    0.0000001,      1.5)# (8.33e-5, 0, 8.44e-5)# 0
    initial_guess_dict['Kmb_help_Ferm'] =   (7.6174,   0.005,    20)      # 10
    initial_guess_dict['Fe3'] =             (15.09,      10,    30)    # 15.587,
    initial_guess_dict['M_Ac'] =            (0.002,   1.3e-08,   0.02) # 0.002
    initial_guess_dict['Inhibition_Ferm']=  (9.4,   0.001,    20) # Je niedriger desto hemmung # Diese Boundaries müssen anhander Acetatekurven angepasst werden

    return(initial_guess_dict)


  