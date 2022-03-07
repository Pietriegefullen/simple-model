import CONSTANTS

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
                "Km_help_Ferm",
                "Km_Ac_Acetate",
                "Km_Homo_CO2",
                "Km_Homo_H2",
                "Km_Hydro_CO2",
                "Km_Hydro_H2",
                "Km_Fe3_Fe3",
                "Km_Fe3_Acetate",
                "Km_Ferm",
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
    initial_guess_dict['Vmax_help_Ferm'] =  (0.347557284254151,   0.0006,     1.)  # 0.05
    initial_guess_dict['Vmax_Ferm'] =       (2.1826192982834263,   0.001,      5.)  # 0.011       # Vmax = 0.5e6 / SOIL_DENSITY # 0.5 from Song
    initial_guess_dict['Vmax_Fe3'] =        (0.0053274970925869525,   0.002,      3.)  # 0.8         # Vprod_max = 0.3* 10**6/ SOIL_DENSITY    # geschätzt
    initial_guess_dict['Vmax_Homo'] =       (0.3842183674485621,    0.005,      1.)   # 0.869       # 0.15 from Song, Laut Ye13 3 bis 6 mal schneller als Hydro
    initial_guess_dict['Vmax_Hydro'] =      (0.6787923256204944,   0.003,      1.)   # 0.182 1.8   # 0.15 mikromol pro cm^3 from Song
    initial_guess_dict['Vmax_Ac'] =         (0.6185701320273594,  0.005,      1.0)  # 0.99           # Vprod_max_Ac = 0.5/ SOIL_DENSITY # 0.5 from song
    initial_guess_dict['Sensenmann'] =      (8.33e-5,    0.0000001,      0.001)# (8.33e-5, 0, 8.44e-5)# 0
    initial_guess_dict['Kmb_help_Ferm'] =   (1322.0602068699052,  0.0005,    2000)      # 10
    initial_guess_dict["Km_help_Ferm"] =    (5188.51170010838/CONSTANTS.SOIL_DENSITY,0.0005, 10000)
    initial_guess_dict["Km_Ac_Acetate"] =    (7.030348185396872/ CONSTANTS.SOIL_DENSITY, 0.0005, 1000)
    initial_guess_dict["Km_Homo_CO2"] =     (848.9503500185294 / CONSTANTS.SOIL_DENSITY,0.00005, 1000)
    initial_guess_dict["Km_Homo_H2"] =      (584.5061484463681 / CONSTANTS.SOIL_DENSITY ,0.00005, 1000)
    initial_guess_dict["Km_Hydro_CO2"] =    (869.0991529683715/CONSTANTS.SOIL_DENSITY ,0.0005, 1000)
    initial_guess_dict["Km_Hydro_H2"] =     (489.1731054166966 / CONSTANTS.SOIL_DENSITY,0.0005, 1000)
    initial_guess_dict["Km_Fe3_Fe3"] =      (533.3685302693583/CONSTANTS.SOIL_DENSITY,0.0005, 1000)
    initial_guess_dict["Km_Fe3_Acetate"] =  (818.2462137022104/CONSTANTS.SOIL_DENSITY,0.00005, 1000)
    initial_guess_dict["Km_Ferm"] =         (630.4916492204161 / CONSTANTS.SOIL_DENSITY,0.0005, 1000)  
    initial_guess_dict['Fe3'] =             (117.3757437908092,      0,    300)    # 15.587,
    initial_guess_dict['Acetate'] =         (1, 0, 100)
    initial_guess_dict['M_Ac'] =            (0.010109432913795091,   1.3e-08,   0.02) # 0.002
    initial_guess_dict['Inhibition_Ferm']=  (10.118875996129882,   0.0001,    20) # Je niedriger desto hemmung # Diese Boundaries müssen anhander Acetatekurven angepasst werden
    initial_guess_dict['CUE_Ferm']=         (0.3698723525408035,   0.,    1)
    initial_guess_dict['CUE_Fe3']=          (0.5572385007984572,   0.,    1)
    initial_guess_dict['CUE_Ac']=           (0.08453242432990166,   0.,    1)
    initial_guess_dict['CUE_Homo']=         (0.3602422303950234,   0.,    1)
    initial_guess_dict['CUE_Hydro']=        (0.7433100141164364,   0.,    1)
    return(initial_guess_dict)


  