

POOL_ORDER = ['C',
            'DOC',
            'CH4',
            'CO2',
            'Fe3',
            'M_Ac',
            'M_Ferm',
            'M_Fe3',
            'M_Hydro',
            'M_Homo',
            'Acetate',
            'H2',
            'Fe2',
            'pH',
            'H2O',
            #'HCO3',
            'weight',
            'water'
            # 'CO2_Ac',
            # 'CO2_Hydro',
            # 'CH4_Hydro',
            # 'CO2_Ferm',
            # 'CO2_Fe3',
            # 'CO2_Homo',
            # 'H2_Homo',
            # 'H2_Hydro',
            # 'DGr_Fe3',
            # 'DGr_Hydro',
            # 'DGr_Homo',
            # 'DGr_Ac'
            ]


def pool_index(pool_name):
    return POOL_ORDER.index(pool_name)
