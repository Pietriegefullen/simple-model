

import CONSTANTS

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
