
from ORDER import POOL_ORDER, pool_index
import CONSTANTS
import chemistry
import numpy as np

DEBUG = False

def builder(defined_pathways, environment, extended_output = None):

    if DEBUG:
        print('pathway parameters:')
        [pathway_formatter(microbe, educts, products) 
         for microbe, educts, products in defined_pathways]

    built_pathways = [pathway_builder(microbe, educts, products, environment, extended_output) 
                      for microbe, educts, products in defined_pathways]

    def right_hand_side(t, system_state):
        """
        This is the function given to the IVP solver.
        """
        pathway_changes = [pathway(t, system_state) for pathway in built_pathways]
        changes = np.sum(np.stack(pathway_changes, axis = 0), axis = 0)
        changes = np.clip(changes, -system_state, np.inf)
        return changes
    
    def extended(t, system_state):
        extended_values = {pathway.__name__ + '_' + name: value
                          for pathway in built_pathways 
                          for name, value in zip(extended_output,pathway(t, system_state))}
        return extended_values
    
    if not extended_output is None:
        return extended
        
    return right_hand_side

def pathway_builder(microbe, educts, products, environment, extended_output = None):
    """
    This wrapper around the pathway(...) function precomputes all vectors that
    are constant during time stepping with the IVP solver.
    It is only called once per pathway, by the builder.
    It provides the necessary vectors to the pathway(...) function.
    """

    # initialize the vectors needed at each time step
    pathway_vector = np.zeros((len(POOL_ORDER),))
    stoich_vector = np.zeros((len(POOL_ORDER),))
    Km_vector = np.zeros((len(POOL_ORDER),))
    inhibition_vector = np.ones((len(POOL_ORDER),))*np.inf
    henrys_law_vector = np.ones((len(POOL_ORDER),))
    deltaG_f = np.zeros((len(POOL_ORDER),))

    microbe_index = pool_index(microbe['name'])
    
    use_thermodynamics = not 'thermodynamics' in microbe or microbe['thermodynamics'] is True

    for educt in educts:
        educt_index = pool_index(educt['name'])
        stoich_vector[educt_index] = -educt['stoich']
        henrys_law_vector[educt_index] = chemistry.henrys_law(educt['name'])
        if use_thermodynamics:
            deltaG_f[educt_index] = chemistry.GIBBS_FORMATION[educt['name']]

        if 'Km' in educt:
            Km_vector[educt_index] = educt['Km']

        if 'C_source' in microbe and 'C_atoms' in educt:
            c_atoms = educt['C_atoms']
            CUE = microbe['CUE']
            pathway_vector[microbe_index] = educt['stoich']*CUE/(1-CUE)*c_atoms*CONSTANTS.MOLAR_MASS_C
            pathway_vector[educt_index] = -educt['stoich']*CUE/(1-CUE)

    for product in products:
        product_index = pool_index(product['name'])
        stoich_vector[product_index] = product['stoich']
        henrys_law_vector[product_index] = chemistry.henrys_law(product['name'])
        if use_thermodynamics:
            deltaG_f[product_index] = chemistry.GIBBS_FORMATION[product['name']]
        
        if 'inhibition' in product:
            inhibition_vector[product_index] = product['inhibition']

    # pathway vector includes microbes (which should not contribute to thermodynamics!)
    pathway_vector += stoich_vector

    # for inverse MM
    if 'Kmb' in microbe:
        Km_vector[microbe_index] = microbe['Kmb']

    v_max = microbe['vmax']
    
    # for thermodynamics
    deltaG_s = np.sum(stoich_vector*deltaG_f) -600000# TODO: normalized stoichimetry? HOW? die -500000 sind um Hydro und Homo zu zwingen

    # handle microbe death
    death_rate_vector = np.zeros((len(POOL_ORDER),))
    death_rate_vector[microbe_index] = microbe['death_rate'] 

    # environment
    T = environment['temperature'] + CONSTANTS.KELVIN

    if DEBUG:
        print_matrix = np.concatenate([np.reshape(henrys_law_vector, (-1, 1)),
                                       np.reshape(Km_vector, (-1, 1)),
                                       np.reshape(inhibition_vector, (-1, 1)),
                                       np.reshape(death_rate_vector, (-1,1))], axis = 1)

        print('')
        print('building: ', microbe['name'], 'pathway')
        print('===========' + '='*len(microbe['name']) + '========')
        print_array(print_matrix,columns = ['henry', 'Km', 'inhib', 'grow'])

    def pathway(t, system_state):
        """
        This is the actual model.
        """

        # compute the available fraction for all pools
        dissolved_system_state = henrys_law_vector * system_state

        # compute the Michalis-Menten factors given the current pools
        eps = np.where(dissolved_system_state == 0, 1e-8, 0) # no effect, only to suppress warning of invalid value
        MM = np.where(Km_vector + dissolved_system_state == 0, 1,
                      np.where(dissolved_system_state == 0, 0,
                               
                               dissolved_system_state/(Km_vector + dissolved_system_state + eps)))
        
        # compute the total MM factor
        total_MM_factor = np.prod(MM)

        # compute the inverse Michaelis-Menten factors
        invMM = np.where(dissolved_system_state == 0, 1,
                         1 - dissolved_system_state/(inhibition_vector + dissolved_system_state))
        invMM = np.where(inhibition_vector == np.inf, 1, invMM)
        # compute the total inhibition factor
        total_inibition_factor = np.prod(invMM)

        # compute the thermodynamic factor
        thermodynamic_factor = 1.
        if use_thermodynamics:
            R = CONSTANTS.GAS_CONSTANT
            log_Q = np.zeros_like(stoich_vector)
            contributes = np.logical_and(stoich_vector != 0, system_state > 0)
            log_Q[contributes] = np.log(1e-6*system_state[contributes])
            log_Q = stoich_vector*log_Q
            deltaG_r = deltaG_s + R*T*np.sum(log_Q)
            
            deltaG_rmin = chemistry.GIBBS_MINIMUM
            thermodynamic_factor = 1 - np.exp(np.minimum(0.,deltaG_r - deltaG_rmin)/(R*T))

        # compute the actual reaction rate
        v = v_max * total_MM_factor * total_inibition_factor * thermodynamic_factor

        # compute the changes for all pools given the current biomass and reaction rate
        biomass = system_state[microbe_index]
        death_rate_vector_biomass = death_rate_vector * biomass
        system_state_changes = biomass * v * pathway_vector - death_rate_vector_biomass

        if DEBUG:
            print_matrix = np.concatenate([np.reshape(dissolved_system_state, (-1, 1)),
                                           np.reshape(MM, (-1, 1)),
                                           np.reshape(invMM, (-1, 1)),
                                           np.reshape(system_state_changes, (-1,1))], axis = 1)

            print('calling: ', microbe['name'])
            print('==========' + '='*len(microbe['name']))
            print_array(print_matrix,columns = ['diss', 'MM', 'inh', 'changes'])

            print('total MM' , total_MM_factor)
            print('total inh', total_inibition_factor)


        if not extended_output is None:
            extended_dict = {'MM': total_MM_factor,
                             'thermo': thermodynamic_factor,
                             'deltaCO2': system_state_changes[pool_index('CO2')],
                             'deltaCH4': system_state_changes[pool_index('CH4')],
                             'deltaH2': system_state_changes[pool_index('H2')],
                             'v':v,
                             'inhibition': total_inibition_factor}
            if use_thermodynamics:
                extended_dict.update({'deltaGr': deltaG_r,
                                      'deltaGs': deltaG_s,
                                      'logQ':log_Q})
            extended_system_state = np.array([extended_dict[k] if k in extended_dict else np.nan 
                                              for k in extended_output ])
            
            return extended_system_state

        return system_state_changes
    
    pathway.__name__ = microbe['name']
    return pathway


def pathway_formatter(microbe, educts, products):
    print('')
    print('parameters for', microbe['name'], 'pathway')
    print('===============' + '='*len(microbe['name'])+ '=========')
    for k,v in microbe.items():
        if k == 'name': continue
        print(f'   {k:10} {v}')

    print('')
    print('educts:')
    for i,educt in enumerate(educts):
        cnt = f'{i+1:2d})'
        for k,v in educt.items():
            if not k == 'name':
                cnt = '   '
            print(f'{cnt}   {k:10} {v}')

    print('')
    print('products:')
    for i,product in enumerate(products):
        cnt = f'{i+1:2d})'
        for k,v in product.items():
            if not k == 'name':
                cnt = '   '
            print(f'{cnt}   {k:10} {v}')

    print('')

def print_array(arr, title = '', columns = None):
    np.set_printoptions(precision = 2,
                        suppress=True)
    if not len(arr.shape) == 2:
        arr = np.reshape(arr, (-1,1))
    arr_str = np.array2string(arr)
    spl = arr_str.split('\n')

    label_width = 10
    labeled_str = '\n'.join(list([f'{p:<{label_width}} {s}' for p, s in zip(POOL_ORDER, spl)]))

    print('')
    if not title == '':
        print(title)
        print('='*len(title))

    if columns is None:
        columns = list()
    print(' '*label_width + ' ' + '     '.join([f'{c:5}' for c in columns]))

    print(labeled_str)
