
from ORDER import POOL_ORDER, pool_index
import CONSTANTS
import chemistry
import numpy as np

DEBUG = False

def builder(defined_pathways):

    if DEBUG:
        print('pathway parameters:')
        [pathway_formatter(microbe, educts, products) for microbe, educts, products in defined_pathways]

    built_pathways = [pathway_builder(microbe, educts, products) for microbe, educts, products in defined_pathways]

    def right_hand_side(t, system_state):
        """
        This is the function given to the IVP solver.
        """
        pathway_changes = [pathway(t, system_state) for pathway in built_pathways]
        changes = np.sum(np.stack(pathway_changes, axis = 0), axis = 0)
        changes = np.clip(changes, -system_state, np.inf)
        return changes

    return right_hand_side

def pathway_builder(microbe, educts, products):
    """
    This wrapper around the pathway(...) function precomputes all vectors that
    are constant during time stepping with the IVP solver.
    It is only called once per pathway, by the builder.
    It provides the necessary vectors to the pathway(...) function.
    """

    # initialize the vectors needed at each time step
    pathway_vector = np.zeros((len(POOL_ORDER),))
    Km_vector = np.zeros((len(POOL_ORDER),))
    inhibition_vector = np.ones((len(POOL_ORDER),))*np.inf
    henrys_law_vector = np.ones((len(POOL_ORDER),))

    microbe_index = pool_index(microbe['name'])

    # TODO: so far, MM is only possible for educts, and inhibition is only possible by products

    for educt in educts:
        educt_index = pool_index(educt['name'])
        pathway_vector[educt_index] = -educt['stoich']
        henrys_law_vector[educt_index] = chemistry.henrys_law(educt['name'])

        if 'Km' in educt:
            Km_vector[educt_index] = educt['Km']

        if 'C_source' in educt:
            c_atoms = educt['C_atoms']
            CUE = microbe['CUE']
            pathway_vector[microbe_index] = educt['stoich']*CUE/(1-CUE)*c_atoms*CONSTANTS.MOLAR_MASS_C
            pathway_vector[educt_index] = -stoich*1/(1-CUE)

    for product in products:
        product_index = pool_index(product['name'])
        pathway_vector[product_index] = product['stoich']
        henrys_law_vector[product_index] = chemistry.henrys_law(product['name'])
        if 'inhibition' in product:
            inhibition_vector[product_index] = product['inhibition']

    # for inverse MM
    if 'Kmb' in microbe:
        Km_vector[microbe_index] = microbe['Kmb']

    v_max = microbe['vmax']

    # handle microbe death
    growth_rate_vector = np.zeros((len(POOL_ORDER),))
    growth_rate_vector[microbe_index] -= microbe['death_rate']

    if DEBUG:
        print_matrix = np.concatenate([np.reshape(henrys_law_vector, (-1, 1)),
                                       np.reshape(Km_vector, (-1, 1)),
                                       np.reshape(inhibition_vector, (-1, 1)),
                                       np.reshape(growth_rate_vector, (-1,1))], axis = 1)

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
        MM = np.where(dissolved_system_state == 0, 0,
                      dissolved_system_state/(Km_vector + dissolved_system_state))
        MM = np.where(Km_vector == 0, 1, MM)
        # compute the total MM factor
        total_MM_factor = np.prod(MM)

        # compute the inverse Michaelis-Menten factors
        invMM = np.where(dissolved_system_state == 0, 1,
                         1 - dissolved_system_state/(inhibition_vector + dissolved_system_state))
        invMM = np.where(inhibition_vector == np.inf, 1, invMM)
        # compute the total inhibition factor
        total_inibition_factor = np.prod(invMM)

        # compute the actual reaction rate
        v = v_max * total_MM_factor * total_inibition_factor

        # compute the changes for all pools given the current biomass and reaction rate
        biomass = system_state[microbe_index]
        system_state_changes = biomass * v * pathway_vector + growth_rate_vector

        if DEBUG:
            print_matrix = np.concatenate([np.reshape(dissolved_system_state, (-1, 1)),
                                           np.reshape(MM_factors, (-1, 1)),
                                           np.reshape(inhibition_factors, (-1, 1)),
                                           np.reshape(system_state_changes, (-1,1))], axis = 1)

            print('calling: ', microbe['name'])
            print('==========' + '='*len(microbe['name']))
            print_array(print_matrix,columns = ['diss', 'MM', 'inh', 'changes'])

            print('total MM' , total_MM_factor)
            print('total inh', total_inibition_factor)

        # TODO: required for solver stability?
        # system_state_changes[system_state_changes < 1e-30] = 0

        return system_state_changes

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
