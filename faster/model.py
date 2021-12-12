
from ORDER import POOL_ORDER, pool_index
import CONSTANTS
import chemistry
import numpy as np

def builder(model_parameters, all_pathways):

    built_pathways = [pathway_builder(model_parameters, *p) for p in all_pathways]

    def right_hand_side(t, system_state):
        pathway_changes = [pathway(t, system_state) for pathway in built_pathways]
        changes = np.sum(pathway_changes, axis = -1)
        changes = np.clip(changes, -system_state, np.inf)

        return changes

    return right_hand_side

def pathway_builder(model_parameters, microbe, educts, products):

    # rows correspond to products, columns correspond to educts
    pathway_vector = np.zeros((len(POOL_ORDER),))
    Km_vector = np.zeros(len(POOL_ORDER))
    henrys_law_vector = np.ones((len(POOL_ORDER),))

    microbe_index = pool_index(microbe['name'])

    for product in products:
        product_index = pool_index(product['name'])
        pathway_vector[product_index] = product['stoich']
        henrys_law_vector[product_index] = chemistry.henrys_law(product['name'])

    for educt in educts:
        educt_index = pool_index(educt['name'])
        pathway_vector[educt_index] = -educt['stoich']
        if 'Km' in educt:
            Km_vector[educt_index] = educt['Km']

        if 'C_source' in educt:
            c_atoms = educt['C_atoms']
            CUE = microbe['CUE']
            pathway_vector[microbe_index] = educt['stoich']*CUE/(1-CUE)*c_atoms*CONSTANTS.MOLAR_MASS_C
            pathway_vector[educt_index] = -stoich*1/(1-CUE)

    # for inverse MM
    if 'Kmb' in microbe:
        Km_vector[microbe_index] = educt['Km']

    v_max = microbe['vmax']

    # handle microbe death
    growth_rate_vector = np.zeros((len(POOL_ORDER),))
    growth_rate_vector[microbe_index] -= microbe['death_rate']

    def pathway(t, system_state):

        dissolved_system_state = henrys_law_vector * system_state
        MM_factors = dissolved_system_state/(Km_vector + dissolved_system_state)
        MM_factors[np.isnan(MM_factors)] = 0
        factor = np.prod(MM_factors)

        v = v_max * factor

        biomass = system_state[microbe_index]
        system_state_changes = biomass * v * pathway_vector + growth_rate_vector*dissolved_system_state

        # TODO: required for solver stability?
        system_state_changes[system_state_changes < 1e-30] = 0

        return system_state_changes

    return pathway
