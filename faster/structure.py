import numpy as np

import system

HENRYS_LAW = system.henrys_law()

# TODO: register model parameters as global variables/singleton?
#       use factory method to put parameters in place? -> module acts as singleton?
#       https://python-patterns.guide/gang-of-four/singleton/

# TODO: check microbe, substance names when building model

# TODO: where to get the initial system state from?
#       C, DOC depend on replica.
#       pH depends on sample

class ObjectiveFunction():
    def __init__(self, replicas):
        pass
        # 1) build a model for each replica? or provide model?
        #    => extract changeable parameters????
        #    => which parameters are specific to and different for each replica?
        #    => sufficient to say that initial C pool (and DOC) is different
        #    => hence, run a separate model for both replicas!!!
        #       but, since the replica-specific parameters are computed (all?),
        #       it should not be a different set of changeable parameters for each!
        #    TODO: check whether only the initial system state differs!!!
        #       should this be the case, the model need be built only once.
    
    def __call__(self, changeable_parameters):
        # 1) setup model with provided values for changeable model parameters
        #   try to avoid another model setup, better to just set the parameter
        #   values for the existing model
        # 2) integrate model (get initial system state, from model?)
        # 3) compute losses for each fitted replica, CO2 and CH4
        # 4) sum the losses for replicas
        return loss

class ModelParameter():
    def __init__(self):
        self._changeable = True


class Model():
    def __init__(self):
        self.contributing_pathways = []
        self.initial_system_state = None
        
    def __call__(self, t, Si):
        dSij_dt = np.stack([pathway_j(t, Si)
                            for pathway_j in self.contributing_pathways], axis = -1)
        dSi_dt = np.sum(dSij_dt, axis = -1)
        dSi_dt = np.clip(dSi_dt, -Si, np.inf) # don't let pools become negative
        return dSi_dt 
    
class Pathway():
    def __init__(self, microbe, educts, products):
        self.Km = np.sum(np.stack([system.vector(educt, educt['Km'])
                                   for educt in educts], axis = -1), axis = -1)
        self.stoichiometry = system.pathway_stoichiometry(educts, products)
        self.inhibition = np.sum(np.stack([system.vector(product,product['Km'])
                                   for product in products], axis = -1), axis = -1)
        self.v_max = microbe['v_max']
        self.death_rate = system.vector(microbe)*microbe['death_rate']
        
        CUE = microbe['CUE']
        C_atoms = microbe['C_source']['C_atoms']
        growth = system.vector(microbe, CUE/(1-CUE)*C_atoms*CONSTANTS.MOLAR_MASS_C)
        self.anabolism = growth - system.vector(microbe['C_source'],CUE/(1-CUE))
        self.microbe_index = system.index(microbe)
        self.use_thermodynamics = microbe['thermodynamics']
        
        
    def thermodynamics(self, Si):
        if not self.use_thermodynamics:
            return 1.
        R = CONSTANTS.GAS_CONSTANT
        log_Q = np.zeros_like(self.stoichiometry)
        
        contributes = np.logical_and(self.stoichiometry != 0, Si > 0)
        log_Q[contributes] = np.log(1e-6*Si[contributes])

        deltaG_r = deltaG_s + R*T*np.sum(self.stoichiometry*log_Q)
        deltaG_rmin = chemistry.GIBBS_MINIMUM
        
        return 1 - np.exp(np.minimum(0.,deltaG_r - deltaG_rmin)/(R*T))
    
    def __call__(self, t, Si):
        
        biomass = Si[self.microbe_index]
        
        pathway_vector = self.stoichiometry + self.anabolism
    
        dissolved_Si = HENRYS_LAW*Si
        
        MM = dissolved_Si/(self.Km + dissolved_Si)
        inhib = 1 - dissolved_S/(self.inhibition + dissolved_Si)
        thermodynamic_factor = self.thermodynamics()
        
        v = self.v_max * np.prod(MM) * np.prod(inhib) * thermodynamic_factor
        
        dSi_dt = biomass * v * pathway_vector - biomass * self.death_rate
        dSi_dt = np.clip(dSi_dt, -Si, np.inf)
        return np.reshape(dSi_dt, (-1,))

    
class Microbe():
    def __init__(self, name):
        self._name = name
        
    def __str__(self):
        return self._name

class Substance():
    def __init__(self, name):
        self._name = name
    
    def __str__(self):
        return self._name

class Hydrolysis(Pathway):
    def __init__(self):
        educts = [Substance('C', ...)]
        products = [Substance('DOC')]
        microbe = Microbe()
        super().__init__(microbe, educts, products)
    
class Fermentation(Pathway):
    pass 

    