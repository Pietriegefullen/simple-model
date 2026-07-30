import numpy as np

import system
import chemistry
import CONSTANTS

HENRYS_LAW = system.henrys_law()

def pathway_by_name(name):
    pwys = {'Homo': Homo,
            'Aceto': Aceto,
            'Hydrolysis': Hydrolysis,
            'Fe3': Fe3,
            'Hydro': Hydro,
            'Fermentation': Fermentation}
    return pwys[name]

class Pathway():
    def __init__(self, microbe, educts, products):        
        self.educts = educts
        self.products = products
        self.microbe = microbe
        
        self.Km = np.sum(np.stack([system.vector(0, educt, educt['Km'])
                                   for educt in educts], axis = -1), axis = -1)
        self.Km += system.vector(0, microbe, microbe['Kmb'])
        
        stoich_vector = np.sum(np.stack([system.vector(0, subst, -subst['stoichiometry'])
                                         for subst in educts], axis = -1), axis = -1)
        stoich_vector += np.sum(np.stack([system.vector(0, subst, subst['stoichiometry'])
                                         for subst in products], axis = -1), axis = -1)
        self.stoichiometry = stoich_vector
        
        self.inhibition = system.vector(np.inf)
        for product in products:
            self.inhibition[system.index(product)] = product['inhibition']
            
        self.v_max = microbe['v_max']
        self.death_rate = system.vector(0, microbe, microbe['death_rate'])
        
        C_source = microbe['C_source']
        self.anabolism = system.vector(0)
        if not C_source is None:
            C_atoms = chemistry.C_atoms(C_source)
            C_source_stoich = [educt for educt in educts 
                               if str(educt) == C_source][0]['stoichiometry']
            CUE = microbe['CUE']
            anabolism_fraction = C_source_stoich*CUE/(1.-CUE + 1e-7)
            microbe_growth = system.vector(0, microbe, anabolism_fraction*C_atoms*CONSTANTS.MOLAR_MASS_C)
            C_source_reduction = system.vector(0, microbe['C_source'], anabolism_fraction)
            self.anabolism = microbe_growth - C_source_reduction
            
        self.pathway_vector = self.stoichiometry + self.anabolism
        
        self.microbe_index = system.index(microbe)
        self.use_thermodynamics = microbe['use_thermodynamics']
        
        if self.use_thermodynamics:
            self.deltaG_f = np.sum(np.stack(
                    [system.vector(0, str(subst), chemistry.GIBBS_FORMATION[str(subst)])
                     for subst in (educts + products)],
                                             axis = -1), axis = -1)
            self.deltaG_s = np.sum(self.stoichiometry*self.deltaG_f)
        
        self.state_logger = None
       
    def system(self):
        syst = [self.microbe.name] + [e.name for e in self.educts] + [p.name for p in self.products]
        return list(set(syst))
        
    def inject_logger(self, state_logger):
        self.state_logger = state_logger
    
    def log(self, name, t, value):
        if not self.state_logger is None:
            self.state_logger.log_snap(self.__class__.__name__ + '_' + name, t, value)
 
     
 
    def thermodynamics(self, t, S): # hier berechnet sich der thermodynamische faktor
        dissolved_S = HENRYS_LAW*S
        thermodynamic_factor = 1.
        if self.use_thermodynamics:
            R = CONSTANTS.GAS_CONSTANT
            T = 4. + CONSTANTS.KELVIN
            log_Q = system.vector(0) # generates a zero-filled vector of system shape
            
            concentrations = system.vector(0)
            contributes = self.stoichiometry != 0
            contributes[system.index('H2O')] = False # water never contributes in thermodynamics
            denom = np.sum(dissolved_S[contributes])
        
            if denom <= 0:
                concentrations[contributes] = 0
            else:
                concentrations[contributes] = dissolved_S[contributes]/denom
            
            contributing_educts = np.logical_and(self.stoichiometry < 0, contributes)
            contributing_products = np.logical_and(self.stoichiometry > 0, contributes)
            educt_concentrations = concentrations[contributing_educts]
            prod_concentrations = concentrations[contributing_products]
            
            if np.any(educt_concentrations == 0):
                deltaG_r = np.inf
                thermodynamic_factor = 0.
                
            elif np.any(prod_concentrations == 0):
                deltaG_r = -np.inf
                thermodynamic_factor = 1.
                
            else:
                log_Q = np.sum(self.stoichiometry[contributes]*np.log(concentrations[contributes]))
                deltaG_r = self.deltaG_s + R*T*log_Q
                deltaG_rmin = chemistry.GIBBS_MINIMUM
                thermodynamic_factor = 1 - np.exp(np.minimum(0.,deltaG_r - deltaG_rmin)/(R*T))
            
            self.log('deltaG_r', t, deltaG_r)
        self.log('thermodynamic_factor', t, thermodynamic_factor)
        return thermodynamic_factor
    
    def __call__(self, t, S): # hier rechnen wir die MM Faktoren
        biomass = S[self.microbe_index]
        biomass = np.clip(biomass, 1e-8, np.inf)
            
        dissolved_S = HENRYS_LAW*S

        # setting eps > 0 where dissolved_S == 0 has no effect
        # because MM will be 0 anyway.
        # only to suppress warnings 
        eps = np.where(dissolved_S == 0, 1e-8, 0) 

        MM = np.where((self.Km + dissolved_S) == 0, 
                      1,
                      np.where(dissolved_S == 0,
                               0,
                               dissolved_S/(self.Km + dissolved_S + eps)))
        
        inhib = 1 - np.where(self.inhibition + dissolved_S == 0, 0, 
                             dissolved_S/(self.inhibition + dissolved_S))
        thermodynamic_factor = self.thermodynamics(t, S)
        
        
        MM_factor = np.prod(MM)
        inhib_factor = np.prod(inhib)
        v = self.v_max * MM_factor * inhib_factor * thermodynamic_factor # hier rechnen wir die rate aus

        dS_dt = biomass * v * self.pathway_vector - biomass * self.death_rate
        dS_dt = np.clip(dS_dt, -S, np.inf)
       
        self.log('MM', t, MM_factor)
        self.log('inhib', t, inhib_factor)
        self.log('v', t, v)
        
        if 'CH4' in self.products:
            ch4_prod = dS_dt[system.index('CH4')]
            self.log('CH4 from ' + self.microbe.name,t,ch4_prod)
            
        
        return np.reshape(dS_dt, (-1,))

    def __str__(self): # schreibt nur den pathway auf, tut aber nix
        educts = ' + '.join([str(s.stoichiometry) + ' ' + str(s) for s in self.educts])
        products = ' + '.join([str(s.stoichiometry) + ' ' + str(s) for s in self.products])
        pwy_string = f'{self.__class__.__name__: <12s}: {educts} -> {products}'    
        return pwy_string



class Microbe():
    def __init__(self, name, v_max, 
                 death_rate = 0,
                 Kmb = 0, 
                 CUE = 0, 
                 C_source = None, 
                 use_thermodynamics = True):
        
        if not CUE == 0:
            assert C_source is not None
            
        self.name = name
        self.death_rate = death_rate
        self.v_max = v_max
        self.Kmb = Kmb
        self.CUE = CUE
        self.C_source = C_source
        self.use_thermodynamics = use_thermodynamics
        
    def __getitem__(self, key):
        return getattr(self,key)
        
    def __str__(self):
        return self.name
    
    def get_config(self):
        return {'name': self.name,
                'death_rate': self.death_rate,
                'v_max':self.v_max,
                'Kmb': self.Kmb,
                'CUE': self.CUE,
                'C_source': self.C_source,
                'use_thermodynamics': self.use_thermodynamics}

class Substance():
    def __init__(self, stoichiometry, name, 
                 Km = None, 
                 inhibition = None):
        assert not (Km is not None and inhibition is not None)
        
        self.name = name
        self.stoichiometry = stoichiometry
        self.Km = Km if not Km is None else 0
        self.inhibition = np.inf if inhibition is None else inhibition
    
    def __getitem__(self, key):
        return getattr(self,key)
    
    def __str__(self):
        return self.name
    
    def __eq__(self, other):
        if isinstance(other, str):
            return self.name == other
        raise NotImplementedError()

    def get_config(self):
        return {'name': self.name,
                'stoichiometry': self.stoichiometry,
                'Km': self.Km,
                'inhibition': self.inhibition}


class Hydrolysis(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'TOC')]
        products = [Substance(1, 'DOC')]
        microbe = Microbe(name = 'M_Ferm',
                          v_max = model_parameters['Hydrolysis_v_max'],
                          Kmb = model_parameters['Hydrolysis_Kmb'],
                          use_thermodynamics = False)
        super().__init__(microbe, educts, products)
    
class Fermentation(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(10, 'DOC', 
                            Km = model_parameters['Ferm_Km'])]
        products = [Substance(3.5, 'Acetate', 
                              inhibition = model_parameters['Ferm_inhibition']),
                    Substance(3, 'CO2'),
                    Substance(6, 'H2')]
        microbe = Microbe(name = 'M_Ferm',
                          v_max = model_parameters['Ferm_v_max'],
                          CUE = model_parameters['Ferm_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'DOC',
                          use_thermodynamics = False)
        super().__init__(microbe, educts, products)

class Hydro(Pathway):
     def __init__(self, model_parameters):
        educts = [Substance(4, 'H2', 
                            Km = model_parameters['Hydro_Km_H2']),
                  Substance(1, 'CO2', 
                            Km = model_parameters['Hydro_Km_CO2'])]
        products = [Substance(1, 'CH4'),
                    Substance(2, 'H2O')]
        microbe = Microbe(name = 'M_Hydro',
                          v_max = model_parameters['Hydro_v_max'],
                          CUE = model_parameters['Hydro_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'CO2')
        super().__init__(microbe, educts, products)

class Homo(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(4, 'H2', 
                            Km = model_parameters['Homo_Km_H2']),
                  Substance(2, 'CO2', 
                            Km = model_parameters['Homo_Km_CO2'])]
        products = [Substance(1, 'Acetate'),
                    Substance(2, 'H2O')]
        microbe = Microbe(name = 'M_Homo',
                          v_max = model_parameters['Homo_v_max'],
                          CUE = model_parameters['Homo_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'CO2')
        super().__init__(microbe, educts, products)

class Aceto(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'Acetate', 
                            Km = model_parameters['Aceto_Km_Ac'])]
        products = [Substance(1, 'CH4'),
                    Substance(1, 'CO2')]
        microbe = Microbe(name = 'M_Ac',
                          v_max = model_parameters['Ac_v_max'],
                          CUE = model_parameters['Ac_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'Acetate')
        super().__init__(microbe, educts, products)

class Fe3(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'Acetate', 
                            Km = model_parameters['Fe3_Km_Ac']),
                  Substance(1, 'H2O'),
                  Substance(8, 'Fe3', 
                            Km = model_parameters['Fe3_Km_Fe3'])]
        products = [Substance(8, 'Fe2'),
                    Substance(2, 'CO2')]
        microbe = Microbe(name = 'M_Fe3',
                          v_max = model_parameters['Fe3_v_max'],
                          CUE = model_parameters['Fe3_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'Acetate')
        super().__init__(microbe, educts, products)
