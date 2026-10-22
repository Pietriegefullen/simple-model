import numpy as np

import system
import chemistry
import CONSTANTS


from parameters import flatten as flat


from CONSTANTS import EPS

def pathway_by_name(name):
    pwys = {'Homo': Homo,
            'Aceto': Aceto,
            'Hydrolysis': Hydrolysis,
            'Fe3': Fe3,
            'Hydro': Hydro,
            'Fermentation': Fermentation}
    return pwys[name]

PATHWAY_COLORS = {
    'Homo': 'tab:pink',
    'Fe3': 'tab:red',
    'Fermentation': 'tab:brown',
    'Hydro': 'tab:cyan',
    'Hydrolysis': 'tab:olive',
    'Aceto': 'tab:purple',
}


def pathway_color(name):
    """Return the project-wide color assigned to a biological pathway."""
    return PATHWAY_COLORS[name]

class Pathway():
    def __init__(self, microbe, educts, products, use_thermodynamics):
        self.educts = educts
        self.products = products
        self.microbe = microbe

        self.use_thermodynamics = use_thermodynamics

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

        self.use_inhib = np.any(np.isfinite(self.inhibition))

        self.v_max = microbe['v_max']
        self.death_rate = system.vector(0, microbe, microbe['death_rate'])

        C_source = microbe['C_source']
        self.anabolism = system.vector(0)
        if not C_source is None:
            C_atoms = chemistry.C_atoms(C_source)
            C_source_stoich = [educt for educt in educts
                               if str(educt) == C_source][0]['stoichiometry']
            CUE = microbe['CUE']
            anabolism_fraction = C_source_stoich*CUE/(1.-CUE + EPS)
            microbe_growth = system.vector(0, microbe, anabolism_fraction*C_atoms*CONSTANTS.MOLAR_MASS_C)
            C_source_reduction = system.vector(0, microbe['C_source'], anabolism_fraction)
            self.anabolism = microbe_growth - C_source_reduction

        self.pathway_vector = self.stoichiometry + self.anabolism

        self.microbe_index = system.index(microbe)

        if float(self.use_thermodynamics) == 1:
            self.deltaG_f = np.sum(np.stack(
                    [system.vector(0, str(subst), chemistry.GIBBS_FORMATION[str(subst)])
                     for subst in (educts + products)],
                                             axis = -1), axis = -1)
            self.deltaG_s = np.sum(self.stoichiometry*self.deltaG_f)

        self.state_logger = None

    def get_config(self, only_structure = False):
        cfg = {'microbe': self.microbe.get_config(only_structure),
               'educts': {e.name: e.get_config(only_structure) for e in sorted(self.educts)},
               'products': {e.name: e.get_config(only_structure) for e in sorted(self.products)},
               'use_thermodynamics': flat(self.use_thermodynamics, only_structure)
               }
        return cfg

    def system(self):
        syst = [self.microbe.name] + [e.name for e in self.educts] + [p.name for p in self.products]
        return list(set(syst))

    def inject_logger(self, state_logger):
        self.state_logger = state_logger

    def log(self, name, t, value):
        if not self.state_logger is None:
            if np.ndim(value) == 0:
                self.state_logger.log_snap(self.__class__.__name__ + '_' + name, t, value)
            else:
                self.state_logger.log(self.__class__.__name__ + '_' + name, t, value)

    def deltaG_r(self, S):
        """Return reaction Gibbs energy for one state or an accepted state series.
        ``S`` is either a single system-state vector or an array with system
        pools along axis zero and time along axis one.  Keeping this calculation
        free of logging allows diagnostics to be evaluated from the accepted
        solver output rather than from LSODA's intermediate trial states.
        """

        R = CONSTANTS.GAS_CONSTANT
        T = 4. + CONSTANTS.KELVIN

        contributes = self.stoichiometry != 0
        contributes[system.index('H2O')] = False # water never contributes in thermodynamics

        # substances are  in micromole per g dry weight of soil.
        # this is not a concentration but a normalized amount.
        # instead, use activities of dissolved species in mol/L
        rho_H2O = CONSTANTS.DENSITY_WATER_4_DEGC # g/L
        m_H2O = CONSTANTS.MOLAR_MASS_H2O # g/mol

        # calculate in micro L/gdw, to match dissolved quantities' units: micro mol/gdw
        water_content = S[system.index('H2O')]*m_H2O/rho_H2O #  micro L/gdw
        activities = S/water_content # mol/L

        stoichiometry = self.stoichiometry[contributes]
        with np.errstate(divide = 'ignore'):
            log_activities = np.log(activities[contributes])


        if log_activities.ndim == 1:
            log_Qj = stoichiometry*log_activities
        else:
            log_Qj = stoichiometry[:, np.newaxis]*log_activities

        # NOTE: numpy returns nan if sum() contains both +inf and -inf.
        # If any +inf is present (any reactant is zero), set log(Q) to inf
        if np.any(log_Qj == np.inf):
            log_Q = np.inf
        else:
            log_Q = np.sum(log_Qj, axis = 0)
    
        return self.deltaG_s + R*T*log_Q

    def thermodynamics(self, t, S): # hier berechnet sich der thermodynamische faktor
        """
        S is the dissolved species vector
        """

        R = CONSTANTS.GAS_CONSTANT
        T = 4. + CONSTANTS.KELVIN
        deltaG_r = self.deltaG_r(S)
        deltaG_rmin = chemistry.GIBBS_MINIMUM
        thermodynamic_factor = 1 - np.exp(np.minimum(0.,deltaG_r - deltaG_rmin)/(R*T))

        self.log('thermodynamic_factor', t, thermodynamic_factor)
        return thermodynamic_factor

    def kinetics(self, S):
        if S.ndim == 1:
            Km = self.Km
            ax = None
        else:
            Km = self.Km[:, np.newaxis]
            ax = 0

        MM = np.where(Km == 0,
                           1.,
                      S/(Km + S + EPS))
        MM_factor = np.prod(MM, axis = ax)
        return MM_factor

    def kinetic_inhibition(self, S):
        if S.ndim == 1:
            inhibition = self.inhibition
            ax = None
        else:
            inhibition = self.inhibition[:, np.newaxis]
            ax = 0

        inhib = 1 - np.where(np.isinf(inhibition), 0,
                             S/(inhibition + S + EPS))
        inhib_factor = np.prod(inhib, axis = ax)
        return inhib_factor

    def compute_v(self, t, S):#ö, f_MM, f_thermo, f_inhib):
        f_MM = self.kinetics(S)

        v = float(self.v_max) * f_MM

        f_inhib = None
        if self.use_inhib:
            f_inhib = self.kinetic_inhibition(S)
            v *= f_inhib
            self.log('inhib', t, f_inhib)

        f_thermo = None
        if float(self.use_thermodynamics) == 1:
            f_thermo = self.thermodynamics(t, S)
            v *= f_thermo
            self.log('thermodynamic_factor', t, f_thermo)

        self.log('MM', t, f_MM)
        self.log('v', t, v)
        return v

    def __call__(self, t, S): # hier rechnen wir die MM Faktoren
        """
        S is the dissolved species vector
        """
        biomass = S[self.microbe_index]

        v = self.compute_v(t, S)
        dS_dt = biomass * v * self.pathway_vector - biomass * self.death_rate

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
                 C_source = None):

        if not CUE == 0:
            assert C_source is not None

        self.name = name
        self.death_rate = death_rate
        self.v_max = v_max
        self.Kmb = Kmb
        self.CUE = CUE
        self.C_source = C_source

    def __getitem__(self, key):
        return getattr(self,key)

    def __str__(self):
        return self.name

    def get_config(self, only_structure = False):
        return {'name': self.name,
                'death_rate': flat(self.death_rate, only_structure),
                'v_max':flat(self.v_max, only_structure),
                'Kmb': flat(self.Kmb, only_structure),
                'CUE': flat(self.CUE, only_structure),
                'C_source': flat(self.C_source, only_structure)}

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

    def get_config(self, only_structure):
        return {'name': self.name,
                'stoichiometry': flat(self.stoichiometry, only_structure),
                'Km': flat(self.Km, only_structure),
                'inhibition': flat(self.inhibition, only_structure)}

    def __lt__(self, other):
        return self.name < other.name



class Hydrolysis(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'TOC')]
        products = [Substance(1, 'DOC')]
        microbe = Microbe(name = 'M_Ferm',
                          v_max = model_parameters['Hydrolysis_v_max'],
                          Kmb = model_parameters['Hydrolysis_Kmb'])
        super().__init__(microbe, educts, products, False)

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
                          C_source = 'DOC')
        super().__init__(microbe, educts, products, False)

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
        super().__init__(microbe, educts, products, model_parameters['Hydro_thermodynamics'])

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
        super().__init__(microbe, educts, products, model_parameters['Homo_thermodynamics'])

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
        super().__init__(microbe, educts, products, model_parameters['Aceto_thermodynamics'])

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
        super().__init__(microbe, educts, products, model_parameters['Fe3_thermodynamics'])
