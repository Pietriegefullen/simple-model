import numpy as np

import model

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

    def get_config(self):
        return {'name': self.name,
                'stoichiometry': self.stoichiometry,
                'Km': self.Km,
                'inhibition': self.inhibition}


class Hydrolysis(model.Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'C')]
        products = [Substance(1, 'DOC')]
        microbe = Microbe(name = 'M_Ferm',
                          v_max = model_parameters['Hydrolysis_v_max'],
                          Kmb = model_parameters['Hydrolysis_Kmb'],
                          use_thermodynamics = False)
        super().__init__(microbe, educts, products)
    
class Fermentation(model.Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(6, 'DOC', 
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
        

class Hydro(model.Pathway):
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

class Homo(model.Pathway):
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

class Aceto(model.Pathway):
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

class Fe3(model.Pathway):
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
