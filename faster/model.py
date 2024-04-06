import os

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import json

import system
import optimizer
import pathways
import parameters

OPTIMIZATION_ALGORITHM = 'PSO' #'dual_annealing' #'differential_evolution' #'direct' # 'gradient' # 'PSO'


def get_pathways(model_type):
    basic = ['Hydrolysis',
             'Fermentation',
             'Hydro',
             'Aceto']
    if model_type == 'complex':
        return basic + ['Homo',
                        'Fe3']
    elif model_type == 'simple':
        return basic
    else:
        raise NotImplementedError()


class Model():
    def __init__(self, pwys):
        self.system_state_log = ModelRun()
        self.model_parameters = parameters.ModelParameters()
        
        pathway_classes = [pathways.pathway_by_name(p) if isinstance(p, str) else p 
                           for p in pwys]
        self._unbuilt_contributing_pathways = pathway_classes
        self.contributing_pathways = None
        self.build(quiet = True)
        
    def build(self, quiet = False):
        self.contributing_pathways = [p(self.model_parameters) 
                                      for p in self._unbuilt_contributing_pathways]
        for p in self.contributing_pathways:
            p.inject_logger(self.system_state_log)
        
        # to initialize model parameters used in initial state
        _ = system.initial_state(None, self.model_parameters)
        
    def __call__(self, t, S):
        S = np.where(S < 1e-40, 0, S)
        
        dSj_dt = np.stack([pathway_j(t, S)
                            for pathway_j in self.contributing_pathways], axis = -1)
        
        dS_dt = np.sum(dSj_dt, axis = -1)
        dS_dt = np.clip(dS_dt, -S, np.inf) # don't let pools become negative
        return dS_dt
    
    def fit(self, replicas, algorithm = OPTIMIZATION_ALGORITHM):
        if not isinstance(replicas, list):
            replicas = [replicas]
            
        algo = optimizer.Algorithm(algorithm, 
                                   **optimizer.algo_kwargs(OPTIMIZATION_ALGORITHM))
        return algo.minimize(self, replicas)
        
        
    def predict(self, replica, t = None, quiet = False):
        if t is None:
            t = replica['days']
        self.build(quiet = quiet)
        S0 = system.initial_state(replica, self.parameters())
        self.parameters().check()
        self.system_state_log.reset()
        
        solver_result = scipy.integrate.solve_ivp(self, (0, max(t)),
                                                  S0, 
                                                  t_eval = t,
                                                  method = 'LSODA',
                                                  max_step = 10,
                                                  first_step = 1e-6)
        
        for t, S in zip(t, np.transpose(solver_result.y)):
            for Si, pool_name in zip(S, system.SYSTEM):
                self.system_state_log.log(pool_name, t, Si)
        
        return self.system_state_log
    
    def parameters(self):
        return self.model_parameters
    
    def __str__(self):
        model_string = f'Model with {len(self.contributing_pathways)} Pathways:\n'
        model_string += len(model_string)*'=' + '\n'
        model_string += '\n'.join([str(p) for p in self.contributing_pathways])
        model_string += '\n'
        model_string += str(self.model_parameters)
        model_string += '\n'.join( [str(p) for p in self.contributing_pathways])
        return model_string

    def save(self, target_directory, file_name):
        cfg = {'pathways': [p.__class__.__name__ 
                            for p in self.contributing_pathways],
               'parameters': self.parameters().get_config()}
        if not os.path.isdir(target_directory):
            os.makedirs(target_directory)
        with open(os.path.join(target_directory, file_name + '.json'), 'w') as df:
            json.dump(cfg, df, indent = 4)
            
    def load(self, file):
        with open(file, 'r') as df:
            cfg = json.load(df)
        
        self.__init__(cfg['pathways'])
        self.model_parameters.set(cfg['parameters'])

class ModelRun():
    def __init__(self):
        self._log = {}
        
    def __eq__(self, other):
        return self._log == other._log
    
    def __getitem__(self, key):
        return self._log[key]
        
    def log(self, name, t, value):
        if not name in self._log:
            self._log[name] = []
        
        self._log[name].append((t,value))
        
    def reset(self):
        self._log = {}
        
    def plot(self, name = None, newfigure = True):
        if name is None:
            name = list(self._log.keys())
            
        if not isinstance(name, list):
            name = [name]
        
        for n in name:
            if not n in self._log:
                print(n + ' not logged')
            if newfigure:
                plt.figure()
            x, y = zip(*self._log[n])
            plt.plot(x, y, '-', label = n)
            plt.title(n)
        
    def __str__(self):
        run_string = 'Model run:'
        run_string += '\n' + '='*len(run_string) + '\n'
        run_string += '\n'.join([name + ' ' + str(self._log[name])
                          for name in sorted(self._log.keys())])
        return run_string
    

if __name__ == '__main__':
    import data
    d = data.get_data_before_carex()
    
    model = Model(get_pathways('simple'))
    
    p = {
    "death_rate": 8.33e-05,
    "Acetate": 1,
#    "temperature": 4.0,
#    "C": 2546.5533333333337,
#    "DOC": 50.93106666666667,
#    "pH": 3.95,
#    "weight": 11.82,
#    "water": 4.0,
#    "H2O": 222033.74024716797,
#    "M_Fe3": 0.15377556552732402,
    "M_Ferm": 0.29386173195040044,
    "M_Hydro": 0.43281383907529236,
    "M_Homo": 0.2574949776769526,
    "Hydrolysis_v_max": 0.6166340111649226,
    "Ferm_v_max": 1.070276371909682,
#    "Vmax_Fe3": 1.3804244562902253,
#    "Vmax_Homo": 0.9318093189492231,
    "Hydro_v_max": 0.7064582813317815,
    "Ac_v_max": 0.4599047701351146,
    "Hydrolysis_Kmb": 288.99678942466437,
    "Aceto_Km_Ac": 145.99705636830586,
#    "Km_Homo_CO2": 376.8013896720814,
#    "Km_Homo_H2": 688.3240608121672,
    "Hydro_Km_CO2": 661.7562751340953,
    "Hydro_Km_H2": 497.8934720994153,
#    "Km_Fe3_Fe3": 173.34626557916957,
#    "Km_Fe3_Acetate": 637.4030208609411,
    "Ferm_Km": 160.15461453008587,
    "Ferm_inhibition": 4.643075236732733,
    "Fe3": 81.99097055433658,
    "M_Ac": 0.014042559314258995,
    "Ferm_CUE": 0.30944032735284144,
#    "CUE_Fe3": 0.012291327263939777,
    "Ac_CUE": 0.5724913805190271,
#    "CUE_Homo": 0.4988932700684491,
    "Hydro_CUE": 0.5054549655151662
}
 
    model.parameters().set('default')
    model.parameters().set(p)
    #results = model.fit([d['13514'], d['13515']], 'PSO')
    model_run = model.predict(d['13514'])
    model_run.plot(['CO2', 'CH4'], newfigure = False)
    d['13514'].plot()
    plt.figure()
    model_run.plot()
    plt.show()