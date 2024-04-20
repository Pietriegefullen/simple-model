import os

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import json

import system
import optimizer
import pathways
import parameters

import USER_VARIABLES

OPTIMIZATION_ALGORITHM = 'differential_evolution' #'dual_annealing' #'differential_evolution' #'direct' # 'gradient' # 'PSO'


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

def r2(predicted, measured):
    measured_mean = np.mean(measured)
    SS_res = np.sum((predicted - measured)**2)
    SS_total = np.sum((measured - measured_mean)**2)
    r2_value = 1 - SS_res/SS_total
    return r2_value

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
    
    def fit(self, replicas, algorithm = OPTIMIZATION_ALGORITHM, log = False):
        if not isinstance(replicas, list):
            replicas = [replicas]
            
        algo = optimizer.Algorithm(algorithm, 
                                   **optimizer.algo_kwargs(OPTIMIZATION_ALGORITHM))
        return algo.minimize(self, replicas, log = log)
        
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

        _, predicted_CO2 = zip(*self.system_state_log['CO2'])
        measured_CO2 = replica['CO2']
        co2_r2 = r2(predicted_CO2, measured_CO2)

        _, predicted_CH4 = zip(*self.system_state_log['CH4'])
        measured_CH4 = replica['CH4']
        ch4_r2 = r2(predicted_CH4, measured_CH4)
         
        self.system_state_log._log['R2'] = {'CO2': co2_r2,
                                            'CH4': ch4_r2}
        return self.system_state_log
    
    def parameters(self):
        return self.model_parameters

    def __str__(self):
        model_string = f'Model with {len(self.contributing_pathways)} Pathways:\n'
        model_string += len(model_string)*'=' + '\n'
        model_string += '\n'.join([str(p) for p in self.contributing_pathways])
        model_string += '\n'*2
        model_string += str(self.model_parameters)
        model_string += '\n'
        return model_string

    def save(self, target_directory, file_name):
        cfg = {'pathways': [p.__class__.__name__ 
                            for p in self.contributing_pathways],
               'parameters': self.parameters().get_config()}
        if not os.path.isdir(target_directory):
            os.makedirs(target_directory)
        with open(os.path.join(target_directory, file_name + '.json'), 'w') as df:
            json.dump(cfg, df, indent = 4)

    def model_type(self):
        simple = sorted(get_pathways('simple'))
        complex = sorted(get_pathways('complex'))
        pwys = sorted([p.__class__.__name__
                      for p in self.contributing_pathways])
        if len(simple) == len(pwys):
            for s,p in zip(simple, pwys):
                if not s == p:
                    raise Exception('Unknown model type')
            return 'simple'
        elif len(complex) == len(pwys):
            for c, p in zip(complex, pwys):
                if not c == p:
                    raise Exception('Unknown model type')
            return 'complex'
        raise Exception('Unknown model type')

    def load(self, file):
        with open(file, 'r') as df:
            cfg = json.load(df)
        
        self.__init__(cfg['pathways'])
        self.model_parameters.set(cfg['parameters'])

class ModelRun():
    def __init__(self):
        self._log = {}
        
    def keys(self):
        return self._log.keys()
    
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
            if not isinstance(self._log[n], list):
                continue
            if newfigure:
                plt.figure()
            x, y = zip(*self._log[n])
            label = n
            if 'R2' in self._log and n in self._log['R2']:
                value = self._log['R2'][n]
                label += ' ' + f'R² = {value:4.2f}'
            plt.plot(x, y, '-', label = label)
            plt.title(n)
            plt.legend()

            if 'MM' in n:
                plt.ylim([0,1])
        
    def __str__(self):
        run_string = 'Model run:'
        run_string += '\n' + '='*len(run_string) + '\n'
        run_string += '\n'.join([name + ' ' + str(self._log[name])
                          for name in sorted(self._log.keys())])
        return run_string
    

def get_best_loss_parameters(parameter_source):
    all_files = []
    for f in os.listdir(parameter_source):
        if 'loss_' in f:
            loss = float(f.split('loss_')[-1])
            file = os.path.join(parameter_source, f)
            all_files.append((loss, file))
    if len(all_files) == 0:
        raise Exception('loading parameters failed')
    best_loss, best_loss_file = list(sorted(all_files))[0]
    with open(best_loss_file, 'r') as pf:
        best_parameters = json.load(pf)
    return best_loss, best_parameters

if __name__ == '__main__':
    import data
    d = data.get_data_before_day()
    #
    model = Model(get_pathways('simple'))
    #print(model)
    #model.parameters().set('default')
    #print(model)
    #replica = d['13546']
    #run = model.predict(replica)
    #run.plot(['CO2', 'CH4'], newfigure = False)
    #replica.plot()
    #plt.figure()
    #run.plot(['DOC'])
    #1/0
    
    plot_log = False
    result_sample = '1375'
    folders = []
    for _d in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
        if result_sample in _d:
            folders.append(_d)
    #results_folder = 'fit_13544_2024-04-18--09-14-54'
    
    for results_folder in folders:
        parameter_source = os.path.join(USER_VARIABLES.LOG_DIRECTORY, results_folder)
        best_loss, p = get_best_loss_parameters(parameter_source)
        print('best loss', best_loss)
        
        #model.parameters().set('default')

        #del p['M_Homo']
        model.parameters().set(p)
        print(model)

        fit_replicas = [s for s in results_folder.replace('simple', '').replace('complex','').replace('fit_', '').replace('log', '').split('_2024')[0].replace(' ', '_').split('_') if not s == '']

        for repl in fit_replicas:
            replica = d[repl]

            model_run = model.predict(replica)

            plt.figure()
            model_run.plot(['CO2', 'CH4'], newfigure = False)
            replica.plot(log = plot_log, newfigure = False)
            ax = plt.gca()
            t = ax.get_title()
            plt.title(t + results_folder)
            
            #plt.figure()
            #model_run.plot(['CO2', 'CH4'], newfigure = False)
            #replica.plot(newfigure = False)
            #ax = plt.gca()
            #t = ax.get_title()
            #plt.title(t + results_folder)
            

        #model_run.plot(['Fermentation_MM'])
        #model_run.plot(['Hydrolysis_MM'])
        #model_run.plot(['Hydro_MM'])
        #model_run.plot(['Aceto_MM'])

    plt.show()
