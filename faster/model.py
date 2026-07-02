import os
import multiprocessing

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import json

import system
import optimizer
from optimizer import r2
import pathways
import parameters

import USER_VARIABLES

OPTIMIZATION_ALGORITHM = 'differential_evolution' #'dual_annealing' #'differential_evolution' #'direct' # 'gradient' # 'PSO'


def integrate(f, t, S0, solver_result, reset_Fe3):
    #print('solving IVP')
    
    t_after = None
    if isinstance(reset_Fe3, int):
        if reset_Fe3 >= max(t):
            raise Exception('Trying to reset iron after run end.')
        t = np.sort(np.unique(np.concatenate([t, (reset_Fe3,)])))
        t_before = t[t <= reset_Fe3]
        t_after = t[t >= reset_Fe3]
        if not reset_Fe3 in t:
            t_before = np.concatenate([t_before, [reset_Fe3]])
            t_after = np.concatenate([[reset_Fe3], t_after], axis = 0)
    else:
        t_before = t

    result = scipy.integrate.solve_ivp(f, (0, max(t_before)),
                                                  S0, 
                                                  t_eval = t_before,#np.arange(0, max(t_before)+1),
                                                  method = 'LSODA',
                                                  max_step = 10,
                                                  first_step = 1e-6, 
                                                  #min_step = 1e-4
                                                  )

    solver_result.append(result)

    if not t_after is None:
        S1 = result.y[:,-1] # system state on last day before resetting
        fe3_index = system.index('Fe3')
        print('resetting Fe3 on day', min(t_after), 'from', S1[fe3_index], 'to', S0[fe3_index])
        S1[fe3_index] = S0[fe3_index]
        print('solving IVP after reset')
        t_eval = np.sort(np.unique(np.concatenate([t_after,
                                                   np.arange(min(t_after), max(t_after), 10) ])))
        result_after = scipy.integrate.solve_ivp(f, (min(t_after), max(t_after)),
                                                 S1,
                                                 t_eval = t_eval,
                                                 method = 'LSODA',
                                                 max_step = 10,
                                                 first_step = 1e-6)

        solver_result[0].y = np.concatenate([solver_result[0].y,
                                             result_after.y], axis = -1)
        
        solver_result[0].t = np.concatenate([solver_result[0].t,
                                             result_after.t])
    
    return solver_result
        


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
        print('model type:', model_type)
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
    
    def fit(self, replicas, algorithm = OPTIMIZATION_ALGORITHM, log_co2 = True, log_ch4 = True,
            fit_from = 0, fit_to = None):
        if not isinstance(replicas, list):
            replicas = [replicas]
            
        algo = optimizer.Algorithm(algorithm, 
                                   **optimizer.algo_kwargs(OPTIMIZATION_ALGORITHM))
        return algo.minimize(self, replicas, log_co2 = log_co2, log_ch4 = log_ch4,
                             fit_from = fit_from, fit_to = fit_to)
        
    def predict(self, replica, t = None, quiet = False, parallel = False, 
                reset_Fe3 = None, days_beyond_reset = 1000):
        measured_days = replica['days']

        t_eval = measured_days
        if not t is None:
            t_eval = np.array(t)
            
        if not reset_Fe3 is None:
            last_day = reset_Fe3 + days_beyond_reset
            t = np.concatenate([t, [last_day]], axis = 0)
            
        # prepare for solving
        self.build(quiet = quiet)
        S0 = system.initial_state(replica, self.parameters())
        self.parameters().check()
        self.system_state_log.reset()
      
        # solve initial value problem
        if parallel:
            manager = multiprocessing.Manager()
            solver_result = manager.list()
            p = multiprocessing.Process(target = integrate,
                                        args = (self, t_eval, S0, solver_result, reset_Fe3))
            p.daemon = True
            p.start()
            p.join(timeout = 10.)
            if p.is_alive():
                p.terminate()
                p.join()
        else:
            solver_result = []
            integrate(self, t_eval, S0, solver_result, reset_Fe3)

        if not len(solver_result) == 1:
            print('o', end = '', flush = True)
            raise Exception('timeout')
            
        solver_result = solver_result[0]

        # add pool values to log
        for Si, pool_name in zip(solver_result.y, system.SYSTEM):
            self.system_state_log.log(pool_name, solver_result.t, Si)
   
        # compute R2 values
        used_measured_indices = np.array([np.nonzero(measured_days == t)[0] for t in t_eval])
        
        _, predicted_CO2 = self.system_state_log['CO2']
        predicted_CO2_on_measured = predicted_CO2
        measured_CO2 = replica['CO2'][used_measured_indices]
        co2_r2 = r2(predicted_CO2_on_measured, measured_CO2, log = True)
        
        _, predicted_CH4 = self.system_state_log['CH4']
        predicted_CH4_on_measured = predicted_CH4
        measured_CH4 = replica['CH4'][used_measured_indices]
        ch4_r2 = r2(predicted_CH4_on_measured, measured_CH4, log = True)
        
        self.system_state_log._log['CO2_on_measured'] = t_eval, predicted_CO2_on_measured
        self.system_state_log._log['CH4_on_measured'] = t_eval, predicted_CH4_on_measured
        self.system_state_log._log['R2'] = {'CO2': co2_r2,
                                            'CH4': ch4_r2}
        
        add_to_log = []
        for k, val in self.system_state_log._log.items():
            if 'CH4 from' in k:
                t, v = val
                delta_t = np.diff(t)
                integral = delta_t*(v[:-1]+v[1:])*0.5 # trapezoidal rule
                n = k + ' (integrated)'
                add_to_log.append((n, t[1:], np.cumsum(integral)))
                
        for name, ts, vs in add_to_log:
            self.system_state_log.log(name, ts, vs)
            
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
       
    def log_snap(self, name, t, value):
        if not name in self._log:
            ts = np.empty((0,))
            vs = np.empty((0,))
            self._log[name] = (ts, vs)
        ts, vs = self._log[name]
        ts = np.concatenate([ts, np.reshape(t,(1,))],
                            axis = 0)
        vs = np.concatenate([vs, np.reshape(value,(1,))],
                            axis = 0)
        self.log(name, ts, vs)
        
    def log_snap_2(self, name, t, value):
        if not name in self._log:
            self._log[name] = [[],[]]
        self._log[name][0].append(t)
        self._log[name][1].append(value)

    def log(self, name, ts, values):
        #if not name in self._log:
        #    self._log[name] = []
        self._log[name] = (ts, values)

    def reset(self):
        self._log.clear()
        
    def plot(self, name = None, newfigure = True, log = False):
        if name is None:
            name = list(self._log.keys())
            
        if not isinstance(name, list):
            name = [name]
        
        for n in name:
            if not n in self._log:
                print(n + ' not logged')
            if newfigure:
                fig, ax = plt.subplots()
            else:
                fig = plt.gcf()
                
            x, y = self._log[n]
            label = n
            mark = '-'
            #if 'R2' in self._log and n in self._log['R2']:
             #   value = self._log['R2'][n]
              #  label += ' ' + f'R² = {value:4.2f}'
               # mark = 'x'
            if n == 'CH4':
                ax = fig.axes
                if isinstance(ax, list):
                    ax = ax[0] 
                    ch4_ax = ax.twinx()
                else:
                    ch4_ax = ax
                ch4_ax.plot(x, y, mark, label = label, color = 'orange')
            else:
                ax = plt.gca()
                ax.plot(x, y, mark, label = label)
                
            plt.title(n)
            plt.legend()
            
            if log:
                plt.yscale('log')
                plt.title(n + ' (log)')
                
            elif n in system.SYSTEM:
                #plt.yscale('log')
                pass
        
            elif 'MM' in n:
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


