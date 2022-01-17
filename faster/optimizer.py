import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import multiprocessing

from predict import predictor
from ORDER import POOL_ORDER, pool_index
import data
import OPTIMIZATION_PARAMETERS
from OPTIMIZATION_PARAMETERS import CHANGEABLES


def fit_specimen(specimen_index, site, pathways, fixed_parameters, algo): 
    # falls run not fit sin fixed_parameters die model_parameters
    for name in CHANGEABLES:
        if name in fixed_parameters.keys():
            print(f'Deleting {name} from fixed parameters')
            del fixed_parameters[name]
            #raise Exception('A changeable model parameter is also in fixed_parameters!')

    initial_guess_dict = OPTIMIZATION_PARAMETERS.get_initial_guesses()
    lower_bounds = [initial_guess_dict[key][1] for key in CHANGEABLES]
    upper_bounds = [initial_guess_dict[key][2] for key in CHANGEABLES]

    initial_guess_bounds = list(zip(lower_bounds, upper_bounds))

    sample_list = [(specimen_index, site)]
    objectives = [SpecimenObjective(pathways,
                                    fixed_parameters,
                                    data.specimen_data(sample, site))
                  for sample, site in sample_list]

    changeables_string = '\nchangeables:\n' + '='*11 + '\n   ' + '\n   '.join(CHANGEABLES) + '\n'
    print(changeables_string)
    workers = OPTIMIZATION_PARAMETERS.WORKERS
    if workers == -1:
        workers = multiprocessing.cpu_count()
        
    if algo == 'differential_evolution':
        
        opti_string = f'fitting {len(objectives):d} samples with DIFFERENTIAL EVOLUTION on {workers:d} workers'
        print('\n' + '#'*len(opti_string) + '\n' + opti_string + '\n' + '#'*len(opti_string) + '\n')

        optimization_result  = scipy.optimize.differential_evolution(ObjectiveFunction(),
                                                                     bounds = initial_guess_bounds,
                                                                     args = (objectives,),
                                                                     disp = True,
                                                                     workers = workers,
                                                                     **OPTIMIZATION_PARAMETERS.DIFF_EVOL_PARAMETERS)
        changeables_optimal_array = optimization_result.x

    elif algo == 'gradient':
        method_name = OPTIMIZATION_PARAMETERS.GRADIENT_PARAMETERS['method']
        opti_string = f'fitting {len(objectives):d} samples with {method_name:}'
        print('\n' + '#'*len(opti_string) + '\n' + opti_string + '\n' + '#'*len(opti_string) + '\n')

        initial_guess_array = np.array([initial_guess_dict[k][0] for k in CHANGEABLES])
        optimization_result = scipy.optimize.minimize(ObjectiveFunction(),
                                                      initial_guess_array,
                                                      args = (objectives,),
                                                      bounds = initial_guess_bounds,
                                                      **OPTIMIZATION_PARAMETERS.GRADIENT_PARAMETERS)
        changeables_optimal_array = optimization_result.x

    elif algo == 'PSO':
        import pyswarms as ps
        options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}

        # Call instance of GlobalBestPSO
        optimizer = ps.single.GlobalBestPSO(n_particles=20, 
                                            dimensions=len(OPTIMIZATION_PARAMETERS.CHANGEABLES),
                                            options=options,
                                            bounds = (np.array(lower_bounds), np.array(upper_bounds)))
        _, changeables_optimal_array = optimizer.optimize(ParticleObjective(),
                                                           iters = 1000,
                                                           specimen_objectives=objectives,
                                                           n_processes = workers)

    changeables_optimal_dict = dict(zip(CHANGEABLES, changeables_optimal_array))

    print('optimal parameters:')
    for k, v in changeables_optimal_dict.items():
        print(f'   {k[:15]:15} {v}')
    print('')

    return changeables_optimal_dict


class ParticleObjective():
    
    def __call__(self, changeable_parameters, specimen_objectives):
        particle_fitnesses = []
        for particle_changeables in changeable_parameters:
            losses = [sample_loss(particle_changeables) for sample_loss in specimen_objectives]
            total_loss = np.sum(losses)
            particle_fitnesses.append(total_loss)
        
        return np.array(particle_fitnesses)


class ObjectiveFunction:
    def __init__(self):
        self.best = np.inf

    def __call__(self,changeable_parameters, specimen_objectives):
        losses = [sample_loss(changeable_parameters) for sample_loss in specimen_objectives]
        total_loss = np.sum(losses)
        print(f'{total_loss:.5e}')

        if total_loss < self.best:
            self.best = total_loss
            specimen_objectives[0].plot(changeable_parameters)

        return total_loss

class SpecimenObjective:
    def __init__(self,pathways, fixed_parameters, measured_data_dict):
        self.pathways = pathways
        self.fixed_parameters = fixed_parameters
        self.measured_data_dict = measured_data_dict

        if OPTIMIZATION_PARAMETERS.PLOT_LIVE_FIT:
            plt.ion()
            self.fig, (self.ax0, self.ax1) = plt.subplots(2,1)

    def __call__(self, changeable_parameters):

        self.fixed_parameters.update({k:v for k,v in zip(CHANGEABLES, changeable_parameters)})

        measure_days = self.measured_data_dict['measured_time']
        y_predicted_dict = predictor(t_eval = measure_days,
                                    model_parameters = self.fixed_parameters,
                                    chosen_pathways = self.pathways,
                                    verbose = False,
                                    mark = CHANGEABLES)

        CO2_predicted = y_predicted_dict['CO2']
        CH4_predicted = y_predicted_dict['CH4']

        CO2_measured = self.measured_data_dict['CO2']
        CH4_measured = self.measured_data_dict['CH4']

        measure_day_weight = np.ones((len(CO2_measured),))
        if OPTIMIZATION_PARAMETERS.MEASURE_DAYS_WEIGHTING:
            measure_day_weight = np.concatenate([np.array([1]), np.diff(measure_days)])

        # Die Berechnung der Abweichung zwischen gemessenem und vorhergesagtem Wert
        error_CO2 = measure_day_weight * (CO2_predicted - CO2_measured)
        error_CH4 = measure_day_weight * (CH4_predicted - CH4_measured)

        # ist es wichtiger an CO2 oder an CH4 gut zu fitten. (je höher desto wichtiger)
        weight_CO2 = 1.
        weight_CH4 = 1.

        sum_of_squared_residuals = weight_CO2*np.sum(error_CO2**2) + weight_CH4*np.sum(error_CH4**2)

        return sum_of_squared_residuals

    def plot(self,changeable_parameters):
        if not OPTIMIZATION_PARAMETERS.PLOT_LIVE_FIT:
            return

        self.fixed_parameters.update({k:v for k,v in zip(CHANGEABLES, changeable_parameters)})

        measure_days = self.measured_data_dict['measured_time']
        y_predicted_dict = predictor(t_eval = measure_days,
                                    model_parameters = self.fixed_parameters,
                                    chosen_pathways = self.pathways,
                                    verbose = False,
                                    mark = CHANGEABLES)

        CO2_predicted = y_predicted_dict['CO2']
        CH4_predicted = y_predicted_dict['CH4']

        CO2_measured = self.measured_data_dict['CO2']
        CH4_measured = self.measured_data_dict['CH4']

        self.ax0.clear()
        self.ax1.clear()

        plt.figure(self.fig.number)
        self.ax0.plot(measure_days, CO2_predicted)
        self.ax0.plot(measure_days, CO2_measured, 'x')
        self.ax0.set_ylim([np.min(CO2_measured), np.max(CO2_measured)])
        self.ax0.set_title('CO2')

        self.ax1.plot(measure_days, CH4_predicted)
        self.ax1.plot(measure_days, CH4_measured, 'x')
        self.ax1 .set_ylim([np.min(CH4_measured), np.max(CH4_measured)])
        self.ax1.set_title('CH4')

        loss = self.__call__(changeable_parameters)
        self.fig.suptitle(f'{loss:>.2f}')

        plt.draw()
        plt.pause(0.0001)
