import numpy as np
import scipy.optimize

from predict import predictor
from ORDER import POOL_ORDER, pool_index
import data
import OPTIMIZATION_PARAMETERS
from OPTIMIZATION_PARAMETERS import CHANGEABLES

def fit_specimen(specimen_index, site, pathways, fixed_parameters, algo, verbose = True):

    for name in CHANGEABLES:
        if name in fixed_parameters.keys():
            print(f'Deleting {name} from fixed parameters')
            del fixed_parameters[name]
            #raise Exception('A changeable model parameter is also in fixed_parameters!')

    initial_guess_dict = OPTIMIZATION_PARAMETERS.get_initial_guesses()
    lower_bounds = [initial_guess_dict[key][1] for key in CHANGEABLES]
    upper_bounds = [initial_guess_dict[key][2] for key in CHANGEABLES]

    objective_function = objective_builder([(specimen_index, site)], pathways, fixed_parameters)

    initial_guess_bounds = list(zip(lower_bounds, upper_bounds))

    if algo == 'differential_evolution':
        optimization_result  = scipy.optimize.differential_evolution(objective_function,
                                                                     bounds = initial_guess_bounds,
                                                                     disp = True,
                                                                     **OPTIMIZATION_PARAMETERS.DIFF_EVOL_PARAMETERS) #welche Größe ist angemessen?

    elif algo == 'gradient':
        initial_guess_array = np.array([initial_guess_dict[k][0] for k in CHANGEABLES])
        optimization_result = scipy.optimize.minimize(objective_function,
                                                      initial_guess_array,
                                                      bounds = initial_guess_bounds,
                                                      **OPTIMIZATION_PARAMETERS.GRADIENT_PARAMETERS)

    changeables_optimal_array = optimization_result.x
    changeables_optimal_dict = dict(zip(CHANGEABLES, changeables_optimal_array))

    print('optimal parameters:')
    for k, v in changeables_optimal_dict.items():
        print(f'   {k[:15]:15} {v}')
    print('')

    return changeables_optimal_dict

def objective_builder(sample_list, pathways, fixed_parameters):

    objectives = [specimen_objective(pathways,
                                     fixed_parameters,
                                     data.specimen_data(sample, site))
                  for sample, site in sample_list]

    def objective_function(changeable_parameters):
        losses = [sample_loss(changeable_parameters) for sample_loss in objectives]
        total_loss = np.sum(losses)
        print(f'{total_loss:.5e}')
        return total_loss

    return objective_function


def specimen_objective(pathways, fixed_parameters, measured_data_dict):
    if OPTIMIZATION_PARAMETERS.PLOT_LIVE_FIT:
        import matplotlib.pyplot as plt
        plt.ion()
        fig, (ax0, ax1) = plt.subplots(2,1)

    def objective_function(changeable_parameters):

        fixed_parameters.update({k:v for k,v in zip(CHANGEABLES, changeable_parameters)})

        measure_days = measured_data_dict['measured_time']
        y_predicted_dict = predictor(t_eval = measure_days,
                                    model_parameters = fixed_parameters,
                                    chosen_pathways = pathways,
                                    verbose = False,
                                    mark = CHANGEABLES)

        CO2_predicted = y_predicted_dict['CO2']
        CH4_predicted = y_predicted_dict['CH4']

        CO2_measured = measured_data_dict['CO2']
        CH4_measured = measured_data_dict['CH4']

        if OPTIMIZATION_PARAMETERS.PLOT_LIVE_FIT:
            ax0.plot(measure_days, CO2_predicted)
            ax0.plot(measure_days, CO2_measured, 'x')
            ax0.set_ylim([np.min(CO2_measured), np.max(CO2_measured)])
            ax0.set_title('CO2')

            ax1.plot(measure_days, CH4_predicted)
            ax1.plot(measure_days, CH4_measured, 'x')
            ax1 .set_ylim([np.min(CH4_measured), np.max(CH4_measured)])
            ax1.set_title('CH4')

            plt.draw()
            plt.pause(0.0001)
            ax0.clear()
            ax1.clear()

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

    return objective_function
