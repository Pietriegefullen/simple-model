import numpy as np
import scipy.optimize

from predict import predictor
from ORDER import POOL_ORDER, pool_index
import data
import OPTIMIZATION_PARAMETERS
from OPTIMIZATION_PARAMETERS import CHANGEABLES

def fit_specimen(specimen_index, site, pathways, parameters, algo):

    model_parameters = data.model_parameters_from_data(specimen_index, site)
    fixed_parameters = {k:v for k,v in model_parameters.items()
                        if not k in CHANGEABLES}
    fixed_parameters.update(parameters)

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

    return changeables_optimal_dict

def objective_builder(sample_list, pathways, fixed_parameters):

    objectives = [specimen_objective(sample, site, pathways, fixed_parameters) for sample, site in sample_list]

    def objective_function(changeable_parameters):
        losses = [sample_loss(changeable_parameters) for sample_loss in objectives]
        return np.sum(losses)

    return objective_function


def specimen_objective(specimen_index, site, pathways, fixed_parameters):

    measured_data_dict = data.specimen_data(specimen_index, site)

    def objective_function(changeable_parameters):

        # get any fixed or changeable parameters that initialize pools
        initial_system_state = np.zeros((len(POOL_ORDER),))
        model_parameters = {}
        for name, value in zip(CHANGEABLES, changeable_parameters):
            if name in POOL_ORDER:
                initial_system_state[pool_index(name)] = value

            else:
                model_parameters[name] = value

        for name, value in fixed_parameters.items():
            if name in POOL_ORDER:
                initial_system_state[pool_index(name)] = value

            else:
                model_parameters[name] = value

        y_predicted_dict = predictor(t_eval = measured_data_dict['measured_time'],
                                    model_parameters = model_parameters,
                                    all_pathways = pathways,
                                    initial_system_state = initial_system_state)

        CO2_predicted = y_predicted_dict['CO2']
        CH4_predicted = y_predicted_dict['CH4']

        CO2_measured = measured_data_dict['CO2']
        CH4_measured = measured_data_dict['CH4']

        # Die Berechnung der Abweichung zwischen gemessenem und vorhergesagtem Wert
        error_CO2 = CO2_predicted - CO2_measured
        error_CH4 = CH4_predicted - CH4_measured

        # ist es wichtiger an CO2 oder an CH4 gut zu fitten. (je höher desto wichtiger)
        weight_CO2 = 1.
        weight_CH4 = 1.

        # TODO: weighting for point density

        sum_of_squared_residuals = weight_CO2*np.sum(error_CO2**2) + weight_CH4*np.sum(error_CH4**2)

        return sum_of_squared_residuals

    return objective_function
