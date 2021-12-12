import numpy as np
from predict import predictor

def fit_single_sample(sample, pathways):

    def objective_function(parameters):

        # TODO: parameters must override any model parameters!

        model_parameters =
        initial_system_state =

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
