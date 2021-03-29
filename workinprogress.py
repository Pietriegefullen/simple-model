# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 13:27:41 2021

@author: Lara
"""

import numpy as np

parameter_order = ['vmax_ferm',
                   'vmax_aceto']


def array_to_dict(array_of_parameters):
    
    dict_of_parameters = dict()
    for i in range(len(parameter_order)):
        parameter_name = parameter_order[i]
        parameter_value = array_of_parameters[i]
        
        dict_of_parameters[parameter_name] = parameter_value
    
    return dict_of_parameters


def dict_to_array(dict_of_parameters):
    
    array_of_parameters = list()
    for parameter_name in parameter_order:
        parameter_value = dict_of_parameters[parameter_name]
    
        array_of_parameters.append(parameter_value)
    
    return np.array(array_of_parameters)































def sort_my_parameters(parameter_as_dict):
    
    parameters_as_array = list()
    
    for key in parameter_order:
        parameters_as_array.append(parameter_as_dict[key])
    
    return np.array(parameters_as_array)


















def predictor(t_out, initial_values, parameters):
    
    
    
    
    return y


































def least_squares_error(all_model_parameters):
    
    initial_values = dict()
    parameters = dict()
    
    # TODO: was sind die pools und deren initialwerte?
    
    y_predicted = predictor(days_with_measurements,
                            initial_values, 
                            parameters)
    
    CO2_predicted = y_predicted['co2']
    CH4_predicted = y_predicted['ch4']
    
    error_CO2 = CO2_predicted - CO2_measured
    error_CH4 = CH4_predicted - CH4_measured
    
    sum_of_squared_residuals = np.sum(error_CO2**2 + error_CH4**2)
    
    

    return sum_of_squared_residuals

























def fit_my_model():
    
    # load data
    
    fitted_parameters_as_dict['vmax_ferm'] = 10
    fitted_parameters_as_dict['vmax_aceto'] = 5
    
    
    fitted_parameters_as_array = dict_to_array(fitted_parameters_as_dict)
    
    
    optimal_parameters = opti_algo(least_squares_error, initial_model_parameters)
    
    vmax_ferm_opti = 
    
    
x0 = [5, 4]
result = scipy.optimize.minimize(least_squares_error, x0)
    
xopt = result.x

parameters_as_dict = helper(xopt)

y_best_fit = predictor(range(1,1000), parameters_as_dict)

for pool_name, pool_curve in y_best_fit.items():
    
    


# xopt = [3.456, 7.45]
