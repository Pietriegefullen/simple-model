# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 14:55:24 2021

@author: Lara
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:45:44 2020

@author: Lara
"""

import copy
import os

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

import scipy
from scipy import optimize
import scipy.signal
from scipy.optimize import basinhopping
from scipy.optimize import differential_evolution
from scipy.optimize import curve_fit
import scipy.integrate as integ
from scipy.integrate import odeint, solve_ivp

import numpy as np

import Thermomodel
from order import pool_order, changeables_order, parameter_units, get_fixed_quantities, get_initial_guesses
from scaling import scale, rescale


#####

from USER_VARIABLES import ROOT_DIRECTORY
from data.load import load_matlab

def plot_the_ratio():

    superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()

     # plot forthe Ratio between CH4 and CO2 for all sites.
    for i in replica_list_superdata_2021_all:

        plt.figure()
        np.seterr(divide='ignore') # ignoriert durch 0 teile


        # smoothen der CH4 und Co2 messung
        CH4_smooth = scipy.signal.savgol_filter(superdata_2021_all[i]['CH4'][:],window_length=3, polyorder=2) # window size 51, polynomial order 3
        CO2_smooth = scipy.signal.savgol_filter(superdata_2021_all[i]['CO2'][:],window_length=3, polyorder=2) # window size 51, polynomial order 3
        # die differenz zwischen den gesmoothen orginaldaten, gibt die produktion zwischen zwei tagen
        CH4_smooth_diff= CH4_smooth[1:]-CH4_smooth[:-1]
        CO2_smooth_diff= CO2_smooth[1:]-CO2_smooth[:-1]
        #ratio der gesmoothen orginaldaten
        Smoothed_CH4_CO2 = CH4_smooth_diff/CO2_smooth_diff


        # differenz zwischen den orginaldaten, gibt die produktion zwischen zwei tagen (ungesmoothed)
        CH4_diff = superdata_2021_all[i]['CH4'][1:]-superdata_2021_all[i]['CH4'][:-1]
        CO2_diff = superdata_2021_all[i]['CO2'][1:]-superdata_2021_all[i]['CO2'][:-1]
        # ratio der ungesmoothed daten
        Ratio = CH4_diff/CO2_diff
        #print(np.round(Ratio))
        #smoothen des Ratios
        Smoothed_ratio = scipy.signal.savgol_filter(Ratio,window_length=7, polyorder=3)#
        #print(np.round(Smoothed_ratio), 3)

        #Anteil CH4 am Gesamtkohlenstoff
        Anteil_gesamt = 100*superdata_2021_all[i]['CH4'][1:]/(superdata_2021_all[i]['CO2'][1:]+ superdata_2021_all[i]['CH4'][1:])

        fig, ax1 = plt.subplots()

        plt.title('CH4/CO2')
        #ax1.plot(superdata_2021_all[i]['measured_time'][1:],Smoothed_CH4_CO2)
        ax1.plot(superdata_2021_all[i]['measured_time'][1:],Smoothed_ratio, 'magenta')
        ax1.plot(superdata_2021_all[i]['measured_time'][1:],Ratio, 'g--')
        ax1.plot(superdata_2021_all[i]['measured_time'][1:],Anteil_gesamt, 'gold')
        ax1.set_ylim([-10,10])
        #ax1.ylabel('Ratio')
        ax1.plot(superdata_2021_all[i]['measured_time'][1:], np.repeat(1, len(superdata_2021_all[i]['measured_time'][1:])), 'k--')
        ax2 = ax1.twinx()
        ax2.plot(superdata_2021_all[i]['measured_time'],superdata_2021_all[i]['CH4'],'ro')
        ax2.plot(superdata_2021_all[i]['measured_time'],superdata_2021_all[i]['CO2'], 'bo')

        plt.legend([i])

def curve_merger_wrapper(fixed_quantities):
    """
    läd die fixed_quantities einmal rein, für alle durchläufe von merged_curves_pred_funciton

    merged_curves_pred_function ruft den predictor auf und fügt die vorhergesagten
    verläufe von CO2 und CH4 aneinander, damit curve_fit sie gegen die gemessenen
    daten von CO2 und CH4 vergleichen kann.

    changeables: die veränderlichen Parameter, die von curve_fit angepasst werden können.
    im ersten schritt des fit-algorithmus sind diese die initiual_guesses

    """
    def merged_curves_pred_function(t, *changeables):

        changeables_dict = dict(zip(changeables_order, changeables))

        # es gibt changeales die Parameter sind (Vmax), oder die Pool-startwerte sind (Fe3), hier werden sie sortiert
        # und die initialwerte zusammengefügt
        initial_pool_values = dict()
        for pool_name in pool_order:
            if pool_name in fixed_quantities:
                initial_pool_values[pool_name] = fixed_quantities[pool_name]
            elif pool_name in changeables_dict:
                initial_pool_values[pool_name] = changeables_dict[pool_name]
            else:
                initial_pool_values[pool_name] = 0

        model_parameters = dict()
        for key, value in changeables_dict.items():
            if not key in pool_order:
                model_parameters[key] = value
        for key, value in fixed_quantities.items():
            if not key in pool_order:
                model_parameters[key] = value

        # rescale the parameters nachdem scaling
        initial_pool_values = rescale(initial_pool_values)
        model_parameters    = rescale(model_parameters)

        pool_dict = predictor(t,initial_pool_values, model_parameters)

        # nimm CO2 und CH4 und pack sie zusammen zu einer einzigen kurve.
        predicted_CO2_and_CH4 = np.concatenate((pool_dict['CO2'],pool_dict['CH4']),axis = 0)

        return predicted_CO2_and_CH4

    return merged_curves_pred_function


best = np.inf
NFeval = 0
save_file = None

def least_squares_error_wrapper(fixed_quantities_dict, measured_data_dict):
    # Hilfsfunktion für den Baisinhopper optimierer zur übergabe von fixed_quantities_dict und Measured_data_dict
    def least_squares_error_hopper(changeables_array):
        global best
        global NFeval
        global save_file

        NFeval += 1

        changeables_dict = dict(zip(changeables_order, changeables_array))

        # Sortieren der changeables und der fixed quantities
        initial_pool_dict = dict()
        model_parameters_dict = dict()
        for key in changeables_dict:
            if key in pool_order:
                initial_pool_dict[key] = changeables_dict[key]
            else:
                model_parameters_dict[key] = changeables_dict[key]

        for key in fixed_quantities_dict:
            if key in pool_order:
                initial_pool_dict[key] = fixed_quantities_dict[key]
            else:
                model_parameters_dict[key] = fixed_quantities_dict[key]


        # rescale the parameters
        initial_pool_dict     = rescale(initial_pool_dict)
        model_parameters_dict = rescale(model_parameters_dict)

        # die Berechneten Poolwerte für den Vergleich
        y_predicted_dict = predictor(measured_data_dict['measured_time'],
                                    initial_pool_dict,
                                    model_parameters_dict)
       # print('y predicted', y_predicted_dict)

        CO2_predicted = y_predicted_dict['CO2']
        CH4_predicted = y_predicted_dict['CH4']

        CO2_measured = measured_data_dict['CO2']
        CH4_measured = measured_data_dict['CH4']

        # Die Berechnung der Abweichung zwischen gemessenem und vorhergesagtem Wert
       # print('CO2predicted',CO2_predicted )
       # print('CO2_measured',CO2_measured )
        error_CO2 = CO2_predicted - CO2_measured
        error_CH4 = CH4_predicted - CH4_measured
        # predicted und measured sind jetzt immer gleich lang (siehe predictor)

        # ist es wichtiger an CO2 oder an CH4 gut zu fitten. (je höher desto wichtiger)
        weight_CO2 = 1.
        weight_CH4 = 1.
        sum_of_squared_residuals = weight_CO2*np.nanmean(error_CO2**2) + weight_CH4*np.nanmean(error_CH4**2)

        if sum_of_squared_residuals < best:
            save_string = ''
            if np.isinf(best):
                save_string += 'fixed:\n=====\n'
                for k,v in fixed_quantities_dict.items():
                    save_string += '{:<20s}'.format(k) + '{:f}'.format(v) + '\n'
                save_string += '\n\n'

            best = sum_of_squared_residuals
            # beste changeales bisher speichern
            save_string += 'SSR {:10.3f}\n'.format(sum_of_squared_residuals)
            save_string += '==============\n'
            for k,v in dict(zip(changeables_order, changeables_array)).items():
                save_string += '{:<20s}'.format(k) + '{:f}'.format(v) + '\n'
            save_string += '\n'
            with open(save_file, 'a+') as sf:
                sf.write(save_string)

        print('{:5d}: SSR {:10.3f}   best {:10.3f}'.format(NFeval,
                                                            sum_of_squared_residuals,
                                                            best))
        return sum_of_squared_residuals

    return least_squares_error_hopper

def least_squares_error(changeables_array, fixed_quantities_dict, measured_data_dict):
    """
    Für jede Parameterkombination wird der least squares error berechnet

    in den folgenden beiden blöcken werden die variablen und fixen
    parameter neu geordnet, weil predictor() die parameter getrennt haben
    will, je nach dem, ob sie initiale pool-pegel sind oder andere modell-parameter
    dazu gehen wir alle pools durch, und prüfen, ob sich der entsprechende
    startwert im dict changeables_array, oder im dict fixed_quantities_dict befindet.

    """

    changeables_dict = dict(zip(changeables_order, changeables_array))

    # rescale the parameters
    changeables_dict = rescale(changeables_dict)

    # Sortieren der changeables und der fixed quantities
    initial_pool_dict = dict()
    model_parameters_dict = dict()
    for key in changeables_dict:
        if key in pool_order:
            initial_pool_dict[key] = changeables_dict[key]
        else:
            model_parameters_dict[key] = changeables_dict[key]

    for key in fixed_quantities_dict:
        if key in pool_order:
            initial_pool_dict[key] = fixed_quantities_dict[key]
        else:
            model_parameters_dict[key] = fixed_quantities_dict[key]


    # die Berechneten Poolwerte für den Vergleich
    y_predicted_dict = predictor(measured_data_dict['measured_time'],
                                initial_pool_dict,
                                model_parameters_dict)

    CO2_predicted = y_predicted_dict['CO2']
    CH4_predicted = y_predicted_dict['CH4']

    CO2_measured = measured_data_dict['CO2']
    CH4_measured = measured_data_dict['CH4']

    # Die Berechnung der Abweichung zwischen gemessenem und vorhergesagtem Wert
    error_CO2 = CO2_predicted - CO2_measured
    error_CH4 = CH4_predicted - CH4_measured

    # ist es wichtiger an CO2 oder an CH4 gut zu fitten. (je höher desto wichtiger)
    weight_CO2 = 100.
    weight_CH4 = 100.
    sum_of_squared_residuals = weight_CO2*np.sum(error_CO2**2) + weight_CH4*np.sum(error_CH4**2)

    print('SSR',sum_of_squared_residuals)

    return sum_of_squared_residuals

def predictor(t, initial_pool_values, model_parameters):
    """
    Der SOLVER
    Hier werden die Verläufe der Pools berechnet aus den Änderungen die von Cdec
    zurückgegeben werden.

    """
    print_on_call = False
    print_only_changeables = True

    # Ausgabe meiner optimierten Werte
    if print_on_call == True:
        print('')
        print('calling predictor:')
        print('initial pool values:')
        for k, v in initial_pool_values.items():
            changeable = '*' if k in changeables_order else ''
            if not print_only_changeables or k in changeables_order:
                print('   {:20s}'.format(k), '{:8.4f}'.format(v), '{}'.format(changeable))
        print('model parameters:')
        for k, v in model_parameters.items():
            changeable = '*' if k in changeables_order else ''
            if not print_only_changeables or k in changeables_order:
                print('   {:20s}'.format(k), '{:8.4f}'.format(v), '{}'.format(changeable))
        print('')


    # alle pools, für die kein initialer zustand explizit definiert ist,
    # werden zum startzeitpunkt (t=0) auf 0  gesetzt.
    initial_pool_list = list()
    name_list = list()
    for pool_name in pool_order:
        if pool_name in initial_pool_values:
            initial_pool_list.append(initial_pool_values[pool_name])
            name_list.append(pool_name)
        else:
            initial_pool_list.append(0)
            name_list.append(pool_name)
    initial_system_state = np.array(initial_pool_list)

    #
    #print(list(zip(name_list, initial_system_state)))

    cdec = Thermomodel.Cdec_wrapper(model_parameters)


    method = 'Radau'#'LSODA' #'BDF' # 'Radau'
    try:
        solver_result = integ.solve_ivp(cdec,
                                      (0, np.max(t)),
                                      initial_system_state,
                                      method=method,
                                      t_eval = t,
                                      atol = 1e-100,
                                      rtol = 1e-1,
                                      max_step = 10)#,
                                      #first_step = 0.01)

    except Exception as ex:
        print('exception in odeint:')
        print(type(ex),':', str(ex))


    # falls solver abbricht, fülle restliche tage mit nan auf
    n_pools, n_days = solver_result.y.shape
    nan_array = np.empty((n_pools, t.size - n_days))
    nan_array[:] = np.nan
    padded_pools = np.concatenate([solver_result.y,
                                   nan_array], axis = -1)
   # print( 'results',solver_result.y)
    pool_dict = dict(zip(pool_order, padded_pools))

    return pool_dict

def compute_extra_info(pool_value_dict, model_parameters_dict):
    """
    Für unsere Plots brauchen wir einen anderen Aufruf von Cdec, weil der solver
    nicht klarkommt mit Ausgaben, die er nicht zur Berechnung berücksichtigen soll

    """

    cdec_fcn = Thermomodel.Cdec_wrapper(model_parameters_dict,
                                        return_thermodynamics=True)

    extra_curves = dict()

    first_key = list(pool_value_dict.keys())[0]
    number_of_steps = len(pool_value_dict[first_key])
    for t in range(number_of_steps):
        # cdec will den system_state als array
        system_state = np.array([pool_value_dict[pool][t] for pool in pool_order])

        dummy_time = 0
        _, extra_dict = cdec_fcn(dummy_time,system_state)

        for key in extra_dict.keys():
            if not key in extra_curves:
                extra_curves[key] = list()

            if 'DGr' in key:
                extra_curves[key].append(extra_dict[key])
            else:
                previous_value = extra_curves[key][-1] if len(extra_curves[key])>0 else 0.0
                # previous_value = 0 # crazy hack für nicht-kumulative plots
                extra_curves[key].append(previous_value + extra_dict[key])

    return extra_curves

def plot_my_data(Realdata, days_for_plot, pool_value_dict, specimen_index):
    """
    wird nur von fit_my_model, nicht von run my model aufgerufen
    """
    measurement_days = Realdata['measured_time'].astype(int)


    # plots der CH4 und CO2 kurven und  den Messwerten
    for a, col in zip(["CH4","CO2"], ["r", "b"]):

        plt.figure()
        plt.plot(measurement_days, Realdata[a], col+"o", label = "Observed")
        plt.plot(pool_value_dict[a], "k-",label = "Predicted")
        plt.ylabel(a)
        plt.legend(Realdata['Probe'])



# --------------Berechne R^2---------------------------------------------------
        measured_values = Realdata[a]
        predicted_values = pool_value_dict[a][measurement_days]
        mean_measured_value = np.mean(measured_values)
        residuals = measured_values - predicted_values
        residual_sum_of_squares = np.sum(residuals**2)
        total_sum_of_squares = np.sum((measured_values - mean_measured_value)**2)
        r_squared = 1 - residual_sum_of_squares/total_sum_of_squares

        # adjust R^2 for number of degrees of freedom
        n = measurement_days.size # number of measurements
        p = len(changeables_order) # number of explanatory terms
        adjusted_r_squared = 1 - (1-r_squared)*(n-1)/(n-p-1)

        print("r2 for "+a+" is", r_squared)
        print("r2adj  "+a+" is", adjusted_r_squared)


    plt.figure()
    plt.plot(measurement_days, Realdata['CH4'], 'r'+"o", label = "Observed")
    plt.plot(measurement_days, Realdata['CO2'], 'b'+"o", label = "Observed")
    plt.plot(pool_value_dict['CH4'], "k-",label = "Predicted")
    plt.plot(pool_value_dict['CO2'], "k-",label = "Predicted")
    #plt.ylabel()
    plt.legend(Realdata['Probe'])


# =============================================================================
    #plots aller pools
    for key in pool_value_dict.keys():
         plt.figure()
         plt.plot(range(0, max(measurement_days)+1), pool_value_dict[key], label = [key])
         plt.ylabel(key)
         plt.legend(Realdata['Probe'])
 #=============================================================================


    save_path = os.path.join(ROOT_DIRECTORY,'Figs', a +'_fit_' +str(specimen_index)+'.png')
    plt.savefig(save_path)

#---------Alle Kurven, die CO2 Produzieren oder Verbrauchen -------------------
    all_CO2_contributers = dict()
    for pool_name, pool_curve in pool_value_dict.items():
        if "CO2" in pool_name:
            all_CO2_contributers[pool_name] = pool_curve
    plt.figure()
    for pool_name, pool_curve in all_CO2_contributers.items():
        plt.plot(days_for_plot, pool_curve, label= pool_name)
    plt.legend()

#---------Alle Kurven, die CH4 Produzieren oder Verbrauchen -------------------
    all_CH4_contributers = dict()
    for pool_name, pool_curve in pool_value_dict.items():
        if "CH4" in pool_name:
            all_CH4_contributers[pool_name] = pool_curve
    plt.figure()
    for pool_name, pool_curve in all_CH4_contributers.items():
        plt.plot(days_for_plot, pool_curve, label= pool_name)
    plt.legend()

#---------CH4 und CO2 Predicted und Observed -------------------
    plt.figure()
    plt.plot( Realdata['measured_time'],Realdata['CH4'],'r.')
    plt.plot(pool_value_dict['CH4'], 'r-')
    plt.plot( Realdata['measured_time'], Realdata['CO2'],'b.')
    plt.plot(pool_value_dict['CO2'], 'b-')

    plt.show()

def callback_wrapper(objective_function):
    def opti_callback(xk, convergence):
        print('')
        print('current x:')
        for k,v in dict(zip(xk,changeables_order)):
            print(k,v)
        print('====')
        print('SSR: ', objective_function(xk))
        print('')
        return False

def fit_my_model(specimens, Site, opt):
    """
    Hier werden die Parameter des Models optimiert

    1)  die Daten werden geladen und vorbereitet
    2)  einer von drei möglichen Optimierern wird optimiert (wahl fällt im funktionsaufruf)
    3)   Mit den optimierten parametern wird das Model laufengelassen
    4)   die optimierten Kurven werden geplottet

    """
    global save_file
    date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_dir = os.path.join(ROOT_DIRECTORY, 'optimization_log')
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)
    save_file = os.path.join(log_dir, date+'_diff_evol.txt')

    #plt.close('all')
#-----------------------------laden der Daten ---------------------------------

    superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()

    for specimen_index in specimens:
        if Site == "S":
            Realdata = superdata_Sam[replica_list_Sam[specimen_index]]
        elif Site == "K":
            Realdata = superdata_Kuru[replica_list_Kuru[specimen_index]]
        elif Site == "all":
            Realdata = superdata_2021_all[replica_list_superdata_2021_all[specimen_index]]

#-----------------------------Vorbereiten der Optimierung---------------------------

        # define the initial pool values,all pools, for which no value is speicified, will be initialized as empty
        fixed_quantities_dict = get_fixed_quantities()
        m_gluc = 180                                                            # molar mass of glucose, g dw pro mol
        TOC =  Realdata['Corg (%)']/100.0                                       # Knoblauchs Daten , g dw
        specimen_mass = Realdata['weight']                                      # g (Knoblauch proben)
        Cpool_init = (10**6)* specimen_mass * TOC / m_gluc

        fixed_quantities_dict['C'] = float(Cpool_init)
        fixed_quantities_dict['DOC'] = float(Cpool_init)*0.02                   #TODO:check this 0.02 ratio in song

        initial_guess_dict = get_initial_guesses()                              # get_initial_guesses() ist eine Funktion in order


        # scale parameters to allow better optimization
        initial_guess_dict = scale(initial_guess_dict)
        #fixed_quantities_dict = scale(fixed_quantities_dict)

        # vorbereitung des startpunktes für die optimierung und der bounds für curve_fit
        initial_guess_array = [initial_guess_dict[key][0] for key in changeables_order]
        lower_bounds = [initial_guess_dict[key][1] for key in changeables_order]
        upper_bounds = [initial_guess_dict[key][2] for key in changeables_order]


        # schreibe die werte, die zwar im initial_guess_dict definiert sind,
        # aber trotzdem während der optimierung Festgehalten werden sollen,
        # in das dict für die nicht mit-optimierten werte.
        fixed_quantities_dict.update({k:v[0] for k,v in initial_guess_dict.items() if not k in changeables_order})


        # gemessene daten aus Realdata laden und zusammenhängen, um CO2 und CH4 gleichzeitig zu fitten.
        measurement_days = Realdata['measured_time']
        merged_measurements = np.concatenate((Realdata['CO2'],Realdata['CH4']), axis = 0)

#----------------------------- OPTIMIERUNG-------------------------------------

        # je nach gewählter opt wird der optimierer genutzt
        if opt =='Curve':
           # print('using curve_fit')
            print('Check if bounds are feasable')
            curve_merger = curve_merger_wrapper(fixed_quantities_dict)
            optimal_parameters , _ = curve_fit(curve_merger,
                                               measurement_days,
                                               merged_measurements, #method="dogbox",
                                               p0 = initial_guess_array,
                                               bounds=(lower_bounds, upper_bounds))

            changeables_optimal_dict = dict(zip(changeables_order, optimal_parameters))
            print(optimal_parameters)

        elif opt == "Min":

            # print('using minimize')
             initial_guess_bounds = list(zip(lower_bounds, upper_bounds))
             optimization_result = scipy.optimize.minimize(least_squares_error,
                                                      initial_guess_array,
                                                      args = (fixed_quantities_dict, Realdata),
                                                      bounds= initial_guess_bounds,
                                                      method = 'L-BFGS-B',
                                                      options = {'maxiter':20,
                                                                'disp':True})
             print(optimization_result)


             changeables_optimal_array = optimization_result.x
             changeables_optimal_dict = dict(zip(changeables_order,changeables_optimal_array))

        elif opt =='Hopper':
            #ein globaler optimierer, der größere Sprünge hüpft zwischen optimierungen
           #print('using Hopper')
            initial_guess_bounds = list(zip(lower_bounds, upper_bounds))
            least_squares_error_hopper = least_squares_error_wrapper(fixed_quantities_dict, Realdata )
            optimization_result  = scipy.optimize.basinhopping(least_squares_error_hopper,
                                                               x0 = initial_guess_array,
                                                               stepsize= 10)   #welche Größe ist angemessen?

            #print(optimization_result)


            changeables_optimal_array = optimization_result.x
            changeables_optimal_dict = dict(zip(changeables_order,changeables_optimal_array))



        elif opt =='Evolve':
            #Differential Evolution is stochastic in nature (does not use gradient methods) to find the minimum, and can search large areas of candidate space, but often requires larger numbers of function evaluations than conventional gradient-based techniques.
            #print('using Evolution')
            initial_guess_bounds = list(zip(lower_bounds, upper_bounds))
            least_squares_error_evolve = least_squares_error_wrapper(fixed_quantities_dict, Realdata )
            optimization_result  = scipy.optimize.differential_evolution(least_squares_error_evolve,
                                                                         bounds = initial_guess_bounds,
                                                                strategy= 'best1bin',
                                                                disp = True,
                                                                updating = 'immediate',
                                                                callback = callback_wrapper(least_squares_error_evolve)) #welche Größe ist angemessen?

            #print(optimization_result)


            changeables_optimal_array = optimization_result.x
            changeables_optimal_dict = dict(zip(changeables_order,changeables_optimal_array))


        changeables_optimal_dict = rescale(changeables_optimal_dict)
        #fixed_quantities_dict = rescale(fixed_quantities_dict)

        #Printing the optimal Parameters and their value
        for parameter_name in changeables_optimal_dict.keys():
            p = changeables_optimal_dict[parameter_name]
            u = parameter_units[parameter_name]
            print("{:<18} {:6.3f} {:<10}".format(parameter_name,p,u))

#-------------------Calculating the model output with optimal parameters:-------
        initial_pool_dict = dict()
        optimal_model_parameters_dict = dict()

        for key in changeables_optimal_dict:
            if key in pool_order:
                initial_pool_dict[key] = changeables_optimal_dict[key]
            else:
                optimal_model_parameters_dict[key] = changeables_optimal_dict[key]

        for key in fixed_quantities_dict:
            if key in pool_order:
                initial_pool_dict[key] = fixed_quantities_dict[key]
            else:
                optimal_model_parameters_dict[key] = fixed_quantities_dict[key]

        # predict for optimal parameters
        all_days = np.arange(max(Realdata['measured_time'])+1)
        pool_value_dict = predictor(all_days,
                                    initial_pool_dict,
                                    optimal_model_parameters_dict)

#-------------------------------Plotting------------------------------------------

        print('calling extra info')
        extra_curves = compute_extra_info( pool_value_dict, optimal_model_parameters_dict)
        pool_value_dict.update(extra_curves)

        plot_my_data(Realdata, all_days, pool_value_dict, specimen_index)

            # if input('repeat fit with new initial parameters or quit (type "q")').lower() == 'q':
            #     break
            # initial_guess_array = [changeables_optimal_dict[key] for key in changeables_order]


def run_my_model(specimens, Site = "all", plotting = "all"):
    """
    Das Model wird mit den vorgebenen Startwerten  OHNE OPTIMIERUNG laufen gelassen

    """
    #plt.close('all')

#-----------------------------laden der Daten ---------------------------------

    superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()

    for specimen_index in specimens:
        if Site == "S":
            Realdata = superdata_Sam[replica_list_Sam[specimen_index]]
        elif Site == "K":
            Realdata = superdata_Kuru[replica_list_Kuru[specimen_index]]
        elif Site == "all":
            Realdata = superdata_2021_all[replica_list_superdata_2021_all[specimen_index]]

#-----------------------------Vorbereiten für den Model run---------------------------------

    # all pools, for which no value is speicified, will be initialized as empty


    fixed_quantities_dict = get_fixed_quantities()
    m_gluc = 180                                                               # molar mass of glucose, g dw pro mol
    TOC =  Realdata['Corg (%)']/100.0                                          # Knoblauchs Daten , g dw
    specimen_mass = Realdata['weight']                                         # g (Knoblauch proben)
    Cpool_init = (10**6)* specimen_mass * TOC / m_gluc

    fixed_quantities_dict['C'] = float(Cpool_init)
    fixed_quantities_dict['DOC'] = float(Cpool_init)*0.02                      # 0.02 ratio in song, 2% sind DOC laut den Christians.

    #print('ph ganze reihe in run my model', Realdata['pH (H2O)'])
    fixed_quantities_dict['pH'] = Realdata['pH (H2O)']#[0]
    fixed_quantities_dict['weight'] =  float(specimen_mass)
    fixed_quantities_dict['water'] = float(Realdata['water'])
    # fixed_quantities_dict['Acetate'] = 1e-30                                    # to fill the initial Acetate pool manually
    # fixed_quantities_dict['H2'] = 1e-30                                    # to fill the initial Acetate pool manually
    #fixed_quantities_dict['H2O'] = 55508.43506179199 * fixed_quantities_dict['water']    # mikromol H2O in 1 ml Wasser * amount of water
    # specify Starting values
    initial_guess_dict = get_initial_guesses()

    changeables_optimal_dict = dict()

    # schreibe die werte, die zwar im initial_guess_dict definiert sind,
    # aber trotzdem während der optimierung Festgehalten werden sollen,
    # in das dict für die nicht mit-optimierten werte.
    fixed_quantities_dict.update({k:v[0] for k,v in initial_guess_dict.items() if not k in changeables_order})
    changeables_optimal_dict.update({k:v[0] for k,v in initial_guess_dict.items() if k in changeables_order})

    #Printing the Parameters and their values
    for parameter_name in changeables_optimal_dict.keys():
        p = changeables_optimal_dict[parameter_name]
        u = parameter_units[parameter_name]
        print("{:<18} {:6.3f} {:<10}".format(parameter_name,p,u))

    #Calculating the model output with choosen parameters:
    initial_pool_dict = dict()
    optimal_model_parameters_dict = dict()

    for key in changeables_optimal_dict:
        if key in pool_order:
            initial_pool_dict[key] = changeables_optimal_dict[key]
        else:
            optimal_model_parameters_dict[key] = changeables_optimal_dict[key]

    for key in fixed_quantities_dict:
        if key in pool_order:
            initial_pool_dict[key] = fixed_quantities_dict[key]
        else:
            optimal_model_parameters_dict[key] = fixed_quantities_dict[key]

    # predict for choosen Parameters
    all_days = np.arange(4500)                                                 # Days to make predictions for
    pool_value_dict = predictor(all_days,                                      # vorhersagen erstellen
                                initial_pool_dict,
                                optimal_model_parameters_dict)

    #for key in pool_value_dict:
       # print(key ,pool_value_dict[key])
#=============================================================================
#=====================================PLOTTING=================================


    print('calling extra info')
    extra_curves = compute_extra_info( pool_value_dict, optimal_model_parameters_dict)
    pool_value_dict.update(extra_curves)


    # print('gibbs acetate',pool_value_dict['DGr_Ac kJ/mol'])
    if plotting == 'all':
        plt.plot(all_days[1:],pool_value_dict['Acetate_used'][1:], 'chocolate')
    #Plots of all pools individually
        for k,v in pool_value_dict.items():
            if k == 'Acetate_used':
                continue

            plt.figure()

            plt.plot(all_days[1:],v[1:],'b-', linewidth = .5)
            plt.plot([0, max(all_days)], [0,0], 'k--', linewidth = .1)
            plt.ylabel(k)

            if k=='Acetate':
                plt.plot(all_days[1:],pool_value_dict['Acetate_used'][1:], 'chocolate')

    #Plots of all CO2 contributors in one plot
        all_CO2_contributers = dict()
        for pool_name, pool_curve in pool_value_dict.items():
            if "CO2" in pool_name:
                all_CO2_contributers[pool_name] = pool_curve
        plt.figure()
        for pool_name, pool_curve in all_CO2_contributers.items():
            plt.plot( pool_curve, label= pool_name)
        plt.legend()

    #Plots of all CH4 contributors in one plot
        all_CH4_contributers = dict()
        for pool_name, pool_curve in pool_value_dict.items():
            if "CH4" in pool_name:
                all_CH4_contributers[pool_name] = pool_curve
        plt.figure()
        for pool_name, pool_curve in all_CH4_contributers.items():
            plt.plot( pool_curve, label= pool_name)
        plt.legend()
    #=============================================================================

    #plot of measured and predicted CO2 and CH4 in one plot
        plt.figure()
        plt.plot( Realdata['measured_time'],Realdata['CH4'],'ro')
        plt.plot(pool_value_dict['CH4'])
        plt.figure()
        plt.plot( Realdata['measured_time'], Realdata['CO2'],'bo')
        plt.plot(pool_value_dict['CO2'])

    #=============================================================================
    #=============================================================================
    #=============================================================================

    #  Plot for important pools and Gibbs reacktion values on two axis

        fig, ax1 = plt.subplots()

        ax1.plot(all_days[1:], pool_value_dict['CO2_Fe3'][1:], label = 'CO2_Fe3')
        ax1.plot(all_days[1:], pool_value_dict['CO2_Ac'][1:], label = 'CO2_Ac')
        ax1.plot(all_days[1:], pool_value_dict['Acetate'][1:], label = 'Acetate')
        ax1.plot(all_days[1:], pool_value_dict['DOC'][1:], label = 'DOC')
        ax1.plot(all_days, np.repeat(0, len(all_days)), 'b--')
        ax2 = ax1.twinx()
       # ax2.plot(all_days[1:], pool_value_dict['DGr_Ac kJ/mol'][1:], 'magenta', label = 'DGr_Ac kJ/mol')
        #ax2.plot(all_days[1:], pool_value_dict['DGr_Fe3 kJ/mol'][1:], 'lime', label = 'DGr_Fe3 kJ/mol')
        ax2.plot(all_days, np.repeat(0, len(all_days)), 'b--')
        fig.tight_layout()
        ax1.legend()
        ax2.legend()

    # plots für die Gibbs reactionsenergien
        plt.figure()
        plt.plot(all_days[1:], pool_value_dict['DGr_Homo kJ/mol'][1:], 'lime', label = 'DGr_Homo kJ/mol')
        plt.plot(all_days, np.repeat(0, len(all_days)), 'k--')
        plt.legend()

        plt.figure()
        plt.plot(all_days[1:], pool_value_dict['DGr_Hydro kJ/mol'][1:], 'lime', label = 'DGr_Hydro kJ/mol')
        plt.plot(all_days, np.repeat(0, len(all_days)), 'k--')
        plt.legend()

        plt.figure()
        plt.plot(all_days[1:], pool_value_dict['DGr_Fe3 kJ/mol'][1:], 'lime', label = 'DGr_Fe3 kJ/mol')
        plt.plot(all_days[1:], pool_value_dict['DGr_Ac kJ/mol'][1:], 'magenta', label = 'DGr_Ac kJ/mol')
        plt.plot(all_days, np.repeat(0, len(all_days)), 'k--')
        plt.legend()
    ##=============================================================================

    #=============================================================================

        plt.figure()
        plt.plot( Realdata['measured_time'],Realdata['CH4'],'r.')
        plt.plot(pool_value_dict['CH4'], 'r-')
        plt.plot( Realdata['measured_time'], Realdata['CO2'],'b.')
        plt.plot(pool_value_dict['CO2'], 'b-')
    else:

        #=============================================================================
        print('me is plotting')
        plt.figure()
        plt.plot( Realdata['measured_time'],Realdata['CH4'],'r.')
        plt.plot(pool_value_dict['CH4'], 'r-')
        plt.plot( Realdata['measured_time'], Realdata['CO2'],'b.')
        plt.plot(pool_value_dict['CO2'], 'b-')
        plt.title("optimized")


if __name__ == '__main__':

# =============================================================================
#     Wahlmöglichkeiten um das Model laufen zu lassen
#
#     Site = K , alle Kurunak daten
#     Site = S, alle Samoylov daten
#     Site = all, alle Daten
#
#     opt = Curve , benutzt Curve fit
#     opt = Min, benutzt Minimize
#     opt = Hopper, benutzt den baisinhopper
#     opt = Evolve, differental evolution
#
#     specimens= der Datensatz der benutzt wird
#
# =============================================================================

    # plottet die Ratios CH4 zu CO2 für alle Proben
    #plot_the_ratio()

    #specimenlistmax = [*range(0, 35, 1)]
    #specimenlist_all= list(range(34))  # Site= "all"
    #specimenlist_Sam= [i for i in range(18)]  # Site = "S"
    #specimenlist_Kuru= [i for i in range(16)] # Site = "K"
    #specimens = [specimenlist_Sam][0]#,1,2,3,4,5,6,7]
    #specimens = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]

    specimens = [9]
    run_my_model(specimens, Site = "all", plotting = "all")
    #fit_my_model(specimens, Site = "all", opt = 'Evolve')

    # for i in list(range(3)):
    #       specimens = [i]
    #       print(specimens)
    #       run_my_model(specimens, Site = "all", plotting = "all")
    #       break

# #%%