# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:45:44 2020

@author: Lara
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import numpy as np
#from SimpleIN import Cdec # importiert das eigentliche mathematische Model 
from scipy.integrate import odeint
import SimpleIN

#Importieren benutzter Funktionen und Daten
from Babypascal import Mol_nach_Pa


pool_order = ['C',
            'Fe',
            'M_Ac',
            'M_Ferm',
            'M_Fe',
            'M_Hydro',
            'M_Homo',
            'CH4',
            'CO2',
            'CO2_Ac',
            'Acetate',
            'H2',
            'CO2_Hydro',
            'CH4_Hydro',
            'H2_Ferm2',
            'M_Ferm2']

changeables_order = ['Vmax_Ferm',
                     'Vprod_max_AltE',
                     'Vprod_max_Homo',
                     'Vprod_max_Hydro',
                     'Vprod_max_Ace',
                     'w_Ferm',
                     'w_AltE',
                     'w_Hydro', 
                     'w_Homo',
                     'w_Ace',
                     'Sensenmann',
                     'Stoch_ALtE',
                     'Kmb_Ferm',
                     'Kmh_Ferm',
                     'Kmb_AltE',
                     'Kmb_Auto',
                     'Kmb_Hydro',
                     'Fe']

parameter_units = {'Vmax_Ferm':'μmol/mg',         
                   'Vprod_max_AltE':'μmol/mg',       
                   'Vprod_max_Homo':'μmol/mg',      
                   'Vprod_max_Hydro':'μmol/mg',          
                   'Vprod_max_Ace':'μmol/mg',      
                   'w_Ferm':'mg/μmol',
                   'w_AltE':'mg/μmol',
                   'w_Hydro':'mg/μmol',
                   'w_Homo':'mg/μmol',
                   'w_Ace':'mg/μmol',
                   'Sensenmann':'-', 
                   'Stoch_ALtE':'-',
                   'Kmb_Ferm':'mg' ,
                   'Kmh_Ferm':'μmol', 
                   'Kmb_AltE':'mg',
                   'Kmb_Auto':'mg',
                   'Kmb_Hydro':'mg',
                   'Fe':'μmol',
                   'C':'μmol',
                   'M_Ac':'mg',
                   'M_Ferm':'mg',
                   'M_Fe':'mg',
                   'M_Hydro':'mg',
                   'M_Homo':'mg',
                   'CH4':'μmol',
                   'CO2':'μmol',
                   'CO2_Ac':'μmol',
                   'Acetate':'μmol',
                   'H2':'μmol',
                   'CO2_Hydro':'μmol',
                   'CH4_Hydro':'μmol',
                   'H2_Ferm2':'μmol',
                   'M_Ferm2':'mg'}

def load_data(specimen_index):
    # TODO load Corg content of samples Datei Corg_Probennummern
    
    day_column_index = specimen_index*3
    CH4_column_index = day_column_index + 1
    CO2_column_index = day_column_index + 2
    
    
    
    data_df = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Trainingdatafull.xlsx',
                       engine = 'openpyxl')

    data_df = data_df.apply(pd.to_numeric, errors='coerce') # Macht " nicht zahlen" zu Nan
    data_df = data_df.dropna() #  löscht Zeilen mit NaN
    
    data_array = data_df.values
        
    data_array[:,day_column_index] = np.around(data_array[:,day_column_index]) # rundet die Tage auf ganze tage (trotzdem floats)
    
    #erstellt Liste wo negative Werte vorhanden sind die ersetzt werden müssen
    indexListCH4 = np.where(np.isin(data_array[:,CH4_column_index],data_array[:,CH4_column_index][np.where(data_array[:,CH4_column_index]<0)]))
    indexListCO2 = np.where(np.isin(data_array[:,CO2_column_index],data_array[:,CO2_column_index][np.where(data_array[:,CO2_column_index]<0)]))

    # ersetzt Messfehler mit neg. werten mit dem vorrangegangenen wert 
    for p in range(len(indexListCH4[0])):   
        data_array[:,CH4_column_index][indexListCH4[0][p]] = data_array[:,CH4_column_index][indexListCH4[0][p]-1]
        
    for r in range(len(indexListCO2)):   
        data_array[:,CH4_column_index][indexListCO2[r]] = data_array[:,CH4_column_index][indexListCO2[r]-1] 
    
    #umwandeln der daten in dict für "fit my model"                        
    Realdata = {'measured_time': data_array[:,day_column_index],
                'CH4':data_array[:,CH4_column_index],
                'CO2':data_array[:,CO2_column_index]}

    return Realdata


def curve_wrapper(fixed_quantities):
    def merged_curves(t, *changeables):
        
        changeables_dict = dict(zip(changeables_order, changeables))
        
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
                
        pool_dict = predictor(t,initial_pool_values, model_parameters)
        
        merged = np.concatenate((pool_dict['CO2'],pool_dict['CH4']),axis = 0)
        return merged
    return merged_curves

def predictor(t, initial_pool_values, model_parameters):
    
    initial_system_state = np.array([initial_pool_values[pool_name] if pool_name in initial_pool_values else 0 
                                     for pool_name in pool_order])

    fitters = [model_parameters[parameter_name] for parameter_name in changeables_order if parameter_name in model_parameters]
    
    pool_curves = odeint(SimpleIN.Cdec, initial_system_state, t, args = (fitters,))
    pool_curves = pool_curves.transpose()
    
    pool_dict = dict(zip(pool_order, pool_curves))

    return pool_dict



def plot_my_data(Realdata, days_for_plot, pool_value_dict, specimen_index): 
    plt.close('all')

    measurement_days = Realdata['measured_time'].astype(int)

    for a, col in zip(["CH4","CO2"], ["r", "b"]):
        
        plt.figure()
        plt.plot(measurement_days, Realdata[a], col+"o", label = "Observed")
        plt.plot(pool_value_dict[a], "k-",label = "Predicted")
        plt.ylabel(a)
        plt.legend()
        
        save_path = os.path.join('C:/Users/Lara/Desktop/simple model/Figs', a +'_fit_' +str(specimen_index)+'.png')
        plt.savefig(save_path)
                
        # berechne R^2
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

    for pool_name, pool_curve in pool_value_dict.items(): 
        plt.figure()
        plt.plot(days_for_plot, pool_curve, label= pool_name)
        plt.title(pool_name)
        plt.ylabel(parameter_units[pool_name])
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/'+ pool_name +'_'+str(specimen_index)+ '.png') 


    all_CO2_contributers = dict()
    for pool_name, pool_curve in pool_value_dict.items():
        if "CO2" in pool_name:
            all_CO2_contributers[pool_name] = pool_curve
    plt.figure()
    for pool_name, pool_curve in all_CO2_contributers.items():
        plt.plot(days_for_plot, pool_curve, label= pool_name)
    plt.legend()

def fit_my_model(specimens):
    plt.close('all')
    
    for specimen_index in specimens:#and2and3and4and5and6and7and8and9:
        Realdata = load_data(specimen_index)
        # Fitting the parameters:    
        
        xlist = [int(measurement_day) for measurement_day in Realdata['measured_time']] # int der Tage an denen wir Messwerte haben
        xdata = xlist + xlist  # aneinandergehängt, weil wir die werte sowohl für CH4 als auch CO2 brauchen
        ydata = list(Realdata['CH4']) + list(Realdata['CO2']) # meine Realdata an die gefittet werden soll.
        

        measurement_days = Realdata['measured_time']
        merged_measurements = np.concatenate((Realdata['CO2'],Realdata['CH4']), axis = 0)

        # define the initial pool values
        # all pools, for which no value is speicified, will be initialized as empty
        fixed_quantities_dict = dict()        
        m_gluc = 180   # molar mass of glucose, g dw pro mol
        m_cell = 162   # molar mass of cellulose, g dw pro mol
        TOC = 0.04     # Knoblauchs Daten , g dw
        labile = 0.005 # Knoblauchs Daten, anteil von TOC einheitslos
        Cpool_init = (TOC*labile/m_gluc + TOC*(1-labile)/m_cell) * (10**6) # 0.00024679012345679013 * (10**6) mikromol pro g
        fixed_quantities_dict['C'] = Cpool_init
        fixed_quantities_dict['M_Ac'] = 0.001
        fixed_quantities_dict['M_Ferm'] = 0.2
        fixed_quantities_dict['M_Fe'] = 0.2
        fixed_quantities_dict['M_Ferm2'] = 0.2 

        # specify initial guesses and bounds for the parameters to be optimized
        initial_guess_dict = dict()         #   init    lower upper
        initial_guess_dict['Vmax_Ferm'] =       (0.1,   0.01, 0.11)
        initial_guess_dict['Vprod_max_AltE'] =  (0.3,   0.029, 1.9)
        initial_guess_dict['Vprod_max_Homo'] =  (0.133, 0.005, 1.)
        initial_guess_dict['Vprod_max_Hydro'] = (0.086, 0.03, 0.2)
        initial_guess_dict['Vprod_max_Ace'] =   (0.207, 0.05, 0.39)
        initial_guess_dict['w_Ferm'] =          (0.05,  0.03, 0.05)
        initial_guess_dict['w_AltE'] =          (0.013, 0.01, 0.05)
        initial_guess_dict['w_Hydro'] =         (0.024, 0.01, 0.05)
        initial_guess_dict['w_Homo'] =          (0.049, 0.01, 0.05)
        initial_guess_dict['w_Ace'] =           (0.04,  0.01, 0.05)
        initial_guess_dict['Sensenmann'] =      (8.33e-5, 0, 8.44e-5)
        initial_guess_dict['Stoch_ALtE'] =      (4,     1,  8)
        initial_guess_dict['Kmb_Ferm'] =        (10,    1,  10)
        initial_guess_dict['Kmh_Ferm'] =        (10,    1,  10)
        initial_guess_dict['Kmb_AltE'] =        (10,    1,  10)
        initial_guess_dict['Kmb_Auto'] =        (10,    1,  10)
        initial_guess_dict['Kmb_Hydro'] =       (10,    1,  10)
        initial_guess_dict['Fe'] =              (5.75,  2,  10)


        initial_guess = [initial_guess_dict[key][0] for key in changeables_order]
        lower_bounds = [initial_guess_dict[key][1] for key in changeables_order]
        upper_bounds = [initial_guess_dict[key][2] for key in changeables_order]
        
        merged_curves = curve_wrapper(fixed_quantities_dict)
        optimal_parameters , _ = curve_fit(merged_curves,
                                           measurement_days, 
                                           merged_measurements, #method="dogbox",
                                           p0 = initial_guess, 
                                           bounds=(lower_bounds, upper_bounds))
        optimal_parameters , _ = curve_fit(merged_curves, 
                                           measurement_days, 
                                           merged_measurements, #method="dogbox",
                                           p0 = optimal_parameters, 
                                           bounds=(lower_bounds, upper_bounds))
        changeables_optimal_dict = dict(zip(changeables_order, optimal_parameters))

        #Printing the Parameter and its value
        for parameter_name in changeables_optimal_dict.keys():
            p = changeables_optimal_dict[parameter_name]
            u = parameter_units[parameter_name]
            print("{:<18} {:6.3f} {:<10}".format(parameter_name,p,u))
        
        #Calculating the model output with optimal parameters:
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

        plot_my_data(Realdata, all_days, pool_value_dict, specimen_index)
        

if __name__ == '__main__':
    specimens = [0]#,1,2,3,4,5,6,7]
    fit_my_model(specimens)
    
