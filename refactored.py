# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:45:44 2020

@author: Lara
"""

import copy
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from scipy.integrate import odeint, solve_ivp
from scipy import stats
from scipy.optimize import curve_fit
import scipy.io as sio
import numpy as np
#from SimpleIN import Cdec # importiert das eigentliche mathematische Model 


import model
from order import pool_order, changeables_order, parameter_units
from scaling import scale, rescale

os.chdir('C:/Users/Lara/Desktop/simple model')

def load_matlab():
    
    Data = sio.loadmat('ActivityData_04062016', appendmat=True)
    
    anaerobic_indices = np.nonzero(Data['anaerobic'])[1] # die indizes im Matlab für alle anaeroben proben
    prob_array = Data['prob'] # das array mit den Probennummern
    measurement_time = Data['duration'] # das array mit allen Messtagen
    
    # erstelle die äußerste struktur, ein leeres dict für jede probe hat
    superdata = dict()
    for index_for_anaerobic_sample in anaerobic_indices:
        prob_name = prob_array[0,index_for_anaerobic_sample] #holt aus prob_array den probennamen an der indexstelle von index_for_anaerobic _sample.
        
        keys_so_far = [ k[:4]  for k in superdata.keys()] #schreibt die ersten vier zifFe3rn aller keys die wir schon haben in eine liste
        replica_count = keys_so_far.count(str(prob_name)) # zählt wie viele keys aus der liste die gleichen 4 anfangszifFe3rn haben wie die aktuelle probname
        
        replica_name = str(prob_name) + str(replica_count)
        superdata[replica_name] = dict()
        
        # erstellt das dict für die anaeroben proben
        superdata[replica_name]['measured_time'] = measurement_time[:,index_for_anaerobic_sample]
        superdata[replica_name]['CH4'] = Data['ch4'][:,index_for_anaerobic_sample]
        superdata[replica_name]['CO2'] = Data['co2'][:,index_for_anaerobic_sample]
    
        
    # fügt den Proben jeweils die relevanten Metadaten hinzu    
    Metadata = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Metadaten_all_86_cleanandneat.xlsx',engine = 'openpyxl')
    
    for meta_prob_names in Metadata['Probe']:
        for replica_name_iter in superdata.keys():
            if str(meta_prob_names) in replica_name_iter:
            #if str(meta_prob_names) == replica_name_iter[:4]:
                index = np.where(Metadata['Probe'] == meta_prob_names)[0]
                for key in Metadata.keys():
                    superdata[replica_name_iter][key] = np.array(Metadata[key][index])
     
    replica_list = list(superdata.keys())    
    
    
    for key in superdata.keys():
        #ersetzt negative messwerte mit den jeweils vorherigen messwerten, und nan mit 0 
        for index in range(len(superdata[key]['CH4'])):
            if superdata[key]['CH4'][index] < 0 :
                superdata[key]['CH4'][index] = superdata[key]['CH4'][index-1]
            if np.isnan(superdata[key]['CH4'][index]):
                superdata[key]['CH4'][index] = 0    
    
            if superdata[key]['CO2'][index] < 0 :
                superdata[key]['CO2'][index] = superdata[key]['CO2'][index-1]         
            if  np.isnan(superdata[key]['CO2'][index]):
                superdata[key]['CO2'][index] = 0
                

    # findet den ersten nan wert in measured_time und überträgt alles danach in superdata_carex        
    superdata_carex = copy.deepcopy(superdata)
    superdata_all = copy.deepcopy(superdata)
    for key in superdata.keys():
        FirstNan = np.where(np.isnan(superdata[key]['measured_time']))[0][0]
        for column_name in ['CO2','CH4','measured_time']:
            superdata_carex[key][column_name] = superdata_carex[key][column_name][(FirstNan+1):]
            superdata[key][column_name] = superdata[key][column_name][:FirstNan]
              # key entry ab wann Carex zugegeben wurde
            superdata[key]['First_Carex_index'] =     len(superdata[key]['measured_time'])
            superdata[key]['Last_non_Carex_day'] =    max(superdata[key]['measured_time'])
            
            superdata_all[key][column_name] = superdata_all[key][column_name][:]
 
    # for key in superdata:
    #     plt.figure()  
    #     plt.plot(superdata[key]['measured_time'], superdata[key]['CH4'], "r", label= "CH4")       
    #     plt.plot(superdata[key]['measured_time'], superdata[key]['CO2'],"b",  label= "CO2")
    #     plt.title(str(superdata[key]['Probe'])  + "_____" + superdata[key]['Site'] + "_____" + superdata[key]['Location']+ "_____" + str(superdata[key]['depth']))
     
    
  


# Ersetzen der alten Carexdaten mit den neuen Daten von Knoblauch bis 2021
    Carex_addition = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Carex_addition_2021_clean.xlsx',engine = 'openpyxl')
    Carex = np.array(Carex_addition)     
    #Proben = np.unique(Carex[:,0])     
   # Probenfree = np.unique(Proben[~np.isnan(Proben)])
    
   # vergleich der keys mit den Probennamen des array und ermitteln der Indizes
    for key in superdata_carex.keys():
        indices = []
        
        for i in range(len([Carex[:,0]][0])):
            if str(int([Carex[i,0]][0])) == str(key):
                indices.append(i)
        #hier werden die neuen Werte eingefügt
        superdata_carex[key]['measured_time'] = Carex[indices,1]
        superdata_carex[key]['CO2'] = Carex[indices,2]
        superdata_carex[key]['CH4'] = Carex[indices,3]
        
        superdata_carex[key]['measured_time']=[int(i) for i in superdata_carex[key]['measured_time']]# von float to int
        
    superdata_2021_all = copy.deepcopy(superdata)
    replica_list_superdata_2021_all = list(superdata_2021_all.keys()) 
    
    for key in  superdata_2021_all.keys():
        #anfügen der Carextage an die nicht carextage
          Time_to_append = [i + max(superdata_2021_all[key]['measured_time'] )  for i in superdata_carex[key]['measured_time']]

          existing_time_values = superdata_2021_all[key]['measured_time']
          time_all_values = np.concatenate([existing_time_values,Time_to_append], axis = 0)
          superdata_2021_all[key]['measured_time'] = time_all_values
          
         #anfügen der Carex CH4 messwerte an die nicht carex CH4 Messwerte
         
          CH4_to_append = [i + superdata_2021_all[key]['CH4'][-1]  for i in superdata_carex[key]['CH4']]

          existing_CH4_values = superdata_2021_all[key]['CH4']
          CH4_all_values = np.concatenate([existing_CH4_values,CH4_to_append], axis = 0)
          superdata_2021_all[key]['CH4'] = CH4_all_values
          
          #anfügen der Carex CH4 messwerte an die nicht carex CH4 Messwerte
         
          CO2_to_append = [i + superdata_2021_all[key]['CO2'][-1]  for i in superdata_carex[key]['CO2']]

          existing_CO2_values = superdata_2021_all[key]['CO2']
          CO2_all_values = np.concatenate([existing_CO2_values,CO2_to_append], axis = 0)
          superdata_2021_all[key]['CO2'] = CO2_all_values
     
          
    #erstellt einen Datensatz nur für die Kurunak daten
    superdata_Kuru = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key]['Site']=='K':
            del superdata_Kuru[key]        
    replica_list_Kuru = list(superdata_Kuru.keys()) 
    
    #erstellt einen datensatz nur für die Samoylov daten    
    superdata_Sam = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key]['Site']=='S':
            del superdata_Sam[key]        
    replica_list_Sam = list(superdata_Sam.keys())    
          
     # erstellen der Datenreihen, die vermutlich ohne Fe3 Pool sind
    Rep_ohne_Fe3 = [13690, 13531,13530,13521,13520,13742, 13741, 13740,13732,13731, 13730,13721,13720]

    superdata_ohne_Fe3 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key] in Rep_ohne_Fe3:
            del superdata_ohne_Fe3[key]        
    
     
     # erstellen der Datenreihen, die vermutlich MIT Fe3 Pool sind

    Rep_mit_Fe3 = [13510,13511,13512,13670,13671,13672,13691,13692,13700,13722, 13731, 13750,13751,13752]
    
    superdata_mit_Fe3 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key] in Rep_mit_Fe3:
            del superdata_mit_Fe3[key]        
    
    print(superdata_Kuru['13510']['First_Carex_index'])
    
                
    return superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3

# superdata sind alle Datensätze VOR dem Carexexperiment,
# superdata_carex sind die Datensätze von Christian NACH dem Carex experiment ( eigentlich alt )
# superdata_all sind alle Datensätze von Christian vor und nach Carex (auch alt)
# Superdata_2021_all sind alle Datensätze vor und Nach Carex Zugabe von Knoblauch, darauf bauen auf: 
# Superdata Kuru sind die Datensätze von Kurunak vor und nach Carex
# Superdata_sam sind die Datensätze von Samoylov vor und nach Carex
# superdata_ohne_Fe3 sind alle Datensätze deren exp  plot des CH4 vermuten lässt dass kein Fe3 vorhanden ist
# superdata_mit_Fe3 sind alle Datensätze deren exp  plot des CH4 vermuten lässt dass Fe3 vorhanden ist


# Alte Funktion für die alten Daten
# =============================================================================
# def load_data(specimen_index):
#     # TODO load Corg content of samples Datei Corg_Probennummern
#     # TODO load pH values
#     
#     day_column_index = specimen_index*3
#     CH4_column_index = day_column_index + 1
#     CO2_column_index = day_column_index + 2
#     
#     
#     
#     data_df = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Trainingdatafull.xlsx',
#                        engine = 'openpyxl')
# 
#     data_df = data_df.apply(pd.to_numeric, errors='coerce') # Macht " nicht zahlen" zu Nan
#     data_df = data_df.dropna() #  löscht Zeilen mit NaN
#     
#     data_array = data_df.values
#         
#     data_array[:,day_column_index] = np.around(data_array[:,day_column_index]) # rundet die Tage auf ganze tage (trotzdem floats)
#     
#     #erstellt Liste wo negative Werte vorhanden sind die ersetzt werden müssen
#     indexListCH4 = np.where(np.isin(data_array[:,CH4_column_index],data_array[:,CH4_column_index][np.where(data_array[:,CH4_column_index]<0)]))
#     indexListCO2 = np.where(np.isin(data_array[:,CO2_column_index],data_array[:,CO2_column_index][np.where(data_array[:,CO2_column_index]<0)]))
# 
#     # ersetzt MessFe3hler mit neg. werten mit dem vorrangegangenen wert 
#     for p in range(len(indexListCH4[0])):   
#         data_array[:,CH4_column_index][indexListCH4[0][p]] = data_array[:,CH4_column_index][indexListCH4[0][p]-1]
#         
#     for r in range(len(indexListCO2)):   
#         data_array[:,CH4_column_index][indexListCO2[r]] = data_array[:,CH4_column_index][indexListCO2[r]-1] 
#     
#     #umwandeln der daten in dict für "fit my model"                        
#     Realdata = {'measured_time': data_array[:,day_column_index],
#                 'CH4':data_array[:,CH4_column_index],
#                 'CO2':data_array[:,CO2_column_index]}
# 
#     return Realdata
# =============================================================================


# diese Funktion erstellt für curve_fit eine funktion, die die zusammengefügten
# werte von CO2 und CH4 zurückgibt.
def curve_merger_wrapper(fixed_quantities):
    
    # curve_fit ruf merged_curves_filled auf.
    #
    # merged_curves_pred funktion ruft den predictor auf und fügt die vorhergesagten verläuFe3
    # von CO2 und CH4 aneinander, damit curve_fit sie gegen die gemessenen 
    # daten von CO2 und CH4 vergleichen kann.
    #
    # changeables sind die veränderlichen parameter, die von curve_fit 
    # angepasst werden können. im ersten schritt des fit-algorithmus sind 
    # diese die initiual_guesses
    def merged_curves_pred(t, *changeables):
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
        for key, value in fixed_quantities.items():
            if not key in pool_order:
                model_parameters[key] = value
        
        # rescale the parameters
        initial_pool_values = rescale(initial_pool_values)
        model_parameters = rescale(model_parameters)
        
        pool_dict = predictor(t,initial_pool_values, model_parameters)
        
        # nimm CO2 und CH4 und pack sie zusammen zu einer einzigen kurve.
        predicted_CO2_and_CH4 = np.concatenate((pool_dict['CO2'],pool_dict['CH4']),axis = 0)
        return predicted_CO2_and_CH4
    
    # diese ZEile gehört zu curve_wrapper. zurückgegeben wird hier die funktion
    # merged_curves
    return merged_curves_pred

def least_squares_error(changeables_array, fixed_quantities_dict, measured_data_dict):  
   # print('')
   # print('ls call')
    # in den folgenden beiden blöcken werden die variablen und fixen
    # parameter neu geordnet, weil predictor() die parameter getrennt haben 
    # will, je nach dem, ob sie initiale pool-pegel sind oder andere modell-parameter
    #
    # dazu gehen wir alle pools durch, und prüFe3n, ob sich der entsprechende 
    # startwert im dict changeables_array, oder im dict fixed_quantities_dict befindet.

    changeables_dict = dict(zip(changeables_order, changeables_array))
    
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
    initial_pool_dict = rescale(initial_pool_dict)
    model_parameters_dict = rescale(model_parameters_dict)
    
    #predictor_start = time.time()
    y_predicted_dict = predictor(measured_data_dict['measured_time'],
                                initial_pool_dict, 
                                model_parameters_dict)
   # predictor_time = time.time() - predictor_start
    #print('predictor time', predictor_time)
    CO2_predicted = y_predicted_dict['CO2']
    CH4_predicted = y_predicted_dict['CH4']
    
    CO2_measured = measured_data_dict['CO2']
    CH4_measured = measured_data_dict['CH4']
    
    error_CO2 = CO2_predicted - CO2_measured
    error_CH4 = CH4_predicted - CH4_measured
    
    # ist es wichtiger an CO2 oder an CH4 gut zu fitten. (je höher desto wichtiger)
    weight_CO2 = 1.
    weight_CH4 = 1.
    sum_of_squared_residuals = np.sum(weight_CO2*error_CO2**2 + weight_CH4*error_CH4**2) 

    return sum_of_squared_residuals 


def predictor(t, initial_pool_values, model_parameters):
    print_on_call = True
    print_only_changeables = True
    
    # gibt mir meine optimierten werte zurück ? 
    if print_on_call:
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
    # werden zum startzeitpunkt (t=0) zu null gesetzt.
    initial_pool_list = list()
    for pool_name in pool_order:
        if pool_name in initial_pool_values:
            initial_pool_list.append(initial_pool_values[pool_name])
        else:
            initial_pool_list.append(0)    
    
    initial_system_state = np.array(initial_pool_list)
    
    cdec = model.Cdec_wrapper(model_parameters) 
    pool_curves = odeint(cdec, 
                         initial_system_state, 
                         t)
    pool_curves = pool_curves.transpose()
    
    pool_dict = dict(zip(pool_order, pool_curves))

    return pool_dict


def plot_my_data(Realdata, days_for_plot, pool_value_dict, specimen_index): 
    #plt.close('all')

    measurement_days = Realdata['measured_time'].astype(int)
    print(pool_value_dict.keys())
    
    for a, col in zip(["CH4","CO2"], ["r", "b"]):
        
        plt.figure()
        plt.plot(measurement_days, Realdata[a], col+"o", label = "Observed")
        plt.plot(pool_value_dict[a], "k-",label = "Predicted")
        plt.ylabel(a)
        plt.legend(Realdata['Probe'])
        
        
# =============================================================================
         #plots aller pools
    for key in pool_value_dict.keys():
         plt.figure()
         plt.plot(range(0, max(measurement_days)+1), pool_value_dict[key], label = [key])
         plt.ylabel([key])
         plt.legend(Realdata['Probe'])
 #=============================================================================
        
    
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

    # for pool_name, pool_curve in pool_value_dict.items(): 
    #     plt.figure()
    #     plt.plot(days_for_plot, pool_curve, label= pool_name)
    #     plt.title(pool_name)
    #     try:
    #         plt.ylabel(parameter_units[pool_name])
    #     except:
    #         pass
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/'+ pool_name +'_'+str(specimen_index)+ '.png') 


    all_CO2_contributers = dict()
    for pool_name, pool_curve in pool_value_dict.items():
        if "CO2" in pool_name:
            all_CO2_contributers[pool_name] = pool_curve
    plt.figure()
    for pool_name, pool_curve in all_CO2_contributers.items():
        plt.plot(days_for_plot, pool_curve, label= pool_name)
    plt.legend()


def fit_my_model(specimens, Site, opt):
    plt.close('all')
    
    superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()
    
    for specimen_index in specimens: 
        if Site == "S":
            Realdata = superdata_Sam[replica_list_Sam[specimen_index]]
            print(replica_list_Sam[specimen_index])
        elif Site == "K":
            Realdata = superdata_Kuru[replica_list_Kuru[specimen_index]]   
            print(replica_list_Kuru[specimen_index])
        elif Site == "all":
            Realdata = superdata_2021_all[replica_list_superdata_2021_all[specimen_index]]
            print(replica_list_superdata_2021_all[specimen_index])
            
        
        # define the initial pool values
        # all pools, for which no value is speicified, will be initialized as empty
        fixed_quantities_dict = dict()        
        m_gluc = 180   # molar mass of glucose, g dw pro mol
        m_cell = 162   # molar mass of cellulose, g dw pro mol
        TOC = 0.04     # Knoblauchs Daten , g dw
        labile = 0.005 # Knoblauchs Daten, anteil von TOC einheitslos
        Cpool_init = (TOC*labile/m_gluc + TOC*(1-labile)/m_cell) * (10**6) # 0.00024679012345679013 * (10**6) mikromol pro g
        fixed_quantities_dict['C'] = Cpool_init
        #fixed_quantities_dict['M_Ac'] = 0.2 # superdata_Kuru = 0.001
        fixed_quantities_dict['M_Ferm'] = 0.2
        fixed_quantities_dict['M_Fe3'] = 0.2
        fixed_quantities_dict['M_Ferm2'] = 0.2 

        # specify initial guesses and bounds for the parameters to be optimized
        initial_guess_dict = dict()         #   init    lower upper
        initial_guess_dict['Vmax_Ferm'] =       (0.07,   0.01,0.11)  
        initial_guess_dict['Vmax_Fe3'] =         (0.3,   0.029, 1.9)
        initial_guess_dict['Vmax_Homo'] =       (0.133, 0.005, 1.)
        initial_guess_dict['Vmax_Hydro'] =      (0.086, 0.03, 0.2)
        initial_guess_dict['Vmax_Ac'] =         (0.207, 0.05, 0.39)
        initial_guess_dict['w_Ferm'] =          (0.05,  0.03, 0.05)
        initial_guess_dict['w_Fe3'] =            (0.013, 0.01, 0.05)
        initial_guess_dict['w_Hydro'] =         (0.024, 0.01, 0.05)
        initial_guess_dict['w_Homo'] =          (0.049, 0.01, 0.05)
        initial_guess_dict['w_Ac'] =            (0.04,  0.01, 0.05)
        initial_guess_dict['Sensenmann'] =      (8.33e-5, 0, 8.44e-5)
        initial_guess_dict['Stoch_Fe3'] =        (4,     1,  8)
        initial_guess_dict['Kmb_Ferm'] =        (10,    1,  10)
        initial_guess_dict['Kmh_Ferm'] =        (10,    1,  10)
        initial_guess_dict['Kmb_Fe3'] =          (10,    1,  10)
        initial_guess_dict['Kmb_Ac'] =          (10,    1,  10)
        initial_guess_dict['Kmb_Hydro'] =       (10,    1,  10)
        initial_guess_dict['Fe3'] =              (10.75,  2,  100)
        initial_guess_dict['M_Ac'] =            (0.05,  0.001,  0.2)
        initial_guess_dict['KmA_Ferm']=         (0.05, 0.05, 20)   # Diese Boundaries müssen anhander Acetatekurven angepasst werden


        # scale parameters to allow better optimization
        initial_guess_dict = scale(initial_guess_dict)
        
        
        # vorbereitung des startpunktes für die optimierung und der bounds 
        # für curve_fit
        initial_guess_array = [initial_guess_dict[key][0] for key in changeables_order]
        lower_bounds = [initial_guess_dict[key][1] for key in changeables_order]
        upper_bounds = [initial_guess_dict[key][2] for key in changeables_order]

        # schreibe die werte, die zwar im initial_guess_dict definiert sind, 
        # aber trotzdem während der optimierung Fe3stgehalten werden sollen, 
        # in das dict für die nicht mit-optimierten werte.
        fixed_quantities_dict.update({k:v[0] for k,v in initial_guess_dict.items() if not k in changeables_order})

        # gemessene daten aus Realdata laden und zusammenhängen, um CO2 und CH4
        # gleichzeitig zu fitten.
        measurement_days = Realdata['measured_time']
        merged_measurements = np.concatenate((Realdata['CO2'],Realdata['CH4']), axis = 0)

        # curve_wrapper hängt die vorhersagen für CO2 und CH4 hintereinander
        # damit curve_fit gleichzeitig gegen die gemessenen daten für CH4 und 
        # CO2 fitten kann.
        #print(merged_measurements)
        #print(type(measurement_days))
        #print(initial_guess_array)
        #print(lower_bounds)
        #print(upper_bounds)
        if opt =='Curve':
            print('using curve_fit')
            
            curve_merger = curve_merger_wrapper(fixed_quantities_dict)
            optimal_parameters , _ = curve_fit(curve_merger,
                                               measurement_days, 
                                               merged_measurements, #method="dogbox",
                                               p0 = initial_guess_array, 
                                               bounds=(lower_bounds, upper_bounds))
            # optimal_parameters , _ = curve_fit(curve_merger, 
            #                                    measurement_days, 
            #                                    merged_measurements, #method="dogbox",
            #                                    p0 = optimal_parameters, 
            #                                    bounds=(lower_bounds, upper_bounds))
            changeables_optimal_dict = dict(zip(changeables_order, optimal_parameters))
            changeables_optimal_dict = rescale(changeables_optimal_dict)
            print(optimal_parameters)
        elif opt == "Min":
             print('using minimize')

             initial_guess_bounds = list(zip(lower_bounds, upper_bounds))
             optimization_result = scipy.optimize.minimize(least_squares_error, 
                                                      initial_guess_array,
                                                      args = (fixed_quantities_dict, Realdata),
                                                      bounds= initial_guess_bounds,
                                                      options = {'maxiter':200000,
                                                                'disp':True})
             print(optimization_result)
             changeables_optimal_array = optimization_result.x
             changeables_optimal_dict = dict(zip(changeables_order,changeables_optimal_array))
             changeables_optimal_dict = rescale(changeables_optimal_dict)

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
        
        print( pool_value_dict.keys())
        

if __name__ == '__main__':

# =============================================================================
#     Wahlmöglichkeiten um das Model lauFe3n zu lassen 
#
#     Site = K , alle Kurunak daten
#     Site = S, alle Samoylov daten
#     Site = all, alle Daten
#     
#     opt = Curve , benutzt Curve fit
#     opt = Min, benutzt Minimize
#
#     specimens muss die Länge des Datensatzes sein
#     
# =============================================================================

#%%
    #specimenlistmax = [*range(0, 35, 1)]
    #specimenlist_all= [i for i in range(34)]
    #specimenlist_Sam= [i for i in range(18)]
    #specimenlist_Kuru= [i for i in range(16)]
    #specimens = [specimenlist_Sam][0]#,1,2,3,4,5,6,7]
    #specimens = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    specimens = [10]
    
    fit_my_model(specimens, Site = "all", opt = 'Min')
    




