# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 13:27:41 2021

@author: Lara
"""

import numpy as np
import scipy
import pandas as pd
import matplotlib as plt
# variablen die optimiert werden
changeables_order = ['Vmax_Ferm', 
                   'w_Ferm',
                   'Kmb_Ferm',
                   'Kmh_Ferm',
                   'Vmax_Fe',
                   'w_Fe',
                   'Stoch_Fe',
                   'Kmb_Fe',                  
                   'Vmax_Ac',
                   'w_Ac',
                   'Kmb_Ac',
                   'Vmax_Homo',
                   'w_Homo',
                   'Vmax_Hydro',
                   'w_Hydro',
                   'Kmb_Hydro',
                   'Sensenmann',
                   'Fe_pool'
                   ]

#pools
pool_order = ['C_pool',
              'Fe_pool',
              'Acetate_pool',
              'M_Ac_pool',
              'M_Ferm_pool',
              'M_Fe_pool',
              'M_Hydro_pool',
              'M_Homo_pool',
              'M_Ferm2_pool',
              'CH4_pool',
              'CH4_Hydro_pool',
              'CO2_pool',
              'CO2_Ac_pool',
              'CO2_Hydro_pool',
              'CO2_Ferm_pool',
              'CO2_Fe_pool',
              'CO2_Homo_pool',
              'H2_pool',
              'H2_Homo_pool',
              'H2_Ferm2_pool',
              'H2_Hydro_pool'
              ]



def load_data(specimen_index):
    day_column_index = specimen_index*3
    CH4_column_index = day_column_index + 1
    CO2_column_index = day_column_index + 2
    
    
    data_df = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Trainingdatafull.xlsx',
                       engine = 'openpyxl')

    data_df = data_df.apply(pd.to_numeric, errors='coerce') # Macht " nicht zahlen" zu Nan
    data_df = data_df.dropna() #  löscht Zeilen mit NaN
    
    data_array = data_df.values
        
    data_array[:,m] = np.around(data_array[:,day_column_index]) # rundet die Tage auf ganze tage (trotzdem floats)
    
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
                'CO2':data_array[:,CH4_column_index],
                'CH4':data_array[:,CO2_column_index]}
   

def Cdec(time, system_state_array, model_parameter_dict):# fehlen hier die parameter? 
    system_dict = dict(zip(pool_order, system_state_array))
    
    C_pool = system_dict['C_pool']
    Fe_pool = system_dict['Fe_pool']
    Acetate_pool = system_dict['Acetate_pool']
    M_Ac_pool = system_dict['M_Ac_pool']
    M_Ferm_pool = system_dict['M_Ferm_pool']
    M_Fe_pool = system_dict['M_Fe_pool']
    M_Hydro_pool = system_dict['M_Hydro_pool']
    M_Homo_pool = system_dict['M_Homo_pool']
    M_Ferm2_pool = system_dict['M_Ferm2_pool']
    CH4_pool = system_dict['CH4_pool']
    CH4_Hydro_pool = system_dict['Ch4_Hydro_pool']
    CO2_pool = system_dict['CO2_pool']
    CO2_Ac_pool = system_dict['CO2_Ac_pool']
    CO2_Hydro_pool = system_dict['CO2_Hydro_pool']
    CO2_Ferm_pool = system_dict['CO2_Ferm_pool']
    CO2_Fe_pool = system_dict['CO2_Fe_pool']
    CO2_Homo_pool = system_dict['CO2_Homo_pool']   
    H2_pool = system_dict['H2_pool']
    H2_Homo_pool = system_dict['H2_Homo_pool']
    H2_Ferm2_pool = system_dict['H2_Ferm2_pool']
    H2_Hydro_pool = system_dict['H2_Hydro_pool']

    


    # FERM FERM FERM 
    deltaM_Ferm, deltaCpool, _, Tot_Ferm =   Fermenters(M_Ferm, C_pool, Acetate, Vmax_Ferm, w_Ferm, Sensenmann, Kmb_Ferm, Kmh_Ferm)
    deltaAcetate_Ferm = - (deltaCpool)
    # Produziert:
    deltaCO2_Ferm = deltaAcetate_Ferm * 0.5 
    deltaH2_Ferm = deltaAcetate_Ferm *(1/6) 
    
    #FERM2 FERM2 FERM2
    #deltaM_Ferm2, deltaCpool2, _, Tot_Ferm2 =   Fermenters(M_Ferm2, Cpool, Acetate, Vmax_Ferm, w_Ferm, Sensenmann, Kmb_Ferm, Kmh_Ferm)
    # Produziert:
    #deltaH2_Ferm2 =  -deltaCpool2 * 1.2 # 1.2 aus cabrol2017microbial
    deltaH2_Ferm2 = 0
    deltaM_Ferm2 = 0
    
    # FE FE FE FE
    deltaM_Fe, deltaAcetate_Fe, deltaFe_pool, Tot_Fe =   Fe(M_Fe, Acetate, Fe_pool, Stoch_Fe, Vmax_Fe, w_FE, Sensenmann, Kmb_Fe)
    deltaCO2_Fe = - deltaAcetate_Fe * 2 # pro 1 Acetate entstehen zwei CO2
    deltaH2_Fe = - deltaAcetate_Fe #* 0 # was ist der wirkliche Wert? 
    
    # ACETO ACETO ACETO   
    deltaM_Ac, deltaAcetate_Ac,Tot_Ac =   Acetoclast(M_Ac, Acetate, w_Ac, Vmax_Ac, Sensenmann, Kmb_Ac)
    deltaCH4_Ac = - deltaAcetate_Ac * 0.5 # pro mol Acetate entsteht 0.5 Mol CH4
    deltaCO2_Ac = - deltaAcetate_Ac * 0.5 # pro mol Acetate entsteht 0.5 Mol CO2
    
    # HYDRO HYDRO HYDRO
    deltaM_Hydro, deltaCO2_Hydro, deltaH2_Hydro, Tot_Hydro =   Hydrotrophes(M_Hydro, CO2, H2, w_Hydro, Vmax_Hydro, Sensenmann, Kmb_Hydro)
    deltaCH4_Hydro = - deltaCO2_Hydro # pro mol CO2 entsteht 1 mol CH4 
    
    # HOMO HOMO HOMO
    deltaM_Homo, deltaCO2_Homo ,deltaH2_Homo, Tot_Homo =   Homo(M_Homo, CO2, H2, w_Homo, Vmax_Homo, Sensenmann)
    deltaAcetate_Homo  = - deltaH2_Homo * 4 # aus 4 mol H2 wird ein mol Acetate

    
    m_C = 12.01*1e-3 # molar mass of carbon
    changes_dict                    = dict()
    changes_dict['C_pool']          = -min(Cpool,-deltaCpool) +  (Tot_Hydro + Tot_Ac + Tot_Fe + Tot_Ferm + Tot_Homo)/m_C
    changes_dict['Fe_pool']         = deltaFe_pool
    changes_dict['Acetate_pool']    = -min(Acetate, -(deltaAcetate_Ferm + deltaAcetate_Fe + deltaAcetate_Ac + deltaAcetate_Homo ))  
    changes_dict['M_Ac_pool']       = deltaM_Ac
    changes_dict['M_Ferm_pool']     = deltaM_Ferm
    changes_dict['M_Fe_pool']       = deltaM_Fe
    changes_dict['M_Hydro_pool']    = deltaM_Hydro
    changes_dict['M_Homo_pool']     = deltaM_Homo
    changes_dict['M_Ferm2_pool']    = deltaM_Ferm2
    changes_dict['CH4_pool']        = deltaCH4_Ac + deltaCH4_Hydro 
    changes_dict['Ch4_Hydro_pool']  = deltaCH4_Hydro
    changes_dict['CO2_pool']        = -min(CO2, -(deltaCO2_Ferm + deltaCO2_Fe + deltaCO2_Ac + deltaCO2_Hydro + deltaCO2_Homo))  
    changes_dict['CO2_Ac_pool']     = deltaCO2_Ac
    changes_dict['CO2_Hydro_pool']  = deltaCO2_Hydro
    changes_dict['CO2_Ferm_pool']   = deltaCO2_Ferm
    changes_dict['CO2_Fe_pool']     = deltaCO2_Fe
    changes_dict['CO2_Homo_pool']   = deltaCO2_Homo
    changes_dict['H2_pool']         = -min(H2, -(deltaH2_Ferm + deltaH2_Fe + deltaH2_Hydro + deltaH2_Ferm2  + deltaH2_Homo ))  
    changes_dict['H2_Homo_pool']    = deltaH2_Homo
    changes_dict['H2_Ferm2_pool']   = deltaH2_Ferm2
    changes_dict['H2_Hydro_pool']   = deltaH2_Hydro
    
    
    changes_to_system_list = [changes_dict[key] for key in pool_order]
    changes_to_system_state_array = np.array(changes_to_system_list)
    return changes_to_system_state_array


def predictor(t_out, initial_pool_values_dict, model_parameters_dict):
    
    # dict to array with entries ordered as defined in pool_order
    initial_pool_values_list = [initial_pool_values_dict[key] for key in pool_order] # diese zeile sorgt für die richtige reihenfolge für Cdec
    initial_pool_values_array = np.array(initial_pool_values_list)
    
    pool_values_at_time_points = odeint(Cdec, 
                                        initial_pool_values_array, 
                                        t_out, 
                                        args = model_parameters_dict,
                                        tfirst = True)
    pool_values_at_time_points = np.transpose(pool_values_at_time_points)
    # VOR dem transponieren ist der output von odeint ein 2D array mit
    # einer Zeile für jeden Zeitpunkt,
    # und einer Spalte für jeden pool
    # NACH dem transponieren umgekehrt, 
    # eine Zeile für jeden Pool, und eine Spalte für jeden Zeitpunkt

    pool_value_dict = dict(zip(pool_order, pool_values_at_time_points))
    
    return pool_value_dict


def least_squares_error(changeables_array, fixed_quantities_dict, measured_data_dict):  
    
    # in den folgenden beiden blöcken werden die variablen und fixen
    # parameter neu geordnet, weil predictor() die parameter getrennt haben 
    # will, je nach dem, ob sie initiale pool-pegel sind oder andere modell-parameter
    #
    # dazu gehen wir alle pools durch, und prüfen, ob sich der entsprechende 
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
    
    y_predicted_dict = predictor(measured_data_dict['measured_time'],
                                initial_pool_dict, 
                                model_parameters_dict)
    
    CO2_predicted = y_predicted_dict['CO2']
    CH4_predicted = y_predicted_dict['CH4']
    
    CO2_measured = measured_data_dict['CO2']
    CH4_measured = measured_data_dict['CH4']
    
    error_CO2 = CO2_predicted - CO2_measured
    error_CH4 = CH4_predicted - CH4_measured
    
    sum_of_squared_residuals = np.sum(error_CO2**2 + error_CH4**2) 

    return sum_of_squared_residuals 


def fit_my_model(Realdata):
      
    
    # Festgelegte Initialwerte
    # Die Berechnung hier für Cpool ist noch für zwei Pools gedacht. Momentan nicht relevant
    m_gluc = 180   # molar mass of glucose, g dw pro mol
    m_cell = 162   # molar mass of cellulose, g dw pro mol
    TOC = 0.04     # Knoblauchs Daten , g dw
    labile = 0.005 # Knoblauchs Daten, anteil von TOC einheitslos
    # die Werte sind alle in mikroMol pro gram Trockengewicht
    
    # orange (fixe parameter)
    fixed_quantities_dict  = dict()
    fixed_quantities_dict['C_pool'] = (TOC*labile/m_gluc + TOC*(1-labile)/m_cell) * (10**6) # 0.00024679012345679013 * (10**6) mikromol pro g
    fixed_quantities_dict['M_Ac_pool'] = 0.001     # Monteux 2020, mg Mikrobielles C pro g dw
    fixed_quantities_dict['M_Ferm_pool'] = 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    fixed_quantities_dict['M_Hydro_pool'] = 0.2     # Monteux 2020, mg Mikrobielles C pro g dw
    fixed_quantities_dict['M_Homo_pool'] = 0.2        # Monteux 2020, mg Mikrobielles C pro g dw
    fixed_quantities_dict['M_Fe_pool']= 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    fixed_quantities_dict['M_Ferm2_pool'] = 0.2
    fixed_quantities_dict['CH4_pool'] = 0
    fixed_quantities_dict['CH4_Hydro_pool'] = 0
    fixed_quantities_dict['CO2_pool'] = 0
    fixed_quantities_dict['CO2_Hydro_pool']= 0
    fixed_quantities_dict['CO2_Ferm_pool'] = 0
    fixed_quantities_dict['CO2_Alte_pool'] = 0
    fixed_quantities_dict['CO2_Homo_pool'] = 0
    fixed_quantities_dict['Acetate_pool'] = 0        # Philben wert ca 3, in mikromol pro g 
    fixed_quantities_dict['AcCO2_pool'] = 0
    fixed_quantities_dict['H2_pool'] = 0
    fixed_quantities_dict['H2_Homo_pool'] = 0
    fixed_quantities_dict['H2_Ferm2_pool'] = 0
    fixed_quantities_dict['H2_Hydro_pool'] = 0        # für plot in multifit

    
    # unterkringelt (zu optimierende parameter)
    changeables_initial_guess_dict = dict()
    changeables_initial_guess_dict['Vmax_Ferm'] = 0.1  # Vmax Ferm 0.01  roden2003competition wert ist 17 !
    changeables_initial_guess_dict['w_Ferm']    = 0.05    # w Ferm 0.04
    changeables_initial_guess_dict['Kmb_Ferm']  = 10     # Kmb ferm, for inverse M-M Biomass 10
    changeables_initial_guess_dict['Kmh_Ferm']  = 10      # Kmh ferm, for hemmung of fermenters by acetate 10
    changeables_initial_guess_dict['Vmax_Fe']   = 0.3    # Vmax Fe 0.9
    changeables_initial_guess_dict['w_Fe']      = 0.013    # w Fe 0.02
    changeables_initial_guess_dict['Stoch_Fe']  = 4        # Stoch AltE 7  #  philben2020anaerobic nehmen einen Wert von 4 für Fe3 an
    changeables_initial_guess_dict['Kmb_Fe' ]   = 10      # Kmb Fe 10            
    changeables_initial_guess_dict['Vmax_Ac']   = 0.207    # Vmax Ac 0.15   roden2003competition wert ist 15 !
    changeables_initial_guess_dict['w_Ac']      = 0.04    # w Ac 0.03
    changeables_initial_guess_dict['Kmb_Ac']    = 10     # Kmb Ac 10
    changeables_initial_guess_dict['Vmax_Homo'] = 0.133   # Vmax Homo 0.05
    changeables_initial_guess_dict['w_Homo']    = 0.049    # w Homo 0.05
    changeables_initial_guess_dict['Vmax_Hydro']= 0.086    # Vmax Hydro 0.05
    changeables_initial_guess_dict['w_Hydro']   = 0.024    # w Hydro 0.02
    changeables_initial_guess_dict['Kmb_Hydro'] = 10
    changeables_initial_guess_dict['Sensenmann']= 0.0000833   # Sensenmann 0.0001   delattre2020thermodynamic nimmt 8.33*10^-4 h^-1für alle mikroben
    changeables_initial_guess_dict['Fe_pool']   = 5.75
    
    
    changeables_initial_guess_list = [changeables_initial_guess_dict[key] for key in changeables_order]
    changeables_initial_guess_array = np.array(changeables_initial_guess_list)
    changeables_optimal_array = scipy.optimize.minimize(least_squares_error, 
                                                 changeables_initial_guess_array,
                                                 args = (fixed_quantities_dict, Realdata))
    
    changeables_optimal_dict = dict(zip(changeables_order,changeables_optimal_array))
    
    
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
            
    return initial_pool_dict, optimal_model_parameters_dict
    

# plots
def plot_my_data(specimen_index): 
    
    Realdata = load_data(specimen_index)
    
    optimal_model_parameters_dict, initial_pool_dict =  fit_my_model(Realdata)

    pool_value_dict = predictor(np.linespace(0,max(Realdata['measured_time'])), initial_pool_dict, optimal_model_parameters_dict) 

    data_len = len(Realdata)
    for x, a, col in zip([0,data_len],["CH4","CO2"], ["r", "b"]):
        plt.figure()
        plt.plot(xlist,[ydata[i] for i in range(x,data_len+x)],col+"o", label = "Observed")
        plt.plot(xlist,[CCH4CO2optList[i] for i in range(x,data_len+x)],"k-",label = "Predicted")
        plt.ylabel(a)
        plt.legend()
        save_path = os.path.join('C:/Users/Lara/Desktop/simple model/Figs', a +'_'+str(m)+'.png')
        plt.savefig(save_path)
    _ ,_ , R2, _ , _ = stats.linregress(CCH4CO2optList, y=ydata)
    print("R2 is", R2)



    for pool_name in pool_value_dict: 
        plt.figure()
        plt.plot(np.linespace(0,max(Realdata['measured_time'])), pool_value_dict[pool_name], label= pool_name)
        plt.title(pool_name)
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/'+ pool_name + '.png') 
        #plt.ylabel('mg Mikrobielles C pro g dw')







# xopt = [3.456, 7.45]

if __name__ == '__main__':
    
    Data1 = [0]
    Data1and2 = [0,1]
    Data1and2and3 =[0,1,2]
    Data1and2and3and4and5and6 =[0,1,2,3,4,5]
    Data1and2and3and4and5and6and7and8and9 =[0,1,2,3,4,5,6,7,8]
    
    
    for m in Data1:#and2and3and4and5and6and7and8and9:
    
        Realdata = load_data(m)    
        
        fit_my_model(Realdata)

