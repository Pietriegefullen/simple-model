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


cdec_pool_order = ['C',
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

fitters_order = ['Vmax_Ferm',
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

parameter_units = {'Vmax_Ferm':"μmol/mg",         
                   'Vprod_max_AltE':"μmol/mg",       
                   'Vprod_max_Homo':"μmol/mg",      
                   'Vprod_max_Hydro':"μmol/mg",          
                   'Vprod_max_Ace':"μmol/mg",      
                   'w_Ferm':"mg/μmol",
                   'w_AltE':"mg/μmol",
                   'w_Hydro':"mg/μmol",
                   'w_Homo':"mg/μmol",
                   'w_Ace':"mg/μmol",
                   'Sensenmann':"-", 
                   'Stoch_ALtE':"-",
                   'Kmb_Ferm':"mg" ,
                   'Kmh_Ferm':"μmol", 
                   'Kmb_AltE':"mg",
                   'Kmb_Auto':"mg",
                   'Kmb_Hydro':"mg",
                   'Fe':"μmol"}

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
        
        changeables_dict = dict(zip(fitters_order, changeables))
        
        initial_pool_values = dict()
        for pool_name in cdec_pool_order:
            if pool_name in fixed_quantities:
                initial_pool_values[pool_name] = fixed_quantities[pool_name]
            elif pool_name in changeables_dict:
                initial_pool_values[pool_name] = changeables_dict[pool_name]
            else:
                initial_pool_values[pool_name] = 0
        
        model_parameters = dict()
        for key, value in changeables_dict.items():
            if not key in cdec_pool_order:
                model_parameters[key] = value
                
        pool_dict = predictor(t,initial_pool_values, model_parameters)
        
        merged = np.concatenate((pool_dict['CO2'],pool_dict['CH4']),axis = 0)
        return merged
    return merged_curves

def predictor(t, initial_pool_values, model_parameters):
    
    initial_system_state = np.array([initial_pool_values[pool_name] if pool_name in initial_pool_values else 0 
                                     for pool_name in cdec_pool_order])

    fitters = [model_parameters[parameter_name] for parameter_name in fitters_order if parameter_name in model_parameters]
    
    pool_curves = odeint(SimpleIN.Cdec, initial_system_state, t, args = (fitters,))
    pool_curves = pool_curves.transpose()
    
    pool_dict = dict(zip(cdec_pool_order, pool_curves))

    return pool_dict


def simplesolve(xtage,*Fitters ):
    
    # Festgelegte Initialwerte
    # Die Berechnung hier für Cpool ist noch für zwei Pools gedacht. Momentan nicht relevant
    m_gluc = 180   # molar mass of glucose, g dw pro mol
    m_cell = 162   # molar mass of cellulose, g dw pro mol
    TOC = 0.04     # Knoblauchs Daten , g dw
    labile = 0.005 # Knoblauchs Daten, anteil von TOC einheitslos
    Cpool_init = (TOC*labile/m_gluc + TOC*(1-labile)/m_cell) * (10**6) # 0.00024679012345679013 * (10**6) mikromol pro g
    
    # die Werte sind alle in mikroMol pro gram Trockengewicht
    AltE_init = Fitters[-1] # 7 #146       # 146 cf. Yao, Conrad 1999 (Hyun2017 sagt werte um 100 sind extrem hoch), Philben: 6 jeweils Mikromol pro g dw
    Fitters = Fitters[:-1]
    M_A_CH4_init = 0.001     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Ferm_init = 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    M_AltE_init= 0.2         # Monteux 2020, mg Mikrobielles C pro g dw
    deltaH2_Hydro = 0        # für plot in multifit
    M_H_CH4_init = 0.2     # Monteux 2020, mg Mikrobielles C pro g dw
    M_Homo_init = 0.2        # Monteux 2020, mg Mikrobielles C pro g dw
    CH4_init = 0
    CO2_init = 0
    Acetate_init = 0        # Philben wert ca 3, in mikromol pro g 
    AceCO2_init = 0
    H2_init = 0
    deltaH2_Homo_init = 0
    CO2_Hydro_init= 0
    CH4_Hydro_init = 0
    H2_Ferm2_init = 0
    M_Ferm2_init = 0.2
    
    xs = np.linspace(0,int(max(xtage)), int(max(xtage)+1))
    y0 = Cpool_init,  AltE_init,  M_A_CH4_init,  M_Ferm_init,  M_AltE_init, M_H_CH4_init,  M_Homo_init,  CH4_init,  CO2_init,  AceCO2_init, Acetate_init,  H2_init,  CO2_Hydro_init,  CH4_Hydro_init, H2_Ferm2_init,  M_Ferm2_init
    
    ys = odeint(SimpleIN.Cdec, y0, xs, args = (Fitters,))
    
    Cpool, AltEpool, M_A_CH4, M_Ferm, M_AltE, M_H_CH4, M_Homo, CH4, CO2, AceCO2, Acetate, H2,  CO2_Hydro, CH4_Hydro, H2_Ferm2, M_Ferm2 = zip(*ys)
    
   # AltEpool[0::int(1/step)]

   
    CH4 = [CH4[int(i)] for i in xtage]
    CO2 = [CO2[int(i)] for i in xtage]


    return CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_Ferm, M_AltE, H2, M_H_CH4, M_Homo, AltE_init, deltaH2_Hydro,[0 for _ in range(len(CO2))], CO2_Hydro, CH4_Hydro,H2_Ferm2
    

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
        p = len(fitters_order) # number of explanatory terms
        adjusted_r_squared = 1 - (1-r_squared)*(n-1)/(n-p-1)
        
#        _ ,_ , R2, _ , _ = stats.linregress(y = Realdata[a], x=pool_value_dict[a+'_pool'][measurement_days])
#        print("R2 for "+a+" is", R2)
        print("r2 for "+a+" is", r_squared)
        print("r2adj  "+a+" is", adjusted_r_squared)

    # TODO: measured_time counts from 0 or 1?



    for pool_name, pool_curve in pool_value_dict.items(): 
        plt.figure()
        plt.plot(days_for_plot, pool_curve, label= pool_name)
        plt.title(pool_name)
        plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/'+ pool_name +'_'+str(specimen_index)+ '.png') 
        #plt.ylabel('mg Mikrobielles C pro g dw')

    
    # for _ in pool_value_dict.keys():
    #     if "CO2" in  pool_value_dict.keys():
    #         all_CO2_contributers = dict()
    #         all_CO2_contributers= pool_value_dict.items()
    #         wanted_keys = pool_value_dict.keys()
            
    # all_CO2_contributers2 = dict((k, pool_value_dict[k]) for k in wanted_keys if k in pool_value_dict)
            
    
    all_CO2_contributers = dict()
    for pool_name, pool_curve in pool_value_dict.items():
        if "CO2" in pool_name:
            all_CO2_contributers[pool_name] = pool_curve
    plt.figure()
    for pool_name, pool_curve in all_CO2_contributers.items():
        plt.plot(days_for_plot, pool_curve, label= pool_name)
    plt.legend()

if __name__ == '__main__':
    plt.close('all')

    specimens = [0]#,1,2,3,4,5,6,7]
    
    for specimen_index in specimens:#and2and3and4and5and6and7and8and9:
        Realdata = load_data(specimen_index)
        # Fitting the parameters:    
        
        xlist = [int(measurement_day) for measurement_day in Realdata['measured_time']] # int der Tage an denen wir Messwerte haben
        xdata = xlist + xlist  # aneinandergehängt, weil wir die werte sowohl für CH4 als auch CO2 brauchen
        ydata = list(Realdata['CH4']) + list(Realdata['CO2']) # meine Realdata an die gefittet werden soll.
        

        measurement_days = Realdata['measured_time']
        merged_measurements = np.concatenate((Realdata['CO2'],Realdata['CH4']), axis = 0)

        fixed_quantities_dict = dict()
        fixed_quantities_dict['M_Ac'] = 0.001
        fixed_quantities_dict['M_Ferm'] = 0.2
        fixed_quantities_dict['M_Fe'] = 0.2
        fixed_quantities_dict['M_Ferm2'] = 0.2    

        initial_guess_dict = dict()
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


        initial_guess = [initial_guess_dict[key][0] for key in fitters_order]
        lower_bounds = [initial_guess_dict[key][1] for key in fitters_order]
        upper_bounds = [initial_guess_dict[key][2] for key in fitters_order]
        
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
        changeables_optimal_dict = dict(zip(fitters_order, optimal_parameters))
        
        # TODO: print song for comparison
        song =  ["0.5",                "",            "0.15",           "0.15",            "0.5",         "",         "",     "",     "",         "",    "",  "",     "" ,"", "", "" , "", "", ""]
        # SOIL_DENSITY = 1.3        
               # print("{:<18} {:6.3f} {:<10} ({:})".format(n,p*SOIL_DENSITY,"mol/m^3", s))

        #Printing the Parameter and its value
        for parameter_name in changeables_optimal_dict.keys():
            p = changeables_optimal_dict[parameter_name]
            u = parameter_units[parameter_name]
            print("{:<18} {:6.3f} {:<10}".format(parameter_name,p,u))
        

        #Calculating the model output with optimal parameters:
        initial_pool_dict = dict()
        optimal_model_parameters_dict = dict()
        
        for key in changeables_optimal_dict:
            if key in cdec_pool_order:
                initial_pool_dict[key] = changeables_optimal_dict[key]
            else:
                optimal_model_parameters_dict[key] = changeables_optimal_dict[key]
                
        for key in fixed_quantities_dict:
            if key in cdec_pool_order:
                initial_pool_dict[key] = fixed_quantities_dict[key]
            else:
                optimal_model_parameters_dict[key] = fixed_quantities_dict[key]

        all_days = np.arange(max(Realdata['measured_time'])+1)
        pool_value_dict = predictor(all_days, initial_pool_dict, optimal_model_parameters_dict)
            
        # for key, curve in pool_value_dict.items():
        #     plt.figure()
        #     plt.plot(all_days, curve)
        #     if key in Realdata:
        #         plt.plot(Realdata['measured_time'], Realdata[key], 'o')
            
        
        # plot_my_data(Realdata, all_days, pool_value_dict, specimen_index)
        print(xlist)
    #   #  Fitters_opt = optimal_parameters
    #     #CCH4CO2opt = simplefun(xlist,  *optimal_parameters)
    #     CCH4CO2opt = simplesolve(xlist,  *optimal_parameters)
    #     CCH4CO2optList = list(CCH4CO2opt[0]) + list(CCH4CO2opt[1])
          
    #     Mol_nach_Pa(n= max(CCH4CO2opt[9]))
        
    # #### Observed and Predicted für CH4 und CO2 geplotted
    # #%%
    #     data_len = len(xlist)
    #     for x, a, col in zip([0,data_len],["CH4","CO2"], ["r", "b"]):
    #         plt.figure()
    #         plt.plot(xlist,[ydata[i] for i in range(x,data_len+x)],col+"o", label = "Observed")
    #         plt.plot(xlist,[CCH4CO2optList[i] for i in range(x,data_len+x)],"k-",label = "Predicted")
    #         plt.ylabel(a)
    #         plt.legend()
    #         save_path = os.path.join('C:/Users/Lara/Desktop/simple model/Figs', a +'_'+str(specimen_index)+'.png')
    #         plt.savefig(save_path)
    #     _ ,_ , R2, _ , _ = stats.linregress(CCH4CO2optList, y=ydata)
    #     print("R2 is", R2)
        
         
    # #%%      
    #     # return order CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_CH4,    M_CO2, M_AltE, M_Hyd
    #     #              CH4, CO2, AltEpool, AceCO2, Acetate, Cpool, M_A_CH4, M_Ferm, M_AltE, H2, M_H_CH4, M_Homo, AltE_init, deltaCO2_Hydro
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[2], label='AltEpool')
    #     plt.title('AltEpool')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/AltEpool.png')    
    #     plt.ylabel('mg Mikrobielles C pro g dw')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[3], label='CO2 aus Acetate')
    #     plt.title('CO2 aus Acetate')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/CO2ausAcetate.png')
    #     plt.ylabel('mg Mikrobielles C pro g dw')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[4], label='Acetate')
    #     plt.title('Acetate')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Acetate.png')
    #     plt.ylabel('μmol')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[5], label='Cpool')
    #     plt.title('Cpool')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Cpool.png')
    #     plt.ylabel('mikromol pro g dw')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[6], label='Acetoclastische_Mikroben')
    #     plt.title('Acetoclastische_Mikroben')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Acetoclastic_Mikroben.png')  
    #     plt.ylabel('mg Mikrobielles C pro g dw')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[7], label='Fermentierer')
    #     plt.title('Fermentierer')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Fermentation_Mikroben.png')  
    #     plt.ylabel('mg Mikrobielles C pro g dw')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[8], label=' AltE_Microben')
    #     plt.title(' AltE_Microben')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/AltE_Mikroben.png') 
    #     plt.ylabel('mg Mikrobielles C pro g dw')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[9], label=' H2') #typically in the nanomolar range (10 to 180 nM, correponding to about 8 to 140 ppmv in the gas phase [Table conrad
    #     plt.title('H2')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/H2.png') 
    #     plt.ylabel('μmol')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[10], label=' Hydro_Microben')
    #     plt.title(' Hydro_Microben')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Hydro_Mikroben.png')    
    #     plt.ylabel('mg Mikrobielles C pro g dw')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[11], label=' Homo_Microben')
    #     plt.title(' Homo_Microben')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/Homo_Mikroben.png') 
    #     plt.ylabel('mg Mikrobielles C pro g dw')
        
    #     #plt.figure()
    #     #plt.plot( CCH4CO2opt[12], label=' AltE_init')
    #     #plt.title(' AltE_init')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[13], label='deltaH2_Hydro')
    #     plt.title('deltaH2_Hydro')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/deltaH2_Hydro.png')
    #     plt.ylabel('μmol')        
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[14], label='deltaH2_Homo')
    #     plt.title('deltaH2_Homo')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/deltaH2_Homo.png')
    #     plt.ylabel('μmol')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[15], label='CO2 aus Hydro')
    #     plt.title('CO2 aus Hydro')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/CO2_aus_Hydro.png')
    #     plt.ylabel('μmol')    
    
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[16], label='CH4 aus Hydro')
    #     plt.title('CH4 aus Hydro')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/CH4_aus_Hydro.png') 
    #     plt.ylabel('μmol')
        
    #     plt.figure()
    #     plt.plot( CCH4CO2opt[17], label='H2 aus Ferm2')
    #     plt.title('H2 aus Ferm2')
    #     plt.savefig('C:/Users/Lara/Desktop/simple model/Figs/H2_aus_Ferm.png')  
    #     plt.ylabel('μmol')
    
    
        
        # Goodness of fit with R^2 automatic
            
    
    ################ ab hier ist es uninteressant #################################   
    #    #--calculating R^2 by hand-----------------------------------------------------
    #    # Aufteilen von CCH4CO2opt in die jeweiligen pools und Berechnung des mean
    #    
    #    
    #    CH4obsmean = np.mean(CCH4CO2optList[0:int(len(ydata)/2)])
    #    CO2obsmean = np.mean(CCH4CO2optList[int(len(ydata)/2):len(ydata)])
    #    allobsmean = np.mean(CCH4CO2optList)
    #    
    #    #---- calculating total sum of squares for each pool and global----------------
    #    
    #    SStotCH4 = []
    #    SStotCO2 = []
    #    SStotall = []
    #    
    #    for i in range(len(xlist)):
    #        SStotCH4.append((ydata[i]- CH4obsmean)**2)
    #        SStotCO2.append((ydata[i+len(xlist)]- CO2obsmean)**2)
    #        
    #    # list comprehension für den globalen datensatz    
    #    SStotall.append([((ydata[i] - allobsmean)**2) for i in range(len(ydata))]) 
    #    
    #    
    #    SStotCH4 = np.sum(SStotCH4)
    #    SStotCO2 = np.sum(SStotCO2)
    #    SStotall = np.sum(SStotall)
    #    
    #    #%%
    #    #---- calculating total sum of residuals for each pool and global-------------- 
    #    
    #    SSresCH4 = []
    #    SSresCO2 = []
    #    SSresall = []
    #    
    #    
    #    for i in range(len(xlist)):
    #        SSresCH4.append((ydata[i]- CCH4CO2optList[i])**2)
    #        SSresCO2.append((ydata[i+len(xlist)]- CCH4CO2optList[i+len(xlist)])**2)
    #        
    #    # list comprehension für den globalen datensatz        
    #    SSresall.append([((ydata[i]- CCH4CO2optList[i])**2) for i in range(len(ydata))])    
    #       
    #    
    #    SSresCH4 = np.sum(SSresCH4)
    #    SSresCO2 = np.sum(SSresCO2)
    #    SSresall = np.sum(SSresall)
    #    
    #    
    #    R2CH4 = 1 - (SStotCH4 / SSresCH4)
    #    print("R2CH4 is", R2CH4)
    #    R2CO2 = 1 - (SStotCO2 / SSresCO2)
    #    print("R2CO2 is", R2CO2)
    #    R2all = 1 - (SStotall / SSresall)
    #    print("R2all is", R2all)
    #    
    plt.show()   
        
        
        
    
    
