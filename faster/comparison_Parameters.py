import os
import numpy as np
import matplotlib.pyplot as plt

import argparse
from datetime import datetime
import json

import itertools

from data import load_matlab


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates
from scipy.stats import pearsonr

import predict
import pathways
import plot
import data
import optimizer
from ORDER import POOL_ORDER
import OPTIMIZATION_PARAMETERS
import USER_VARIABLES



import seaborn as sns

Initials = OPTIMIZATION_PARAMETERS.get_initial_guesses()

Bounds= {}

for key in Initials.keys():
    Bounds[key] = Initials[key][1:]

print(Bounds)

    






replica_list_No_CH4, superdata_No_CH4_vor_Impfung,superdata_No_CH4_nach_Impfung,superdata_after_No_CH4,superdata_bevor_No_CH4, superdata_No_CH4,  superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()




#holt meine ergebnisse aus dem ordner
def load_model_parameters(file):

    file_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file)
    
    with open(file_path, 'r') as pf:
        model_parameters = json.load(pf)
          
   
    return model_parameters







def load_singel_file(filename):

    all_losses =  os.path.join(USER_VARIABLES.LOG_DIRECTORY, filename + '.txt')
    
    
    with open(all_losses, 'r') as p:
        all_losses = json.load(p)
       
        
    return all_losses




    


def corrfunc(x, y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    r, _ = pearsonr(x, y)
    ax = ax or plt.gca()
    ax.annotate(f'ρ = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)







# holt aus den datensätzen die Parameterwerte und die parameternamen
def load_and_compare(file, extended_output = None):

    splitparts = file.split('_')
    specimen_index = str(splitparts[splitparts.index('specimen')+1])
    
    model_parameters = load_model_parameters(file)
        
    specimen_dict = {specimen_index: model_parameters}
    
    for k,v in specimen_dict.items():
        columns = list(v.keys())

    
    return specimen_dict, columns








def plot_preperations():
    
    # läd die losses  
    all_losses = load_singel_file('all_losses')
    
    Microbe_groups = ['Ferm', 'Fe3', 'Homo', 'Hydro', '_Ac']

    # ruft load and compare auf     
    all_speciemen_dict = {}
    all_parameter_list = []


    #läd alle parametersätze rein
    for file in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
         
         if file.endswith(".json") : 
             
             specimen_dict, columns = load_and_compare(file)
      
             all_speciemen_dict.update(specimen_dict)
    
             
    # schmeist die nicht parameter raus
    get_out_indices =[]
    for specimen_index, optimal_parameter_dict in all_speciemen_dict.items():
        row = [specimen_index] + list(optimal_parameter_dict.values())
        all_parameter_list.append(row)
        
    columns_list = ['specimen_index'] + columns + ['site']
    get_out = ['temperature', 'C', 'DOC', 'pH', 'weight', 'H2O','water']
    get_out_indices = list(reversed(sorted([columns_list.index(name) for name in get_out])))
        
    for column_index in get_out_indices:
        del columns_list[column_index]
        for sp in range(len(all_parameter_list)):
            del all_parameter_list[sp][column_index]
    
    
    # damit ich nach Site aufteilen kann 
    for key in superdata_2021_all.keys():
        for r, row in enumerate(all_parameter_list):
            specimen = row[0]
            if str(key) == str(specimen):
                site = superdata_2021_all[key]['Site']
                site_int = 0 if site == 'K' else 1
                all_parameter_list[r].append(site_int)
        
    
    # wir brauchen die datan als dataframe für den parallel plot         
    specimen_df = pd.DataFrame(all_parameter_list, 
                                columns = columns_list)

    specimen_df_Kuru = specimen_df[specimen_df['site'] == 0]
    
    specimen_df_Sam  = specimen_df[specimen_df['site'] == 1]
    
 
    
    ## PARAMETER COMPARISM nach Gruppen 
    # alle Km Werte
    Km_Values = [col for col in specimen_df.columns if 'Km' in col or 'Inh' in col] 
    
    # alle Mikroben gruppen
    Microbes_Values = [col for col in specimen_df.columns if 'M_' in col ] 
    
    # alle Vmax werte
    Vmax_Values = [col for col in specimen_df.columns if 'Vmax_' in col ] 
    
    # alle CUE werte
    CUE_Values = [col for col in specimen_df.columns if 'CUE' in col ] 
    
    # Km_Values.append('Fe3')
    # Km_Values.append( 'Sensenmann')
    # Km_Values.append('Acetate')
    # Km_Values.append('site')
    # Para_Values_all = [col for col in specimen_df.columns ] 
    # Para_Values_no_Km =  [p for p in Para_Values_all if not p in Km_Values]

    # Para_df = specimen_df[Para_Values_no_Km]
    # Para_df_Kuru = specimen_df_Kuru[Para_Values_no_Km]
    # Para_df_Sam = specimen_df_Sam[Para_Values_no_Km]
    
    # Para_df = specimen_df[Para_Values_no_Km]
    
    
    
    
    
    
    return  specimen_df, specimen_df_Kuru, specimen_df_Sam, Km_Values, Microbes_Values, Vmax_Values, CUE_Values, columns_list
    
def boxplot():   
   
    import itertools 
   
    (specimen_df,
     specimen_df_Kuru,
     specimen_df_Sam,
     Km_Values,
     Microbes_Values,
     Vmax_Values,
     CUE_Values,
     columns_list) = plot_preperations()
   
    sns.set(style="whitegrid") 
    

    Selected_Params = [ Km_Values, Microbes_Values, Vmax_Values, CUE_Values]


    specimen_df_Kuru = specimen_df_Kuru.rename({column: column + '_K' for column in specimen_df_Kuru.columns}, axis = 1)
    specimen_df_Sam = specimen_df_Sam.rename({column: column + '_S' for column in specimen_df_Sam.columns}, axis = 1)
    
    
    Datasets = [specimen_df, specimen_df_Kuru, specimen_df_Sam]
            
    plot_data = pd.concat([specimen_df_Kuru,specimen_df_Sam], axis = 1)
    for z in Selected_Params:
        base_colors = itertools.cycle(sns.color_palette('deep'))
        modified_colors = itertools.cycle(sns.color_palette('pastel'))        

        plot_columns = [col + site for col in z for site in ('_K', '_S')]
        plot_colors = []
        for column in plot_columns:
            plot_colors.append(next(base_colors) if column.endswith('_K') else next(modified_colors))
            
        ax = sns.boxplot( data=plot_data[plot_columns], showfliers = False, palette = plot_colors)
        plt.xticks(rotation=90)
        
        plt.update_layout(boxmode='group')
        plt.show()      
  
    
  
    
     
    
    
    
    
def PAIRWISE_Correlation_plot():
    # plots für die verschiedenen Mikrobengruppen    
    (specimen_df,
     specimen_df_Kuru,
     specimen_df_Sam,
     Km_Values,
     Microbes_Values,
     Vmax_Values,
     CUE_Values,
     columns_list) = plot_preperations()
    

    # for microbe in Microbe_groups:
    #     print('')
    #     print(microbe)
        
    #     plot_list = [s for s in columns_list if microbe in s]  
    #     print(plot_list)
        
    #     Microbes_sam ={}
    #     Microbes_kuru = {}
                    
    #     Microbes_sam[microbe] = specimen_df_Sam[plot_list]
    #     Microbes_kuru[microbe] = specimen_df_Kuru[plot_list]
    
    #     plt.figure()
    #     sns.pairplot(Microbes_sam[microbe],kind="reg")
    #     plo_sam = sns.pairplot(Microbes_sam[microbe],kind="reg")
    #     plo_sam.fig.suptitle(microbe + ' Sam')
    #     plo_sam.map_lower(corrfunc)
    #     print('plotting Sam ', microbe)
        
    #     plt.figure()
    #     plo_kuru = sns.pairplot(Microbes_kuru[microbe],kind="reg")
    #     plo_kuru.fig.suptitle(microbe + ' Kuru')
    #     plo_kuru.map_lower(corrfunc)
    #     print('plotting Kuru', microbe)
    #     print('')
        
        
    # plt.show()

 
    
    
    
    # normalisieren der Daten für PARALLEL PLOTTING
    minmaxdict = {}
    for col in columns_list[1:-1]:
        col_min = min(specimen_df[col])
        col_max = max(specimen_df[col])
        minmaxdict.update({col: [col_min, col_max]})
        specimen_df[col] = (specimen_df[col] - col_min) /(col_max - col_min)
        
    minmaxdict_Kuru = {}
    for col in columns_list[1:-1]:
        col_min = min(specimen_df_Kuru[col])
        col_max = max(specimen_df_Kuru[col])
        minmaxdict_Kuru.update({col: [col_min, col_max]})
        specimen_df_Kuru[col] = (specimen_df_Kuru[col] - col_min) /(col_max - col_min)    
        
    minmaxdict_Sam = {}
    for col in columns_list[1:-1]:
        col_min = min(specimen_df_Sam[col])
        col_max = max(specimen_df_Sam[col])
        minmaxdict_Sam.update({col: [col_min, col_max]})
        specimen_df_Sam[col] = (specimen_df_Sam[col] - col_min) /(col_max - col_min)    
    
    
    
    # # plottet alle parameter
    # columns_list.remove('specimen_index')
    # ax = parallel_coordinates(specimen_df,"site", cols= columns_list, color = ['r','b'])
    # #ax.get_legend().remove()
    # plt.xticks(rotation=90)
    # plt.show()
    
    # # plottet Kuru 
    # ax = parallel_coordinates(specimen_df_Kuru,"site", cols= columns_list, color = ['r'])
    # #ax.get_legend().remove()
    # plt.xticks(rotation=90)
    # plt.show()
    
    # # plottet Sam        
    # ax = parallel_coordinates(specimen_df_Sam,"site", cols= columns_list, color = ['b'])
    # #ax.get_legend().remove()
    # plt.xticks(rotation=90)
    # plt.show()
    
     
    ## plots für die verschiedenen Mikrobengruppen   
     
    # Fe3_columns_list = [s for s in columns_list if "Fe3" in s]
    
    # #erstellt alle möglichen kombinationen für den plot
    # # permut_Fe3 = itertools.permutations(Fe3_columns_list)
    # # for j in list(permut_Fe3):
    # #     j = list(j)
    
    # Fe3_columns_list = ['M_Fe3', 'Vmax_Fe3', 'Km_Fe3_Fe3', 'Km_Fe3_Acetate', 'Fe3', 'CUE_Fe3']

    
    
    # plt.figure()
    # # axs[0].plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))
    
    # # axs[1].plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))
    # # axs[0].plot(parallel_coordinates(specimen_df_Sam,"site", cols= Fe3_columns_list,color = ['b']))
    
    # plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))
    
    # plt.figure()
    
    # plot(parallel_coordinates(specimen_df_Sam,"site", cols= Fe3_columns_list,color = ['b']))
    
    
    

    
    
    # #ax.get_legend().remove()
    # plt.xticks(rotation=90)
    # plt.show() 
                
    
    # Ferm_columns_list = [s for s in columns_list if "Ferm" in s]
    # permut_Ferm = itertools.permutations(Ferm_columns_list)
    # for j in list(permut_Ferm):
    #     j = list(j)
    
    #     #Ferm_columns_list = ['M_Ferm', 'Vmax_help_Ferm', 'Vmax_Ferm', 'Kmb_help_Ferm', 'Km_Ferm', 'Inhibition_Ferm', 'CUE_Ferm']
    #     ax = parallel_coordinates(specimen_df,"site", cols= j,color = ['r','b'])
    #     ax.get_legend().remove()
    #     plt.xticks(rotation=90)
    #     plt.show() 
            
    
    # Homo_columns_list = [s for s in columns_list if "Homo" in s]
    # ax = parallel_coordinates(specimen_df,"site", cols= Homo_columns_list,color = ['r','b'])
    # ax.get_legend().remove()
    # plt.xticks(rotation=90)
    # plt.show() 
            
    
    
    # Hydro_columns_list = [s for s in columns_list if "Hydro" in s]
    # ax = parallel_coordinates(specimen_df,"site", cols= Hydro_columns_list,color = ['r','b'])
    # ax.get_legend().remove()
    # plt.xticks(rotation=90)
    # plt.show() 
            
          
    # Ac_columns_list = [s for s in columns_list if "Ac" in s]
    # Ac_columns_list.remove('')
    # ax = parallel_coordinates(specimen_df,"site", cols= Ac_columns_list,color = ['r','b'])
    # ax.get_legend().remove()
    # plt.xticks(rotation=90)
    # plt.show() 













































                   
    
# def parallel_plot():
    
#     # läd die losses
    
#     all_losses = load_singel_file('all_losses')
    
    

#     # ruft load and compare auf     
#     all_speciemen_dict = {}
#     all_parameter_list = []


#     #läd alle parametersätze rein
#     for file in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
         
#          if file.endswith(".json") : 
             
#              specimen_dict, columns = load_and_compare(file)
      
#              all_speciemen_dict.update(specimen_dict)
    
             
#     # schmeist die nicht parameter raus
#     get_out_indices =[]
#     for specimen_index, optimal_parameter_dict in all_speciemen_dict.items():
#         row = [specimen_index] + list(optimal_parameter_dict.values())
#         all_parameter_list.append(row)
        
#     columns_list = ['specimen_index'] + columns + ['site']
#     get_out = ['temperature', 'C', 'DOC', 'pH', 'weight', 'H2O','water']
#     get_out_indices = list(reversed(sorted([columns_list.index(name) for name in get_out])))
        
#     for column_index in get_out_indices:
#         del columns_list[column_index]
#         for sp in range(len(all_parameter_list)):
#             del all_parameter_list[sp][column_index]
    
    
    
#     # def extract(lst):
#     #     return [item[0] for item in lst]
    
#     # damit ich nach Site aufteilen kann 
#     for key in superdata_2021_all.keys():
#         for r, row in enumerate(all_parameter_list):
#             specimen = row[0]
#             if str(key) == str(specimen):
#                 site = superdata_2021_all[key]['Site']
#                 site_int = 0 if site == 'K' else 1
#                 all_parameter_list[r].append(site_int)
        
    
#     # wir brauchen die datan als dataframe für den parallel plot         
#     specimen_df = pd.DataFrame(all_parameter_list, 
#                                 columns = columns_list)

#     specimen_df_Kuru = specimen_df[specimen_df['site'] == 0]
    
#     specimen_df_Sam  = specimen_df[specimen_df['site'] == 1]
    
 
    
#     ## PARAMETER COMPARISM PLOT
#     Km_Values = [col for col in specimen_df.columns if 'Km' in col or 'Inh' in col] 
#     K_df = specimen_df [Km_Values]
   
#     Km_Values_Kuru = [col for col in specimen_df_Kuru.columns if 'Km' in col or 'Inh' in col] 
#     K_df_Kuru = specimen_df_Kuru [Km_Values_Kuru]
    
#     Km_Values_Sam = [col for col in specimen_df_Sam.columns if 'Km' in col or 'Inh' in col] 
#     K_df_Sam = specimen_df_Sam [Km_Values_Sam]
   
   
#     sns.set(style="whitegrid") 
    
#     #ax = sns.boxplot( data=K_df, showfliers = False)

#     ax = sns.swarmplot( data=K_df_Kuru, color="r", s=2)
#     ax = sns.swarmplot( data=K_df_Sam, color="b", s=2)
#     plt.xticks(rotation=90)
#     plt.show()
        
  
#     Km_Values.append('Fe3')
#     Km_Values.append( 'Sensenmann')
#     Km_Values.append('Acetate')
#     Km_Values.append('site')
#     Para_Values_all = [col for col in specimen_df.columns ] 
#     Para_Values_no_Km =  [p for p in Para_Values_all if not p in Km_Values]
#     print('Km_Values_no_LM',Para_Values_no_Km )
    

#     Para_df = specimen_df[Para_Values_no_Km]
#     Para_df_Kuru = specimen_df_Kuru[Para_Values_no_Km]
#     Para_df_Sam = specimen_df_Sam[Para_Values_no_Km]
   
    
#     ax = sns.swarmplot( data=Para_df_Kuru, color="r" , s=2)
#     plt.xticks(rotation=90)
    
    
#     plt.figure()
#     ax = sns.swarmplot( data=Para_df_Sam, color="b", s=2)
#     plt.xticks(rotation=90)
    
#     plt.figure()
#     ax = sns.boxplot( data=Para_df_Kuru, showfliers = False)
#     ax.set_ylim([0, 5])
#     ax.set_title('Kuru')
#     plt.xticks(rotation=90)
    
#     plt.figure()
#     ax = sns.boxplot( data=Para_df_Sam, showfliers = False)
#     ax.set_ylim([0, 5])
#     ax.set_title('Sam')
#     plt.xticks(rotation=90)
    
#     plt.show()
    
    
#     plt.figure()
    
    
    
    
    
    
#     # PAIRWISE CORRELATION PLOTS
#     # plots für die verschiedenen Mikrobengruppen    
#     Microbe_groups = ['Ferm', 'Fe3', 'Homo', 'Hydro', '_Ac']
    

#     # for microbe in Microbe_groups:
#     #     print('')
#     #     print(microbe)
        
#     #     plot_list = [s for s in columns_list if microbe in s]  
#     #     print(plot_list)
        
#     #     Microbes_sam ={}
#     #     Microbes_kuru = {}
                    
#     #     Microbes_sam[microbe] = specimen_df_Sam[plot_list]
#     #     Microbes_kuru[microbe] = specimen_df_Kuru[plot_list]
    
#     #     plt.figure()
#     #     sns.pairplot(Microbes_sam[microbe],kind="reg")
#     #     plo_sam = sns.pairplot(Microbes_sam[microbe],kind="reg")
#     #     plo_sam.fig.suptitle(microbe + ' Sam')
#     #     plo_sam.map_lower(corrfunc)
#     #     print('plotting Sam ', microbe)
        
#     #     plt.figure()
#     #     plo_kuru = sns.pairplot(Microbes_kuru[microbe],kind="reg")
#     #     plo_kuru.fig.suptitle(microbe + ' Kuru')
#     #     plo_kuru.map_lower(corrfunc)
#     #     print('plotting Kuru', microbe)
#     #     print('')
        
        
#     # plt.show()

 
    
    
    
#     # normalisieren der Daten für PARALLEL PLOTTING
#     minmaxdict = {}
#     for col in columns_list[1:-1]:
#         col_min = min(specimen_df[col])
#         col_max = max(specimen_df[col])
#         minmaxdict.update({col: [col_min, col_max]})
#         specimen_df[col] = (specimen_df[col] - col_min) /(col_max - col_min)
        
#     minmaxdict_Kuru = {}
#     for col in columns_list[1:-1]:
#         col_min = min(specimen_df_Kuru[col])
#         col_max = max(specimen_df_Kuru[col])
#         minmaxdict_Kuru.update({col: [col_min, col_max]})
#         specimen_df_Kuru[col] = (specimen_df_Kuru[col] - col_min) /(col_max - col_min)    
        
#     minmaxdict_Sam = {}
#     for col in columns_list[1:-1]:
#         col_min = min(specimen_df_Sam[col])
#         col_max = max(specimen_df_Sam[col])
#         minmaxdict_Sam.update({col: [col_min, col_max]})
#         specimen_df_Sam[col] = (specimen_df_Sam[col] - col_min) /(col_max - col_min)    
    
    
    
#     # # plottet alle parameter
#     # columns_list.remove('specimen_index')
#     # ax = parallel_coordinates(specimen_df,"site", cols= columns_list, color = ['r','b'])
#     # #ax.get_legend().remove()
#     # plt.xticks(rotation=90)
#     # plt.show()
    
#     # # plottet Kuru 
#     # ax = parallel_coordinates(specimen_df_Kuru,"site", cols= columns_list, color = ['r'])
#     # #ax.get_legend().remove()
#     # plt.xticks(rotation=90)
#     # plt.show()
    
#     # # plottet Sam        
#     # ax = parallel_coordinates(specimen_df_Sam,"site", cols= columns_list, color = ['b'])
#     # #ax.get_legend().remove()
#     # plt.xticks(rotation=90)
#     # plt.show()
    
     
#     ## plots für die verschiedenen Mikrobengruppen   
     
#     # Fe3_columns_list = [s for s in columns_list if "Fe3" in s]
    
#     # #erstellt alle möglichen kombinationen für den plot
#     # # permut_Fe3 = itertools.permutations(Fe3_columns_list)
#     # # for j in list(permut_Fe3):
#     # #     j = list(j)
    
#     # Fe3_columns_list = ['M_Fe3', 'Vmax_Fe3', 'Km_Fe3_Fe3', 'Km_Fe3_Acetate', 'Fe3', 'CUE_Fe3']

    
    
#     # plt.figure()
#     # # axs[0].plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))
    
#     # # axs[1].plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))
#     # # axs[0].plot(parallel_coordinates(specimen_df_Sam,"site", cols= Fe3_columns_list,color = ['b']))
    
#     # plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))
    
#     # plt.figure()
    
#     # plot(parallel_coordinates(specimen_df_Sam,"site", cols= Fe3_columns_list,color = ['b']))
    
    
    

    
    
#     # #ax.get_legend().remove()
#     # plt.xticks(rotation=90)
#     # plt.show() 
                
    
#     # Ferm_columns_list = [s for s in columns_list if "Ferm" in s]
#     # permut_Ferm = itertools.permutations(Ferm_columns_list)
#     # for j in list(permut_Ferm):
#     #     j = list(j)
    
#     #     #Ferm_columns_list = ['M_Ferm', 'Vmax_help_Ferm', 'Vmax_Ferm', 'Kmb_help_Ferm', 'Km_Ferm', 'Inhibition_Ferm', 'CUE_Ferm']
#     #     ax = parallel_coordinates(specimen_df,"site", cols= j,color = ['r','b'])
#     #     ax.get_legend().remove()
#     #     plt.xticks(rotation=90)
#     #     plt.show() 
            
    
#     # Homo_columns_list = [s for s in columns_list if "Homo" in s]
#     # ax = parallel_coordinates(specimen_df,"site", cols= Homo_columns_list,color = ['r','b'])
#     # ax.get_legend().remove()
#     # plt.xticks(rotation=90)
#     # plt.show() 
            
    
    
#     # Hydro_columns_list = [s for s in columns_list if "Hydro" in s]
#     # ax = parallel_coordinates(specimen_df,"site", cols= Hydro_columns_list,color = ['r','b'])
#     # ax.get_legend().remove()
#     # plt.xticks(rotation=90)
#     # plt.show() 
            
          
#     # Ac_columns_list = [s for s in columns_list if "Ac" in s]
#     # Ac_columns_list.remove('')
#     # ax = parallel_coordinates(specimen_df,"site", cols= Ac_columns_list,color = ['r','b'])
#     # ax.get_legend().remove()
#     # plt.xticks(rotation=90)
#     # plt.show() 
            
          
if __name__ == '__main__':
    
    
    
    
    
    #parallel_plot()
    boxplot()  

    PAIRWISE_Correlation_plot()      
   
    all_losses = load_singel_file('all_losses')
        
    # all_losses_hist = [item[1] for item in all_losses]
    
    # all_losses_hist2 = sorted(np.log10(all_losses_hist))
    
    
    # plt.hist(all_losses_hist2, bins = 50)
    # plt.figure()
    # plt.plot(all_losses_hist2)    
    # plt.show()
    
    
            
      
        
      
        
      
        
      
        
      
        
      
        
      
 
         
         
         
         
         
         
         
         
         
         
         