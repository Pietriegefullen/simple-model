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

import predict
import pathways
import plot
import data
import optimizer
from ORDER import POOL_ORDER
import OPTIMIZATION_PARAMETERS
import USER_VARIABLES

from sklearn.datasets import load_iris

iris = load_iris()
#print(iris)


replica_list_No_CH4, superdata_No_CH4_vor_Impfung,superdata_No_CH4_nach_Impfung,superdata_after_No_CH4,superdata_bevor_No_CH4, superdata_No_CH4,  superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()





directory = os.fsencode(USER_VARIABLES.LOG_DIRECTORY)

#holt meine ergebnisse aus dem ordner
def load_model_parameters(file):

    file_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file)

    with open(file_path, 'r') as pf:
        model_parameters = json.load(pf)

    return model_parameters


# holt aus den datensätzen die Parameterwerte und die parameternamen
def load_and_compare(file, extended_output = None):

    splitparts = file.split('_')
    specimen_index = str(splitparts[splitparts.index('specimen')+1])
    
    model_parameters = load_model_parameters(file)
        
    specimen_dict = {specimen_index: model_parameters}
    
    for k,v in specimen_dict.items():
        columns = list(v.keys())

    
    return specimen_dict, columns
                   
    
    
# ruft load and compare auf     
all_speciemen_dict = {}
all_parameter_list = []
normalization_dict = {}
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     
     if filename.endswith(".json") : 
         
         specimen_dict, columns = load_and_compare(filename)

         all_speciemen_dict.update(specimen_dict)

         
# schmeist die nicht parameter raus
get_out_indices =[]
for specimen_index, optimal_parameter_dict in all_speciemen_dict.items():
    row = [specimen_index] + list(optimal_parameter_dict.values())
    all_parameter_list.append(row)
    #print(row[0])

columns_list = ['specimen_index'] + columns + ['site']
get_out = ['temperature', 'C', 'DOC', 'pH', 'weight', 'H2O','water']
get_out_indices = list(reversed(sorted([columns_list.index(name) for name in get_out])))
    
for column_index in get_out_indices:
    del columns_list[column_index]
    for sp in range(len(all_parameter_list)):
        del all_parameter_list[sp][column_index]



def Extract(lst):
    return [item[0] for item in lst]

ext = Extract(all_parameter_list)

for key in superdata_2021_all.keys():
    for r, row in enumerate(all_parameter_list):
        #print(type(specimen))
        specimen = row[0]
        if str(key) == str(specimen):
            site = superdata_2021_all[key]['Site']
            #print(site)
            site_int = 0 if site == 'K' else 1
            all_parameter_list[r].append(site_int)
    

# wir brauchen die datan als dataframe           
specimen_df = pd.DataFrame(all_parameter_list, 
                            columns = columns_list)

specimen_df_Kuru = specimen_df[specimen_df['site'] == 0]

specimen_df_Sam  = specimen_df[specimen_df['site'] == 1]




# normalisieren der Daten
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


# ax = parallel_coordinates(specimen_df_Kuru,"site", cols= columns_list, color = ['r'])
# #ax.get_legend().remove()
# plt.xticks(rotation=90)
# plt.show()

         
# ax = parallel_coordinates(specimen_df_Sam,"site", cols= columns_list, color = ['b'])
# #ax.get_legend().remove()
# plt.xticks(rotation=90)
# plt.show()

 
# plots für die verschiedenen Mikrobengruppen        
Fe3_columns_list = [s for s in columns_list if "Fe3" in s]

#erstellt alle möglichen kombinationen für den plot
permut_Fe3 = itertools.permutations(Fe3_columns_list)
for j in list(permut_Fe3):
    j = list(j)

Fe3_columns_list = ['M_Fe3', 'Vmax_Fe3', 'Km_Fe3_Fe3', 'Km_Fe3_Acetate', 'Fe3', 'CUE_Fe3']

# fig, axs = plt.subplots(2)

# x = [1,2,3,4]
# y = [1,2,3,4]

# #axs[0].plot(x,y, 'r')
# #axs[1].plot(x,y, 'b')



# axs[0].plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))

# axs[1].plot(parallel_coordinates(specimen_df_Kuru,"site", cols= Fe3_columns_list,color = ['r']))
# axs[0].plot(parallel_coordinates(specimen_df_Sam,"site", cols= Fe3_columns_list,color = ['b']))


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
        
      
        
      
        
      
        
      
        
      
        
      
        
      
        
      
        
      
 
         
         
         
         
         
         
         
         
         
         
         