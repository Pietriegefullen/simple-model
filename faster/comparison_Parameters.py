import os
import numpy as np
import matplotlib.pyplot as plt

import argparse
from datetime import datetime
import json


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


directory = os.fsencode(USER_VARIABLES.LOG_DIRECTORY)

def load_model_parameters(file):

    file_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file)

    with open(file_path, 'r') as pf:
        model_parameters = json.load(pf)

    return model_parameters



def load_and_compare(file, extended_output = None):

    splitparts = file.split('_')
    specimen_index = str(splitparts[splitparts.index('specimen')+1])
   # site = splitparts[splitparts.index('site')+1]
    
    model_parameters = load_model_parameters(file)
    #print(model_parameters)
    
    # for k,v in model_parameters.items():
    #     new_dict = {specimen_index:{k,v}}
        
    specimen_dict = {specimen_index: model_parameters}
    
    for k,v in specimen_dict.items():
        columns = list(v.keys())
       # print(columns)
    
  
    
    #print('specimen_dict', specimen_dict)
    
    return specimen_dict, columns
                   
    
    
    
all_speciemen_dict = {}
all_parameter_list = []
normalization_dict = {}
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     
     if filename.endswith(".json") : 
         #print('filenmae', filename)
         
         specimen_dict, columns = load_and_compare(filename)

         all_speciemen_dict.update(specimen_dict)
         #print(all_speciemen_dict.keys())

         

         get_out_indices =[]
         for specimen_index, optimal_parameter_dict in all_speciemen_dict.items():
             row = [specimen_index] + list(optimal_parameter_dict.values())
             all_parameter_list.append(row)

columns_list = ['specimen_index'] + columns
get_out = ['temperature', 'C', 'DOC', 'pH', 'weight', 'H2O','water']
get_out_indices = list(reversed(sorted([columns_list.index(name) for name in get_out])))
    
for column_index in get_out_indices:
    del all_parameter_list[column_index]
    del columns_list[column_index]

# for name in get_out:
#     columns_list.remove(name)


# for c, name in enumerate(columns_list):
#     print(name, all_parameter_list[0][c])
    
  
            
         
            
specimen_df = pd.DataFrame(all_parameter_list, 
                           columns = ['specimen_index'] + columns)


ax = parallel_coordinates(specimen_df,"specimen_index", cols= columns)
ax.get_legend().remove()
plt.xticks(rotation=90)
plt.show()
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         