import os
import numpy as np
import matplotlib.pyplot as plt

import argparse
from datetime import datetime
import json

import yellowbrick
from yellowbrick.features import ParallelCoordinates
from yellowbrick.datasets import load_occupancy

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
    del columns_list[column_index]
    for sp in range(len(all_parameter_list)):
        del all_parameter_list[sp][column_index]
    



    
  
#print(len(all_parameter_list[0])) 
#print(len(columns_list)) 


            
specimen_df = pd.DataFrame(all_parameter_list, 
                            columns = columns_list)

minmaxdict = {}
for col in columns_list[1:]:
    col_min = min(specimen_df[col])
    col_max = max(specimen_df[col])
    minmaxdict.update({col: [col_min, col_max]})
    specimen_df[col] = (specimen_df[col] - col_min) /(col_max - col_min)
# Calculate the limits on the data




columns_list.remove('specimen_index')
ax = parallel_coordinates(specimen_df,"specimen_index", cols= columns_list)
ax.get_legend().remove()
plt.xticks(rotation=90)
plt.show()
         
 
# plots für die verschiedenen Mikrobengruppen        
Fe3_columns_list = [s for s in columns_list if "Fe3" in s]
ax = parallel_coordinates(specimen_df,"specimen_index", cols= Fe3_columns_list)
ax.get_legend().remove()
plt.xticks(rotation=90)
plt.show() 
        

Ferm_columns_list = [s for s in columns_list if "Ferm" in s]
ax = parallel_coordinates(specimen_df,"specimen_index", cols= Ferm_columns_list)
ax.get_legend().remove()
plt.xticks(rotation=90)
plt.show() 
        

Homo_columns_list = [s for s in columns_list if "Homo" in s]
ax = parallel_coordinates(specimen_df,"specimen_index", cols= Homo_columns_list)
ax.get_legend().remove()
plt.xticks(rotation=90)
plt.show() 
        


Hydro_columns_list = [s for s in columns_list if "Hydro" in s]
ax = parallel_coordinates(specimen_df,"specimen_index", cols= Hydro_columns_list)
ax.get_legend().remove()
plt.xticks(rotation=90)
plt.show() 
        
      
Ac_columns_list = [s for s in columns_list if "Ac" in s]
Ac_columns_list.remove('Acetate')
ax = parallel_coordinates(specimen_df,"specimen_index", cols= Ac_columns_list)
ax.get_legend().remove()
plt.xticks(rotation=90)
plt.show() 
        
      
        
      
        
      
        
      
        
      
        
      
        
      
        
      
        
      
# anderer ansatz der dumm ist
         

# # Load the classification data set
# X= specimen_df
# y1 = np.ones(len(specimen_df['specimen_index']))
# print(y1)

# #X,y = load_occupancy()

# print(X)
# #print(y)
# # Specify the features of interest and the classes of the target
# features = columns_list
# #classes = ["unoccupied", "occupied"]

# # Instantiate the visualizer
# visualizer = ParallelCoordinates(
#       features=features,
#     normalize='standard', sample=0.05, shuffle=True,
# )

# # Fit the visualizer and display it
# visualizer.fit_transform(specimen_df[features], y1)
# #ax = plt.gca()
# plt.xticks(rotation=90)
# visualizer.show() 

         
         
         
         
         
         
         
         
         
         
         
         
         
         
         