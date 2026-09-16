import sys
import os
import traceback
import matplotlib.pyplot as plt
from datetime import datetime
import json
import numpy as np

import model
import data
import USER_VARIABLES

# TODO: consistent suffix handling!
# => only as support, don't rely on it!

# TODO: list included pathways, ditch 'simple'/'complex', but maintain backwards compatibility?

# optionally load initial parameters from different source
#   => checkpoint ID? (-> hash), 


# SET DEFAULTS
chosen = {
            'sample':                   1351,
            'validation_replica':       4, 
            'fit_mode':                 'single',
            'model_type':               'simple',
            'normalized_parameters':    True,
            'algorithm':                'differential_evolution',
            
            }

algo_config = {'differential_evolution': {},
               'powell':                 {}
               }

initial_parameters = {}


hasargs = len(sys.argv) > 1
if hasargs:
    chosen['sample'] = sys.argv[1]

# TODO: initialise parameters
#       load parameters from chosen or default source
#       store initial parameters (and parameter range) in config

# PREPARE OPTIMISER
# choose algorithm
# set checkpoint path
# configure algorithm
# get variables
# set bounds
# configure replica objectives and total objective 

# save hyperparameters
#   timestamp
#   model version
#   fit/val replicas
#   variables, initial parameter values, parameter ranges
#   optimiser and objective configuration

# for hashing, make sure to unify datatypes! e.g. sample number as int/str

# save intermediate results (handled by Objective? -> or callback to Optimiser?)
# catch interrupt


#TODO: set model variable/constant parameters , initial values LATER

# save hyperparameters with every plot (how?) -> maintain origin: model version, ...

# get sample from dataset
dataset = data.get_data_before_day()
sample = dataset[chosen['sample']]
split = sample.get_split(chosen['validation_replica'], 
                         chosen['fit_mode'])

# build and configure model
chosen_pathways = model.get_pathways(chosen['model_type'])
pathway_model = model.Model(chosen_pathways)
pathway_model.parameters().set('default', normalized = chosen['normalized_parameters'])
pathway_model.parameters().set(initial_parameters)


# select optimiser
algo = optimizer.get(chosen['algorithm'])
algo.configure(algo_config[chosen['algorithm']])

# build objective:
    # keep optimizer.Objective and optimizer.ReplicaObjective for now? but simplify!
objective = optimizer.Objective(pathway_model)

for replica in fit_replicas:
    objective.add_loss(CO2_loss_function(replica))
    objective.add_loss(CH4_loss_function(replica))
    
# for each replica, 
# define handling of CO2 and CH4 individually, possibly Ac, Fe...
# allow t_start to t_end fitting
# => using run log, compute loss
# consider normalization, weights, log transform, 

# what is input to loss function?
# run log, data, i.e. pred, true. anything else?
# -> transformations, normalizations, ...


# define checkpoint callback


algo.minimize(objective)
