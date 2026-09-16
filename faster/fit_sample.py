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

# MODEL CONFIGURATION:
# TODO: enable/disable pathways
# TODO: enable/disable thermodynamics


# SET DEFAULTS
chosen = {
            'sample':                   1351,
            'validation_replica':       4, 
            'fit_mode':                 'single',
            'model_type':               'simple',
            'normalized_parameters':    True,
            'algorithm':                'differential_evolution',
            'loss_weight':              {'CO2': 1.,
                                         'CH4': 1.}
            'reduction':                {'CO2': 'mse',
                                         'CH4': 'mse'},
            'transform':                {'CO2': ['normalize'],
                                         'CH4': ['log', 'normalize']}
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


class Objective():
    
    def __init__(self, model):
        self._model = model
    
    def add_loss(self, replica, loss_function, weight = 1.0)

objective = optimizer.Objective(pathway_model)


# TODO: move elsewhere
def mse(true, pred):
    return np.sqrt(np.sum((true - pred)**2))

loss_functions = {'mse': mse}

def loss_function(pool, reduction = 'mse', transform = None):
    reduction_function = loss_functions['mse']
    def loss(run_log):
        t_pred, pool_pred = run_log[pool]
        t_true, pool_true = replica[pool]
        
        if callable(transform):
            pool_pred = transform(pool_pred)
            pool_true = transform(pool_true)
    
        loss_value = reduction_function(pool_true, pool_pred)
    
        return loss_value
    return loss

transforms = {'log': parameters.LogTransform,
              'normalize': parameters.MinMaxNormalization}

def get_transform(function_names):
    if not isinstance(function_names, list):
        function_names = [function_names]
    
    transform_function = None
    for f in function_names:
        transform_function = transforms[f](transform_function)
        
    return transform_function


for replica in fit_replicas:
    for pool in ['CO2', 'CH4']:
        transform = get_transform(chosen['transform'][pool])
        loss = loss_function(pool, chosen['reduction'][pool]], transform)
        objective.add_loss(replica, loss, chosen['loss_weight'][pool])
    
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
