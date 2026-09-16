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
            'model_type':               'complex',
            
            'normalized_parameters':    True,
            'algorithm':                'differential_evolution',
            'loss_weight':              {'CO2': 1.,
                                         'CH4': 1.}
            'reduction':                {'CO2': 'mse',
                                         'CH4': 'mse'},
            'transform':                {'CO2': ['normalize'],
                                         'CH4': ['log', 'normalize']}
            't_start':                  None,
            't_end':                    None
            }

algo_config = {'differential_evolution': {}, # empty dict uses default
               'powell':                 {}
               }

initial_parameters = {}


hasargs = len(sys.argv) > 1
if hasargs:
    chosen['sample'] = sys.argv[1]

# TODO: initialise parameters / set bounds
#       load parameters from chosen or default source
#       store initial parameters (and parameter range) in config

#TODO: set model variable/constant parameters , initial values LATER

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

# build objective function
replica_objectives = []
for replica in fit_replicas:
    replica_objective = Objective(pathway_model, replica)
    for pool in ['CO2', 'CH4']:
        transform = parameters.get_transform(chosen['transform'][pool])
        pool_loss = optimizer.loss_function(pool, chosen['reduction'][pool]], transform,
                                            t_start = chosen['t_start'],
                                            t_end = chosen['t_end'])
        replica_objective.add_loss(pool_loss, chosen['loss_weight'][pool])
total_objective = sum(replica_objectives)

# print callback?
# store checkpoint callback!
#   => TODO: design a sensible structure!
# if using hashes for model version, keep a lookup table to describe models!
total_objective.add_callback(CheckpointCallback())

algo.minimize(objective, initial_parameters)
