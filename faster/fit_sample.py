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

# SET DEFAULTS
default = {
            'sample':       '1351',
            'model_type':   'complex',
            
            }
# query and change as you go
hasargs = len(sys.argv) > 1

if hasargs:
    default['sample'] = sys.argv[1]


# PREPARE DATA
# load data
# get fit and validation replicas

# CONFIGURE MODEL
# choose pathways
# build model
# initialise parameters

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

# run optimiser (until stopping criterion?)
# save intermediate results (handled by Objective? -> or callback to Optimiser?)
# catch interrupt

d = data.get_data_before_day()
sample = d[default['sample']]

if fit_mode == 'split':
    splits = sample.leave_one_out_split()
elif fit_mode == 'single':
    splits = [{'fit':replica, 'val': replica}
              for replica in sample.replicas]
else:
    raise NotImplementedError()


if not local_search:
    algo = optimizer.DifferentialEvolution()
    
    
else:
    algo = optimizer.Powell()

try:
    algo.fit(pathway_model, fit_replicas, val_replicas)

except KeyboardInterrupt:
    quit()