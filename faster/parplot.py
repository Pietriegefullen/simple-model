# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 20:35:56 2026

@author: Lara
"""

from parallel_coordinates import pcp as parallel_coordinates

from main_file import load_fitted_parameters
import parameters

if __name__ == '__main__':
    
    replicas = [
                13514,
                13515
        ]
    
    all_parameters = {}
    
    for rep in replicas:
        sample_name = str(rep)[:4]
        val_replica = str(rep)[-1]
    
        try:
            loaded_parameters = load_fitted_parameters(sample_name, 
                                                       val_replica,
                                                       model_type = 'complex',
                                                       best = True)
        except Exception as ex:
            if 'No best result' in str(ex):
                continue
            raise ex
    
        all_parameters[str(rep)] = loaded_parameters
    
    parameter_ordering = sorted(loaded_parameters.keys())
    results = [([str(r)] + [all_parameters[str(r)][p] for p in parameter_ordering]) for r in replicas]
    ytype = ['log' if par.log else 'linear' for par in parameters.default_model_parameters()]
    parallel_coordinates(results, ['replica']+ parameter_ordering, ytype = None)
    