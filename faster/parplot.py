# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 20:35:56 2026

@author: Lara
"""

from parallel_coordinates import pcp as parallel_coordinates

from main_file import load_fitted_parameters
import parameters
import model
import pathways

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
    
    groups = []
    for pathway_name in model.get_pathways('complex'):
        pathway = pathways.pathway_by_name(pathway_name)
        parameter_registry = parameters.ModelParameters()
        pathway(parameter_registry)
        pathway_parameters = parameter_registry.as_dict()
        groups.append(list(pathway_parameters.keys()))
    
    dp = parameters.ModelParameters({p.name:p for p in parameters.default_model_parameters()})
    for parameter_group in groups:
        results = [([str(r)] + [all_parameters[str(r)][p] 
                                for p in parameter_group])
                   for r in replicas]
        ytype = [par.scale for par in dp]
        ylims = [[]] + [[dp[par].low, dp[par].high] if dp[par].is_variable() else [] 
                        for par in parameter_group]
        parallel_coordinates(results, ['replica'] + parameter_group, 
                             ytype = None,
                             ylim = ylims,
                             curves = False)
    