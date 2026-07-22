# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 20:35:56 2026

@author: Lara
"""

from parallel_coordinates import pcp as parallel_coordinates
import matplotlib.pyplot as plt
import data
from main_file import load_fitted_parameters
import parameters
import model
import pathways
import numpy as np

if __name__ == '__main__':
    
    d = data.get_data_before_carex()
    
    log_CO2 = False
    log_CH4 = True
    
    model_type = 'complex'
    
    replicas = [str(r) for r in d.replicas()]
    
    replicas = [
                '13514', 
                '13515', 
                '13516', 
                '13674', 
                '13675', 
                '13676', 
                '13694', 
                '13695', 
                '13696', 
                '13704', 
                '13705', 
                '13706', 
                '13724', 
                '13725', 
                '13726', 
                '13734', 
                '13735', 
                '13736', 
                '13744', 
                '13745', 
                '13754', 
                '13755', 
                '13756', 
                '13764', 
                '13765', 
                '13766'
                ]
    
    all_parameters = {}
    all_replicas = []
    for rep in replicas:
        sample_name = str(rep)[:4]
        val_replica = str(rep)[-1]
    
        try:
            loaded_parameters, loaded_loss = load_fitted_parameters(sample_name, 
                                                       val_replica,
                                                       model_type = model_type,
                                                       best = True, 
                                                       return_loss = True)
        except Exception as ex:
            if 'No best result' in str(ex):
                continue
            raise ex
    
        all_replicas.append(rep)
        all_parameters[str(rep)] = loaded_parameters
        all_parameters[str(rep)]['loss'] = float(loaded_loss)
        
        selected_pathways = model.get_pathways(model_type)
        pathway_model = model.Model(selected_pathways)
        pathway_model.parameters().set(loaded_parameters)
        log = pathway_model.predict(d[rep])
        val_r2_CO2 = log.R2('CO2', log_fit = log_CO2)
        val_r2_CH4 = log.R2('CH4', log_fit = log_CO2)
        
        all_parameters[str(rep)]['R2_CO2'] = val_r2_CO2
        all_parameters[str(rep)]['R2_CH4'] = val_r2_CH4
    
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
        results = [([int(r), 
                     all_parameters[str(r)]['loss'],
                     all_parameters[str(r)]['R2_CO2'], all_parameters[str(r)]['R2_CH4']] + [all_parameters[str(r)][p] 
                                for p in parameter_group])
                   for r in all_replicas]
        ytype = ['linear', 'log', 'linear', 'linear'] + [dp[p].scale for p in parameter_group]
        ylims = [[13510, 13806], [], [0,1], [0,1]] + [[dp[par].low, dp[par].high] if dp[par].is_variable() else [] 
                        for par in parameter_group] 
        
        labels =  ['replica', 'loss', 'R2 CO2', 'R2 CH4'] + parameter_group
    
        fig = parallel_coordinates(results, labels, 
                             ytype = ytype,
                             ylim = ylims,
                             curves = False,
                             colorbar = False,
                             #alpha = np.maximum(0,results[2])
                             )
    plt.show()
    