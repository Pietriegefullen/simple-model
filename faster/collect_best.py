# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 14:41:30 2026

@author: Lara
"""
import os
import shutil
import USER_VARIABLES

after = '2026-07-08--08-30'

source = USER_VARIABLES.LOG_DIRECTORY
target = os.path.join(USER_VARIABLES.simple_model_dir, 'best')

for f in os.listdir(source):
    folder = os.path.join(source,f)
    if not os.path.isdir(folder): continue
    
    model_type = 'complex' if 'complex' in f else 'simple'
    fit = f.split()[0].replace('fit_','')
    fit_replicas = fit.split('_')[:2]
    sample_name = fit_replicas[0][:-1]
    fit_replicas = ''.join(sorted([fr.replace(sample_name, '') for fr in fit_replicas]))
    
    date = f.split('_')[-1]
    if date < after:
        continue
    
    replica_target = os.path.join(target, sample_name, fit_replicas)
    
    if not os.path.isdir(replica_target):
        os.makedirs(replica_target)
    
    current_best = None
    for pf in os.listdir(folder):
        
        best_files = os.listdir(replica_target)
        if len(best_files) == 1:
            current_best = float(best_files[0].split('loss_')[-1])
        elif len(best_files) > 1:
            raise Exception()
        
        file = os.path.join(folder, pf)
        if not os.path.isfile(file): continue
        
        loss = float(pf.split('loss_')[-1])
        
        best_file_name = f'{sample_name}_{fit_replicas}_' + pf
        
        if current_best is None or (len(best_files) > 0 and loss < current_best):
            shutil.copy2(file, os.path.join(replica_target, best_file_name))
            current_best = loss
        
            if len(best_files) > 0:
                os.remove(os.path.join(replica_target, best_files[0]))
        
for sample in os.listdir(target):
    print(sample)
    for replica in os.listdir(os.path.join(target, sample)):
        results = os.listdir(os.path.join(target, sample, replica))
        if len(results) == 1:
            loss = float(results[0].split('loss_')[-1])
            print('   '+ replica, f'{loss:.2f}')
    print()
