# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 14:41:30 2026

@author: Lara
"""
import os
import shutil
import stat
import USER_VARIABLES

after = '2026-07-08--08-30'

plot = False #['1351']
plot_only_missing = False

log_co2 = False
log_ch4 = True

source = USER_VARIABLES.LOG_DIRECTORY
source_suffix = ''#'0-400' # '0-200'

if not source_suffix == '' and not source_suffix[0] =='_':
    source_suffix = '_' + source_suffix
    
source = source + source_suffix
target = os.path.join(USER_VARIABLES.simple_model_dir, 'best' + source_suffix)

for f in os.listdir(source):
    folder = os.path.join(source,f)
    if not os.path.isdir(folder): continue
    
    model_type = 'complex' if 'complex' in f else 'simple'
    if model_type == 'simple':
        raise NotImplementedError()
        
    fit = f.split()[0].replace('fit_','')
    fit_replicas = [k for k in fit.split('complex')[0].replace('fit_','').split('_')
                    if not k == '']
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
        
        best_files =  [f for f in os.listdir(replica_target)
                       if not os.path.isdir(os.path.join(replica_target, f))]
        if len(best_files) == 1:
            current_best = float(best_files[0].split('loss_')[-1])
            
        elif len(best_files) > 1:
            print(best_files)
            raise Exception()
        
        file = os.path.join(folder, pf)
        if not os.path.isfile(file): continue
        
        loss = float(pf.split('loss_')[-1])
        
        best_file_name = f'{sample_name}_{fit_replicas}_' + pf
        
        if current_best is None or (len(best_files) > 0 and loss < current_best):
            current_best = loss
            
            plot_target = os.path.join(replica_target, 'plot')
            if os.path.isdir(plot_target):
                os.chmod(plot_target, stat.S_IWRITE)
                shutil.rmtree(plot_target)
                
            shutil.copy2(file, os.path.join(replica_target, best_file_name))
    
            if len(best_files) > 0:
                os.remove(os.path.join(replica_target, best_files[0]))
    
if not plot is None:
    import data
    dataset = data.get_data_before_carex()
    from main_file import load_fitted_parameters, plot_fit
    import model
    import matplotlib.pyplot as plt
    
    day_limits = None
    if not source_suffix == '':
        day_limits = [int(i) for i in source_suffix.replace('_','').split('-')]
    
for sample_name in os.listdir(target):
    print(sample_name)        
    for fit_replicas in os.listdir(os.path.join(target, sample_name)):
        results = [f for f in os.listdir(os.path.join(target, sample_name, fit_replicas))
                   if os.path.isfile(os.path.join(target, sample_name, fit_replicas, f))]
        if len(results) == 0:
            print('   '+fit_replicas, 'no results')
        elif len(results) > 1:
            print('   '+ fit_replicas, 'more than 1 result!')
        elif len(results) == 1:
            loss = float(results[   0].split('loss_')[-1])
            print('   '+ fit_replicas, f'{loss:.2f}')

            if isinstance(plot, list) and not sample_name in plot:
                continue
        
            if not plot:
                continue

            sample = dataset[sample_name]
            plot_target = os.path.join(target, sample_name, fit_replicas, 'plot')
            
            print(sample_name, plot_target)
            if plot_only_missing and os.path.isdir(plot_target) and len(os.listdir(plot_target)) > 0:
                continue
            
            plt.close('all')

            val_replica = ''.join([str(r.replica_number) 
                                    for r in sample.replicas])
            for r in str(fit_replicas):
                val_replica = val_replica.replace(str(r), '')
            
            try:
                print('loading', val_replica)
                loaded_parameters = load_fitted_parameters(sample_name, 
                                                           val_replica,
                                                           model_type = 'complex',
                                                           best = 'best' + source_suffix)
            except Exception as ex:
                if 'No best result' in str(ex):
                    continue
                print(ex)
                raise ex

            selected_pathways = model.get_pathways(model_type)
            pathway_model = model.Model(selected_pathways)
            pathway_model.parameters().set(loaded_parameters)
            
            log = pathway_model.predict(sample[val_replica])
            
            if not os.path.isdir(plot_target):
                os.makedirs(plot_target)
                
            for m in ['CO2', 'CH4']:
                plot_fit(sample[val_replica], log, m, log_co2)
                
                if not day_limits is None:
                    plt.gca().set_xlim(day_limits)
                
                file_name = '_'.join(['00', sample_name, val_replica, m, 'fit'])
                plt.savefig(os.path.join(plot_target,file_name + '.png') , dpi = 300)
            
            log.plot(save_target = plot_target, 
                     xlim = day_limits)
            
    print()
