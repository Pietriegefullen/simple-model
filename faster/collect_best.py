
import os
import shutil
import stat
from itertools import cycle
from USER_VARIABLES import RESULTS_DIRECTORY, PROJECT_DIRECTORY
import parameters
from fit_sample import run
import data
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plot
import model
import optimizer

# TODO: 
# validation replica is specified.
# => fit replicas are determined automatically.
# => in split: other(s), in single, same.
# in split mode, 
# load checkpoint from split run
#   if sample has only two replicas, load from single
# 
# TODO: in plots, show which replicas are fit/validation.
# in single mode, both are same, i.e. ONLY fit, no validation!

parser = argparse.ArgumentParser()

parser.add_argument('input', default = None, nargs = '+')
args = parser.parse_args()

samples = None
validation_replica = None
fit_modes = None
if not args.input is None:
    for i in args.input:
        try:
            int(i)
            if len(str(i)) == 1:
                validation_replica = str(i)
            elif len(str(i)) == 4:
                samples = str(i)
            else:
                raise NotImplementedError()
        except:
            if i == 'split' or i == 'single':
                fit_modes = i
            else:
                raise NotImplementedError()

dataset = data.get_data_before_carex()
if samples is None:
    samples = [s.sample_name for s in dataset.samples]
elif not isinstance(samples, list):
    samples = [samples]
    
if fit_modes is None:
    fit_modes = ['split', 'single']
elif not isinstance(fit_modes, list):
    fit_modes = [fit_modes]

all_replicas = validation_replica is None

def get_best_checkpoint(criteria = None, exclude = None):
    candidate = None
    cp_id = None
    for root, dirs, files in os.walk(RESULTS_DIRECTORY):
        for file in files:
            if file.startswith('.'): continue
            file_path = os.path.join(root, file)
            
            pos = criteria is None or all([c in file_path for c in criteria])
            neg = exclude is None or not any([c in file_path for c in exclude])
            if not pos or not  neg:
                continue
            
            loaded_parameters, loss, run_config = parameters.load_parameter_file(file_path)
            cp_id = [f for f in file.split('_') if 'cp-' in f][0]
            run_id = [f for f in file.split('_') if 'run-' in f][0]
            if candidate is None or (loss < candidate[1]):
                candidate = (loaded_parameters, loss, run_config, cp_id)
    if candidate is None:
        print('found no candidates')
        print('criteria:', criteria, 'exclude:', exclude)
        return None, None
    return candidate[:-1], candidate[-1]

for fit_mode in fit_modes:
    for sample_number in samples:
        if all_replicas:
            validation_replica = [r.replica_number for r in dataset[sample_number].replicas]
        elif not isinstance(validation_replica, list):
            validation_replica = [validation_replica]
        for replica in validation_replica:
            replica_name = str(sample_number) + str(replica)
            fit_criteria = 'single' if len(dataset[sample_number].replicas) <= 2 else fit_mode
            cp, cp_id = get_best_checkpoint(criteria = ['_'.join([str(sample_number), 
                                                            str(replica)]), 
                                                  fit_criteria])
            if cp is None: continue
            cp_parameters, _, run_config = cp
            pathway_model, objective, run_log = run(run_config, cp_parameters)

            val_objective = optimizer.build_objective_function(pathway_model,
                                                               dataset[replica_name],
                                                               run_config['objective'],
                                                               t_start = 0, t_end = None)
            fit_replicas = dataset[sample_number].get_split(replica, fit_mode)['fit']
            axs = None
            markersize = 3
            markers = cycle(['^', 's', 'D', 'v', 'P', 'X'])
            handles = []
            for fit_rep in fit_replicas:
                marker = next(markers)
                fig, axs = plot.plot_data(fit_rep, ax = axs,
                                          marker = marker,
                                          ms = markersize,
                                          mfc = 'none')
                handles.append(Line2D(
                    [], [], color = 'k', marker = marker, linestyle = 'None',
                    markerfacecolor = 'none', markersize = markersize,
                    label = f'fit replica {fit_rep.replica_number}'))
            
            r = dataset[replica_name]
            if fit_mode == 'split': 
                fig, axs = plot.plot_data(r, ax = axs, marker = 'o', ms = markersize)
                handles.append(Line2D([], [], color='k', marker='o', linestyle='None',
                                   markersize=markersize,
                                   label=f'validation replica {r.replica_number}'))
            fig, axs = plot.plot_fit(run_log, ax = axs)

            r2_fit = objective.R2(run_log)
            r2_val = val_objective.R2(run_log)

            for loss_f, loss_w in objective.loss_contributions():
                print(loss_f._replica, loss_f.pool, loss_f.R2(run_log))

            print(r2_fit, r2_val)

            handles.append( Line2D([], [], color='k', linestyle='-', label='model'))

            axs['CH4'].legend(handles=handles,
                               fancybox = False,
                               edgecolor = 'k', 
                               loc = 'lower right')
                   
            for pool, ax in axs.items():
                plot.format_ax(ax,
                              log_scale = 'log' in run_config['objective']['transform'][pool])

            # store with cp_id
            plot_target = os.path.join(PROJECT_DIRECTORY, 'best_' + fit_mode, 
                                       sample_number, replica)

            plt.show()

            # in provided axes.
            # do a def plot_fit() wrapper around this.
            # store in plot target using 'fit', sample, replica, fit_mode, loss, cp_id
            # plot CO2 and CH4 side-by-side???
            # using the log-axis as transform dictates in run_config.
            
# store plots using cp_id

quit()


after = '2026-08-04--08-30'

fit_mode = 'split' # 'single' OR 'split

initial_mean_days = 0

plot = True #['1351']
plot_only_missing = True

log_co2 = False
log_ch4 = True

source = USER_VARIABLES.LOG_DIRECTORY
source_suffix = ''#'0-400' # '0-200'

if not source_suffix == '' and not source_suffix[0] =='_':
    source_suffix = '_' + source_suffix

if not fit_mode == 'split':
    source_suffix = '_'.join([source_suffix, fit_mode])
    
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
    print('replica target', replica_target)
    
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
    from main_file import load_fitted_parameters, plot_fit2 as plot_fit
    import model
    import matplotlib.pyplot as plt
    
    day_limits = None
    if not source_suffix == '' and '-' in source_suffix:
        day_limits = [int(i) for i in source_suffix.replace('_','').split('-')]
    
for sample_name in os.listdir(target):
    print(sample_name)        
    sample_dir = os.path.join(target, sample_name)
    if not os.path.isdir(sample_dir):
        continue
    for fit_replicas in os.listdir(sample_dir):
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
            if fit_mode == 'single':
                val_replica = fit_replicas
            
            elif fit_mode == 'split':
                for r in str(fit_replicas):
                    val_replica = val_replica.replace(str(r), '')
            else:
                raise NotImplementedError()
            
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
            
            log = pathway_model.predict(sample[val_replica],
                                        initial_mean_days = initial_mean_days)
            
            if not os.path.isdir(plot_target):
                os.makedirs(plot_target)
                
            for m, log_plot in zip(['CO2', 'CH4'], [log_co2, log_ch4]):
                plot_fit(sample[val_replica], log, m, log_plot,
                         single_fit = fit_mode == 'single')
                
                if not day_limits is None:
                    plt.gca().set_xlim(day_limits)
                
                file_name = '_'.join(['00', sample_name, val_replica, m, 'fit'])
                plt.savefig(os.path.join(plot_target,file_name + '.png') , dpi = 300)
            
            log.plot(save_target = plot_target, 
                     xlim = day_limits)
            
    print()
