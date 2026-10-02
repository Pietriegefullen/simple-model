
import os
import shutil
import stat
from itertools import cycle
from USER_VARIABLES import RESULTS_DIRECTORY, PROJECT_DIRECTORY
import parameters
from fit_sample import run
import data
import argparse
import hashing
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plot
import model
import optimizer

# TODO: get_best_... MUST respect model id, loss id, ...
#       => otherwise, checkpoints will be mixed!!!
#       => get all checkpoints first, separated
#       => or return all matching, separated by model id, ...

# TODO: improve plotting interface!!!
# in provided axes.
# do a def plot_fit() wrapper around this.
# store in plot target using 'fit', sample, replica, fit_mode, loss, cp_id
# plot CO2 and CH4 side-by-side???
# using the log-axis as transform dictates in run_config.

# TODO: val objective should use same normalization as fit?
#       => does normalization make ANY sense for R2?

parser = argparse.ArgumentParser()

parser.add_argument('input', default = None, nargs = '?')
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
        #print('found no candidates')
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

            # store with cp_id
            if cp is None: continue
            cp_parameters, _, run_config = cp
            model_id = hashing.build_model_id(run_config['model'])
            print(model_id,  sample_number, replica, fit_mode)
            plot_target = os.path.join(PROJECT_DIRECTORY, 
                                       '_'.join(['best', model_id, fit_mode]),
                                       '-'.join([sample_number, replica]))
            file_name = '_'.join(['00',
                                   '-'.join([sample_number, replica]),
                                   model_id, fit_mode, cp_id])

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
                print(loss_f._replica, loss_f.pool, f'{loss_f.R2(run_log):.2f}')

            r2str = ', '.join([rf'$R^2_{{\mathrm{{{m}}}}} = {r:.2f} $' 
                              for r, m in zip([r2_fit, r2_val], ['fit', 'val']) if not (m == 'val' and fit_mode == 'single')])

            handles.append( Line2D([], [], color='k', linestyle='-',
                                  label=' '.join([r'$\mathrm{model}$', r2str])))

            axs['CH4'].legend(handles=handles,
                               fancybox = False,
                               edgecolor = 'k', 
                               loc = 'lower right')
                   
            for pool, ax in axs.items():
                plot.format_ax(ax,
                              log_scale = 'log' in run_config['objective']['transform'][pool])

            if not os.path.isdir(plot_target):
                os.makedirs(plot_target)
            plt.savefig(os.path.join(plot_target, file_name), dpi = 300)
            #plt.show()
