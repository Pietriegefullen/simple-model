import os

import numpy as np
import json
import matplotlib.pyplot as plt

import parameters
import model
import data
from USER_VARIABLES import LOG_DIRECTORY

def save_all_plots():
    PARENT_DIR, _ = os.path.split(LOG_DIRECTORY)
    plot_target = os.path.join(PARENT_DIR, 'plots')
    print('plotting to', plot_target)
    if os.path.isdir(plot_target):
        raise Exception('Target directory for plots already exists.')

    os.makedirs(plot_target)

    d = data.get_data_before_day()

    plt.close('all')
    all_parameters, found = parameters.load_best()
    parameters.boxplots(all_parameters, save_target = plot_target)

    for sample in d.samples:
        for model_type in ['simple', 'complex']:
            try:
                splits = sample.leave_one_out_split()
            except Exception as ex:
                print('While getting splits:')
                print(ex)
                continue

            for split_number, split in enumerate(splits):
                fit_replica = splits[split_number]['fit']
                val_replica = splits[split_number]['val']

                name = 'fit_' + '_'.join([str(r) for r in sorted(fit_replica, key = lambda r: r.replica_number)])
               
                found = []
                for f in os.listdir(LOG_DIRECTORY):
                    if name in f and model_type in f and 'log' in f:
                        found.append(f)

                if not found:
                    s = name + ' ' + model_type
                    print(f'Found no results for {s}')
                    break

                elif len(found) > 1:
                    dates = []
                    for f in found:
                        timestamp = f.split('_')[-1]
                        dates.append(timestamp)
                    latest_index = dates.index(max(dates))
                    latest_found = [found[latest_index]]
                    losses = []
                    for f in found:
                        pf = os.path.join(LOG_DIRECTORY, f)
                        try:
                            best_loss, _ = model.get_best_loss_parameters(pf)
                        except Exception as ex:
                            print('Failed to load best parameters', name)
                            print(ex)
                            best_loss = 1e12

                        losses.append(best_loss)
                    best_index = losses.index(min(losses))
                    found = [found[best_index]]
                    if not latest_index == best_index: 
                        print('Found several results, using best. (not equal latest!)')
                    else:
                        print('Found several results, using best.')
                   
                if not len(found) == 1:
                    raise Exception('found')

                parameter_source = os.path.join(LOG_DIRECTORY, found[0])
                best_loss, best_parameters = model.get_best_loss_parameters(parameter_source)

                loaded_model = model.Model(model.get_pathways(model_type))
                loaded_model.parameters().set(best_parameters)
                #val_run = loaded_model.predict(val_replica)
                val_run = loaded_model.predict(val_replica, t = range(int(val_replica.last_day)))

                figure_postfix = '_'.join([
                                        sample.sample_name + val_replica.replica_number,
                                        model_type,
                                        'val'
                                        ])

                figure_target = os.path.join(plot_target, model_type, sample.sample_name + val_replica.replica_number)
                if not os.path.isdir(figure_target):
                    os.makedirs(figure_target)

                parameter_file = os.path.join(figure_target, f'parameters_loss_{best_loss:.4f}.json')
                with open(parameter_file, 'w') as pf:
                    json.dump(best_parameters, pf, indent = 4)
                

                plt.close('all')
                _, measured_CO2 = val_replica.CO2()
                _, measured_CH4 = val_replica.CH4()
                _, predicted_CO2 = val_run['CO2_on_measured']
                _, predicted_CH4 = val_run['CH4_on_measured']
                max_val = max([np.max(measured_CO2), np.max(measured_CH4), np.max(predicted_CO2), np.max(predicted_CH4)])
                plt.plot(measured_CO2, predicted_CO2, 'rx')
                plt.plot(measured_CH4, predicted_CH4, 'bx')
                plt.plot([0,max_val], [0, max_val], 'k-')
                figure_name = 'measured_vs_predicted'
                plt.savefig(os.path.join(figure_target, figure_name + figure_postfix)+'.svg')

                plt.close('all')
                val_run.plot(['CO2', 'CH4'], newfigure = False)
                for plot_log in [False, True]:
                    figure_name = 'CO2-CH4' + ('log' if plot_log else 'lin')
                    val_replica.plot(log = plot_log, newfigure = False)
                    plt.savefig(os.path.join(figure_target, figure_name + figure_postfix)+'.svg')

                for k in val_run._log.keys():
                    plt.close('all')
                    val_run.plot(k)
                    figure_name = k
                    plt.savefig(os.path.join(figure_target, figure_name + figure_postfix)+ '.svg')

if __name__ == '__main__':
    save_all_plots()
