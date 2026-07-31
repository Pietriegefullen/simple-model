import matplotlib.pyplot as plt
import os
import data
import model
import USER_VARIABLES
from optimizer import r2
from data import all_sample_numbers
import numpy as np
import parameters

def load_parameter_range(sample_name, replica_name, model_type, best_N, after = None, best = 'best'):
    import data
    d = data.get_data_before_carex()
    sample = d[sample_name]
    repls = ''.join([str(r.replica_number) for r in sample.replicas])
    fit_replicas = repls.replace(str(replica_name),'')
    
    
    suffix = best.replace('_', '').replace('best', '')
    sample_name = str(sample_name)
    folders = []
    for _d in os.listdir(USER_VARIABLES.LOG_DIRECTORY + suffix):
        if ((sample_name + str(fit_replicas[0])) in _d or \
        (sample_name + str(fit_replicas[1])) in _d) and \
         not (sample_name + str(replica_name)) in _d and model_type in _d:
            folders.append(_d)
    if len(folders) == 0:
        raise Exception('Found no fit results directory in ' + str(USER_VARIABLES.LOG_DIRECTORY + suffix))
        
    elif len(folders) == 1:
        raise Exception('Found only a single result file')
        
    elif len(folders) > 1:
        print('loading')
        best = []
        for f in folders:
            parameter_source = os.path.join(USER_VARIABLES.LOG_DIRECTORY, f)
            try:
                loaded_loss, loaded_parameters = model.get_best_loss_parameters(parameter_source)#
            except:
                print('empty folder')
                continue
            
            date = f.split('_')[-1]
            if not after is None and date <= after:
                continue
            
            best.append((loaded_loss, loaded_parameters))
        
            
        best = sorted(best, key = lambda tpl: tpl[0])
        best = best[:best_N]
        largest_range = parameters.ModelParameters(best[0][1])
        print('Parameter range:', 'best loss', best[0][0], 'worst loss', best[-1][0])
        for loss, loaded_parameters in best:
            for p, par in loaded_parameters.items():
                largest_p = largest_range[p]
                if largest_p.high is None or par > largest_p.high:
                    largest_p.high = par
                if largest_p.low is None or par < largest_p.low:
                    largest_p.low = par
                
    return largest_range

def load_fitted_parameters(sample_name, replica_name, model_type, after = None, best = 'best', return_loss = False):
    import data
    d = data.get_data_before_day()
    sample = d[sample_name]
    repls = ''.join([str(r.replica_number) for r in sample.replicas])
    fit_replicas = repls.replace(str(replica_name),'')
    
    sample_name = str(sample_name)
    
    folders = []
    if best is None:
        for _d in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
            if ((sample_name + str(fit_replicas[0])) in _d or \
            (sample_name + str(fit_replicas[1])) in _d) and \
             not (sample_name + str(replica_name)) in _d and model_type in _d:
                folders.append(_d)
        
    else:
        source = os.path.join(USER_VARIABLES.simple_model_dir, best)
        sample_source = os.path.join(source, str(sample_name))
        replica_source = os.path.join(sample_source, fit_replicas)
        if not os.path.isdir(replica_source):
            raise Exception('No best result for this replica: ' + replica_source)
        folders.append(replica_source)
            
            
    if len(folders) == 0:
        raise Exception('Found no fit results directory')
        
    elif len(folders) > 0:
        best_loss = None
        best_parameters = None
        par_source = None
        for f in folders:
            parameter_source = os.path.join(USER_VARIABLES.LOG_DIRECTORY, f)
            try:
                loaded_loss, loaded_parameters = model.get_best_loss_parameters(parameter_source)#
            except:
                print('empty folder')
                continue
            
            date = f.split('_')[-1]
            if not after is None and date <= after:
                continue
            if best_parameters is None or loaded_loss < best_loss: 
                best_loss = loaded_loss
                best_parameters = loaded_parameters
                par_source = parameter_source
        
        print('loading parameters from', par_source)
        print('loss', best_loss)
        if not return_loss:
            return best_parameters
        return best_parameters, best_loss
    
def plot_fit(val_replica, run_log, measurement, log_fit):
    val_replica.plot(measurements = [measurement])
    
    sample = val_replica.sample
    sample_name = sample.sample_name
    
    fit_replicas = ''.join([str(r.replica_number) 
                            for r in sample.replicas]).replace(str(val_replica.replica_number), '')

    r2_val = run_log.R2(measurement, log_fit = log_fit)
    run_log.plot([measurement], newfigure = False, label = f'val R² = {r2_val:4.2f}')
    t_pred, pred_val = run_log[measurement]
    for repl in fit_replicas:
        if str(repl) in sample:
            fit_replica = sample[str(repl)]
            if measurement == 'CO2':
                t_meas, meas_val = fit_replica.CO2()

            else:
                t_meas, meas_val = sample[str(repl)].CH4()
            t = np.intersect1d(np.round(t_pred,2), np.round(t_meas,2))
            idx_pred = np.squeeze([np.nonzero(np.round(t_pred,2) == np.round(_t,2))[0] for _t in t])
            idx_meas = np.squeeze([np.nonzero(np.round(t_meas,2) == np.round(_t,2))[0] for _t in t])
            r2_val = r2(pred_val[idx_pred], meas_val[idx_meas], log = log_fit)
            
            fit_replica.plot(measurements = [measurement], 
                                                  label = f'fit R² = {r2_val:.2f}', 
                                                  marker = '.',
                                                  newfigure = False)
    plt.gca().set_title(f'{str(sample_name)} validation: {val_replica.replica_number}')


def plot_fitted_ratio(pathway_model, val_replica):
    # plot ratio of change in CO2 to change in CH4
    derivative_log = pathway_model.system_change_log
    
    ax = val_replica.plot_ratio()

    t, dCO2_dt = derivative_log['CO2']
    t, dCH4_dt = derivative_log['CH4']
    ax.plot(t,dCO2_dt/dCH4_dt)

if __name__ == '__main__':
    #model_type = 'simple' # or 'complex'
  
    model_type= 'complex'
    sample_name =  '1355' #1353 1351, 1367, 1369, 1370, 
    replica_name = 5 # 4?, 5?, 6?
    reset_Fe3 = None #2000 # set the day on which to reset Fe3 to initial value, None to omit reset
    log_co2 = False
    log_ch4 = True
    
    best = 'best' # 'best' or 'best_0-400 or ...
    dataset = data.get_data_before_day()
        
    #for sample_name in all_sample_numbers:
    for _ in range(1):
        
        val_replica = dataset[sample_name + str(replica_name)]

        loaded_parameters = load_fitted_parameters(sample_name, replica_name, model_type,
                                               after = '2026-07-02--09-47',
                                               best = best)
    
    
        selected_pathways = model.get_pathways(model_type)
        pathway_model = model.Model(selected_pathways)
        pathway_model.parameters().set(loaded_parameters)
        print(pathway_model)
    
        print()

        log = pathway_model.predict(val_replica, reset_Fe3 = reset_Fe3)
    
        #plot_fitted_ratio(pathway_model, val_replica)

        log.plot()

        #plot_fit(val_replica, log, 'CO2', log_co2)
        #plot_fit(val_replica, log, 'CH4', log_ch4)
        
        # NOTE: xlim not adjusted to fitting range (e.g. before day 400)

        plt.show()

