import matplotlib.pyplot as plt
import os
import data
import model
import USER_VARIABLES
from optimizer import r2
import numpy as np
import parameters

def load_parameter_range(sample_name, replica_name, model_type, best_N, after = None):
    fit_replicas = str(456).replace(str(replica_name),'')
    
    sample_name = str(sample_name)
    folders = []
    for _d in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
        if ((sample_name + str(fit_replicas[0])) in _d or \
        (sample_name + str(fit_replicas[1])) in _d) and \
         not (sample_name + str(replica_name)) in _d and model_type in _d:
            folders.append(_d)
    if len(folders) == 0:
        raise Exception('Found no fit results directory')
        
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
            
        largest_range = None
        best = sorted(best, key = lambda tpl: tpl[0])
        best = best[:best_N]
        for loss, loaded_parameters in best:
            if largest_range is None:
                largest_range = parameters.ModelParameters(loaded_parameters)
                
            else:
                for p, par in loaded_parameters.items():
                    largest_p = largest_range[p]
                    if largest_p.high is None or par > largest_p.high:
                        largest_p.high = par
                    if largest_p.low is None or par < largest_p.low:
                        largest_p.low = par
                    
    return largest_range

def load_fitted_parameters(sample_name, replica_name, model_type, after = None):
    fit_replicas = str(456).replace(str(replica_name),'')
    
    sample_name = str(sample_name)
    
    folders = []
    for _d in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
        if ((sample_name + str(fit_replicas[0])) in _d or \
        (sample_name + str(fit_replicas[1])) in _d) and \
         not (sample_name + str(replica_name)) in _d and model_type in _d:
            folders.append(_d)
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
        return best_parameters

if __name__ == '__main__':
    #model_type = 'simple' # or 'complex'
    model_type= 'complex'
    sample_name = '1351' # 1351, 1367, 1369, 1370, 
    replica_name = 4 # 4?, 5?, 6?
    reset_Fe3 = None #2000 # set the day on which to reset Fe3 to initial value, None to omit reset
    log_co2 = False
    log_ch4 = True
        

    dataset = data.get_data_before_carex()
    replica = dataset[sample_name + str(replica_name)]

    loaded_parameters = load_fitted_parameters(sample_name, replica_name, model_type,
                                               after = '2026-07-02--09-47')
    
    
    selected_pathways = model.get_pathways(model_type)
    pathway_model = model.Model(selected_pathways)
    pathway_model.parameters().set(loaded_parameters)
    print(pathway_model)
    
    print()

    log = pathway_model.predict(replica, reset_Fe3 = reset_Fe3)

    replicas = '456'.replace(str(replica_name), '')

    #log.plot()

    replica.plot(measurements = ['CO2'])
    log.plot(['CO2'], newfigure = False)
    for repl in replicas:
        if sample_name in dataset and str(repl) in dataset[sample_name]:
            
            t_pred, pred_co2 = log['CO2']
            t_meas, meas_co2 = dataset[sample_name + str(repl)].CO2()
            t = np.intersect1d(np.round(t_pred,2), np.round(t_meas,2))
            idx_pred = np.squeeze([np.nonzero(np.round(t_pred,2) == np.round(_t,2))[0] for _t in t])
            idx_meas = np.squeeze([np.nonzero(np.round(t_meas,2) == np.round(_t,2))[0] for _t in t])
            r2_co2 = r2(pred_co2[idx_pred], meas_co2[idx_meas], log = log_co2)
            
            dataset[sample_name + str(repl)].plot(measurements = ['CO2'], 
                                                  label = f'fit R² = {r2_co2:.2f}', 
                                                  marker = '.',
                                                  newfigure = False)
    sample = dataset[sample_name]
    plt.gca().set_title(f'{str(sample_name)} validation: {replica_name}')

    replica.plot(measurements = ['CH4'])
    log.plot(['CH4'], newfigure = False, log = True)
    for repl in replicas:
        if sample_name in dataset and str(repl) in dataset[sample_name]:
            
            t_pred, pred_ch4 = log['CH4']
            t_meas, meas_ch4 = dataset[sample_name + str(repl)].CH4()
            t = np.intersect1d(np.round(t_pred,2), np.round(t_meas,2))
            idx_pred = np.squeeze([np.nonzero(np.round(t_pred,2) == _t)[0] for _t in t])
            idx_meas = np.squeeze([np.nonzero(np.round(t_meas,2) == _t)[0] for _t in t])
            r2_ch4 = r2(pred_ch4[idx_pred], meas_ch4[idx_meas], log = log_ch4)
            
            dataset[sample_name + str(repl)].plot(measurements = ['CH4'], 
                                              label = f'fit R² = {r2_ch4:.2f}', 
                                              marker = '.',
                                              newfigure = False)
    plt.gca().set_title(f'{str(sample_name)} validation: {replica_name}')

    plt.show()

