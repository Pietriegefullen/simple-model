import matplotlib.pyplot as plt
import os
import data
import model
import USER_VARIABLES

def load_fitted_parameters(sample_name, replica_name, after = None):
    fit_replicas = str(456).replace(str(replica_name),'')
    
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
    sample_name = '1370' # 1351, 1367, 1369, 1370, 
    replica_name = 4 # 4?, 5?, 6?
    reset_Fe3 = None #2000 # set the day on which to reset Fe3 to initial value, None to omit reset
    
    
    loaded_parameters = load_fitted_parameters(sample_name, replica_name, 
                                               after = '2026-07-02--09-47')
    
    dataset = data.get_data_before_carex()
    replica = dataset[sample_name + str(replica_name)]

    selected_pathways = model.get_pathways(model_type)
    pathway_model = model.Model(selected_pathways)
    pathway_model.parameters().set(loaded_parameters)
    print(pathway_model)

    log = pathway_model.predict(replica, reset_Fe3 = reset_Fe3)

    replicas = '456'.replace(str(replica_name), '')

    #log.plot()

    replica.plot(measurements = ['CO2'])
    log.plot(['CO2'], newfigure = False)
    for repl in replicas:
        if sample_name in dataset and str(repl) in dataset[sample_name]:
            dataset[sample_name + str(repl)].plot(measurements = ['CO2'], 
                                                  label = 'fit', 
                                                  marker = '.',
                                                  newfigure = False)
    sample = dataset[sample_name]
    plt.gca().set_title(f'{str(sample)} {sample.site} ({sample.origin})')

    replica.plot(measurements = ['CH4'])
    log.plot(['CH4'], newfigure = False, log = True)
    for repl in replicas:
        if sample_name in dataset and str(repl) in dataset[sample_name]:
            dataset[sample_name + str(repl)].plot(measurements = ['CH4'], 
                                              label = 'fit', 
                                              marker = '.',
                                              newfigure = False)
    plt.gca().set_title(f'{str(sample)} {sample.site} ({sample.origin})')
    
    print('DOC on day', log['DOC'][0][-1], log['DOC'][1][-1])
    
    plt.show()

