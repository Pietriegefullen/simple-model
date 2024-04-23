import os
import sys
import matplotlib.pyplot as plt
import data
import USER_VARIABLES
import model

if __name__ == '__main__':
    d = data.get_data_before_day()
   
    model_type = None
    if 'complex' in sys.argv:
        model_type = 'complex'
    elif 'simple' in sys.argv:
        model_type = 'simple'

    plot_log = False
    if 'log' in sys.argv:
        plot_log = True
    
    result_sample = sys.argv[1]
    folders = []
    for _d in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
        if result_sample in _d:
            folders.append(_d)
    #results_folder = 'fit_13544_2024-04-18--09-14-54'
    
    for results_folder in folders:
        parameter_source = os.path.join(USER_VARIABLES.LOG_DIRECTORY, results_folder)
        best_loss, p = model.get_best_loss_parameters(parameter_source)
        print('best loss', best_loss)
       
        loaded_model_type = 'simple'
        if 'complex' in results_folder:
            loaded_model_type = 'complex'
        if model_type is None:
            model_type = loaded_model_type
        elif not loaded_model_type == model_type:
            continue
        loaded_model = model.Model(model.get_pathways(model_type))
        loaded_model.parameters().set(p)
        print(loaded_model)

        fit_replicas = [s for s in results_folder.replace('simple', '').replace('complex','').replace('fit_', '').replace('log', '').split('_2024')[0].replace(' ', '_').split('_') if not s == '']

        for repl in fit_replicas:
            replica = d[repl]

            model_run = loaded_model.predict(replica)

            plt.figure()
            model_run.plot(['CO2', 'CH4'], newfigure = False)
            replica.plot(log = plot_log, newfigure = False)
            ax = plt.gca()
            t = ax.get_title()
            plt.title(t + results_folder)
            
            #plt.figure()
            #model_run.plot(['CO2', 'CH4'], newfigure = False)
            #replica.plot(newfigure = False)
            #ax = plt.gca()
            #t = ax.get_title()
            #plt.title(t + results_folder)
            

        #model_run.plot(['Fermentation_MM'])
        #model_run.plot(['Hydrolysis_MM'])
        #model_run.plot(['Hydro_MM'])
        #model_run.plot(['Aceto_MM'])

    plt.show()
