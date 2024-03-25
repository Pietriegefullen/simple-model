import os
import matplotlib.pyplot as plt
from datetime import datetime

import model
import pathways

# TODO: parallel coordinates plot for optimization runs
# TODO: throw unused pools out of the system!
#       -> build system as model is configured? add only as needed
#       -> substances and microbes
# TODO: shorter, cleaner print
# TODO: encapsulate optimization.

# TODO: rebuild compare_model functionality:
#       run simle or complex
#       save
#       load/plot
# TODO: check matlab data and compare. why different?
# TODO: global color scheme?
# TODO: compute measures of fit (Whose responsibility?)


def get_pathways(model_type):
    basic = [pathways.Hydrolysis,
             pathways.Fermentation,
             pathways.Hydro,
             pathways.Aceto]
    if model_type == 'complex':
        return basic + [pathways.Homo,
                        pathways.Fe3]
    elif model_type == 'simple':
        return basic
    else:
        raise NotImplementedError()

if __name__ == '__main__':
    import data
    import USER_VARIABLES
    
    data = data.get_data_before_carex()
    target_directory = USER_VARIABLES.LOG_DIRECTORY
    
    for sample in data.samples:
        splits = sample.leave_one_out_split()
        for split in splits:
            fit_replicas = split['fit']
            validation_replica = split['val']
            
            for model_type in ['simple', 'complex']:
                chosen_pathways = get_pathways(model_type)
                pathway_model = model.Model(chosen_pathways)
                pathway_model.parameters().set('default')
            
                best_loss, _ = pathway_model.fit(fit_replicas)
                
                fit_repl = '_'.join([str(r) for r in fit_replicas])
                now = datetime.now()
                timestamp = '_' + now.strftime('%Y-%m-%d_%H-%M-%S')
                file_name = f'fit_{fit_repl}_{model_type}_{best_loss:.3g}' + timestamp
                pathway_model.save(target_directory, file_name)
                                
                results = pathway_model.predict(validation_replica)
                
                plt.close('all')
                results.plot(['CH4', 'CO2'], newfigure = False)
                validation_replica.plot()
                figure_name = file_name + '.svg'
                plt.savefig(os.path.join(target_directory, figure_name))
    
    
    
    