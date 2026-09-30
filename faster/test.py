import os
from fit_sample import fit, run
from USER_VARIABLES import PROJECT_DIRECTORY
import parameters

TMP = os.path.join(PROJECT_DIRECTORY, 'tmp')

def test_round_trip():
    chosen = {
            'sample':                   1351,
            'validation_replica':       4, 
            
            't_start':                  0,
            't_end':                    None,
            
            'fit_mode':                 'split',
            'pathways':                 ['Hydrolysis',
                                         'Fermentation',
                                         'Hydro',
                                         'Aceto',
                                         'Homo',
                                         'Fe3'],
            'parameter_override':       {},  # e.g. 'Homo_thermodynamics': False
            'normalized_parameters':    True,
            'algorithm':                'differential_evolution',
            }

    objective_config = {
            'loss_weight':              {'CO2': 1.,
                                         'CH4': 1.},
            'reduction':                {'CO2': 'mse',
                                         'CH4': 'mse'},
            'transform':                {'CO2': ['normalize'],
                                         'CH4': ['log', 'normalize']},
            }
    
    algo_config = {
            'differential_evolution':   {}, # empty dict uses default
            'powell':                   {}
                }
    
    init_config = {
            'best_N':                   None,
            'sample':                   chosen['sample'],
            'validation_replica':       chosen['validation_replica'],
            'model':                    None,
            'run_ID':                   None,
            'file':                     None,
    }
    
    #fit(chosen, objective_config, algo_config, init_config, cp_target = TMP, minimize = False)
    
    init_config = {
            'best_N':                   1,
            'sample':                   chosen['sample'],
            'validation_replica':       chosen['validation_replica'],
            'model':                    None,
            'run_ID':                   None,
    }
    
    initial_parameters = parameters.load_parameters(init_config, source_directory = None)

    run(run_config = {'chosen': chosen, 
                      'objective': objective_config, 
                      'algo': algo_config,
                      'range': initial_parameters
                          },
        initial_parameters = initial_parameters, cp_target = TMP)
    
if __name__ == '__main__':
    test_round_trip()