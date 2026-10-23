
from fit_sample import run
from plot_checkpoint import resolve_checkpoint, replay_inputs
from USER_VARIABLES import RESULTS_DIRECTORY

if __name__ == '__main__':
    cp_id = 'QKE4'
    source = resolve_checkpoint(cp_id, RESULTS_DIRECTORY)
    _, run_config, initial_parameters, contents = replay_inputs(source)
    model, objective, model_run = run(run_config, initial_parameters)

    for k in model_run.keys():
        print(k)

    print(model_run['CO2'])


