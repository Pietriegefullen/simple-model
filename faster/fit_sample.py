import os

import model
import data
import optimizer
import parameters
import checkpoint
import hashing
import numpy as np


def _requested_checkpoint(path):
    """Load an annotated checkpoint selected explicitly on the command line."""
    try:
        contents = parameters.load_file(path)
    except (OSError, ValueError) as error:
        raise ValueError(f'Could not read checkpoint {path}: {error}') from error
    required = {'parameters', 'total_loss', 'run_config'}
    if not isinstance(contents, dict) or not required.issubset(contents):
        present = set(contents) if isinstance(contents, dict) else set()
        raise ValueError(
            f'{path} is not an annotated checkpoint; missing '
            + ', '.join(sorted(required - present)))
    return contents


def checkpoint_compatibility_warnings(checkpoint_data, chosen, model_config):
    """Describe metadata differences that need explicit operator approval."""
    saved_run = checkpoint_data.get('run_config')
    if not isinstance(saved_run, dict):
        return ['checkpoint has no run configuration to verify']

    warnings = []
    saved_chosen = saved_run.get('chosen', {})
    for key, label in (
            ('sample', 'sample'),
            ('validation_replica', 'validation replica'),
            ('fit_mode', 'fit mode')):
        expected = chosen.get(key)
        actual = saved_chosen.get(key)
        if str(actual) != str(expected):
            warnings.append(f'{label}: checkpoint={actual!r}, requested={expected!r}')

    saved_model = saved_run.get('model')
    try:
        expected_model_id = hashing.build_model_id(model_config)
        actual_model_id = hashing.build_model_id(saved_model)
    except (AttributeError, KeyError, TypeError, ValueError):
        warnings.append('checkpoint has no valid model variant metadata')
    else:
        if actual_model_id != expected_model_id:
            warnings.append(
                f'model variant: checkpoint={actual_model_id}, '
                f'requested={expected_model_id}')
    return warnings


def checkpoint_range_errors(checkpoint_data, fit_parameters):
    """Return values that cannot be used in this fit's declared search range."""
    saved_parameters = parameters.ModelParameters(checkpoint_data['parameters'])
    errors = []
    for saved in saved_parameters:
        # A variant mismatch can contain parameters unused by the requested
        # model.  Compatibility confirmation covers that case; values for
        # parameters the model does use must still be in its own range.
        if saved not in fit_parameters:
            continue
        target = fit_parameters[saved.name]
        if not target.is_variable():
            continue
        if target.options is not None:
            in_range = saved.value in target.options
            range_text = repr(target.options)
        else:
            in_range = target.low <= saved.value <= target.high
            range_text = f'[{target.low:g}, {target.high:g}]'
        if not in_range:
            errors.append(
                f'{saved.name}={saved.value!r} is outside {range_text}')
    return errors


def confirm_checkpoint_compatibility(path, warnings, input_function=None):
    """Require a deliberate acknowledgement for an incompatible restart."""
    if not warnings:
        return
    print(f'WARNING: checkpoint {path} is not compatible with this fit:')
    for warning in warnings:
        print(f'  - {warning}')
    if input_function is None:
        input_function = input
    try:
        confirmed = input_function('Use this checkpoint anyway? [y/N] ')
    except EOFError as error:
        raise RuntimeError('Checkpoint use was not confirmed.') from error
    if confirmed.strip().lower() not in {'y', 'yes'}:
        raise RuntimeError('Checkpoint use was not confirmed.')

# TODO: for samples with two replicas, split is equivalent to single
#       => store only single folder
#       => but plotting differs!

# TODO: when loading candidates for fit_sample.py 1351 6 --best 10, 
#       Hydro.use_thermodynamics is variable and should be bool?
#       => when computing run ID for checkpoint
#       CAUSE: model use_thermodynamics is seen as variable. should not!?
#       SUSPICION: are checkpoints loaded from incompatible model variants!!!!

# TODO: handle few usable sample points!!!
# TODO: rename checkpoints:
#       recompute model_id from config and rename files.
#       => always compute id from config, b/c hashing makes config canonical.

def run(run_config, initial_parameters, **kwargs):
    chosen = run_config['chosen']
    objective_config = run_config['objective']
    algo_config = {chosen['algorithm']: run_config['algo']}
    rng = {k:v for k, v in run_config['range'].items() if k in initial_parameters}
    init_config = {'range': parameters.ModelParameters(rng),
                   'parameters': initial_parameters}
    pathway_model, objective, run_log = fit(chosen, objective_config, algo_config, init_config, 
                                          store_checkpoints = False,
                                          minimize = False,
                                          **kwargs)
    return pathway_model, objective, run_log

def fit(chosen, objective_config, algo_config, init_config, 
        store_checkpoints = True,
        checkpoint_keep_only_n = 10,
        verbose_callback = False,
        minimize = True, 
        cp_target = None,
        checkpoint_source = None,
        overwrite_existing_checkpoints = False,
        checkpoint_confirmation = None):
    # get sample from dataset
    dataset = data.get_data_before_day()
    sample = dataset[chosen['sample']]
    split = sample.get_split(chosen['validation_replica'], 
                             chosen['fit_mode'])

    # build model and configure parameters
    pathway_model = model.build_model(chosen['pathways'], chosen['normalized_parameters'])

    requested_checkpoint_path = init_config.get('file')
    if requested_checkpoint_path is not None:
        requested_checkpoint = _requested_checkpoint(requested_checkpoint_path)
        range_errors = checkpoint_range_errors(
            requested_checkpoint, pathway_model.parameters())
        if range_errors:
            raise ValueError(
                'Requested checkpoint has values outside the fit parameter range:\n  '
                + '\n  '.join(range_errors))
        compatibility_warnings = checkpoint_compatibility_warnings(
            requested_checkpoint, chosen,
            pathway_model.get_config(only_structure = True))
        confirm_checkpoint_compatibility(
            requested_checkpoint_path, compatibility_warnings,
            input_function=checkpoint_confirmation)

    if 'best_N' in init_config and not init_config['best_N'] is None:
        model_id = hashing.build_model_id(pathway_model.get_config(only_structure = True))
        init_config['model'] = model_id

    legacy_path = None
    initial_parameters = parameters.load_parameters(init_config, source_directory=checkpoint_source)
        
    # override model parameters
    for p_name, p_value in chosen['parameter_override'].items():
        initial_parameters[p_name].constant(p_value)

    pathway_model.parameters().set(initial_parameters)

    # select optimiser
    algo = optimizer.get(chosen['algorithm'])
    algo.configure(algo_config[chosen['algorithm']])
    
    total_objective = optimizer.build_objective_function(pathway_model,
                                                         split['fit'],
                                                         objective_config)

    # make replica-provided parameters nan
    for name in ['H2O', 'CH4', 'CO2', 'TOC', 'DOC']:
        if name in initial_parameters:
            del initial_parameters._parameters[name] 

    run_config = {'model': pathway_model.get_config(only_structure = True),
                 'chosen': chosen,
                 'objective': objective_config,
                 'algo': algo.get_config(),
                 'range': initial_parameters.get_config(only_range = True)
                 }
    if not legacy_path is None:
        run_config.update({'legacy': legacy_path, 'legacy_file': legacy_file})
    
    existing_best = None
    existing_best_r2 = None
    if store_checkpoints:
        checkpoint_callback = checkpoint.CheckpointCallback(
            run_config, keep_only_n = checkpoint_keep_only_n,
            verbose = verbose_callback, target = cp_target,
            overwrite_existing = overwrite_existing_checkpoints)
        existing_best = checkpoint_callback.load_best_existing_checkpoint()

        if existing_best is not None:
            checkpoint_file, checkpoint_data = existing_best
            # Retain the range requested for this fit, but start from values
            # in the best compatible checkpoint.  This matters for local
            # minimisers and makes the loaded checkpoint a real incumbent,
            # rather than merely a number used for retention.
            saved_parameters = parameters.ModelParameters(
                checkpoint_data['parameters'])
            for saved_parameter in saved_parameters:
                if saved_parameter in initial_parameters:
                    initial_parameters[saved_parameter.name].set(
                        saved_parameter.value)
            pathway_model.parameters().set(initial_parameters)
            # Recreate the saved candidate's run log without incrementing the
            # fit's call count or running its top-level callbacks.  The
            # resulting R² combines the replica-specific logs of split fits.
            total_objective._call(initial_parameters, transformed = False)
            existing_best_r2 = total_objective.last_R2()
            # Command-line overrides define this fit and must take precedence
            # over the values saved by an earlier checkpoint.
            for p_name, p_value in chosen['parameter_override'].items():
                initial_parameters[p_name].constant(p_value)
            pathway_model.parameters().set(initial_parameters)
            print('Loaded existing best checkpoint '
                  f'{os.path.basename(checkpoint_file)} '
                  f'(loss {float(checkpoint_data["total_loss"]):.6g}, '
                  f'R2 {existing_best_r2:.2f})')

        total_objective.add_callback(checkpoint_callback)
    total_objective.add_callback(checkpoint.PrintCallback(
        run_config,
        initial_best_loss=(None if existing_best is None
                           else float(existing_best[1]['total_loss'])),
        initial_best_r2=existing_best_r2))

    if not minimize:
        total_objective(initial_parameters, transformed = False)
        print()

    else:    
        algo.minimize(total_objective, initial_parameters)
    
    run_log = total_objective.model().system_state_log

    return pathway_model, total_objective, run_log

if __name__ == '__main__':
    import user_input
    
    args = user_input.parse_args_fit()

    chosen = {
            'sample':                   args.sample,
            'validation_replica':       args.validation_replica, 
            
            't_start':                  args.t_start,
            't_end':                    args.t_end,
            
            'fit_mode':                 'single' if args.single else 'split',
            'pathways':                 ['Hydrolysis',
                                         'Fermentation',
                                         'Hydro',
                                         'Aceto',
                                         'Homo',
                                         'Fe3'],
            'parameter_override':       args.override,  # e.g. 'Homo_thermodynamics': False
            'normalized_parameters':    True,
            'algorithm':                'powell' if args.local else 'differential_evolution',
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
            'differential_evolution':   {'workers' : -1}, # empty dict uses default
            'powell':                   {}
                }
    
    if args.checkpoint is not None:
        # Keep the configured default ranges while taking every starting value
        # from the selected checkpoint.
        init_config = {'file': args.checkpoint, 'range': 'default'}
    else:
        init_config = {
                'default':                  args.default,
                'best_N':                   None if args.default else args.best,
                'sample':                   None if args.default else chosen['sample'],
                'validation_replica':       None if args.default else chosen['validation_replica'],
                'model':                    None,
                'run_ID':                   None,
                'file':                     None,
        }
    
    for omitted_pathway in args.omit:
        chosen['pathways'].remove(omitted_pathway)
    
    fit(chosen, objective_config, algo_config, init_config,
        store_checkpoints = not args.dry,
        overwrite_existing_checkpoints = args.overwrite_checkpoints)
    
    
