"""Fit many replicas with shared and replica-local model parameters.

Each outer optimiser candidate fixes the requested shared
parameters.  It then creates one model per selected replica and fits every
remaining (local) parameter independently with an inner differential
evolution or CMA-ES search.  The outer loss is the sum of those inner best losses.

The script deliberately evaluates outer candidates serially.  This is needed
to retain the local best parameter vectors that belong to an outer candidate;
the existing multiprocessing DE interface only publishes aggregate losses.
"""

import argparse
import fnmatch
import math
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import checkpoint
import data
import hashing
import model
import optimizer
import parameters
import user_input


PATHWAYS = ['Hydrolysis', 'Fermentation', 'Hydro', 'Aceto', 'Homo', 'Fe3']
REPLICA_PROVIDED_PARAMETERS = frozenset({'H2O', 'TOC', 'DOC', 'CH4', 'CO2'})


def objective_config():
    """Return the objective definition used by ``fit_sample.py``."""
    return {
        'loss_weight': {'CO2': 1.0, 'CH4': 1.0},
        'reduction': {'CO2': 'mse', 'CH4': 'mse'},
        'transform': {'CO2': ['normalize'], 'CH4': ['log', 'normalize']},
    }


def copy_parameters(source):
    """Copy parameters without retaining mutable ``Parameter`` instances."""
    return parameters.ModelParameters({
        parameter.name: parameter.get_config()
        for parameter in source
    })


ALGORITHM_CHOICES = {'de': 'differential_evolution', 'cma-es': 'cma_es'}


def parameters_x0(saved, parameter_set):
    """Return ``saved`` values in the transformed coordinates of a new fit."""
    values = []
    for parameter in parameter_set.variables():
        if parameter.name not in saved:
            return None
        value = saved[parameter.name].value
        if (not np.isfinite(value) or value < parameter.lower()
                or value > parameter.upper()):
            return None
        values.append(parameter.transform(value))
    return np.asarray(values)


def checkpoint_x0(checkpoint_data, parameter_set):
    """Return a checkpoint candidate in DE's transformed coordinates.

    ``None`` means the checkpoint cannot safely seed this (possibly narrowed)
    search space.
    """
    return parameters_x0(parameters.ModelParameters(checkpoint_data['parameters']),
                         parameter_set)


def select_shared_default_parameters(reference, shared_names):
    """Keep the declared full fitting ranges for the selected shared names."""
    shared = copy_parameters(reference)
    for name in list(shared.as_dict()):
        if name not in shared_names:
            del shared.as_dict()[name]
    return shared


def build_parameter_template(pathways, normalized_parameters=True):
    """Build a model and apply every declared default fitting range.

    ``Model.build`` creates several initial-state pools as constants before
    the defaults are applied.  Fitting must instead start from the declared
    parameter catalogue, just as ``fit_sample`` does; otherwise valid local
    parameters such as biomass and initial Fe3 would silently stay constant.
    """
    pathway_model = model.build_model(pathways, normalized_parameters)
    defaults = parameters.default_model_parameters(
        normalize_parameters=normalized_parameters)
    pathway_model.parameters().set(defaults)
    return pathway_model


def expand_parameter_names(requested, available):
    """Expand exact names and shell-style globs against model parameters."""
    if not requested:
        raise ValueError('Specify at least one parameter with --shared.')
    available = list(available)
    expanded = []
    for pattern in requested:
        matches = ([pattern] if pattern in available else
                   [name for name in available if fnmatch.fnmatchcase(name, pattern)])
        if not matches:
            raise ValueError(f'--shared pattern {pattern!r} matches no model parameter.')
        expanded.extend(matches)
    return list(dict.fromkeys(expanded))


def replica_identifier(replica):
    return f'{replica.sample.sample_name}{replica.replica_number}'


def select_replicas(dataset, samples=None, replicas=None, sites=None, origins=None,
                    depth=None, split=None, random_count=None, seed=None):
    """Select replicas using sample metadata and unambiguous replica IDs.

    A replica filter accepts either a replica number (``4``) or a full ID
    (``13514``).  ``split`` is a sample/replica pair and is convenient when
    reproducing a particular old ``fit_sample`` split.
    """
    sample_patterns = None if not samples else [str(value) for value in samples]
    replica_filter = None if not replicas else {str(value) for value in replicas}
    site_filter = None if not sites else {str(value).casefold() for value in sites}
    origin_filter = None if not origins else {str(value).casefold() for value in origins}
    if depth is not None and depth[0] > depth[1]:
        raise ValueError('--depth requires MIN <= MAX.')

    selected = []
    for sample in dataset.samples:
        if (sample_patterns is not None
                and not any(fnmatch.fnmatchcase(str(sample.sample_name), pattern)
                            for pattern in sample_patterns)):
            continue
        if site_filter is not None and str(sample.site).casefold() not in site_filter:
            continue
        # Historical usage calls Samoylov/Kurugnakh an "origin", while the
        # dataset stores those labels as ``site``.  Match both here; ``--site``
        # remains available when a strict site filter is preferred.
        if origin_filter is not None and not {
                str(sample.origin).casefold(), str(sample.site).casefold()
        } & origin_filter:
            continue
        if depth is not None and (sample.depth is None
                                  or not depth[0] <= float(sample.depth) <= depth[1]):
            continue
        for replica in sample.replicas:
            full_id = replica_identifier(replica)
            if replica_filter is not None and (str(replica.replica_number) not in replica_filter
                                               and full_id not in replica_filter):
                continue
            if split is not None and full_id != f'{split[0]}{split[1]}':
                continue
            selected.append(replica)

    if random_count is not None:
        if random_count > len(selected):
            raise ValueError(
                f'--random requested {random_count} replicas, but only {len(selected)} match.')
        generator = np.random.default_rng(seed)
        indices = sorted(generator.choice(len(selected), size=random_count, replace=False))
        selected = [selected[index] for index in indices]

    if not selected:
        raise ValueError('No replicas match the requested filters.')
    return selected


def _checkpoint_matches_model(run_config, model_id):
    try:
        return hashing.build_model_id(run_config['model']) == model_id
    except (KeyError, TypeError, ValueError):
        return False


def single_replica_checkpoint_candidates(model_id, source_directory=None):
    """Return compatible checkpoints produced by ``fit_sample --single``."""
    candidates = parameters.get_candidates({}, source_directory=source_directory)
    selected = []
    for candidate in candidates:
        _, _, run_config, _ = candidate
        chosen = run_config.get('chosen', {})
        if (str(chosen.get('fit_mode')) == 'single'
                and _checkpoint_matches_model(run_config, model_id)):
            selected.append(candidate)
    return selected


def checkpoint_r2(candidate, pathways, normalized_parameters, objective):
    """Re-evaluate a single-replica checkpoint to obtain its missing R²."""
    candidate_parameters, _, run_config, _ = candidate
    chosen = run_config['chosen']
    dataset = data.get_data_before_day()
    sample = dataset[str(chosen['sample'])]
    replica = sample[str(chosen['validation_replica'])]
    pathway_model = model.build_model(pathways, normalized_parameters)
    pathway_model.parameters().set(candidate_parameters)
    candidate_objective = optimizer.build_objective_function(
        pathway_model, [replica], objective)
    candidate_objective._call(candidate_parameters, transformed=False)
    return candidate_objective.last_R2()


def shared_bounds_from_checkpoints(shared_names, pathways=PATHWAYS,
                                   normalized_parameters=True,
                                   objective=None, r2_threshold=0.9,
                                   source_directory=None):
    """Derive shared search ranges from qualifying single-replica fits.

    Checkpoints predate R² storage, so R² is recomputed using the current
    model and objective.  Bounds are the extrema of every qualifying value.
    """
    if objective is None:
        objective = objective_config()
    reference_model = build_parameter_template(pathways, normalized_parameters)
    model_id = hashing.build_model_id(reference_model.get_config(only_structure=True))
    reference = reference_model.parameters()
    for name in shared_names:
        if name not in reference:
            raise ValueError(f'Unknown shared parameter {name!r}.')

    qualifying = []
    for candidate in single_replica_checkpoint_candidates(model_id, source_directory):
        try:
            r2_value = checkpoint_r2(candidate, pathways, normalized_parameters, objective)
        except Exception as error:
            print(f'Skipping checkpoint {candidate[3]}: could not evaluate R² ({error}).')
            continue
        if r2_value is not None and np.isfinite(r2_value) and r2_value > r2_threshold:
            qualifying.append((candidate, float(r2_value)))

    if not qualifying:
        raise ValueError(
            'No compatible single-replica checkpoints exceed '
            f'--r2-threshold {r2_threshold:g}.')

    shared = copy_parameters(reference)
    for name in list(shared.as_dict()):
        if name not in shared_names:
            del shared.as_dict()[name]

    for name in shared_names:
        values = [candidate[0][name].value for candidate, _ in qualifying
                  if name in candidate[0] and np.isfinite(candidate[0][name].value)]
        if len(values) != len(qualifying):
            raise ValueError(f'Qualifying checkpoints do not all contain {name!r}.')
        low, high = min(values), max(values)
        if np.isclose(low, high):
            raise ValueError(
                f'Checkpoint values for shared parameter {name!r} do not span a range.')
        parameter = shared[name]
        if parameter.options is not None:
            raise ValueError(f'Shared parameter {name!r} is categorical and cannot use DE.')
        parameter.variable(parameter.inverse_transform(
            (parameter.transform(low) + parameter.transform(high)) / 2), (low, high))
        # ``Parameter`` caches its normalization transform.  It must be
        # rebuilt after narrowing the physical range.
        parameter.transformer = None
    return shared, qualifying


def _local_candidates(replica, model_id, count, source_directory=None):
    candidates = single_replica_checkpoint_candidates(model_id, source_directory)
    matching = [candidate for candidate in candidates if (
        str(candidate[2]['chosen'].get('sample')) == str(replica.sample.sample_name)
        and str(candidate[2]['chosen'].get('validation_replica')) == str(replica.replica_number))]
    matching.sort(key=lambda candidate: candidate[1])
    selected = matching[:count]
    if len(selected) < count:
        raise ValueError(
            f'Replica {replica_identifier(replica)} has only {len(selected)} compatible '
            f'single-replica checkpoints; --local-best {count} was requested.')
    return selected


def local_initial_parameters(replica, template, shared_names, fixed_parameters,
                             model_id, best_count=None, source_directory=None):
    """Build the per-replica inner search space."""
    local = copy_parameters(template)
    for name, value in fixed_parameters.items():
        local[name].constant(value)
    for name in shared_names:
        # The outer objective supplies the actual value immediately before
        # each inner optimisation.
        local[name].constant(local[name].value)

    if best_count is None:
        return local

    candidates = _local_candidates(replica, model_id, best_count, source_directory)
    local_names = [parameter.name for parameter in local.variables()]
    for name in local_names:
        values = [candidate[0][name].value for candidate in candidates
                  if name in candidate[0]]
        if len(values) != len(candidates):
            raise ValueError(
                f'Local checkpoint candidates for {replica_identifier(replica)} '
                f'do not all contain {name!r}.')
        low, high = min(values), max(values)
        if np.isclose(low, high):
            local[name].constant(low)
        else:
            local[name].variable(local[name].inverse_transform(
                (local[name].transform(low) + local[name].transform(high)) / 2),
                (low, high))
            local[name].transformer = None
    return local


def prepare_local_parameters(replicas, template, shared_names, fixed_parameters,
                             model_id, best_count=None, source_directory=None):
    """Create and validate each replica's local search space once."""
    return {
        replica_identifier(replica): local_initial_parameters(
            replica, template, shared_names, fixed_parameters, model_id,
            best_count, source_directory)
        for replica in replicas
    }


def _parameter_display(parameter):
    if parameter.is_variable():
        return f'[{parameter.lower():.6g}, {parameter.upper():.6g}]'
    if parameter.options is not None:
        return '{' + ', '.join(map(str, sorted(parameter.options, key=str))) + '}'
    return str(parameter.value)


def print_setup_summary(replicas, pathways, omitted_pathways, shared_names,
                        shared_parameters, local_parameters, fixed_parameters,
                        objective, qualifying_checkpoints, local_best=None,
                        local_call_budget=3000, local_r2_threshold=None,
                        local_loss_threshold=None,
                        local_workers=1, resume=True,
                        outer_algorithm='differential_evolution',
                        local_algorithm='differential_evolution',
                        local_warm_start=True,
                        shared_bounds_source='declared model ranges'):
    """Print the exact model, data, parameter, and nested-objective setup."""
    def algorithm_label(name):
        return {'differential_evolution': 'DE', 'cma_es': 'CMA-ES'}.get(name, name)

    outer_label = algorithm_label(outer_algorithm)
    local_label = algorithm_label(local_algorithm)
    print('\nJOINT REPLICA FIT SETUP')
    print('=' * 24)
    print('Replicas in fit (' + str(len(replicas)) + '): ' + ', '.join(
        replica_identifier(replica) for replica in replicas))
    print('Active pathways: ' + ', '.join(pathways))
    print('Omitted pathways: ' + (', '.join(omitted_pathways) if omitted_pathways else 'none'))
    print('Shared parameter bounds: ' + shared_bounds_source)
    if qualifying_checkpoints:
        print('Shared-bound checkpoints: ' + str(len(qualifying_checkpoints))
              + ' qualifying single-replica fits')
    print('Local initialization: ' + (
        f'best {local_best} compatible single-replica checkpoints per replica'
        if local_best else 'default model bounds'))
    print('Outer algorithm: ' + outer_label)
    print('Inner algorithm: ' + local_label)
    print(f'Inner call budget: {local_call_budget} per replica')
    print(f'Inner workers: {local_workers}')
    print('Checkpoint continuation: ' + ('enabled' if resume else 'disabled'))
    print('Inner warm start from per-replica incumbent: ' +
          ('enabled' if local_warm_start else 'disabled'))
    print('Inner early-stop R²: ' + (
        f'{local_r2_threshold:g}' if local_r2_threshold is not None else 'disabled'))
    print('Inner loss guard / outer-candidate rejection threshold: ' + (
        f'{local_loss_threshold:g}' if local_loss_threshold is not None else 'disabled'))
    if fixed_parameters:
        print('Fixed parameters (--p): ' + ', '.join(
            f'{name}={value}' for name, value in sorted(fixed_parameters.items())))
    else:
        print('Fixed parameters (--p): none')

    reference = next(iter(local_parameters.values()))
    rows = []
    counts = {'global': 0, 'local': 0, 'const': 0}
    for parameter in reference:
        name = parameter.name
        if name in shared_names:
            role = 'global (outer ' + outer_label + ')'
            detail = _parameter_display(shared_parameters[name])
            counts['global'] += 1
        elif name in fixed_parameters:
            role = 'constant'
            detail = str(fixed_parameters[name])
            counts['const'] += 1
        elif name in REPLICA_PROVIDED_PARAMETERS:
            role = 'constant (data)'
            detail = 'derived separately for each replica'
        elif parameter.is_variable():
            role = 'local (inner ' + local_label + ')'
            detail = ('replica-specific bounds from best checkpoints'
                      if local_best else _parameter_display(parameter))
            counts['local'] += 1
        else:
            # These are fixed by the selected model structure (currently the
            # thermodynamics switches and any parameter with no declared
            # fitting range), not by the shared/local fitting policy.
            role = 'constant (model fixed)'
            detail = _parameter_display(parameter)
        rows.append((name, role, detail))

    widths = [max(len(header), *(len(row[index]) for row in rows))
              for index, header in enumerate(('parameter', 'role', 'range / value'))]
    print('\nMODEL PARAMETERS')
    print('  '.join(header.ljust(widths[index]) for index, header in enumerate(
        ('parameter', 'role', 'range / value'))))
    print('  '.join('-' * width for width in widths))
    for row in rows:
        print('  '.join(value.ljust(widths[index]) for index, value in enumerate(row)))

    print(f'global: {counts["global"]} parameters')
    print(f'local:  {counts["local"]} parameters')
    print(f'const:  {counts["const"]} parameters')

    print('\nOBJECTIVE STRUCTURE')
    print('  Per replica r:')
    print('    L_r = MSE(normalize(CO2_model), normalize(CO2_data))')
    print('        + MSE(normalize(log(CH4_model)), normalize(log(CH4_data)))')
    print('  Inner optimization: θ_local,r*(g) = argmin_θ_local,r L_r(g, θ_local,r)')
    print('  Outer optimization: L(g) = Σ_r L_r(g, θ_local,r*(g))')
    print('  Global parameters g are constant during every inner fit; each replica')
    print('  has its own model instance and independent local parameter vector.')
    print()


def confirm_setup(input_function=input):
    """Require an explicit acknowledgement before beginning nested DE."""
    try:
        answer = input_function('Start the optimisation with this setup? [y/N] ')
    except EOFError as error:
        raise RuntimeError('Optimisation was not confirmed.') from error
    if answer.strip().casefold() not in {'y', 'yes'}:
        print('Optimisation cancelled.')
        return False
    return True


def confirm_continue_after_inner_interrupt(input_function=input):
    """Ask whether a Ctrl-C should stop the outer run or skip one replica."""
    try:
        answer = input_function(
            'Inner optimisation interrupted. Continue with the next replica? [Y/n] ')
    except EOFError:
        # Non-interactive jobs retain the documented default: continue.
        return True
    return answer.strip().casefold() in {'', 'y', 'yes'}


class LocalR2ThresholdReached(Exception):
    """A replica fit already meets its requested quality threshold."""


class LocalLossThresholdReached(Exception):
    """A replica evaluation exceeded its requested loss guardrail."""


class StopWhenLocalR2Reached(checkpoint.Callback):
    """Stop the current inner optimisation run at a target R²."""
    def __init__(self, threshold):
        super().__init__(None)
        self.threshold = threshold
        self.objective = None

    def set_objective(self, objective):
        self.objective = objective

    def __call__(self):
        r2_value = self.objective.last_R2()
        if (r2_value is not None and np.isfinite(r2_value)
                and r2_value >= self.threshold):
            raise LocalR2ThresholdReached()


class StopWhenLocalLossExceeds(checkpoint.Callback):
    """Stop an inner search immediately after a prohibitively large loss."""
    def __init__(self, threshold):
        super().__init__(None)
        self.threshold = threshold
        self.objective = None

    def set_objective(self, objective):
        self.objective = objective

    def __call__(self):
        loss = self.objective.last_call()[-1]
        if np.isfinite(loss) and loss > self.threshold:
            raise LocalLossThresholdReached()


class _SharedModel:
    """Small model facade required by the existing optimiser interface."""
    def __init__(self, shared_parameters):
        self._parameters = shared_parameters

    def parameters(self):
        return self._parameters

    def __str__(self):
        return 'Shared parameters\n' + str(self._parameters)


class SharedReplicaObjective(optimizer.Addable):
    """Outer objective that performs a local fit for each selected replica."""
    def __init__(self, replicas, shared_parameters, local_template, shared_names,
                 fixed_parameters, pathways, normalized_parameters, objective,
                 local_algorithm, local_algorithm_name='differential_evolution',
                 local_best=None, checkpoint_source=None,
                 prepared_local_parameters=None, local_r2_threshold=None,
                 local_loss_threshold=None, resume=True, local_warm_start=True):
        super().__init__()
        self._model = _SharedModel(shared_parameters)
        self.replicas = list(replicas)
        self.local_template = local_template
        self.shared_names = tuple(shared_names)
        self.fixed_parameters = dict(fixed_parameters)
        self.pathways = list(pathways)
        self.normalized_parameters = normalized_parameters
        self.objective_config = objective
        self.local_algorithm = dict(local_algorithm)
        self.local_algorithm_name = local_algorithm_name
        self.local_best = local_best
        self.checkpoint_source = checkpoint_source
        self.model_id = hashing.build_model_id(
            build_parameter_template(pathways, normalized_parameters).get_config(only_structure=True))
        self._last_results = []
        self._last_weighted_mse = (None, None)
        self._aborted_by_loss_guard = False
        self.prepared_local_parameters = prepared_local_parameters
        self.local_r2_threshold = local_r2_threshold
        self.local_loss_threshold = local_loss_threshold
        self.resume = resume
        self.local_warm_start = local_warm_start
        self._local_warm_starts = {}
        self._local_pool = None

    def _parallel_resources(self):
        """Create one reusable worker pool for all inner fits."""
        workers = self.local_algorithm.get('workers', 1)
        if workers == 1:
            return None, None
        if self._local_pool is None:
            self._local_pool = ProcessPoolExecutor(
                max_workers=None if workers == -1 else workers)
        return None, self._local_pool

    def close_parallel_resources(self):
        """Release the persistent resources even if outer DE raises."""
        if self._local_pool is not None:
            self._local_pool.shutdown()
            self._local_pool = None

    def reset_parallel_resources(self):
        """Discard stalled inner workers after a user interrupts one fit."""
        if self._local_pool is not None:
            self._local_pool.shutdown(wait=False, cancel_futures=True)
            self._local_pool = None

    def _inner_run_config(self, local_model, local_parameters, replica, shared):
        """Build an ordinary annotated-checkpoint configuration for one inner fit."""
        shared_values = {name: shared[name].value for name in self.shared_names}
        return {
            'model': local_model.get_config(only_structure=True),
            'chosen': {
                'sample': int(replica.sample.sample_name),
                'validation_replica': int(replica.replica_number),
                't_start': None,
                't_end': None,
                'fit_mode': 'single',
                'pathways': list(self.pathways),
                'parameter_override': {**self.fixed_parameters, **shared_values},
                'normalized_parameters': self.normalized_parameters,
                'algorithm': self.local_algorithm_name,
            },
            'objective': self.objective_config,
            'algo': self.local_algorithm,
            'range': local_parameters.get_config(only_range=True),
            'multi_replica': {
                'role': 'inner',
                'outer_shared_values': shared_values,
                'replica': replica_identifier(replica),
                'r2_threshold': self.local_r2_threshold,
                'loss_threshold': self.local_loss_threshold,
            },
        }

    def set_parameters(self, parameter_list, values, transformed=False):
        for parameter, value in zip(parameter_list, np.atleast_1d(values)):
            parameter.set(parameter.inverse_transform(value) if transformed else value)

    def _call(self, values, transformed=True):
        shared = self.model().parameters()
        self.set_parameters(shared.variables(), values, transformed=transformed)
        results = []
        loss_contributions = []
        #total_loss = 0.0
        residual_total = 0.0
        variance_total = 0.0
        metrics_available = True
        self._aborted_by_loss_guard = False

        for replica in self.replicas:
            if self.prepared_local_parameters is None:
                local_parameters = local_initial_parameters(
                    replica, self.local_template, self.shared_names, self.fixed_parameters,
                    self.model_id, self.local_best, self.checkpoint_source)
            else:
                local_parameters = copy_parameters(
                    self.prepared_local_parameters[replica_identifier(replica)])
            for name in self.shared_names:
                local_parameters[name].constant(shared[name].value)
            local_model = model.build_model(self.pathways, self.normalized_parameters)
            local_model.parameters().set(local_parameters)
            local_objective = optimizer.build_objective_function(
                local_model, [replica], self.objective_config)
            inner_run_config = self._inner_run_config(
                local_model, local_parameters, replica, shared)
            # Each run ID encodes its replica and fixed shared values.  Keep
            # its best local result in the normal annotated-checkpoint format
            # so the parent can recover it after parallel worker evaluation.
            inner_run_id = hashing.build_run_id(inner_run_config)
            inner_target = os.path.join(
                checkpoint.CP_ROOT, 'joint_inner', inner_run_id,
                replica_identifier(replica))
            inner_checkpoints = checkpoint.CheckpointCallback(
                inner_run_config, keep_only_n=1, target=inner_target)
            inner_cma_state = checkpoint.CMAStateCheckpoint(
                inner_run_config, target=inner_target)
            local_objective.add_callback(inner_checkpoints)
            # With a parallel pool, callbacks execute in isolated worker
            # copies of this objective.  Their apparent "best" values are
            # only best for one worker, not for the whole local population.
            # Keep live progress for serial fits; parallel fits print the
            # verified best checkpoint below in the parent process instead.
            show_live_local_progress = self.local_algorithm.get('workers', 1) == 1
            if show_live_local_progress:
                local_objective.add_callback(
                    checkpoint.PrintCallback(None, header=' ' * 3 + str(replica)))

            if self.local_r2_threshold is not None:
                local_objective.add_callback(
                    StopWhenLocalR2Reached(self.local_r2_threshold))
            if self.local_loss_threshold is not None:
                local_objective.add_callback(
                    StopWhenLocalLossExceeds(self.local_loss_threshold))
            local_algorithm = dict(self.local_algorithm)
            seed = None
            if self.resume:
                existing_inner = inner_checkpoints.load_best_existing_checkpoint()
                if existing_inner is not None:
                    seed = checkpoint_x0(existing_inner[1], local_parameters)
            if seed is None and self.local_warm_start:
                seed = parameters_x0(
                    self._local_warm_starts.get(replica_identifier(replica), {}),
                    local_parameters)
            if seed is not None:
                local_algorithm['x0'] = seed
            local_optimizer = optimizer.get(self.local_algorithm_name)
            local_optimizer.configure(local_algorithm)
            if self.local_algorithm_name == 'cma_es':
                local_optimizer.set_state_checkpoint(inner_cma_state)
            local_manager, local_pool = self._parallel_resources()
            if local_manager is not None:
                local_optimizer.set_process_manager(local_manager)
            if local_pool is not None:
                local_optimizer.set_worker_map(local_pool.map)
            # Keep the local best vector available in this process.  The
            # existing DE multiprocessing mode shares only aggregate scores.
            threshold_reached = False
            loss_threshold_reached = False
            interrupted = False
            try:
                local_optimizer.minimize(local_objective, local_parameters, verbose = False)
            except LocalR2ThresholdReached:
                threshold_reached = True
            except LocalLossThresholdReached:
                loss_threshold_reached = True
            except KeyboardInterrupt:
                interrupted = True
                # ``Pool.map`` may have queued stalled work.  It cannot be
                # safely reused after Ctrl-C, so terminate it before moving
                # to the next replica; it is recreated lazily on demand.
                self.reset_parallel_resources()
                if not confirm_continue_after_inner_interrupt():
                    raise
            stored_inner = inner_checkpoints.load_best_existing_checkpoint()
            if stored_inner is None:
                if not interrupted:
                    raise RuntimeError(
                        f'Inner DE for replica {replica_identifier(replica)} '
                        'finished without a checkpoint.')
                # No candidate completed before Ctrl-C.  Continue fitting the
                # other replicas, but make this outer candidate uncompetitive
                # rather than silently omitting a replica from its loss.
                checkpoint_path = None
                best_parameters = local_parameters
                best_loss = optimizer.INVALID_OBJECTIVE_PENALTY
                local_r2 = None
                mse_residual, mse_total = None, None
            else:
                checkpoint_path, checkpoint_data = stored_inner
                best_parameters = parameters.ModelParameters(checkpoint_data['parameters'])
                if self.local_warm_start:
                    self._local_warm_starts[replica_identifier(replica)] = copy_parameters(
                        best_parameters)
                best_loss = local_objective._call(best_parameters, transformed=False)
                local_r2 = local_objective.last_R2()
                mse_residual, mse_total = local_objective.last_weighted_mse()
            if mse_residual is None or mse_total is None:
                metrics_available = False
            else:
                residual_total += mse_residual
                variance_total += mse_total
            loss_contributions.append(best_loss)
            #total_loss += best_loss
            results.append({
                'replica': replica_identifier(replica),
                'sample': str(replica.sample.sample_name),
                'replica_number': str(replica.replica_number),
                'loss': float(best_loss),
                'r2': (None if local_r2 is None or not np.isfinite(local_r2)
                       else float(local_r2)),
                'r2_threshold_reached': threshold_reached,
                'loss_threshold_reached': loss_threshold_reached,
                'interrupted': interrupted,
                'checkpoint': checkpoint_path,
                'parameters': {
                    parameter.name: parameter.value
                    for parameter in local_model.parameters()
                    if parameter.name not in self.shared_names
                },
            })
            if show_live_local_progress:
                print()
            else:
                r2_text = 'n/a' if local_r2 is None or not np.isfinite(local_r2) else f'{local_r2:.2f}'
                print(f'   {replica}: local best loss: {best_loss:8.3g}, '
                      f'best R2 = {r2_text}')
            if loss_threshold_reached:
                # This candidate has explicitly failed the local loss guard.
                # Its partial sum cannot stand in for the omitted replicas:
                # that could make a failed candidate look artificially good.
                # Return the same finite rejection penalty used for invalid
                # model evaluations and let the outer optimiser advance.
                self._last_results = results
                self._last_weighted_mse = (None, None)
                self._aborted_by_loss_guard = True
                print('   local loss guard reached; rejecting this outer candidate')
                return optimizer.INVALID_OBJECTIVE_PENALTY

        self._last_results = results
        self._last_weighted_mse = ((residual_total, variance_total)
                                   if metrics_available else (None, None))
        total_loss = np.mean(loss_contributions)
        return total_loss

    def last_weighted_mse(self):
        return self._last_weighted_mse

    def last_R2(self):
        residual, total = self.last_weighted_mse()
        return None if residual is None or total == 0 else 1 - residual / total

    def checkpoint_metadata(self):
        return {'replicas': self._last_results,
                'aborted_by_local_loss_guard': self._aborted_by_loss_guard}


def fit(replicas, shared_parameters, shared_names, local_template,
        fixed_parameters=None, pathways=PATHWAYS, normalized_parameters=True,
        objective=None, outer_algorithm=None, local_algorithm=None,
        outer_algorithm_name='differential_evolution',
        local_algorithm_name='differential_evolution',
        local_best=None, checkpoint_source=None, store_checkpoints=True,
        checkpoint_target=None, prepared_local_parameters=None,
        local_r2_threshold=None, local_loss_threshold=None, resume=True,
        local_warm_start=True):
    """Run the nested shared/local fit and return its outer objective."""
    if objective is None:
        objective = objective_config()
    fixed_parameters = {} if fixed_parameters is None else fixed_parameters
    outer_algorithm = {} if outer_algorithm is None else dict(outer_algorithm)
    local_algorithm = {} if local_algorithm is None else dict(local_algorithm)
    if outer_algorithm.get('workers', 1) != 1:
        raise ValueError(
            'Nested shared fitting requires outer workers=1 so each '
            'candidate retains its replica results.')
    local_workers = local_algorithm.get('workers', 1)
    if (not isinstance(local_workers, (int, np.integer)) or local_workers == 0
            or local_workers < -1):
        raise ValueError('Inner DE workers must be 1, -1, or a positive integer.')
    outer_algorithm['workers'] = 1
    local_algorithm['workers'] = local_workers
    # A nested run can otherwise inherit SciPy's very large generation
    # default for every replica and every outer candidate.
    if local_algorithm.get('max_evaluations') is None:
        local_algorithm['max_evaluations'] = 3000

    model_id = hashing.build_model_id(
        build_parameter_template(pathways, normalized_parameters).get_config(only_structure=True))
    if prepared_local_parameters is None:
        prepared_local_parameters = prepare_local_parameters(
            replicas, local_template, shared_names, fixed_parameters, model_id,
            local_best, checkpoint_source)
    outer_objective = SharedReplicaObjective(
        replicas, shared_parameters, local_template, shared_names, fixed_parameters,
        pathways, normalized_parameters, objective, local_algorithm,
        local_algorithm_name=local_algorithm_name, local_best=local_best,
        checkpoint_source=checkpoint_source,
        prepared_local_parameters=prepared_local_parameters,
        local_r2_threshold=local_r2_threshold,
        local_loss_threshold=local_loss_threshold, resume=resume,
        local_warm_start=local_warm_start)
    run_config = {
        'model': build_parameter_template(pathways, normalized_parameters).get_config(
            only_structure=True),
        # The legacy checkpoint directory naming convention requires these
        # fields.  The real selection is preserved in ``multi_replica``.
        'chosen': {'sample': 0, 'validation_replica': 0, 't_start': None,
                   't_end': None, 'fit_mode': 'multi', 'pathways': list(pathways),
                   'parameter_override': fixed_parameters,
                   'normalized_parameters': normalized_parameters,
                   'algorithm': outer_algorithm_name},
        'objective': objective,
        'algo': dict(outer_algorithm),
        'range': shared_parameters.get_config(only_range=True),
        'multi_replica': {'shared': list(shared_names),
                          'replicas': [replica_identifier(replica) for replica in replicas],
                          'local_initialization': ('best' if local_best else 'default'),
                          'local_best': local_best,
                          'local_algorithm': local_algorithm_name,
                          'local_algo': dict(local_algorithm),
                          'local_r2_threshold': local_r2_threshold,
                          'local_loss_threshold': local_loss_threshold,
                          'local_warm_start': local_warm_start},
    }
    outer_checkpoints = None
    outer_cma_state = None
    if store_checkpoints:
        if checkpoint_target is None:
            # The legacy directory key only includes one sample/replica.  A
            # joint fit needs the full run ID to avoid mixing distinct replica
            # selections or shared-parameter sets during continuation.
            checkpoint_target = os.path.join(
                checkpoint.CP_ROOT, 'joint_multi', hashing.build_run_id(run_config))
        outer_checkpoints = checkpoint.CheckpointCallback(
            run_config, target=checkpoint_target)
        outer_cma_state = checkpoint.CMAStateCheckpoint(
            run_config, target=checkpoint_target)
        outer_objective.add_callback(outer_checkpoints)
        if resume:
            existing_outer = outer_checkpoints.load_best_existing_checkpoint()
            if existing_outer is None:
                # One-time compatibility path for joint fits saved before the
                # dedicated ``joint_multi`` directory existed.  Never use a
                # merely nearby legacy checkpoint: it must have the exact
                # full run identity, including its replica selection.
                legacy_checkpoints = checkpoint.CheckpointCallback(run_config)
                legacy_outer = legacy_checkpoints.load_best_existing_checkpoint()
                if (legacy_outer is not None
                        and hashing.build_run_id(legacy_outer[1]['run_config'])
                        == hashing.build_run_id(run_config)):
                    existing_outer = legacy_outer
            if existing_outer is not None:
                seed = checkpoint_x0(existing_outer[1], shared_parameters)
                if seed is not None:
                    outer_algorithm['x0'] = seed
                    print('Resuming from existing joint checkpoint '
                          + os.path.basename(existing_outer[0]))
    outer_objective.add_callback(checkpoint.PrintCallback(run_config, carriage_return = False, header = '\nTotal: ', footer = '\n\n'))

    outer_optimizer = optimizer.get(outer_algorithm_name)
    outer_optimizer.configure(outer_algorithm)
    if outer_algorithm_name == 'cma_es' and outer_cma_state is not None:
        outer_optimizer.set_state_checkpoint(outer_cma_state)
    try:
        outer_optimizer.minimize(outer_objective, shared_parameters)
    finally:
        outer_objective.close_parallel_resources()
    return outer_objective


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='fit_replicas',
        description='Jointly fit replicas with independently selectable outer and inner DE/CMA-ES searches.')
    parser.add_argument('--shared', nargs='+', required=True,
                        help='Shared parameter names or globs, e.g. Hydro_v_max M_*.')
    parser.add_argument('--p', nargs='+', dest='fixed',
                        help='Fixed parameter/value pairs, as in fit_sample.')
    parser.add_argument('--omit', choices=PATHWAYS, nargs='+', default=[],
                        help='Pathways to omit from the model.')
    parser.add_argument('--samples', nargs='+',
                        help='Restrict to sample IDs or globs, e.g. 13*; combines with --site.')
    parser.add_argument('--replicas', nargs='+', help='Replica number or full sample+replica ID.')
    parser.add_argument('--split', nargs=2, metavar=('SAMPLE', 'REPLICA'),
                        help='Select one particular sample/replica split.')
    parser.add_argument('--site', '--location', dest='sites', nargs='+',
                        help='Filter site, e.g. Samoylov or Kurugnakh.')
    parser.add_argument('--origin', nargs='+', help='Filter origin, e.g. core or cliff.')
    parser.add_argument('--depth', nargs=2, type=float, metavar=('MIN', 'MAX'),
                        help='Keep samples with depth in this inclusive range.')
    parser.add_argument('--random', type=int, metavar='N',
                        help='Uniformly select N replicas after applying filters.')
    parser.add_argument('--seed', type=int, help='Random seed for --random and either optimiser.')
    parser.add_argument('--r2-threshold', type=float, default=0.9,
                        help='Checkpoint-bound filter: require recomputed R² above this (default: 0.9).')
    parser.add_argument('--shared-checkpoint-bounds', action='store_true',
                        help=('Restrict shared ranges to extrema of qualifying '
                              'single-fit checkpoints; default uses full declared model ranges.'))
    local_initialization = parser.add_mutually_exclusive_group()
    local_initialization.add_argument('--local-default', action='store_true',
                                      help='Use default bounds for every local parameter (default).')
    local_initialization.add_argument('--local-best', type=int, metavar='N',
                                      help='Derive each replica local bounds from its best N single-fit checkpoints.')
    parser.add_argument('--checkpoint-source',
                        help='Directory containing existing single-fit checkpoints.')
    parser.add_argument('--outer-max-evaluations', type=int,
                        help='Stop the outer search after this many completed evaluations.')
    parser.add_argument('--outer-algorithm', choices=ALGORITHM_CHOICES,
                        default='de', help='Outer search algorithm (default: de).')
    parser.add_argument('--local-algorithm', choices=ALGORITHM_CHOICES,
                        default='de', help='Per-replica search algorithm (default: de).')
    parser.add_argument('--outer-popsize', type=int,
                        help='Population size for outer DE or CMA-ES.')
    parser.add_argument('--local-popsize', type=int,
                        help='Population size for each local DE or CMA-ES fit.')
    parser.add_argument('--outer-cma-sigma', type=float,
                        help='Initial CMA-ES step size in unit-bound coordinates (default: 0.30).')
    parser.add_argument('--local-cma-sigma', type=float,
                        help='Initial local CMA-ES step size in unit-bound coordinates (default: 0.30).')
    parser.add_argument('--outer-cma-restarts', type=int, default=0,
                        help='Number of CMA-ES restarts after its first outer run (default: 0).')
    parser.add_argument('--local-cma-restarts', type=int, default=0,
                        help='Number of CMA-ES restarts after each first local run (default: 0).')
    parser.add_argument('--local-call-budget', '--local-max-evaluations',
                        dest='local_max_evaluations', type=int, default=3000,
                        help=('Maximum completed inner-optimiser evaluations for each '
                              'individual replica (default: 3000).'))
    parser.add_argument('--local-r2-threshold', type=float,
                        help=('Stop an individual replica inner-DE fit once its '
                              'R² reaches this value, e.g. 0.7.'))
    parser.add_argument('--local-loss-threshold', type=float,
                        help=('Abort an individual inner fit and reject its entire '
                              'outer candidate after an evaluated loss exceeds this value.'))
    parser.add_argument('--local-workers', type=int, default=1,
                        help=('Worker processes for each inner DE; use -1 for '
                              'all available cores (default: 1).'))
    parser.set_defaults(local_warm_start=True)
    parser.add_argument('--no-local-warm-start', dest='local_warm_start',
                        action='store_false',
                        help='Do not seed a replica fit from its best earlier inner result.')
    parser.set_defaults(resume=True)
    parser.add_argument('--resume', dest='resume', action='store_true',
                        help=('Continue compatible CMA-ES outer and inner searches from their '
                              'saved optimiser state; other searches use saved best candidates (default).'))
    parser.add_argument('--no-resume', dest='resume', action='store_false',
                        help='Do not use saved optimiser-state or best-candidate checkpoints.')
    parser.add_argument('--dry', action='store_true',
                        help='Select inputs and derive bounds, but do not run or save a fit.')
    parser.add_argument('--yes', action='store_true',
                        help='Start without the interactive confirmation prompt.')
    args = parser.parse_args(argv)
    if args.random is not None and args.random <= 0:
        parser.error('--random must be positive.')
    if args.local_best is not None and args.local_best <= 0:
        parser.error('--local-best must be positive.')
    if not -math.inf < args.r2_threshold <= 1:
        parser.error('--r2-threshold must be at most 1.')
    for name in ('outer_max_evaluations', 'local_max_evaluations'):
        value = getattr(args, name)
        if value is not None and value <= 0:
            parser.error(f'--{name.replace("_", "-")} must be positive.')
    if (args.local_r2_threshold is not None
            and not -math.inf < args.local_r2_threshold <= 1):
        parser.error('--local-r2-threshold must be at most 1.')
    if (args.local_loss_threshold is not None
            and (not np.isfinite(args.local_loss_threshold)
                 or args.local_loss_threshold < 0)):
        parser.error('--local-loss-threshold must be finite and non-negative.')
    if args.local_workers == 0 or args.local_workers < -1:
        parser.error('--local-workers must be 1, -1, or a positive integer.')
    for name in ('outer_popsize', 'local_popsize'):
        value = getattr(args, name)
        if value is not None and value <= 0:
            parser.error('--' + name.replace('_', '-') + ' must be positive.')
    for name in ('outer_cma_sigma', 'local_cma_sigma'):
        value = getattr(args, name)
        if value is not None and (not np.isfinite(value) or value <= 0):
            parser.error('--' + name.replace('_', '-') + ' must be finite and positive.')
    for name in ('outer_cma_restarts', 'local_cma_restarts'):
        if getattr(args, name) < 0:
            parser.error('--' + name.replace('_', '-') + ' cannot be negative.')
    args.fixed = user_input.parse_override(args.fixed)
    return args


def run(args):
    dataset = data.get_data_before_day()
    replicas = select_replicas(
        dataset, args.samples, args.replicas, args.sites, args.origin, args.depth,
        args.split, args.random, args.seed)
    pathways = [pathway for pathway in PATHWAYS if pathway not in args.omit]
    if not pathways:
        raise ValueError('At least one pathway must remain after --omit.')
    reference = build_parameter_template(pathways, True).parameters()
    shared_names = expand_parameter_names(args.shared, reference.as_dict())
    overlap = set(shared_names) & set(args.fixed)
    if overlap:
        raise ValueError('--shared and --p cannot both name ' + ', '.join(sorted(overlap)))
    unavailable = set(shared_names) & REPLICA_PROVIDED_PARAMETERS
    if unavailable:
        raise ValueError('Replica-provided parameters cannot be shared: '
                         + ', '.join(sorted(unavailable)))
    unknown_fixed = set(args.fixed) - set(reference.as_dict())
    if unknown_fixed:
        raise ValueError('Unknown --p parameter(s): ' + ', '.join(sorted(unknown_fixed)))
    if args.shared_checkpoint_bounds:
        shared_parameters, qualifying = shared_bounds_from_checkpoints(
            shared_names, pathways=pathways, objective=objective_config(),
            r2_threshold=args.r2_threshold,
            source_directory=args.checkpoint_source)
        shared_bounds_source = (
            f'qualifying checkpoint extrema (R² > {args.r2_threshold:g})')
    else:
        shared_parameters = select_shared_default_parameters(reference, shared_names)
        qualifying = []
        shared_bounds_source = 'full declared model ranges (default)'
    template = build_parameter_template(pathways, True).parameters()
    model_id = hashing.build_model_id(
        build_parameter_template(pathways, True).get_config(only_structure=True))
    prepared_local_parameters = prepare_local_parameters(
        replicas, template, shared_names, args.fixed, model_id,
        args.local_best, args.checkpoint_source)
    outer_algorithm_name = ALGORITHM_CHOICES[args.outer_algorithm]
    local_algorithm_name = ALGORITHM_CHOICES[args.local_algorithm]
    print_setup_summary(
        replicas, pathways, args.omit, shared_names, shared_parameters,
        prepared_local_parameters, args.fixed, objective_config(), qualifying,
        args.local_best, args.local_max_evaluations, args.local_r2_threshold,
        args.local_loss_threshold, args.local_workers, args.resume,
        outer_algorithm_name, local_algorithm_name,
        args.local_warm_start, shared_bounds_source)
    if args.dry:
        return replicas, shared_parameters, qualifying
    #if not args.yes and not confirm_setup():
    #    return None

    outer_algorithm = {'workers': 1, 'seed': args.seed}
    local_algorithm = {'workers': args.local_workers, 'seed': args.seed}
    if args.outer_popsize is not None:
        outer_algorithm['popsize'] = args.outer_popsize
    if args.local_popsize is not None:
        local_algorithm['popsize'] = args.local_popsize
    if outer_algorithm_name == 'cma_es':
        if args.outer_cma_sigma is not None:
            outer_algorithm['sigma0'] = args.outer_cma_sigma
        outer_algorithm['restarts'] = args.outer_cma_restarts
    if local_algorithm_name == 'cma_es':
        if args.local_cma_sigma is not None:
            local_algorithm['sigma0'] = args.local_cma_sigma
        local_algorithm['restarts'] = args.local_cma_restarts
    if args.outer_max_evaluations is not None:
        outer_algorithm['max_evaluations'] = args.outer_max_evaluations
    local_algorithm['max_evaluations'] = args.local_max_evaluations
    return fit(
        replicas, shared_parameters, shared_names, template, args.fixed,
        pathways=pathways, objective=objective_config(), outer_algorithm=outer_algorithm,
        local_algorithm=local_algorithm, local_best=args.local_best,
        checkpoint_source=args.checkpoint_source, store_checkpoints=True,
        prepared_local_parameters=prepared_local_parameters,
        local_r2_threshold=args.local_r2_threshold,
        local_loss_threshold=args.local_loss_threshold, resume=args.resume,
        outer_algorithm_name=outer_algorithm_name,
        local_algorithm_name=local_algorithm_name,
        local_warm_start=args.local_warm_start)


if __name__ == '__main__':
    run(parse_args())
