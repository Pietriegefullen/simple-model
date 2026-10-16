"""Fit many replicas with shared and replica-local model parameters.

Each outer differential-evolution candidate fixes the requested shared
parameters.  It then creates one model per selected replica and fits every
remaining (local) parameter independently with an inner differential
evolution search.  The outer loss is the sum of those inner best losses.

The script deliberately evaluates outer candidates serially.  This is needed
to retain the local best parameter vectors that belong to an outer candidate;
the existing multiprocessing DE interface only publishes aggregate losses.
"""

import argparse
import fnmatch
import math
import os

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
    sample_filter = None if not samples else {str(value) for value in samples}
    replica_filter = None if not replicas else {str(value) for value in replicas}
    site_filter = None if not sites else {str(value).casefold() for value in sites}
    origin_filter = None if not origins else {str(value).casefold() for value in origins}
    if depth is not None and depth[0] > depth[1]:
        raise ValueError('--depth requires MIN <= MAX.')

    selected = []
    for sample in dataset.samples:
        if sample_filter is not None and str(sample.sample_name) not in sample_filter:
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
                        local_workers=1):
    """Print the exact model, data, parameter, and nested-objective setup."""
    print('\nJOINT REPLICA FIT SETUP')
    print('=' * 24)
    print('Replicas in fit (' + str(len(replicas)) + '): ' + ', '.join(
        replica_identifier(replica) for replica in replicas))
    print('Active pathways: ' + ', '.join(pathways))
    print('Omitted pathways: ' + (', '.join(omitted_pathways) if omitted_pathways else 'none'))
    print('Shared-bound checkpoints: ' + str(len(qualifying_checkpoints))
          + ' qualifying single-replica fits')
    print('Local initialization: ' + (
        f'best {local_best} compatible single-replica checkpoints per replica'
        if local_best else 'default model bounds'))
    print(f'Inner-DE call budget: {local_call_budget} per replica')
    print(f'Inner-DE workers: {local_workers}')
    print('Inner-DE early-stop R²: ' + (
        f'{local_r2_threshold:g}' if local_r2_threshold is not None else 'disabled'))
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
            role = 'global (outer DE)'
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
            role = 'local (inner DE)'
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


class LocalR2ThresholdReached(Exception):
    """A replica fit already meets its requested quality threshold."""


class StopWhenLocalR2Reached:
    """Stop the current inner differential-evolution run at a target R²."""
    def __init__(self, threshold):
        self.threshold = threshold
        self.objective = None

    def set_objective(self, objective):
        self.objective = objective

    def __call__(self):
        r2_value = self.objective.last_R2()
        if (r2_value is not None and np.isfinite(r2_value)
                and r2_value >= self.threshold):
            raise LocalR2ThresholdReached()


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
                 local_algorithm, local_best=None, checkpoint_source=None,
                 prepared_local_parameters=None, local_r2_threshold=None):
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
        self.local_best = local_best
        self.checkpoint_source = checkpoint_source
        self.model_id = hashing.build_model_id(
            build_parameter_template(pathways, normalized_parameters).get_config(only_structure=True))
        self._last_results = []
        self._last_weighted_mse = (None, None)
        self.prepared_local_parameters = prepared_local_parameters
        self.local_r2_threshold = local_r2_threshold

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
                'algorithm': 'differential_evolution',
            },
            'objective': self.objective_config,
            'algo': self.local_algorithm,
            'range': local_parameters.get_config(only_range=True),
            'multi_replica': {
                'role': 'inner',
                'outer_shared_values': shared_values,
                'replica': replica_identifier(replica),
            },
        }

    def set_parameters(self, parameter_list, values, transformed=False):
        for parameter, value in zip(parameter_list, np.atleast_1d(values)):
            parameter.set(parameter.inverse_transform(value) if transformed else value)

    def _call(self, values, transformed=True):
        shared = self.model().parameters()
        self.set_parameters(shared.variables(), values, transformed=transformed)
        results = []
        total_loss = 0.0
        residual_total = 0.0
        variance_total = 0.0
        metrics_available = True

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
            local_objective.add_callback(inner_checkpoints)
            local_objective.add_callback(checkpoint.PrintCallback(None, header = str(replica)))
            if self.local_r2_threshold is not None:
                local_objective.add_callback(
                    StopWhenLocalR2Reached(self.local_r2_threshold))
            local_optimizer = optimizer.get('differential_evolution')
            local_optimizer.configure(self.local_algorithm)
            # Keep the local best vector available in this process.  The
            # existing DE multiprocessing mode shares only aggregate scores.
            threshold_reached = False
            try:
                local_optimizer._minimize(local_objective, local_parameters)
            except LocalR2ThresholdReached:
                threshold_reached = True
            stored_inner = inner_checkpoints.load_best_existing_checkpoint()
            if stored_inner is None:
                raise RuntimeError(
                    f'Inner DE for replica {replica_identifier(replica)} '
                    'finished without a checkpoint.')
            checkpoint_path, checkpoint_data = stored_inner
            best_parameters = parameters.ModelParameters(checkpoint_data['parameters'])
            best_loss = local_objective._call(best_parameters, transformed=False)
            local_r2 = local_objective.last_R2()
            mse_residual, mse_total = local_objective.last_weighted_mse()
            if mse_residual is None or mse_total is None:
                metrics_available = False
            else:
                residual_total += mse_residual
                variance_total += mse_total
            total_loss += best_loss
            results.append({
                'replica': replica_identifier(replica),
                'sample': str(replica.sample.sample_name),
                'replica_number': str(replica.replica_number),
                'loss': float(best_loss),
                'r2': (None if local_r2 is None or not np.isfinite(local_r2)
                       else float(local_r2)),
                'r2_threshold_reached': threshold_reached,
                'checkpoint': checkpoint_path,
                'parameters': {
                    parameter.name: parameter.value
                    for parameter in local_model.parameters()
                    if parameter.name not in self.shared_names
                },
            })

        self._last_results = results
        self._last_weighted_mse = ((residual_total, variance_total)
                                   if metrics_available else (None, None))
        return total_loss

    def last_weighted_mse(self):
        return self._last_weighted_mse

    def last_R2(self):
        residual, total = self.last_weighted_mse()
        return None if residual is None or total == 0 else 1 - residual / total

    def checkpoint_metadata(self):
        return {'replicas': self._last_results}


def fit(replicas, shared_parameters, shared_names, local_template,
        fixed_parameters=None, pathways=PATHWAYS, normalized_parameters=True,
        objective=None, outer_algorithm=None, local_algorithm=None,
        local_best=None, checkpoint_source=None, store_checkpoints=True,
        checkpoint_target=None, prepared_local_parameters=None,
        local_r2_threshold=None):
    """Run the nested differential-evolution fit and return its objective."""
    if objective is None:
        objective = objective_config()
    fixed_parameters = {} if fixed_parameters is None else fixed_parameters
    outer_algorithm = {} if outer_algorithm is None else dict(outer_algorithm)
    local_algorithm = {} if local_algorithm is None else dict(local_algorithm)
    if outer_algorithm.get('workers', 1) != 1:
        raise ValueError(
            'Nested shared fitting requires outer DE workers=1 so each '
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
        pathways, normalized_parameters, objective, local_algorithm, local_best,
        checkpoint_source, prepared_local_parameters, local_r2_threshold)
    run_config = {
        'model': build_parameter_template(pathways, normalized_parameters).get_config(
            only_structure=True),
        # The legacy checkpoint directory naming convention requires these
        # fields.  The real selection is preserved in ``multi_replica``.
        'chosen': {'sample': 0, 'validation_replica': 0, 't_start': None,
                   't_end': None, 'fit_mode': 'multi', 'pathways': list(pathways),
                   'parameter_override': fixed_parameters,
                   'normalized_parameters': normalized_parameters,
                   'algorithm': 'differential_evolution'},
        'objective': objective,
        'algo': outer_algorithm,
        'range': shared_parameters.get_config(only_range=True),
        'multi_replica': {'shared': list(shared_names),
                          'replicas': [replica_identifier(replica) for replica in replicas],
                          'local_initialization': ('best' if local_best else 'default'),
                          'local_best': local_best},
    }
    if store_checkpoints:
        outer_objective.add_callback(checkpoint.CheckpointCallback(
            run_config, target=checkpoint_target))
    outer_objective.add_callback(checkpoint.PrintCallback(run_config))

    outer_optimizer = optimizer.get('differential_evolution')
    outer_optimizer.configure(outer_algorithm)
    outer_optimizer.minimize(outer_objective, shared_parameters)
    return outer_objective


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='fit_replicas',
        description='Jointly fit replicas with shared outer and local inner DE searches.')
    parser.add_argument('--shared', nargs='+', required=True,
                        help='Shared parameter names or globs, e.g. Hydro_v_max M_*.')
    parser.add_argument('--p', nargs='+', dest='fixed',
                        help='Fixed parameter/value pairs, as in fit_sample.')
    parser.add_argument('--omit', choices=PATHWAYS, nargs='+', default=[],
                        help='Pathways to omit from the model.')
    parser.add_argument('--samples', nargs='+', help='Restrict to sample IDs.')
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
    parser.add_argument('--seed', type=int, help='Random seed for --random and DE.')
    parser.add_argument('--r2-threshold', type=float, default=0.9,
                        help='Use checkpoint values only when recomputed R² is above this (default: 0.9).')
    local_initialization = parser.add_mutually_exclusive_group()
    local_initialization.add_argument('--local-default', action='store_true',
                                      help='Use default bounds for every local parameter (default).')
    local_initialization.add_argument('--local-best', type=int, metavar='N',
                                      help='Derive each replica local bounds from its best N single-fit checkpoints.')
    parser.add_argument('--checkpoint-source',
                        help='Directory containing existing single-fit checkpoints.')
    parser.add_argument('--outer-max-evaluations', type=int,
                        help='Stop outer DE after this many completed candidate evaluations.')
    parser.add_argument('--local-call-budget', '--local-max-evaluations',
                        dest='local_max_evaluations', type=int, default=3000,
                        help=('Maximum completed inner-DE evaluations for each '
                              'individual replica (default: 3000).'))
    parser.add_argument('--local-r2-threshold', type=float,
                        help=('Stop an individual replica inner-DE fit once its '
                              'R² reaches this value, e.g. 0.7.'))
    parser.add_argument('--local-workers', type=int, default=1,
                        help=('Worker processes for each inner DE; use -1 for '
                              'all available cores (default: 1).'))
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
    if args.local_workers == 0 or args.local_workers < -1:
        parser.error('--local-workers must be 1, -1, or a positive integer.')
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
    shared_parameters, qualifying = shared_bounds_from_checkpoints(
        shared_names, pathways=pathways, objective=objective_config(),
        r2_threshold=args.r2_threshold,
        source_directory=args.checkpoint_source)
    template = build_parameter_template(pathways, True).parameters()
    model_id = hashing.build_model_id(
        build_parameter_template(pathways, True).get_config(only_structure=True))
    prepared_local_parameters = prepare_local_parameters(
        replicas, template, shared_names, args.fixed, model_id,
        args.local_best, args.checkpoint_source)
    print_setup_summary(
        replicas, pathways, args.omit, shared_names, shared_parameters,
        prepared_local_parameters, args.fixed, objective_config(), qualifying,
        args.local_best, args.local_max_evaluations, args.local_r2_threshold,
        args.local_workers)
    if args.dry:
        return replicas, shared_parameters, qualifying
    if not args.yes and not confirm_setup():
        return None

    outer_algorithm = {'workers': 1, 'seed': args.seed}
    local_algorithm = {'workers': args.local_workers, 'seed': args.seed}
    if args.outer_max_evaluations is not None:
        outer_algorithm['max_evaluations'] = args.outer_max_evaluations
    local_algorithm['max_evaluations'] = args.local_max_evaluations
    return fit(
        replicas, shared_parameters, shared_names, template, args.fixed,
        pathways=pathways, objective=objective_config(), outer_algorithm=outer_algorithm,
        local_algorithm=local_algorithm, local_best=args.local_best,
        checkpoint_source=args.checkpoint_source, store_checkpoints=True,
        prepared_local_parameters=prepared_local_parameters,
        local_r2_threshold=args.local_r2_threshold)


if __name__ == '__main__':
    run(parse_args())
