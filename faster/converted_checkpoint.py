"""Recognition and replay support for litre-normalised checkpoint archives."""

from __future__ import annotations

import copy

import optimizer
import parameters
from model_factory import build_model


FORMAT = 'litre-normalized-v1'
REPLAY_MODEL_VERSION = '0.4.2'
REPLAY_PATHWAYS = ('Hydrolysis', 'Fermentation', 'Hydro', 'Aceto', 'Homo', 'Fe3')


def metadata_for_document(document):
    """Return converted-checkpoint metadata, including legacy archives.

    Earlier converted files have ``unit_conversion`` provenance but no replay
    declaration.  They remain usable with the historical v0.4.2 target.
    """
    conversion = document.get('unit_conversion') if isinstance(document, dict) else None
    declared = document.get('converted_checkpoint') if isinstance(document, dict) else None
    if not isinstance(conversion, dict):
        return None
    if conversion.get('target_units', {}).get('M_*') != 'micromol/L':
        return None
    if not isinstance(declared, dict):
        declared = {
            'format': FORMAT,
            'source_model_version': document.get('run_config', {}).get('model', {}).get('version'),
            'replay_model_version': REPLAY_MODEL_VERSION,
            'pathways': list(REPLAY_PATHWAYS),
            'normalized_parameters': False,
            'legacy_inferred': True,
        }
    return {'converted_checkpoint': copy.deepcopy(declared),
            'unit_conversion': copy.deepcopy(conversion)}


def replay_spec(metadata):
    """Return the explicit model settings required for a converted replay."""
    declared = (metadata or {}).get('converted_checkpoint')
    if not isinstance(declared, dict) or declared.get('format') != FORMAT:
        return None
    version = str(declared.get('replay_model_version', REPLAY_MODEL_VERSION))
    pathways = tuple(declared.get('pathways', REPLAY_PATHWAYS))
    return version, pathways, bool(declared.get('normalized_parameters', False))


def replica_for_conversion(conversion, dataset):
    """Find the one measured replica named by converted-file provenance."""
    identity = (conversion or {}).get('replica')
    if not isinstance(identity, dict) or not {'sample', 'replica_number'} <= set(identity):
        raise ValueError('Converted checkpoint has no replica provenance.')
    sample_id = str(identity['sample'])
    replica_id = str(identity['replica_number'])
    try:
        sample = dataset[sample_id]
    except (KeyError, TypeError) as error:
        raise ValueError(f'No data sample {sample_id!r} for converted checkpoint.') from error
    for replica in sample.replicas:
        if str(replica.replica_number) == replica_id:
            return replica
    raise ValueError(f'No replica {sample_id}/{replica_id} in the loaded experiment data.')


def replay(checkpoint, dataset, parameter_values=None, *, gas_limits=True):
    """Evaluate a converted checkpoint with its declared target model.

    Converted parameters are already expressed in the v0.4.2 litre-based
    system, so they must not be replayed through ``fit_sample.run`` using the
    source checkpoint's historical configuration.
    """
    metadata = checkpoint.metadata or {}
    spec = replay_spec(metadata)
    if spec is None:
        raise ValueError('Checkpoint is not a recognised converted checkpoint.')
    version, pathways, normalized = spec
    replica = replica_for_conversion(metadata.get('unit_conversion'), dataset)
    pathway_model = build_model(version, pathways, normalized=normalized)
    pathway_model.use_gas_limits = gas_limits
    values = checkpoint.parameters if parameter_values is None else parameter_values
    pathway_model.parameters().set(parameters.ModelParameters(values))
    objective = optimizer.build_objective_function(
        pathway_model, [replica], checkpoint.run_config['objective'])
    run_log = pathway_model.predict_with_logging(replica, reference='gdw', microbe_unit='mg')
    return pathway_model, objective, run_log
