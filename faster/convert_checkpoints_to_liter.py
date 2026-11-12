"""Convert saved checkpoints from per-g-dry-weight to per-litre units.

The historical checkpoints store microbial pools (``M_*``) in mg/gdw and
other numeric pool, kinetic, and rate parameters in micromol/gdw (or the
corresponding per-gdw rate unit).  This command writes *new* checkpoint JSON
files in which those quantities are expressed per litre of pore water.

One conversion factor is needed for each measured replica::

    gdw / L = dry_weight_g / (water_content_ml * 1e-3)

Microbial pools also divide by ``MOLAR_MASS_BIOMASS`` (mg/micromol).  A joint
fit is consequently expanded to one converted file per saved replica: its
shared values are merged with the replica's local result before conversion.
This avoids assigning one litre-based value to a parameter that was shared
only in the old per-gdw representation.

Run from the repository root, for example::

    python faster/convert_checkpoints_to_liter.py model-3T3P

The default output directory is ``converted_checkpoints/<model-id>/``.  It is
outside ``results/`` so converted files are never accidentally reused as
optimisation inputs.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import defaultdict
from pathlib import Path

import CONSTANTS
import data
from converted_checkpoint import FORMAT, REPLAY_MODEL_VERSION, REPLAY_PATHWAYS
from model_v0_4.model import Model
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY


# Dimensionless, categorical, and inverse-time parameters do not acquire a
# volume term.  ``*_v_max`` is a specific rate: its old denominator was mg of
# biomass, while the v0.4 denominator is micromol of biomass.  Every other
# numeric checkpoint parameter is a pool, concentration, or per-gdw rate and
# is converted by the replica factor.
UNCHANGED_NAMES = frozenset({'death_rate'})
BIOMASS_BASIS_NAMES = frozenset({'Hydrolysis_Kmb'})


def is_unchanged(name, value):
    """Return whether *value* is not expressed per gdw."""
    return (isinstance(value, bool) or name in UNCHANGED_NAMES
            or name.endswith('_CUE') or name.endswith('_thermodynamics')
            )


def uses_biomass_basis(name):
    """Whether a parameter is expressed in the biomass basis, not per gdw."""
    # Hydrolysis_Kmb is placed in the M_Ferm row of Pathway.Km.  It is a
    # biomass half-saturation constant despite its historic parameter name.
    return name.startswith('M_') or name in BIOMASS_BASIS_NAMES


def conversion_factor(dry_weight_g, water_content_ml, name):
    """Return the old-unit to litre-unit multiplier for one parameter."""
    if dry_weight_g <= 0 or water_content_ml <= 0:
        raise ValueError('Replica dry weight and water content must both be positive.')
    gdw_per_litre = dry_weight_g / (water_content_ml * 1e-3)
    if uses_biomass_basis(name):
        return gdw_per_litre / CONSTANTS.MOLAR_MASS_BIOMASS
    if name.endswith('_v_max'):
        return CONSTANTS.MOLAR_MASS_BIOMASS
    return gdw_per_litre


def convert_values(values, dry_weight_g, water_content_ml):
    """Convert a parameter-value mapping and return it with conversion metadata."""
    converted = {}
    converted_names = []
    unchanged_names = []
    for name, value in values.items():
        if is_unchanged(name, value):
            converted[name] = value
            unchanged_names.append(name)
        elif isinstance(value, (int, float)):
            converted[name] = value * conversion_factor(
                dry_weight_g, water_content_ml, name)
            converted_names.append(name)
        else:
            raise ValueError(f'Cannot convert non-numeric parameter {name!r}: {value!r}')
    return converted, {
        'dry_weight_g': dry_weight_g,
        'water_content_ml': water_content_ml,
        'gdw_per_litre': dry_weight_g / (water_content_ml * 1e-3),
        'biomass_molar_mass_mg_per_micromol': CONSTANTS.MOLAR_MASS_BIOMASS,
        'converted_parameters': sorted(converted_names),
        'unchanged_parameters': sorted(unchanged_names),
    }


def convert_declared_ranges(range_config, dry_weight_g, water_content_ml):
    """Convert persisted search ranges alongside checkpoint values.

    Checkpoint ``run_config`` stores ranges as serialised ``Parameter``
    configurations.  Keeping those ranges in their old units would make a
    converted checkpoint internally inconsistent, even though the envelope
    report below is intentionally based on the fitted values.
    """
    converted = copy.deepcopy(range_config)
    if not isinstance(converted, dict):
        return converted
    for name, config in converted.items():
        if not isinstance(config, dict):
            continue
        bounds = config.get('range')
        if (is_unchanged(name, config.get('value')) or not isinstance(bounds, list)
                or len(bounds) != 2):
            continue
        config['range'] = [
            bound * conversion_factor(dry_weight_g, water_content_ml, name)
            for bound in bounds
        ]
    return converted


def replica_data_by_id(dataset):
    """Index loaded experiment data by the checkpoint's (sample, replica) ID."""
    replicas = {}
    for sample in dataset.samples:
        for replica in sample.replicas:
            replicas[(str(sample.sample_name), str(replica.replica_number))] = replica
    return replicas


def checkpoint_model_id(document):
    """Read the canonical generated model ID without trusting a filename."""
    # Import lazily to keep the conversion helpers usable in small unit tests.
    import hashing
    return hashing.build_model_id(document['run_config']['model'])


def matching_checkpoints(results_directory, model_id):
    """Yield parseable checkpoints whose complete model ID contains *model_id*."""
    selector = model_id.casefold().removeprefix('model-')
    for path in Path(results_directory).rglob('*'):
        if not path.is_file() or path.name.startswith('.'):
            continue
        try:
            with path.open() as handle:
                document = json.load(handle)
            found_model_id = checkpoint_model_id(document)
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
        if selector in found_model_id.casefold().removeprefix('model-'):
            yield path, document, found_model_id


def best_checkpoint_per_run(checkpoints):
    """Keep only the lowest-loss checkpoint from every complete run.

    Checkpoint files are periodic optimiser snapshots.  Their run ID is based
    on the complete run configuration, whereas the loss is the quantity that
    ranks snapshots within that run.  A lexical path tie-break makes output
    deterministic if two files have exactly the same loss.
    """
    import hashing
    best = {}
    for path, document, model_id in checkpoints:
        try:
            run_id = hashing.build_run_id(document['run_config'])
            loss = float(document['total_loss'])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f'Checkpoint {path} has no usable run ID or loss.') from error
        candidate = (loss, str(path), path, document, model_id)
        previous = best.get(run_id)
        if previous is None or candidate[:2] < previous[:2]:
            best[run_id] = candidate
    return [(candidate[2], candidate[3], candidate[4])
            for candidate in sorted(best.values(), key=lambda item: item[1])]


def expanded_solutions(document):
    """Yield ``(sample, replica, parameters, source_kind)`` for a checkpoint."""
    chosen = document['run_config']['chosen']
    if chosen.get('fit_mode') != 'multi':
        yield (str(chosen['sample']), str(chosen['validation_replica']),
               document['parameters'], 'checkpoint')
        return

    shared = document['parameters']
    records = document.get('replicas', [])
    if not isinstance(records, list) or not records:
        raise ValueError('Joint checkpoint has no saved per-replica results.')
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get('parameters'), dict):
            continue
        sample = str(record.get('sample'))
        replica = str(record.get('replica_number', record.get('replica')))
        if sample == 'None' or replica == 'None':
            raise ValueError('Joint checkpoint replica result is missing its sample or replica ID.')
        values = dict(shared)
        values.update(record['parameters'])
        yield sample, replica, values, 'joint_replica'


def converted_documents(document, source_path, replicas):
    """Yield converted document copies and their sample/replica identities."""
    for sample, replica_number, values, source_kind in expanded_solutions(document):
        replica = replicas.get((sample, replica_number))
        if replica is None:
            raise ValueError(f'No experiment data for checkpoint replica {sample}/{replica_number}.')
        converted_values, metadata = convert_values(
            values, replica.dry_weight, replica.water_content)
        converted = copy.deepcopy(document)
        converted['parameters'] = converted_values
        if 'range' in converted.get('run_config', {}):
            converted['run_config']['range'] = convert_declared_ranges(
                converted['run_config']['range'], replica.dry_weight,
                replica.water_content)
        # A retained joint-results list is misleading now that this document
        # represents exactly one replica-specific, litre-based solution.
        if source_kind == 'joint_replica':
            converted.pop('replicas', None)
        # The retained run_config records the source fit, whose model is often
        # v0.3.  Converted parameter values instead require this explicit
        # v0.4.2 replay contract; consumers must not infer it from run_config.
        converted['converted_checkpoint'] = {
            'format': FORMAT,
            'source_model_version': converted['run_config']['model'].get('version'),
            'replay_model_version': REPLAY_MODEL_VERSION,
            'pathways': list(REPLAY_PATHWAYS),
            'normalized_parameters': False,
        }
        converted['unit_conversion'] = {
            'source_checkpoint': str(source_path),
            'source_kind': source_kind,
            'replica': {'sample': sample, 'replica_number': replica_number},
            'source_units': {
                'M_*': 'mg/gdw',
                'other_numeric_parameters': 'micromol/gdw or compatible',
            },
            'target_units': {
                'M_*': 'micromol/L',
                'other_numeric_parameters': 'micromol/L or compatible',
            },
            **metadata,
        }
        yield converted, sample, replica_number


def default_bounds():
    """Return model-v0.4 physical parameter bounds in the target unit system."""
    return {
        parameter.name: (parameter.low, parameter.high)
        for parameter in Model.default_model_parameters(normalize=False)
        if parameter.low is not None and parameter.high is not None
    }


def envelope_rows(values_by_parameter):
    """Compare converted value envelopes with model-v0.4 default bounds."""
    defaults = default_bounds()
    rows = []
    for name in sorted(values_by_parameter):
        values = values_by_parameter[name]
        low, high = min(values), max(values)
        default_low, default_high = defaults.get(name, (None, None))
        if default_low is None:
            comparison = 'no_default_bounds'
        elif low >= default_low and high <= default_high:
            comparison = 'inside_default_bounds'
        elif high < default_low or low > default_high:
            comparison = 'outside_default_bounds'
        else:
            comparison = 'overlaps_default_bounds'
        rows.append({
            'parameter': name,
            'converted_low': low,
            'converted_high': high,
            'default_low': default_low,
            'default_high': default_high,
            'comparison': comparison,
            'value_count': len(values),
        })
    return rows


def write_envelope_report(rows, output_directory):
    """Write machine-readable and spreadsheet-friendly bounds comparisons."""
    with (output_directory / 'parameter_envelopes.json').open('w') as handle:
        json.dump(rows, handle, indent=2, sort_keys=True)
    fields = ('parameter', 'converted_low', 'converted_high', 'default_low',
              'default_high', 'comparison', 'value_count')
    with (output_directory / 'parameter_envelopes.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model_id', help='Full model ID or an unambiguous model-ID substring.')
    parser.add_argument('--results', type=Path, default=Path(RESULTS_DIRECTORY),
                        help='Checkpoint root to scan (default: results/).')
    parser.add_argument('--output-dir', type=Path,
                        help='Destination (default: converted_checkpoints/<model-id>/).')
    args = parser.parse_args(argv)

    discovered = list(matching_checkpoints(args.results, args.model_id))
    matches = best_checkpoint_per_run(discovered)
    if not matches:
        parser.error(f'No checkpoints match model ID {args.model_id!r} under {args.results}.')
    model_ids = {model_id for _, _, model_id in matches}
    if len(model_ids) != 1:
        parser.error('Model-ID selector is ambiguous: ' + ', '.join(sorted(model_ids)))
    model_id = next(iter(model_ids))
    output_directory = args.output_dir or Path(PROJECT_DIRECTORY) / 'converted_checkpoints' / model_id
    output_directory.mkdir(parents=True, exist_ok=True)

    replicas = replica_data_by_id(data.get_data_before_day())
    values_by_parameter = defaultdict(list)
    written = 0
    for source_path, document, _ in matches:
        for converted, sample, replica_number in converted_documents(document, source_path, replicas):
            suffix = f'__liter__{sample}_{replica_number}'
            destination = output_directory / (source_path.name + suffix)
            if destination.exists():
                raise FileExistsError(f'Refusing to overwrite existing converted checkpoint: {destination}')
            with destination.open('w') as handle:
                json.dump(converted, handle, indent=2, sort_keys=True)
            for name, value in converted['parameters'].items():
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    values_by_parameter[name].append(value)
            written += 1

    rows = envelope_rows(values_by_parameter)
    write_envelope_report(rows, output_directory)
    print(f'Converted {written} replica-specific checkpoint(s) from {len(matches)} best checkpoint(s), '
          f'selected from {len(discovered)} matching source checkpoint(s).')
    print(f'Output: {output_directory}')
    print(f'Envelope comparison: {output_directory / "parameter_envelopes.csv"}')


if __name__ == '__main__':
    main()
