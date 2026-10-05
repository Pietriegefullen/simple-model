import hashlib
import base64
import copy
import json
import math
import numbers


def _scalar(value):
    """Unbox NumPy-style scalar values without taking a NumPy dependency."""
    if not isinstance(value, (str, bytes)) and hasattr(value, 'item'):
        try:
            return value.item()
        except ValueError:
            pass
    return value


def as_bool(value):
    """Convert one of the supported boolean representations to ``bool``.

    ``bool('False')`` is deliberately not used here: Python considers every
    non-empty string true, which would make a malformed configuration hash as
    the opposite value without warning.
    """
    value = _scalar(value)
    if isinstance(value, bool):
        return value
    if isinstance(value, numbers.Real) and not isinstance(value, bool):
        if value == 0:
            return False
        if value == 1:
            return True
    if isinstance(value, str):
        normalized = value.strip().casefold()
        if normalized in {'true', '1', 'yes', 'on'}:
            return True
        if normalized in {'false', '0', 'no', 'off'}:
            return False
    raise ValueError(f'Expected a boolean value, got {value!r}')


def as_int(value):
    """Convert an integral number or integral numeric string to ``int``."""
    value = _scalar(value)
    if isinstance(value, bool):
        raise ValueError(f'Expected an integer value, got {value!r}')
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real) and math.isfinite(value) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        try:
            parsed = float(value.strip())
        except ValueError as error:
            raise ValueError(f'Expected an integer value, got {value!r}') from error
        if math.isfinite(parsed) and parsed.is_integer():
            return int(parsed)
    raise ValueError(f'Expected an integer value, got {value!r}')


def as_float(value):
    """Convert a numeric value or numeric string to ``float``."""
    value = _scalar(value)
    if isinstance(value, bool):
        raise ValueError(f'Expected a numeric value, got {value!r}')
    try:
        return float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f'Expected a numeric value, got {value!r}') from error


def as_float_or_variable(value):
    """Convert a model number while retaining the model's ``'variable'`` tag."""
    if isinstance(value, str) and value == 'variable':
        return value
    return as_float(value)


def as_thermodynamics_override(value, key=None):
    """Cast parameter overrides according to the parameter they address."""
    if key is not None and str(key).endswith('_thermodynamics'):
        return as_bool(value)
    return as_float(value)


def _cast_scalar(value, converter, path, key=None):
    if value is None:
        return None
    if converter is str:
        return str(value)
    try:
        if converter is as_thermodynamics_override:
            return converter(value, key)
        return converter(value)
    except (TypeError, ValueError) as error:
        rendered_path = '.'.join(map(str, path)) or '<root>'
        raise ValueError(f'Invalid value at {rendered_path}: {error}') from error


def cast_config(config, schema, path=(), key=None):
    """Return a deep, schema-cast copy of a configuration.

    A schema mirrors the relevant parts of a configuration.  Leaf values are
    converters (usually ``bool``, ``int``, ``float``, or ``str``); a ``'*'``
    entry applies to every otherwise unspecified key in a mapping.  Lists with
    one schema item describe homogeneous lists.  Values absent from the schema
    are copied unchanged, so new configuration fields remain part of a hash.
    """
    if config is None or schema is None:
        return copy.deepcopy(config)
    if schema is bool:
        return _cast_scalar(config, as_bool, path, key)
    if schema is int:
        return _cast_scalar(config, as_int, path, key)
    if schema is float:
        return _cast_scalar(config, as_float, path, key)
    if callable(schema):
        return _cast_scalar(config, schema, path, key)
    if isinstance(schema, dict):
        if not isinstance(config, dict):
            rendered_path = '.'.join(map(str, path)) or '<root>'
            raise ValueError(f'Expected a mapping at {rendered_path}, got {config!r}')
        wildcard = schema.get('*')
        return {
            item_key: cast_config(
                value, schema.get(item_key, wildcard), path + (item_key,), item_key)
            for item_key, value in config.items()
        }
    if isinstance(schema, list):
        if len(schema) != 1:
            raise ValueError('A list schema must contain exactly one item schema')
        if not isinstance(config, (list, tuple)):
            rendered_path = '.'.join(map(str, path)) or '<root>'
            raise ValueError(f'Expected a list at {rendered_path}, got {config!r}')
        return [cast_config(value, schema[0], path + (index,), index)
                for index, value in enumerate(config)]
    if isinstance(schema, tuple):
        if not isinstance(config, (list, tuple)) or len(config) != len(schema):
            rendered_path = '.'.join(map(str, path)) or '<root>'
            raise ValueError(f'Expected {len(schema)} values at {rendered_path}, got {config!r}')
        return [cast_config(value, item_schema, path + (index,), index)
                for index, (value, item_schema) in enumerate(zip(config, schema))]
    raise TypeError(f'Unsupported type schema {schema!r}')


MODEL_CONFIG_SCHEMA = {
    'version': str,
    '*': {
        'microbe': {
            'name': str,
            'death_rate': as_float,
            'v_max': as_float_or_variable,
            'Kmb': as_float_or_variable,
            'CUE': as_float_or_variable,
            'C_source': str,
        },
        'educts': {'*': {
            'name': str,
            'stoichiometry': as_float,
            'Km': as_float_or_variable,
            'inhibition': as_float_or_variable,
        }},
        'products': {'*': {
            'name': str,
            'stoichiometry': as_float,
            'Km': as_float_or_variable,
            'inhibition': as_float_or_variable,
        }},
        'use_thermodynamics': as_bool,
    },
}

OBJECTIVE_CONFIG_SCHEMA = {
    'loss_weight': {'*': as_float},
    'reduction': {'*': str},
    'transform': {'*': [str]},
}

CHOSEN_CONFIG_SCHEMA = {
    'sample': as_int,
    'validation_replica': as_int,
    't_start': as_float,
    't_end': as_float,
    'fit_mode': str,
    'pathways': [str],
    'parameter_override': {'*': as_thermodynamics_override},
    'normalized_parameters': as_bool,
    'algorithm': str,
}

ALGORITHM_CONFIG_SCHEMA = {
    'strategy': str,
    'updating': str,
    'popsize': as_int,
    'workers': as_int,
    'tol': as_float,
    'init': str,
    'polish': as_bool,
    'recombination': as_float,
    'mutation': [as_float],
}

PARAMETER_CONFIG_SCHEMA = {
    '*': {
        'name': str,
        'value': as_float,
        'range': [as_float],
        'scale': str,
        'normalize': as_bool,
    },
}

RUN_CONFIG_SCHEMA = {
    'model': MODEL_CONFIG_SCHEMA,
    'chosen': CHOSEN_CONFIG_SCHEMA,
    'objective': OBJECTIVE_CONFIG_SCHEMA,
    'algo': ALGORITHM_CONFIG_SCHEMA,
    'range': PARAMETER_CONFIG_SCHEMA,
}

CHECKPOINT_CONFIG_SCHEMA = {
    'parameters': {'*': as_thermodynamics_override},
    'total_loss': as_float,
    'run_config': RUN_CONFIG_SCHEMA,
}

def freeze(obj):
    if isinstance(obj, dict):
        return {str(k): freeze(v) for k, v in sorted(obj.items(), key = lambda kv: str(kv[0]))}
    if isinstance(obj, (list, tuple)):
        return [freeze(v) for v in obj]
    if isinstance(obj, (set, frozenset)):
        values = [freeze(v) for v in obj]
        return {'type': 'set', 'values': sorted(values, key = _canonical_json)}
    if isinstance(obj, float):
        if math.isnan(obj):
            return {'type': 'float', 'value': 'NaN'}
        if math.isinf(obj):
            return {'type': 'float', 'value': 'Infinity' if obj > 0 else '-Infinity'}
    if hasattr(obj, 'item'):
        return freeze(obj.item())
    return obj

def _canonical_json(obj):
    return json.dumps(obj, sort_keys = True, separators = (',', ':'),
                      ensure_ascii = True, allow_nan = False)

def compute_hash(config, schema=None):
    """Hash a configuration after applying its optional type schema."""
    assert isinstance(config, dict)
    payload = _canonical_json(freeze(cast_config(config, schema))).encode('utf-8')
    return hashlib.blake2b(payload, digest_size = 16).digest()

def build_id(config, group_length = 3, length = 3, schema=None):
    digest = compute_hash(config, schema)
    encoded = base64.b32encode(digest).decode('ascii').rstrip('=')
    s_hash = encoded[:length*group_length]
    
    identifier = "-".join(
                        s_hash[i:i+group_length] 
                        for i in range(0, len(s_hash), group_length)
                    )
    return identifier

def build_run_id(config):
    cp = copy.deepcopy(config)
    ignore = ['legacy', 'legacy_file', ]
    _ = [cp.pop(i,None) for i in ignore]
    return 'run-' + build_id(cp, 3, 2, RUN_CONFIG_SCHEMA)

def add_missing_thermodynamics_switch(config):
    import parameters
    pathway_keys = [k for k in config.keys() if not k == 'version']
    default_values = parameters.default_model_parameters()

    for pwy_name in pathway_keys:
        parameter_name = 'use_thermodynamics'
        if not parameter_name in config[pwy_name].keys():
            idx = [p.name for p in default_values].index(pwy_name + '_thermodynamics')
            default_p = default_values[idx]
            config[pwy_name][parameter_name] = default_p.value
        config[pwy_name][parameter_name] = as_bool(config[pwy_name][parameter_name])

def build_model_id(config):
    cp = copy.deepcopy(config)
    add_missing_thermodynamics_switch(cp)
    return 'model-' + build_id(cp, 4, 2, MODEL_CONFIG_SCHEMA)

def build_loss_id(config):
    return 'loss-' + build_id(config, 4, 2, OBJECTIVE_CONFIG_SCHEMA)

def build_checkpoint_id(config):
    return 'cp-' + build_id(config, 4, 2, CHECKPOINT_CONFIG_SCHEMA)

if __name__ == '__main__':
    d = {'b': 456, 'a': 123, 'c': set([1,2,4])}
    print(build_id(d))
