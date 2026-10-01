import hashlib
import base64
import json
import math

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

def compute_hash(config):
    assert isinstance(config, dict)
    payload = _canonical_json(freeze(config)).encode('utf-8')
    return hashlib.blake2b(payload, digest_size = 16).digest()

def build_id(config, group_length = 3, length = 3):
    digest = compute_hash(config)
    encoded = base64.b32encode(digest).decode('ascii').rstrip('=')
    s_hash = encoded[:length*group_length]
    
    identifier = "-".join(
                        s_hash[i:i+group_length] 
                        for i in range(0, len(s_hash), group_length)
                    )
    return identifier

def build_run_id(config):
    cp = config.copy() # shallow copy!
    ignore = ['legacy', 'legacy_file', ]
    _ = [cp.pop(i,None) for i in ignore]
    return 'run-' + build_id(cp, 3,2)

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
        
def build_model_id(config):
    add_missing_thermodynamics_switch(config)
    return 'model-' + build_id(config, 4,2)

def build_loss_id(config):
    # TODO: compatibility layer here
    return 'loss-' + build_id(config, 4,2)

def build_checkpoint_id(config):
    return 'cp-' + build_id(config, 4,2)

if __name__ == '__main__':
    d = {'b': 456, 'a': 123, 'c': set([1,2,4])}
    print(build_id(d))
