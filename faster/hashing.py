import hashlib
ALPHABET = 'abcdefghijklmnopqrstuvwxyz'.upper()

def freeze(obj):
    if isinstance(obj, dict):
        return tuple(sorted((k, freeze(v)) for k, v in obj.items()))
    if isinstance(obj, (list, tuple)):
        return tuple(freeze(v) for v in obj)
    if isinstance(obj, set):
        return tuple(freeze(v) for v in obj)
    return obj

def compute_hash(config):
    assert isinstance(config, dict)
    digest = hash(tuple(sorted(freeze(config))))
    return digest

def build_id(config, group_length = 3, length = 2):
    digest = compute_hash(config)
    
    number = abs(digest)

    result = []
    base = len(ALPHABET)

    while number:
        number, remainder = divmod(number, base)
        result.append(ALPHABET[remainder])

    s_hash = "".join(reversed(result)).zfill(length*group_length)[-length*group_length:]
    
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
    return 'model-' + build_id(config, 3,1)

def build_loss_id(config):
    # TODO: compatibility layer here
    return 'loss-' + build_id(config, 3,1)

def build_checkpoint_id(config):
    return 'cp-' + build_id(config, 3,2)

if __name__ == '__main__':
    d = {'b': 456, 'a': 123, 'c': set([1,2,4])}
    print(build_id(d))
