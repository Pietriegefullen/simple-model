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
    # TODO: compatibility layer here
    return 'run-' + build_id(config, 4,1)

def build_model_id(config):
    # TODO: compatibility layer here
    return 'model-' + build_id(config, 3,1)

if __name__ == '__main__':
    d = {'b': 456, 'a': 123, 'c': set([1,2,4])}
    print(build_id(d))
