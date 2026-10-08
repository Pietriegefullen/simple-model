import argparse

def numeric(value):
    numbers = '0123456789'
    if value.strip().lower() == 'inf':
        return True
    decimal = value.count('.') <= 1
    return decimal and all([v in numbers for v in str(value).replace('.', '')])
    
def parse(value):
    s = value.strip().lower()
    if s == 'true':
        return True
    if s == 'false':
        return False
    if s == 'none':
        return None
    
    if numeric(value):
        if '.' in s or value.strip().lower() == 'inf':
            return float(value.lower())
        return int(value)
    
    return value

def parse_args_fit():
    parser = argparse.ArgumentParser(
                        prog='fit_sample',
                        description='Fits a model to replica data.',
                        epilog='')
    parser.add_argument('sample', type = int)
    parser.add_argument('validation_replica', type = int)
    parser.add_argument('--default', action = 'store_true')
    parser.add_argument('--omit', nargs = '+', default = [])
    parser.add_argument('--single', action = 'store_true', default = None)
    parser.add_argument('--split', action = 'store_true', default = None)
    parser.add_argument('--t', nargs = 2)
    parser.add_argument('--p', nargs = '+')
    parser.add_argument('--local', action = 'store_true')
    parser.add_argument('--best', default = None)
    parser.add_argument('--dry', action = 'store_true')
    parser.add_argument('--overwrite-checkpoints', '--overwrite', action = 'store_true',
                        help = 'replace checkpoints already stored for this run')
    parser.add_argument('--checkpoint', '--initial', dest='checkpoint', metavar='PATH',
                        help='checkpoint JSON to use as the initial parameters')
    
    args = parser.parse_args()

    if args.single is None and args.split is None:
        args.split = True

    if not args.single is None and not args.split is None:
        raise Exception('Specify either "split" or "single", not both.')

    args.single = not args.split if not args.split is None else args.single
    
    t_start, t_end = None, None
    if not args.t is None:
        t_start, t_end = [parse(t_) for t_ in args.t]
        if t_start == 0:
            t_start = None
        elif t_start < 0:
            raise ValueError()
        if not t_start is None and not t_end is None:
            assert t_end > t_start
    setattr(args, 't_start', t_start)
    setattr(args, 't_end', t_end)
   
    setattr(args, 'override', parse_override(args.p))

    return args

def parse_override(inp):
    import parameters
    if inp is None:
        return {}
    keys = inp[::2]
    values = inp[1::2]
    parsed = {}
    parameter_names = [p.name for p in parameters.default_model_parameters()]
    for k,v in zip(keys, values):
        if '*' in k:
            for p_name in parameter_names:
                if k.replace('*', '') in p_name:
                    parsed[p_name] = parse(v)
        else:
            parsed[k] = parse(v)
    return parsed
