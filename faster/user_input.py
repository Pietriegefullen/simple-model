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
    parser.add_argument('--single', action = 'store_true')
    parser.add_argument('--t', nargs = 2)
    parser.add_argument('--p', nargs = '+')
    parser.add_argument('--local', action = 'store_true')
    parser.add_argument('--best', default = None)
    
    args = parser.parse_args()
    
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
    
    setattr(args, 'override', {} if args.p is None else {k: parse(v) 
                              for k,v in list(zip(args.p[::2], args.p[1::2]))})
    
    # init_config 
    # empty uses default
    # normalized: True/False
    
    
    
    return args
