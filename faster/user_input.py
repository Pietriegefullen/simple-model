import argparse
import fnmatch

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

def parse_args_fit(argv=None):
    parser = argparse.ArgumentParser(
                        prog='fit_sample',
                        description='Fits a model to replica data.',
                        epilog='')
    parser.add_argument('sample', type = int)
    parser.add_argument('validation_replica', type = int)
    parameter_source = parser.add_mutually_exclusive_group()
    parameter_source.add_argument(
        '--default', action='store_true',
        help='start a fresh fit from the default model parameters')
    parameter_source.add_argument(
        '--estimate', action='store_true',
        help='start a fresh fit from the estimated model parameters')
    parser.add_argument(
        '--fix-v-max-except', nargs='+', metavar='PARAMETER',
        help='fix every *v_max parameter except the named parameter(s)')
    parser.add_argument(
        '--fix-v-max-value', type=float, metavar='VALUE',
        help='value to assign to each *v_max parameter fixed by --fix-v-max-except')
    parser.add_argument('--omit', nargs = '+', default = [])
    parser.add_argument('--single', action = 'store_true', default = None)
    parser.add_argument('--split', action = 'store_true', default = None)
    parser.add_argument('--t', nargs = 2)
    parser.add_argument('--p', nargs = '+')
    parser.add_argument('--local', action = 'store_true')
    parser.add_argument('--algorithm', choices=('de', 'cma-es'), default='de',
                        help='Global optimiser; DE is the default.')
    parser.add_argument('--popsize', type=int,
                        help='Population size for DE or CMA-ES.')
    parser.add_argument('--cma-sigma', type=float,
                        help='Initial CMA-ES step size in unit-bound coordinates (default: 0.30).')
    parser.add_argument('--cma-restarts', type=int, default=0,
                        help='Number of CMA-ES restarts (default: 0).')
    parser.add_argument('--max-evaluations', type=int,
                        help='Maximum completed DE/CMA-ES evaluations.')
    parser.add_argument('--best', default = None)
    parser.add_argument('--workers', default = -1)
    parser.add_argument('--dry', action = 'store_true')
    parser.add_argument('--overwrite-checkpoints', '--overwrite', action = 'store_true',
                        help = 'replace checkpoints already stored for this run')
    parser.add_argument('--checkpoint', '--initial', dest='checkpoint', metavar='PATH',
                        help='checkpoint JSON to use as the initial parameters')
    
    args = parser.parse_args(argv)

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
    setattr(args, 'workers', int(args.workers))

    if args.local and args.algorithm != 'de':
        parser.error('--local (Powell) cannot be combined with --algorithm.')
    if args.popsize is not None and args.popsize <= 0:
        parser.error('--popsize must be positive.')
    if args.cma_sigma is not None and args.cma_sigma <= 0:
        parser.error('--cma-sigma must be positive.')
    if args.cma_restarts < 0:
        parser.error('--cma-restarts cannot be negative.')
    if args.max_evaluations is not None and args.max_evaluations <= 0:
        parser.error('--max-evaluations must be positive.')
    if args.fix_v_max_value is not None and args.fix_v_max_except is None:
        parser.error('--fix-v-max-value requires --fix-v-max-except.')


    return args

def parse_override(inp, parameter_names=None):
    """Parse ``--p`` pairs, expanding globs against an explicit catalogue.

    Parameter defaults are model-version-specific.  Callers which know the
    model therefore pass its parameter names here; argument parsing itself
    preserves globs until that catalogue is available.
    """
    if inp is None:
        return {}
    if len(inp) % 2:
        raise ValueError('--p requires PARAMETER VALUE pairs.')
    keys = inp[::2]
    values = inp[1::2]
    parsed = {}
    names = None if parameter_names is None else tuple(parameter_names)
    for key, value in zip(keys, values):
        value = parse(value)
        if '*' in key and names is not None:
            matches = [name for name in names if fnmatch.fnmatchcase(name, key)]
            if not matches:
                raise ValueError(f'--p pattern {key!r} matches no model parameter.')
            parsed.update({name: value for name in matches})
        else:
            parsed[key] = value
    return parsed
