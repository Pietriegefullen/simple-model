import argparse
import fnmatch
from pathlib import Path

import numpy as np

from checkpoint_io import load_checkpoint_document, resolve_checkpoint
from converted_checkpoint import metadata_for_document, replay_spec
from USER_VARIABLES import CONVERTED_CHECKPOINTS_DIRECTORY

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
        '--converted-envelope', action='store_true',
        help='start a fresh fit from the converted-envelope model parameters')
    parameter_source.add_argument(
        '--hypothesis', action='store_true',
        help='start a fresh fit from the hypothesis model parameters')
    parameter_source.add_argument(
        '--estimate', action='store_true',
        help='start a fresh fit from the estimated model parameters')
    parameter_source.add_argument(
        '--bounds', metavar='ID',
        help='use parameter bounds stored by parameter_envelope.py (full or partial ID)')
    parser.add_argument(
        '--fix-v-max-except', nargs='+', metavar='PARAMETER',
        help='fix every *v_max parameter except the named parameter(s)')
    parser.add_argument(
        '--fix-v-max-value', type=float, metavar='VALUE',
        help='value to assign to each *v_max parameter fixed by --fix-v-max-except')
    parser.add_argument('--omit', nargs = '+', default = [])
    fit_mode = parser.add_mutually_exclusive_group()
    fit_mode.add_argument('--single', action='store_true')
    fit_mode.add_argument('--split', action='store_true')
    parser.add_argument('--t', nargs = 2)
    parser.add_argument('--p', nargs='+', action='append', metavar='TOKEN',
                        help='Fix NAME VALUE pairs, or use PATTERN [PATTERN ...] CHECKPOINT.')
    parser.add_argument('--local', action = 'store_true')
    parser.add_argument('--algorithm', choices=('de', 'cma-es'), default='de',
                        help='Global optimiser; DE is the default.')
    parser.add_argument('--popsize', type=int,
                        help='Population size for DE or CMA-ES.')
    parser.add_argument('--cma-sigma', type=float,
                        help='Initial CMA-ES step size in unit-bound coordinates (default: 0.30).')
    parser.add_argument('--funtol', type=float, default=1e-4,
                        help='CMA-ES function-value convergence tolerance (default: 1e-4).')
    parser.add_argument('--cma-restarts', type=int, default=0,
                        help='Number of CMA-ES restarts (default: 0).')
    parser.add_argument('--cma-restart-strategy', choices=('ipop', 'bipop'), default='ipop',
                        help='CMA-ES restart schedule: IPOP growth or BIPOP mixed populations (default: ipop).')
    parser.add_argument('--max-evaluations', type=int,
                        help='Maximum completed DE/CMA-ES evaluations.')
    parser.add_argument('--best', default = None)
    parser.add_argument('--workers', default = -1)
    parser.add_argument('--dry', action = 'store_true')
    parser.add_argument('--gas-limits', action='store_true',
                        help=('stop integrations when CO2 or CH4 exceeds a '
                              'replica-specific runaway limit; disabled by default'))
    parser.add_argument('--overwrite-checkpoints', '--overwrite', action = 'store_true',
                        help = 'replace checkpoints already stored for this run')
    parser.add_argument('--checkpoint', '--initial', dest='checkpoint', metavar='PATH',
                        help='checkpoint JSON to use as the initial parameters')
    
    args = parser.parse_args(argv)

    if not (args.single or args.split):
        parser.error('Specify --single, or use --split for a split fit.')
    
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
   
    setattr(args, 'workers', int(args.workers))

    if args.local and args.algorithm != 'de':
        parser.error('--local (Powell) cannot be combined with --algorithm.')
    if args.popsize is not None and args.popsize <= 0:
        parser.error('--popsize must be positive.')
    if args.cma_sigma is not None and args.cma_sigma <= 0:
        parser.error('--cma-sigma must be positive.')
    if not np.isfinite(args.funtol) or args.funtol <= 0:
        parser.error('--funtol must be finite and positive.')
    if args.cma_restarts < 0:
        parser.error('--cma-restarts cannot be negative.')
    if args.max_evaluations is not None and args.max_evaluations <= 0:
        parser.error('--max-evaluations must be positive.')
    if args.fix_v_max_value is not None and args.fix_v_max_except is None:
        parser.error('--fix-v-max-value requires --fix-v-max-except.')
    if args.bounds is not None and args.checkpoint is not None:
        parser.error('--bounds cannot be combined with --checkpoint.')


    return args

def fixed_literal(token):
    """Parse a numeric or boolean literal, leaving checkpoint selectors as text."""
    if token.casefold() in ('true', 'false'):
        return token.casefold() == 'true'
    try:
        return int(token)
    except ValueError:
        try:
            return float(token)
        except ValueError:
            return token


def checkpoint_override_selector(tokens):
    last = fixed_literal(tokens[-1])
    return last if isinstance(last, str) else None


def checkpoint_parameter_values(selector, model_version):
    converted_root = Path(CONVERTED_CHECKPOINTS_DIRECTORY) / f'v{model_version}'
    source = resolve_checkpoint(selector, additional_roots=(converted_root,))
    document = load_checkpoint_document(source)
    metadata = metadata_for_document(document)
    if metadata is None:
        version = str(document['run_config']['model']['version'])
    else:
        spec = replay_spec(metadata)
        if spec is None or spec[2]:
            raise ValueError(f'{source} does not declare physical checkpoint parameter values.')
        version = spec[0]
    if model_version is not None and version != model_version:
        message = f'{source} supplies v{version} values, but the fit uses v{model_version}.'
        if version == '0.3':
            message += ' Run convert_checkpoints_to_liter.py and use the converted checkpoint ID.'
        raise ValueError(message)
    return source, document['parameters']


def parse_override(inp, parameter_names=None, *, model_version=None, exclude_globs=()):
    """Resolve one ``--p`` group to physical values for the target model."""
    if not inp:
        return {}
    selector = checkpoint_override_selector(inp)
    from_checkpoint = selector is not None
    if len(inp) < 2 or (not from_checkpoint and len(inp) % 2):
        raise ValueError('--p requires NAME VALUE pairs or one or more patterns followed by a checkpoint.')
    pairs = ((pattern, None) for pattern in inp[:-1]) if from_checkpoint else (
        (pattern, fixed_literal(value)) for pattern, value in zip(inp[::2], inp[1::2]))
    names = None if parameter_names is None else tuple(parameter_names)
    parsed = {}
    for pattern, value in pairs:
        matches = [pattern]
        if names is not None:
            matches = [name for name in names if fnmatch.fnmatchcase(name, pattern)]
            if not matches:
                raise ValueError(f'--p pattern {pattern!r} matches no model parameter.')
            if any(character in pattern for character in '*?['):
                matches = [name for name in matches if name not in exclude_globs]
        if isinstance(value, str):
            raise ValueError(f'--p expected a numeric or boolean value for {pattern!r}.')
        parsed.update(dict.fromkeys(matches, value))

    if from_checkpoint and parsed:
        source, values = checkpoint_parameter_values(selector, model_version)
        missing = sorted(name for name in parsed if name not in values)
        if missing:
            raise ValueError(f'Checkpoint {source} has no saved value for: ' + ', '.join(missing))
        parsed = {name: values[name] for name in parsed}
        invalid = sorted(name for name, value in parsed.items()
                         if not isinstance(value, (int, float)) or not np.isfinite(value))
        if invalid:
            raise ValueError(f'Checkpoint {source} has non-finite or non-numeric values for: '
                             + ', '.join(invalid))
    return parsed
