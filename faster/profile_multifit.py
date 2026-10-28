#!/usr/bin/env python3
"""Profile one representative nested multi-replica fit.

This wrapper intentionally runs local optimisers in the current process by
default.  ``cProfile`` cannot collect call statistics from a
``ProcessPoolExecutor`` worker, so profiling a normal parallel run mostly
reports waiting in the parent process rather than the actual model work.

Example (use a deliberately small evaluation budget first)::

    python faster/profile_multifit.py --output profiling/quick \
        -- --samples 1367 --replicas 4 5 --outer-max-evaluations 2 \
        --local-max-evaluations 20 --yes

The arguments after ``--`` are passed to ``fit_replicas.py``.  Results are
written as normal checkpoints; use an isolated result directory/configuration
or a small budget for a profiling run.
"""

from __future__ import annotations

import argparse
import cProfile
import os
import pstats
import sys
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Profile the nested multi-replica fitting path.')
    parser.add_argument('--output', type=Path, default=Path('profiling/multi_fit'),
                        help='Output prefix (creates .pstats and .txt files).')
    parser.add_argument('--top', type=int, default=40,
                        help='Number of individual functions to include (default: 40).')
    parser.add_argument('--keep-local-workers', action='store_true',
                        help=('Keep --local-workers from the fit command. Child-process '
                              'work will not appear in the profile.'))
    parser.add_argument('--profile-integration', action='store_true',
                        help=('Also profile successful LSODA watchdog child processes and '
                              'write an aggregated integration report.'))
    parser.add_argument('fit_args', nargs=argparse.REMAINDER,
                        help='Arguments for fit_replicas.py; put them after --.')
    args = parser.parse_args(argv)
    if args.fit_args[:1] == ['--']:
        args.fit_args = args.fit_args[1:]
    if not args.fit_args:
        parser.error('Supply fit_replicas arguments after --.')
    if args.top <= 0:
        parser.error('--top must be positive.')
    return args


def _function_label(key):
    filename, line, name = key
    try:
        filename = str(Path(filename).resolve().relative_to(PROJECT_ROOT.parent))
    except ValueError:
        filename = Path(filename).name
    return f'{filename}:{line}({name})'


def _category(filename):
    """Give a useful high-level view without relying on package internals."""
    # cProfile uses ``~`` as the pseudo-filename for C built-ins. In this
    # project it commonly represents ``select.poll`` while the parent waits
    # for the watchdog integration process, not time spent in project code.
    if filename == '~' or filename.startswith('~'):
        return 'native/system calls'
    path = Path(filename)
    try:
        relative = path.resolve().relative_to(PROJECT_ROOT)
    except ValueError:
        text = str(path)
        if 'scipy' in text:
            return 'SciPy'
        if 'numpy' in text:
            return 'NumPy'
        return 'Python/other libraries'
    return f'project: {relative.parts[0]}' if len(relative.parts) > 1 else f'project: {relative.name}'


def write_report(stats_path, text_path, top):
    stats = pstats.Stats(str(stats_path))
    total_time = stats.total_tt
    categories = defaultdict(lambda: [0.0, 0])
    for (filename, _line, _name), (primitive, total, _tottime, cumtime, _callers) in stats.stats.items():
        # Cumtime is intentionally not summed here: parent and child call
        # times overlap.  Self time is additive and therefore meaningful by
        # category.
        bucket = categories[_category(filename)]
        bucket[0] += _tottime
        bucket[1] += primitive

    with text_path.open('w', encoding='utf-8') as report:
        report.write('Multi-fit profile\n')
        report.write('=================\n\n')
        report.write(f'Profiled elapsed time: {total_time:.3f} s\n')
        report.write('This run used one local worker unless --keep-local-workers was set.\n\n')
        report.write('Self time by code area (additive; best guide to where CPU time goes)\n')
        report.write('-------------------------------------------------------------------\n')
        for name, (self_time, calls) in sorted(categories.items(),
                                               key=lambda item: item[1][0], reverse=True):
            fraction = 0.0 if total_time == 0 else 100 * self_time / total_time
            report.write(f'{self_time:10.3f} s  {fraction:5.1f}%  {calls:10d} calls  {name}\n')

        report.write('\nTop functions by inclusive time\n')
        report.write('-------------------------------\n')
        # ``Stats.print_stats`` provides established, readable cProfile
        # formatting and includes both call count and self/cumulative time.
        stats.sort_stats(pstats.SortKey.CUMULATIVE)
        original_stream = stats.stream
        stats.stream = report
        stats.print_stats(top)
        stats.stream = original_stream

        rows = sorted(stats.stats.items(), key=lambda item: item[1][3], reverse=True)[:top]
        report.write('\nCompact ranking\n')
        report.write('---------------\n')
        for key, (primitive, total, self_time, cumulative, _callers) in rows:
            report.write(f'{cumulative:10.3f} s cumulative  {self_time:10.3f} s self  '
                         f'{primitive:8d}/{total:8d} calls  {_function_label(key)}\n')


def write_merged_integration_report(profile_directory, output_prefix, top):
    """Merge cProfile data emitted by successful integration child processes."""
    profiles = sorted(profile_directory.glob('integration-*.pstats'))
    if not profiles:
        print('No integration child profile was written (all runs may have timed out).')
        return
    merged_path = output_prefix.with_name(output_prefix.name + '-integration.pstats')
    pstats.Stats(*(str(path) for path in profiles)).dump_stats(str(merged_path))
    text_path = output_prefix.with_name(output_prefix.name + '-integration.txt')
    write_report(merged_path, text_path, top)
    print(f'Integration profiles merged from {len(profiles)} successful worker(s): {text_path}')
    print(f'Merged integration data: {merged_path}')


def main(argv=None):
    args = parse_args(argv)
    fit_args = list(args.fit_args)
    if not args.keep_local_workers:
        if '--local-workers' in fit_args:
            index = fit_args.index('--local-workers')
            if index + 1 >= len(fit_args):
                raise SystemExit('--local-workers needs a value.')
            fit_args[index + 1] = '1'
        else:
            fit_args.extend(['--local-workers', '1'])

    # The fit scripts expect their own directory on sys.path when invoked as
    # scripts. Make importing them from this wrapper behave identically.
    script_directory = str(PROJECT_ROOT)
    if script_directory not in sys.path:
        sys.path.insert(0, script_directory)
    import fit_replicas

    output_prefix = args.output
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    stats_path = output_prefix.with_suffix('.pstats')
    text_path = output_prefix.with_suffix('.txt')
    integration_directory = output_prefix.with_name(output_prefix.name + '-integration')
    previous_integration_directory = None
    if args.profile_integration:
        integration_directory.mkdir(parents=True, exist_ok=True)
        previous_integration_directory = os.environ.get(
            'SOIL_INTEGRATION_PROFILE_DIRECTORY')
        os.environ['SOIL_INTEGRATION_PROFILE_DIRECTORY'] = str(integration_directory)

    print('Profiling multi-fit with arguments:', ' '.join(fit_args))
    if not args.keep_local_workers:
        print('Forcing --local-workers 1 so the profile includes local optimisation work.')
    else:
        print('Keeping local worker processes: their CPU work is not included in cProfile.')

    profiler = cProfile.Profile()
    try:
        profiler.enable()
        fit_replicas.run(fit_replicas.parse_args(fit_args))
    finally:
        profiler.disable()
        profiler.dump_stats(str(stats_path))
        write_report(stats_path, text_path, args.top)
        print(f'Profile report: {text_path}')
        print(f'Raw profile data: {stats_path}')
        if args.profile_integration:
            if previous_integration_directory is None:
                os.environ.pop('SOIL_INTEGRATION_PROFILE_DIRECTORY', None)
            else:
                os.environ['SOIL_INTEGRATION_PROFILE_DIRECTORY'] = previous_integration_directory
            write_merged_integration_report(integration_directory, output_prefix, args.top)


if __name__ == '__main__':
    main()
