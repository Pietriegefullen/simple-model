"""Rebuild checkpoint IDs after configuration type canonicalisation.

Run without ``--execute`` first.  The dry run loads every checkpoint, validates
that it can be canonicalised, and reports how many files will move.  An execute
run copies canonical checkpoint JSON to the current results tree and moves every
original to a timestamped sibling archive; no original is deleted.
"""

import argparse
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import shutil
import tempfile

import hashing
from USER_VARIABLES import RESULTS_DIRECTORY


CHECKPOINT_KEYS = {'parameters', 'total_loss', 'run_config'}


@dataclass
class Migration:
    source: Path
    destination: Path
    archived: Path
    contents: dict


def _loss_label(loss):
    """Use exactly the filename representation used by CheckpointCallback."""
    label = f'{loss:.6f}'.replace('.', '')
    return '9' * 8 if len(label) > 8 else label.zfill(8)


def checkpoint_destination(source_root, checkpoint):
    """Return the canonical result-file path for a loaded checkpoint."""
    run_config = checkpoint['run_config']
    chosen = run_config['chosen']
    directory = '_'.join((
        hashing.build_loss_id(run_config['objective']),
        hashing.build_model_id(run_config['model']),
        str(chosen['sample']),
        str(chosen['validation_replica']),
        str(chosen['fit_mode']),
    ))
    filename = '_'.join((
        directory,
        f"loss-{_loss_label(checkpoint['total_loss'])}",
        hashing.build_run_id(run_config),
        hashing.build_checkpoint_id(checkpoint),
    ))
    return source_root / directory / filename


def canonical_checkpoint(checkpoint):
    """Cast known checkpoint configuration fields without dropping unknown ones."""
    return hashing.cast_config(checkpoint, hashing.CHECKPOINT_CONFIG_SCHEMA)


def _same_checkpoint(left, right):
    return hashing._canonical_json(hashing.freeze(left)) == hashing._canonical_json(hashing.freeze(right))


def _checkpoint_files(source_root):
    for path in source_root.rglob('*'):
        if path.is_file() and not path.name.startswith('.'):
            yield path


def plan_migration(source_root, archive_root):
    """Read and validate every checkpoint before any filesystem mutation."""
    migrations = []
    skipped = 0
    errors = []
    destinations = {}

    for source in _checkpoint_files(source_root):
        try:
            with source.open() as handle:
                loaded = json.load(handle)
        except (OSError, json.JSONDecodeError) as error:
            # Results folders can also contain logs and non-checkpoint files.
            skipped += 1
            continue

        if not isinstance(loaded, dict) or not CHECKPOINT_KEYS.issubset(loaded):
            skipped += 1
            continue

        try:
            checkpoint = canonical_checkpoint(loaded)
            destination = checkpoint_destination(source_root, checkpoint)
        except (KeyError, TypeError, ValueError) as error:
            errors.append(f'{source}: {error}')
            continue

        migration = Migration(
            source=source,
            destination=destination,
            archived=archive_root / source.relative_to(source_root),
            contents=checkpoint,
        )
        prior = destinations.get(destination)
        if prior is not None and not _same_checkpoint(prior.contents, checkpoint):
            errors.append(f'{source}: canonical ID collides with different checkpoint {prior.source}')
        else:
            destinations[destination] = migration
            migrations.append(migration)

    for migration in migrations:
        if migration.archived.exists():
            errors.append(f'{migration.source}: archive target already exists: {migration.archived}')
        existing = migration.destination
        if existing.exists() and existing not in {item.source for item in migrations}:
            try:
                with existing.open() as handle:
                    current = canonical_checkpoint(json.load(handle))
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
                errors.append(f'{existing}: cannot validate existing destination: {error}')
            else:
                if not _same_checkpoint(current, migration.contents):
                    errors.append(f'{migration.source}: destination already exists with different contents: {existing}')

    return migrations, skipped, errors


def _write_checkpoint(destination, contents):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=destination.parent,
                                     prefix='.reindex-', delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(contents, handle, indent=4, allow_nan=True)
        handle.write('\n')
    temporary.replace(destination)


def remove_empty_directories(source_root):
    """Remove empty source subdirectories, deepest first, but keep the root."""
    removed = 0
    directories = sorted(
        (path for path in source_root.rglob('*') if path.is_dir()),
        key=lambda path: len(path.parts), reverse=True)
    for directory in directories:
        try:
            directory.rmdir()
        except OSError:
            # A directory containing a checkpoint, log, or active lock remains.
            continue
        removed += 1
    return removed


def execute_migration(migrations, source_root=None):
    """Archive sources, write canonical replacements, then remove empty folders."""
    for migration in migrations:
        migration.archived.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(migration.source), str(migration.archived))

    written = set()
    for migration in migrations:
        if migration.destination in written or migration.destination.exists():
            continue
        _write_checkpoint(migration.destination, migration.contents)
        written.add(migration.destination)

    return 0 if source_root is None else remove_empty_directories(source_root)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(RESULTS_DIRECTORY),
                        help='results directory to reindex (default: %(default)s)')
    parser.add_argument('--archive-dir', type=Path,
                        help='where to retain originals (default: timestamped sibling)')
    parser.add_argument('--execute', action='store_true',
                        help='perform the migration; without this flag, only report the plan')
    parser.add_argument('--verbose', action='store_true',
                        help='print every source-to-destination mapping')
    args = parser.parse_args(argv)

    source_root = args.source.resolve()
    timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    archive_root = (args.archive_dir or source_root.with_name(
        f'{source_root.name}-pre-reindex-{timestamp}')).resolve()
    if not source_root.is_dir():
        parser.error(f'checkpoint source does not exist: {source_root}')
    if archive_root == source_root or source_root in archive_root.parents:
        parser.error('archive directory must be outside the checkpoint source')

    migrations, skipped, errors = plan_migration(source_root, archive_root)
    print(f'Found {len(migrations)} checkpoint(s); skipped {skipped} non-checkpoint file(s).')
    if errors:
        print(f'Validation failed for {len(errors)} file(s); no files were moved:')
        print('\n'.join(errors))
        return 1

    if args.verbose:
        for migration in migrations:
            print(f'{migration.source} -> {migration.destination}')
    if not args.execute:
        print(f'Dry run only. Re-run with --execute; originals will be retained in {archive_root}.')
        return 0

    removed_directories = execute_migration(migrations, source_root)
    print(f'Reindexed {len(migrations)} checkpoint(s); removed '
          f'{removed_directories} empty folder(s). Originals are in {archive_root}.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
