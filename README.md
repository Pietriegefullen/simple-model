simple model

## Replaying a checkpoint

`faster/run_checkpoint.py` evaluates a saved annotated checkpoint once.  It
uses the model structure, fitting split, objective, and parameter ranges stored
in the checkpoint; it does not launch an optimiser or modify the checkpoint.

```sh
python faster/run_checkpoint.py KRQY-LMRT --output-dir replay_results/KRQY-LMRT
python faster/run_checkpoint.py KRQY --output-dir replay_results/KRQY-override \
  --p Hydro_v_max 0.02 Homo_thermodynamics false
python faster/run_checkpoint.py KRQY --output-dir replay_results/KRQY-window --t 5 40
```

`--p` values replace (and fix) the corresponding saved parameter values for
that one replay.  `--t` replaces the checkpoint's observation window; use `0`
as its first value to begin from the start of the experiment.
The checkpoint argument accepts a file path or a full/partial `cp-` identifier;
an ID must identify exactly one file below `results/`.
The target directory receives the standard fit and pathway diagnostic plots
from `collect_best`, along with `metadata.json` containing the source
checkpoint parameters, effective run configuration, requested overrides, and
result metrics.  It must be new or empty.
Relative `--output-dir` paths are rooted at this project directory (not
`faster/`), so the examples above write below `simple_model/replay_results/`.

## Resuming nested CMA-ES fits

`faster/fit_replicas.py` checkpoints both the outer and per-replica CMA-ES
search state.  Pressing Ctrl-C leaves the latest state in the normal joint
checkpoint directory; rerun the identical command (the default `--resume`)
to continue with the saved covariance, evolution paths, random state, and any
already completed candidates in an interrupted population.  Use
`--no-resume` to deliberately begin a fresh search.

## Polishing replicas from a multi fit

`faster/polish_replicas.py` takes an annotated multi-fit checkpoint (a path or
a unique full/partial `cp-` ID), freezes its saved shared parameter values, and
then refines each selected replica's local parameters with CMA-ES. It seeds
each search from the incumbent stored in the multi checkpoint and leaves the
source fit untouched. The default unit-bound `sigma0` is `0.10`, deliberately
smaller than the `0.30` exploratory default used by a fresh CMA-ES fit.

```sh
python faster/polish_replicas.py JI7G
python faster/polish_replicas.py JI7G --replicas 13514 13526 --max-evaluations 2000
python faster/polish_replicas.py JI7G --sigma 0.05 --popsize 16
```

Replica IDs are the full sample+replica IDs recorded by the multi fit (for
example `13514` or `1351/4` means sample `1351`, replica `4`); omitting
`--replicas` polishes every replica in that checkpoint. Results are stored under
`results/joint_polish/cp-<source-id>/`, alongside a `manifest.json` recording
the source checkpoint and CMA-ES settings.

## Finding alternative local minima

`faster/find_minima.py` performs independent, seeded differential-evolution
searches followed by Powell local refinement.  It keeps one final checkpoint
for every restart, rather than retaining only the globally best N candidates.

```sh
python faster/find_minima.py 1366 6 --restarts 40 --global-maxiter 300 \
  --run-name 1366-6-40-restarts
```

The run is stored under `minima_results/<run-name>/`, separately from
`results/`.  It writes the per-restart checkpoints, a candidate table, a basin
summary, a basin-discovery curve, parameter-position heatmap, and held-out
CO2/CH4 trajectory comparison.  Only locally converged candidates within the
configured relative/absolute loss cutoff are clustered; all candidates,
including failed or unconverged runs, remain in the per-restart archive.
The default serial run reports the start, duration, and loss of its first and
every 25th ODE evaluation; adjust this with `--progress-every`.

## Analysing existing checkpoints

`faster/analyze_checkpoint_basins.py` reads the already-retained annotated
checkpoints without changing `results/`, groups compatible fits, and identifies
low-loss parameter basins.  Its reports are written separately below
`checkpoint_basin_results/`.

```sh
python faster/analyze_checkpoint_basins.py --sample 1366 --replica 6 \
  --fit-mode split --run-name 1366-6-existing
```

## Parameter distributions for successful single fits

`faster/plot_single_parameter_distributions.py` retains the best (lowest-loss)
checkpoint for every single-replica fit, replays it, and plots each fitted
parameter grouped by generated model ID.  Only checkpoints whose combined fit
R² is strictly above the threshold are included; the default is `R² > 0.8`.
It writes a paginated PDF plus a CSV containing every plotted value and its
checkpoint metadata.

```sh
python faster/plot_single_parameter_distributions.py
python faster/plot_single_parameter_distributions.py --r2-threshold 0.9
```
