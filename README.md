simple model

## Resuming nested CMA-ES fits

`faster/fit_replicas.py` checkpoints both the outer and per-replica CMA-ES
search state.  Pressing Ctrl-C leaves the latest state in the normal joint
checkpoint directory; rerun the identical command (the default `--resume`)
to continue with the saved covariance, evolution paths, random state, and any
already completed candidates in an interrupted population.  Use
`--no-resume` to deliberately begin a fresh search.

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
