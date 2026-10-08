simple model

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
