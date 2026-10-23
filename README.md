simple model

## Composable checkpoint plots

`faster/plot.py` separates figure layout, model data, and the routines that
draw curves.  Build a `PlotBook` to state which routines belong on each named
figure, then give it a saved checkpoint.  The checkpoint is replayed only
once, and every figure receives the same model trajectory.

```python
import plot

book = plot.PlotBook({
    'gases': plot.FigurePlan(
        plot.FigureSpec(ncols=2, axis_names=('CO2', 'CH4'), figsize=(9, 4)),
        calls=(
            plot.PlotCall(plot.measured_pool, 'CO2', {'pool': 'CO2'}),
            plot.PlotCall(plot.modelled_pool, 'CO2', {'pool': 'CO2'}),
            # A CH4 curve on CO2's right-hand scale:
            plot.PlotCall(plot.modelled_pool, 'CO2', {'pool': 'CH4'},
                          secondary=True),
            plot.PlotCall(plot.measured_pool, 'CH4', {'pool': 'CH4'}),
            plot.PlotCall(plot.modelled_pool, 'CH4', {'pool': 'CH4'}),
        ),
        suptitle='checkpoint comparison',
    ),
})

figures = plot.render_checkpoint('KRQY-LMRT', book)
figures['gases'][0].savefig('gases.svg')
```

A plotting routine has the simple signature
`routine(context, axis, **options)`.  Use `plot.log_quantity` for an arbitrary
logged series, or write a small routine with that signature for a specialised
curve.  Several calls can target the same axis; calls with `secondary=True`
share its named `<axis>_right` secondary axis.  Text is rendered through LaTeX
globally; `plot.label('time', 'd')` produces text-mode `time` with only the
unit typeset in math mode.

## Plotting checkpoints and comparisons

`faster/plot_checkpoint.py` accepts full checkpoint paths or unique partial
checkpoint IDs.  With one checkpoint it writes the complete standard
report.  With two or more checkpoints, matching fit and
diagnostic plots are combined so each model result appears in the same figure.

```sh
python faster/plot_checkpoint.py KRQY --output-dir checkpoint_plots/KRQY
python faster/plot_checkpoint.py KRQY LMRT --output-dir checkpoint_plots/comparison
python faster/plot_checkpoint.py KRQY --compare D "*_thermodynamics" False
python faster/plot_checkpoint.py --best --sample 1351 --fit-mode split
```

`--best` selects the lowest-loss checkpoint for every matching model, fit mode,
sample, and validation replica. It writes one report per checkpoint. The
older `plot_best.py` is only a short compatibility shortcut for
`plot_checkpoint.py --best`; `run_checkpoint.py` is likewise retained only as
a compatibility entry point, so there is one plotting interface to learn.

Without `--output-dir`, reports are written below
`USER_VARIABLES.PLOT_DIRECTORY`: `checkpoint/<fit-mode>/<sample>-<replica>/<id>`
for one checkpoint and
`checkpoint-comparison/<fit-mode>/<sample>-<replica>/<id>__<id>` for a
comparison. The target directory is created when needed.

Existing figure files are skipped, so re-running `--best` only creates missing
figures. Add `--overwrite` to regenerate files already present in the target
directory.

Comparison titles contain checkpoint metadata shared by all inputs.  Legends
always identify the checkpoint and add model, replica, fit mode, or objective
only when that attribute differs between inputs.

`--compare` adds a named counterfactual replay for every requested checkpoint.
It takes a custom name followed by one or more parameter/value pairs; parameter
names accept shell-style wildcards and values are parsed as booleans, numbers,
or text.  For example, model D can reuse a model A checkpoint while disabling
all of its thermodynamic factors:

```sh
python faster/plot_checkpoint.py KRQY --compare D "*_thermodynamics" False
python faster/plot_best.py --sample 1351 --fit-mode split \
  --compare D "*_thermodynamics" False
```

The D replay is labelled `model D`, while the original checkpoint retains its
model label.  The wildcard expands only to parameters present in each source
checkpoint, and `*_thermodynamics=False` also updates the corresponding model
switches.  A named comparison uses the default directory
`checkpoint-comparison/<fit-mode>/<sample>-<replica>/<id>__<name>` below
`PLOT_DIRECTORY` unless
`--output-dir` is supplied to `plot_checkpoint`.  `plot_best` creates one such
comparison report per selected checkpoint; joint multi-replica fits are not
currently supported for named replays.

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

## Local sensitivity around one fitted checkpoint

`faster/sensitivity.py` performs a one-at-a-time, two-sided local sensitivity
analysis around the values in a saved checkpoint. It estimates each centred
finite-difference derivative, `dQ/dp ≈ (Q(p+h)-Q(p-h))/(2h)`, using a default
step of `h = 1e-4 × |p|`, without altering the checkpoint. The influence plot
shows this derivative by default. The CSV also reports the dimensionless
elasticity `d ln(Q)/d ln(p)`; use `--plot elasticity` when comparing effects
across parameters with unlike units. By default it measures the final modelled
CO2; any logged model quantity—or an expression such as `CO2 + CH4` or
`CO2 / CH4`—can be used, reduced to its final value, maximum, or time
integral.

```sh
python faster/sensitivity.py KRQY --quantity CO2 --stat final
python faster/sensitivity.py KRQY --quantity CH4 --stat integral --fraction 0.05
python faster/sensitivity.py KRQY --target 'CO2 + CH4' --stat final
python faster/sensitivity.py KRQY --target 'CO2 / CH4' --stat final
python faster/sensitivity.py KRQY --parameters '*_v_max' 'M_*'
python faster/sensitivity.py KRQY --plot elasticity
```

For a nested joint checkpoint, select the replica whose local solution should
be combined with the checkpoint's shared solution. `--scope shared` and
`--scope local` then test only the corresponding class of parameters.

```sh
python faster/sensitivity.py JOINT --replica 1351 5 --scope shared
```

Results are written by default under
`results/sensitivity/<sample>-<replica>/cp-<checkpoint>_<quantity>_<stat>/`.
Each run contains `metadata.json` (source checkpoint, replica, settings, and
parameter roles), `sensitivity.csv` (the baseline and both perturbation
outcomes), and `parameter_influence.png` (the per-parameter response plot).
Derivative plots use a signed logarithmic (`symlog`) horizontal axis, which
keeps positive and negative central differences while making small and large
magnitudes visible together. Recreate a plot without rerunning any model
simulations with:

```sh
python faster/replot_sensitivity.py results/sensitivity/1352-5/cp-ABCD_CO2_final
python faster/replot_sensitivity.py results/sensitivity/1352-5/cp-ABCD_CO2_final \
  --plot elasticity --output parameter_elasticity.png
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
