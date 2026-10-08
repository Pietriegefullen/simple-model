"""Compare fitted model checkpoints on a common, scale-aware parameter map.

The left panel maps each fitted value to its position in the configured search
range.  The right panel shows the signed change in that position relative to a
reference checkpoint.  This makes a value close to a bound, and a meaningful
change, directly comparable for parameters with different units and log/linear
scales.

Examples
--------
python faster/compare_checkpoints.py 7RV6-HVP2 KCB5 \\
    --output plots/checkpoint-comparison.svg

Checkpoint identifiers are looked up recursively in
``USER_VARIABLES.RESULTS_DIRECTORY``. A partial identifier is accepted when it
resolves to exactly one checkpoint.
"""

import argparse
import contextlib
import io
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch

import parameters
import hashing
import data
import optimizer
from fit_sample import run
from model_variants import model_variant_from_id
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY

LOSS_PATTERN = re.compile(r"loss[-_]\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?\d+)?)")
CHECKPOINT_ID_PATTERN = re.compile(r"(?:^|[_-])cp-([a-z0-9-]+)$", re.IGNORECASE)


def checkpoint_id(path):
    """Extract the short checkpoint identifier from a results filename."""
    match = CHECKPOINT_ID_PATTERN.search(Path(path).name)
    return match.group(1) if match else None


def checkpoint_index(results_directory):
    """Index result files that carry a ``cp-...`` checkpoint identifier."""
    root = Path(results_directory)
    if not root.is_dir():
        raise ValueError(f"Results directory does not exist: {root}")
    return [(checkpoint_id(path), path) for path in root.rglob("*")
            if path.is_file() and checkpoint_id(path) is not None]


def resolve_checkpoints(identifiers, results_directory):
    """Resolve direct paths or unique (possibly partial) checkpoint IDs."""
    indexed = checkpoint_index(results_directory)
    resolved = []
    for identifier in identifiers:
        direct_path = Path(identifier)
        if direct_path.is_file():
            resolved.append(direct_path)
            continue

        query = identifier.casefold().removeprefix("cp-")
        matches = [(checkpoint_id, path) for checkpoint_id, path in indexed
                   if query in checkpoint_id.casefold()]
        if not matches:
            raise ValueError(f"No checkpoint ID matching {identifier!r} in {results_directory}")
        if len(matches) > 1:
            suggestions = "\n".join(
                f"  {checkpoint_id}: {path.relative_to(results_directory)}"
                for checkpoint_id, path in matches[:12]
            )
            remaining = "" if len(matches) <= 12 else f"\n  ... and {len(matches) - 12} more"
            raise ValueError(
                f"Checkpoint ID {identifier!r} is ambiguous; use more characters:\n"
                f"{suggestions}{remaining}"
            )
        resolved.append(matches[0][1])
    return resolved


def read_checkpoint(path):
    """Return the parameter mapping from plain or annotated JSON checkpoints."""
    with Path(path).open() as handle:
        checkpoint = json.load(handle)
    return checkpoint.get("parameters", checkpoint)


def metadata():
    """Parameter scales and bounds declared by the model, keyed by name."""
    return {
        parameter.name: ((parameter.low, parameter.high)
                         if parameter.is_variable() else None, parameter.scale)
        for parameter in parameters.default_model_parameters()
    }


# This mirrors the pool and pathway ordering used by ``collect_best`` plots.
# Keeping a pathway's biomass, kinetic parameters, and CUE together makes it
# much easier to scan a model change than grouping all parameters by type.
PATHWAY_PARAMETER_ORDER = {
    "Initial pools": ("TOC", "DOC", "Acetate", "Fe3", "H2", "H2O", "Fe2"),
    "Hydrolysis": ("Hydrolysis_v_max", "Hydrolysis_Kmb"),
    "Fermentation": ("Ferm_v_max", "Ferm_Km", "Ferm_inhibition", "Ferm_CUE", "M_Ferm"),
    "Fe3": ("Fe3_Km_Ac", "Fe3_Km_Fe3", "Fe3_v_max", "Fe3_CUE", "M_Fe3"),
    "Aceto": ("Aceto_Km_Ac", "Ac_v_max", "Ac_CUE", "M_Ac"),
    "Hydro": ("Hydro_Km_CO2", "Hydro_v_max", "Hydro_CUE", "Hydro_Km_H2", "M_Hydro"),
    "Homo": ("Homo_Km_H2", "Homo_Km_CO2", "Homo_v_max", "Homo_CUE", "M_Homo"),
    "Shared": ("death_rate",),
}
GROUP_ORDER = list(PATHWAY_PARAMETER_ORDER) + ["Other"]
PARAMETER_TO_GROUP = {
    parameter: group
    for group, ordered_parameters in PATHWAY_PARAMETER_ORDER.items()
    for parameter in ordered_parameters
}
PARAMETER_ORDER = {
    parameter: index
    for index, parameter in enumerate(
        parameter for ordered_parameters in PATHWAY_PARAMETER_ORDER.values()
        for parameter in ordered_parameters
    )
}


def parameter_group(name):
    """Associate every known parameter with the pathway drawn in the plots."""
    return PARAMETER_TO_GROUP.get(name, "Other")

def range_position(value, bounds, scale):
    """Map a raw parameter value to [0, 1] using its declared fit scale."""
    if value is None or bounds is None:
        return np.nan
    lower, upper = bounds
    if scale == "log":
        if value <= 0 or lower <= 0 or upper <= 0:
            return np.nan
        value, lower, upper = map(math.log10, (value, lower, upper))
    if upper == lower:
        return np.nan
    return (value - lower) / (upper - lower)


def loss_from_checkpoint(path):
    """Use the loss stored in modern checkpoints, then fall back to filenames."""
    with Path(path).open() as handle:
        checkpoint = json.load(handle)
    if "total_loss" in checkpoint:
        return checkpoint["total_loss"]
    match = LOSS_PATTERN.search(Path(path).name)
    return float(match.group(1)) if match else None


def variant_from_checkpoint(path):
    """Return the mapped model variant, retaining an unknown model ID visibly."""
    with Path(path).open() as handle:
        checkpoint = json.load(handle)
    run_config = checkpoint.get("run_config")
    if not run_config or "model" not in run_config:
        return None
    model_id = hashing.build_model_id(run_config["model"])
    return model_variant_from_id(model_id) or f"unknown ({model_id})"


def r2_from_checkpoint(path):
    """Recreate a checkpoint run and return its fit and validation R² values."""
    with Path(path).open() as handle:
        checkpoint = json.load(handle)
    if "run_config" not in checkpoint or "parameters" not in checkpoint:
        return None, None

    try:
        # ``run`` emits optimisation progress even for a single replay; the
        # comparison command reports one concise warning if a replay fails.
        with contextlib.redirect_stdout(io.StringIO()):
            pathway_model, objective, run_log = run(
                checkpoint["run_config"],
                parameters.ModelParameters(checkpoint["parameters"]),
            )
            fit_r2 = objective.R2(run_log)
            chosen = checkpoint["run_config"]["chosen"]
            validation_r2 = None
            if chosen["fit_mode"] != "single":
                sample = data.get_data_before_carex()[chosen["sample"]]
                validation_replica = sample[chosen["validation_replica"]]
                validation_objective = optimizer.build_objective_function(
                    pathway_model, validation_replica,
                    checkpoint["run_config"]["objective"], t_start=0, t_end=None,
                )
                validation_r2 = validation_objective.R2(run_log)
        return fit_r2, validation_r2
    except Exception as error:
        print(f"Warning: could not calculate R² for {Path(path).name}: {error}")
        return None, None


def display_label(path, loss, variant):
    name = checkpoint_id(path) or Path(path).name
    lines = [name]
    if variant is not None:
        lines.append(f"model {variant}")
    if loss is not None:
        lines.append(f"loss {loss:g}")
    return "\n".join(lines)


def load_rows(paths, show_constant):
    fitted = [read_checkpoint(path) for path in paths]
    model_metadata = metadata()
    names = set().union(*[values.keys() for values in fitted])
    rows = []
    for name in names:
        bounds, scale = model_metadata.get(name, (None, "linear"))
        # Constants and unknown fields do not have a meaningful search position.
        if bounds is None:
            continue
        values = [checkpoint.get(name) for checkpoint in fitted]
        positions = [range_position(value, bounds, scale) for value in values]
        finite = np.asarray([value for value in positions if np.isfinite(value)])
        # Fixed parameters add visual noise.  They remain available on request.
        if (not show_constant and len(finite) == len(paths)
                and np.ptp(finite) < 1e-12):
            continue
        rows.append({
            "name": name,
            "group": parameter_group(name),
            "scale": scale,
            "bounds": bounds,
            "values": values,
            "positions": positions,
        })
    return rows


def order_rows(rows, reference_index, sort):
    for row in rows:
        positions = np.asarray(row["positions"], dtype=float)
        reference = positions[reference_index]
        differences = positions - reference
        row["max_change"] = np.nanmax(np.abs(differences)) if np.isfinite(differences).any() else -1
    if sort == "change":
        return sorted(rows, key=lambda row: (-row["max_change"], row["name"]))
    group_index = {name: index for index, name in enumerate(GROUP_ORDER)}
    return sorted(rows, key=lambda row: (
        group_index.get(row["group"], 99),
        PARAMETER_ORDER.get(row["name"], len(PARAMETER_ORDER)),
        row["name"],
    ))


def add_group_separators(axis, rows, row_offset=0):
    last_group = rows[0]["group"]
    for index, row in enumerate(rows[1:], start=1):
        if row["group"] != last_group:
            axis.axhline(row_offset + index - .5, color="black", linewidth=.7)
            last_group = row["group"]


def annotate_r2_rows(axis, values, signed=False):
    """Print exact R² values in the first rows of a checkpoint-aligned grid."""
    for row_index, row_values in enumerate(values):
        for column_index, value in enumerate(row_values):
            text = "—" if not np.isfinite(value) else f"{value:+.2f}" if signed else f"{value:.2f}"
            axis.text(column_index, row_index, text, ha="center", va="center",
                      fontsize=8, color="black",
                      bbox={"facecolor": "white", "edgecolor": "none", "alpha": .75, "pad": .5})


def plot_comparison(paths, output, reference_index, sort, show_constant, show_r2=True):
    rows = order_rows(load_rows(paths, show_constant), reference_index, sort)
    if not rows:
        raise ValueError("No varying, ranged parameters found. Try --include-constant.")

    position_data = np.asarray([row["positions"] for row in rows], dtype=float)
    reference = position_data[:, [reference_index]]
    difference_data = position_data - reference
    losses = [loss_from_checkpoint(path) for path in paths]
    variants = [variant_from_checkpoint(path) for path in paths]
    labels = [display_label(path, loss, variant)
              for path, loss, variant in zip(paths, losses, variants)]
    r2_values = [r2_from_checkpoint(path) for path in paths] if show_r2 else []
    has_r2 = any(value is not None and np.isfinite(value)
                 for values in r2_values for value in values)

    r2_row_count = 3 if has_r2 else 0  # fit, validation, then a visual divider
    r2_labels = ["R² fit", "R² validation", ""] if has_r2 else []
    if has_r2:
        fit_r2 = np.asarray([values[0] for values in r2_values], dtype=float)
        validation_r2 = np.asarray([values[1] for values in r2_values], dtype=float)
        r2_data = np.vstack((fit_r2, validation_r2, np.full(len(paths), np.nan)))
        r2_delta = r2_data - r2_data[:, [reference_index]]
        position_data = np.vstack((r2_data, position_data))
        difference_data = np.vstack((r2_delta, difference_data))

    height = max(5, .28 * (len(rows) + r2_row_count) + 2.3)
    width = max(10, 4.5 + 1.05 * len(paths))
    figure = plt.figure(figsize=(width, height))
    grid = figure.add_gridspec(1, 2, wspace=.08)
    position_axis = figure.add_subplot(grid[0, 0])
    delta_axis = figure.add_subplot(grid[0, 1], sharey=position_axis)

    position_map = plt.get_cmap("viridis").copy()
    delta_map = plt.get_cmap("coolwarm").copy()
    position_map.set_bad("#d9d9d9")
    delta_map.set_bad("#d9d9d9")
    position_image = position_axis.imshow(position_data, aspect="auto", cmap=position_map,
                                           vmin=0, vmax=1, interpolation="nearest")
    delta_image = delta_axis.imshow(
        difference_data, aspect="auto", cmap=delta_map,
        norm=TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1), interpolation="nearest",
    )

    position_title = ("R² and fitted position in declared search range"
                      if has_r2 else "Fitted position in declared search range")
    for axis, title in (
        (position_axis, position_title),
        (delta_axis, f"Change from {checkpoint_id(paths[reference_index]) or Path(paths[reference_index]).name}"),
    ):
        axis.set_title(title)
        axis.set_xticks(np.arange(len(paths)), labels, rotation=55, ha="right")
        axis.tick_params(axis="x", length=0)
        add_group_separators(axis, rows, row_offset=r2_row_count)
        if has_r2:
            axis.axhline(r2_row_count - .5, color="black", linewidth=1.1)

    position_axis.set_yticks(np.arange(len(rows) + r2_row_count),
                             r2_labels + [row["name"] for row in rows])
    position_axis.tick_params(axis="y", length=0, labelsize=8)
    delta_axis.tick_params(axis="y", left=False, labelleft=False)
    if has_r2:
        annotate_r2_rows(position_axis, r2_data[:2])
        annotate_r2_rows(delta_axis, r2_delta[:2], signed=True)
    figure.colorbar(position_image, ax=position_axis, fraction=.035, pad=.03,
                    label="Parameter rows: 0 = lower bound; 1 = upper bound")
    figure.colorbar(delta_image, ax=delta_axis, fraction=.035, pad=.03,
                    label="change in range position")
    figure.legend(
        handles=[Patch(facecolor="#d9d9d9", label="not present / no declared range")],
        loc="lower center", ncol=1, frameon=False,
    )
    figure.subplots_adjust(left=.28, bottom=.25, top=.92)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint_ids", nargs="+",
                        help="unique full or partial IDs from results/ (or direct checkpoint paths)")
    parser.add_argument("--results-dir", default=RESULTS_DIRECTORY,
                        help="checkpoint search directory (default: USER_VARIABLES.RESULTS_DIRECTORY)")
    parser.add_argument("--output", "-o", default=str(Path(PROJECT_DIRECTORY) / "comparison.png"),
                        help="output figure (SVG, PDF, or PNG)")
    parser.add_argument("--reference", type=int, default=0,
                        help="zero-based checkpoint index used as the delta reference (default: 0)")
    parser.add_argument("--sort", choices=("group", "change"), default="group",
                        help="row order: biological group or largest change")
    parser.add_argument("--include-constant", action="store_true",
                        help="also show parameters identical across the selected checkpoints")
    parser.add_argument("--no-r2", action="store_true",
                        help="skip checkpoint replays and omit the R² summary")
    args = parser.parse_args()
    if not 0 <= args.reference < len(args.checkpoint_ids):
        parser.error("--reference must index one of the supplied checkpoint IDs")
    try:
        paths = resolve_checkpoints(args.checkpoint_ids, args.results_dir)
    except ValueError as error:
        parser.error(str(error))
    rows = plot_comparison(paths, args.output, args.reference,
                           args.sort, args.include_constant, show_r2=not args.no_r2)
    print(f"Wrote {args.output} ({len(rows)} parameters, {len(paths)} checkpoints)")


if __name__ == "__main__":
    main()
