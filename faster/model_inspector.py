"""Interactive checkpoint and parameter inspector for the soil model.

Run from the project root with::

    python faster/model_inspector.py

The inspector is deliberately a local Tk application: it needs no web server,
and it calls the same checkpoint replay and report plotting code used by
``plot_checkpoint.py``.  Checkpoints are read-only; changing a control creates
a clearly labelled manual replay in memory.
"""

from __future__ import annotations

import copy
import math
import traceback
from dataclasses import dataclass

import matplotlib

# Select the interactive backend before importing pyplot or Tk canvas helpers.
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import messagebox, ttk

import data
import parameters
import checkpoint_report as report
import plot_checkpoint
from fit_sample import run
from model_factory import build_model


def checkpoint_label(checkpoint, hidden=(), details=True) -> str:
    """Show useful catalogue metadata, omitting dimensions already filtered."""
    hidden = set(hidden)
    fields = [] if "checkpoint_id" in hidden else [f"{checkpoint.id}"]
    fields.append(f"loss {checkpoint.loss:.7g}")
    if not details:
        return "  |  ".join(fields)
    sample, replica = checkpoint.replica
    details = (
        ("model_version", f"version {checkpoint.run_config['model'].get('version', '?')}"),
        ("sample", f"sample {sample}"),
        ("replica", f"replica {replica}"),
        ("fit_mode", checkpoint.fit_mode),
        ("model_id", checkpoint.model_id),
        ("run_id", checkpoint.run_id),
    )
    fields.extend(value for name, value in details if name not in hidden)
    return "  |  ".join(fields)


def is_log_slider(parameter) -> bool:
    """Only use logarithmic controls for strictly-positive logarithmic ranges."""
    return (parameter.scale == "log" and parameter.is_variable()
            and parameter.lower() > 0 and parameter.upper() > 0)


def slider_to_value(position: float, low: float, high: float, logarithmic: bool) -> float:
    """Map a 0..1000 Tk scale coordinate to a physical parameter value."""
    fraction = min(1.0, max(0.0, float(position) / 1000.0))
    if logarithmic:
        return math.exp(math.log(low) + fraction * (math.log(high) - math.log(low)))
    return low + fraction * (high - low)


def value_to_slider(value: float, low: float, high: float, logarithmic: bool) -> float:
    """Map a physical parameter value to a 0..1000 Tk scale coordinate."""
    if high <= low:
        return 0.0
    if logarithmic:
        value = min(high, max(low, value))
        return 1000.0 * (math.log(value) - math.log(low)) / (math.log(high) - math.log(low))
    return 1000.0 * (min(high, max(low, value)) - low) / (high - low)


@dataclass
class ParameterControl:
    parameter: parameters.Parameter
    frame: ttk.Frame
    slider: ttk.Scale
    value: tk.StringVar
    low: tk.StringVar
    high: tk.StringVar
    scale_label: ttk.Label


class FigurePane:
    """A Matplotlib canvas which can live in the main or a detached window."""

    def __init__(self, parent, figure, title=""):
        self.parent = parent
        self.figure = figure
        self.canvas = FigureCanvasTkAgg(figure, master=parent)
        self.toolbar = NavigationToolbar2Tk(self.canvas, parent, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(fill="x")
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        if title and hasattr(parent, "wm_title"):
            parent.wm_title(title)
        self.canvas.draw_idle()

    def destroy(self):
        self.toolbar.destroy()
        self.canvas.get_tk_widget().destroy()
        plt.close(self.figure)


class ModelInspector(ttk.Frame):
    """Checkpoint catalogue, parameter editor, and reusable report viewer."""

    def __init__(self, root):
        super().__init__(root, padding=8)
        self.root = root
        self.dataset = None
        self.checkpoints = []
        self.filtered = []
        self.active_checkpoint = None
        self.working_parameters = None
        self.controls: dict[str, ParameterControl] = {}
        self.plot_names: list[str] = ["fit"]
        self.current_plot = tk.StringVar(value="fit")
        self.log_y = tk.BooleanVar(value=False)
        # A newly loaded checkpoint is already one of the comparison curves.
        # Turn this on automatically as soon as a value is edited.
        self.include_manual = tk.BooleanVar(value=False)
        self._main_pane = None
        self._detached: list[tuple[tk.Toplevel, FigurePane]] = []
        self._result_cache = {}

        self._build()
        self.reload_catalogue()

    def _build(self):
        self.pack(fill="both", expand=True)
        self.root.title("Soil model inspector")
        self.root.geometry("1450x930")
        self.columnconfigure(2, weight=1)
        self.rowconfigure(0, weight=1)

        self._build_catalogue(self)
        self._build_parameters(self)
        self._build_viewer(self)

    def _build_catalogue(self, parent):
        panel = ttk.LabelFrame(parent, text="Checkpoints", padding=6)
        panel.grid(row=0, column=0, sticky="nsew", padx=(0, 6))
        panel.rowconfigure(2, weight=1)
        panel.columnconfigure(0, weight=1)

        self.filters = {name: tk.StringVar() for name in (
            "model_version", "sample", "replica", "fit_mode", "model_id", "checkpoint_id", "run_id")}
        filter_row = ttk.Frame(panel)
        filter_row.grid(row=0, column=0, sticky="ew")
        filter_row.columnconfigure(1, weight=1)
        for row, (label, variable) in enumerate((
                ("Model version", self.filters["model_version"]),
                ("Sample", self.filters["sample"]),
                ("Replica", self.filters["replica"]),
                ("Fit mode", self.filters["fit_mode"]),
                ("Model ID", self.filters["model_id"]),
                ("Checkpoint ID", self.filters["checkpoint_id"]),
                ("Run ID", self.filters["run_id"]))):
            ttk.Label(filter_row, text=label).grid(row=row, column=0, sticky="w", padx=(0, 4))
            entry = ttk.Entry(filter_row, textvariable=variable, width=20)
            entry.grid(row=row, column=1, sticky="ew", pady=1)
            entry.bind("<KeyRelease>", lambda _: self.apply_filter())
        options = ttk.Frame(panel)
        options.grid(row=1, column=0, sticky="ew", pady=5)
        self.best_per_run = tk.BooleanVar(value=False)
        ttk.Checkbutton(options, text="Best loss per run", variable=self.best_per_run,
                        command=self.apply_filter).pack(side="left")
        ttk.Button(options, text="Reload", command=self.reload_catalogue).pack(side="right")

        self.catalogue = tk.Listbox(panel, selectmode="extended", exportselection=False,
                                    width=21, height=28)
        self.catalogue.grid(row=2, column=0, sticky="nsew")
        self.catalogue.bind("<<ListboxSelect>>", self._checkpoint_selected)
        scroll = ttk.Scrollbar(panel, orient="vertical", command=self.catalogue.yview)
        scroll.grid(row=2, column=1, sticky="ns")
        self.catalogue.configure(yscrollcommand=scroll.set)
        self.catalogue_status = ttk.Label(panel, text="Loading checkpoint catalogue…")
        self.catalogue_status.grid(row=3, column=0, sticky="w", pady=(5, 0))

    def _build_parameters(self, parent):
        panel = ttk.LabelFrame(parent, text="Parameter values", padding=6)
        panel.grid(row=0, column=1, sticky="nsew", padx=(0, 6))
        panel.rowconfigure(1, weight=1)
        panel.columnconfigure(0, weight=1)
        toolbar = ttk.Frame(panel)
        toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 5))
        ttk.Button(toolbar, text="Restore checkpoint", command=self.restore_checkpoint).pack(side="left")
        ttk.Button(toolbar, text="Refresh plot", command=self.refresh_main_plot).pack(side="left", padx=4)
        ttk.Label(toolbar, text="(bounds are editable)").pack(side="left")

        canvas = tk.Canvas(panel, width=510, highlightthickness=0)
        canvas.grid(row=1, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(panel, orient="vertical", command=canvas.yview)
        scroll.grid(row=1, column=1, sticky="ns")
        canvas.configure(yscrollcommand=scroll.set)
        self.parameter_body = ttk.Frame(canvas)
        self.parameter_window = canvas.create_window((0, 0), window=self.parameter_body, anchor="nw")
        self.parameter_body.bind("<Configure>", lambda _: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(self.parameter_window, width=event.width))

    def _build_viewer(self, parent):
        panel = ttk.LabelFrame(parent, text="Plots", padding=6)
        panel.grid(row=0, column=2, sticky="nsew")
        panel.rowconfigure(1, weight=1)
        panel.columnconfigure(0, weight=1)
        toolbar = ttk.Frame(panel)
        toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 5))
        toolbar.columnconfigure(1, weight=1)
        ttk.Label(toolbar, text="Plot").grid(row=0, column=0, sticky="w")
        self.plot_picker = ttk.Combobox(toolbar, state="readonly", textvariable=self.current_plot,
                                        values=self.plot_names)
        self.plot_picker.grid(row=0, column=1, sticky="ew", padx=(5, 0))
        self.plot_picker.bind("<<ComboboxSelected>>", lambda _: self.refresh_main_plot())
        ttk.Checkbutton(toolbar, text="Log y", variable=self.log_y,
                        command=self.refresh_main_plot).grid(row=1, column=0, sticky="w", pady=(3, 0))
        ttk.Checkbutton(toolbar, text="Include manual values", variable=self.include_manual,
                        command=self.refresh_main_plot).grid(row=1, column=1, sticky="w", pady=(3, 0))
        ttk.Button(toolbar, text="Detach / compare", command=self.detach_plot).grid(
            row=2, column=0, columnspan=2, sticky="ew", pady=(3, 0))
        self.viewer = ttk.Frame(panel)
        self.viewer.grid(row=1, column=0, sticky="nsew")

    def reload_catalogue(self):
        try:
            self.checkpoints = report.collect_all()

        except Exception as error:
            traceback.print_exc()
            messagebox.showerror("Could not load checkpoints", str(error))
            self.checkpoints = []
        self._clear_caches()
        self.apply_filter()

    def apply_filter(self):
        filters = {name: variable.get().casefold().strip()
                   for name, variable in self.filters.items()}

        # The common initial view needs no metadata comparisons: the catalogue
        # has already been loaded and loss-sorted in ``reload_catalogue``.
        # Avoid calculating every checkpoint's run hash just to conclude that
        # it belongs in the unfiltered list.
        if not any(filters.values()) and not self.best_per_run.get():
            self.filtered = list(self.checkpoints)
            self.catalogue.delete(0, "end")
            for checkpoint in self.filtered:
                self.catalogue.insert("end", checkpoint_label(checkpoint, details=False))
            self.catalogue_status.configure(
                text=f"{len(self.filtered)} of {len(self.checkpoints)} checkpoint(s)")
            return

        def matches(checkpoint):
            sample, replica = checkpoint.replica
            values = {
                "model_version": checkpoint.run_config['model'].get('version', ''),
                "sample": sample,
                "replica": replica,
                "fit_mode": checkpoint.fit_mode,
                "model_id": checkpoint.model_id,
                "checkpoint_id": checkpoint.id,
                "run_id": checkpoint.run_id,
            }
            return all(not requested or requested in str(values[name]).casefold()
                       for name, requested in filters.items())

        candidates = self.checkpoints
        if self.best_per_run.get():
            best_by_run = {}
            for checkpoint in candidates:
                prior = best_by_run.get(checkpoint.run_id)
                if prior is None or checkpoint.loss < prior.loss:
                    best_by_run[checkpoint.run_id] = checkpoint
            candidates = list(best_by_run.values())
        self.filtered = sorted((checkpoint for checkpoint in candidates if matches(checkpoint)),
                               key=lambda checkpoint: checkpoint.loss)
        self.catalogue.delete(0, "end")
        hidden = {name for name, requested in filters.items() if requested}
        for checkpoint in self.filtered:
            self.catalogue.insert("end", checkpoint_label(checkpoint, hidden))
        self.catalogue_status.configure(text=f"{len(self.filtered)} of {len(self.checkpoints)} checkpoint(s)")

    def selected_checkpoints(self):
        return [self.filtered[index] for index in self.catalogue.curselection()]

    def _checkpoint_selected(self, _event=None):
        selected = self.selected_checkpoints()
        if not selected:
            return
        # The first selected item supplies editable values; all selected items are
        # retained as comparison curves when opening or refreshing a plot.
        checkpoint = selected[0]
        if checkpoint is not self.active_checkpoint:
            self.active_checkpoint = checkpoint
            self.load_parameters(checkpoint)

    def load_parameters(self, checkpoint):
        # Persisted checkpoint values intentionally omit their fitting ranges.
        # Reapply the model defaults so every numerical model parameter remains
        # editable, while preserving its physical/logarithmic scale.
        chosen = checkpoint.run_config['chosen']
        self.working_parameters = parameters.ModelParameters(build_model(
            checkpoint.run_config['model']['version'], chosen['pathways'],
            chosen.get('normalized_parameters', True)).parameters())
        for saved in checkpoint.parameters:
            if saved.name in self.working_parameters:
                self.working_parameters[saved.name].set(saved.value)
                continue
            if isinstance(saved.value, bool):
                continue
            value = float(saved.value)
            if saved.scale == "log" and value > 0:
                bounds = (value / 1_000, value * 1_000)
            else:
                width = max(abs(value), 1.0)
                bounds = (value - width, value + width)
            self.working_parameters[saved.name].variable(value, bounds)
        self.include_manual.set(False)
        self._clear_caches()
        self._build_parameter_controls()
        self.refresh_main_plot()

    def restore_checkpoint(self):
        if self.active_checkpoint is not None:
            self.load_parameters(self.active_checkpoint)

    def _build_parameter_controls(self):
        for child in self.parameter_body.winfo_children():
            child.destroy()
        self.controls.clear()
        if self.working_parameters is None:
            return
        for row, parameter in enumerate(sorted(self.working_parameters, key=lambda item: item.name)):
            if not parameter.is_variable() or parameter.options is not None:
                continue
            self._add_parameter_control(row, parameter)

    def _add_parameter_control(self, row, parameter):
        frame = ttk.Frame(self.parameter_body, padding=(0, 3))
        frame.grid(row=row, column=0, sticky="ew")
        frame.columnconfigure(2, weight=1)
        ttk.Label(frame, text=parameter.name, width=20).grid(row=0, column=0, sticky="w")
        value = tk.StringVar(value=f"{parameter.value:.6g}")
        low = tk.StringVar(value=f"{parameter.lower():.6g}")
        high = tk.StringVar(value=f"{parameter.upper():.6g}")
        logarithmic = is_log_slider(parameter)
        slider = ttk.Scale(frame, from_=0, to=1000, orient="horizontal")
        ttk.Entry(frame, textvariable=value, width=11).grid(row=1, column=0, sticky="w")
        slider.grid(row=0, column=1, columnspan=5, sticky="ew", padx=(4, 4))
        ttk.Label(frame, text="min").grid(row=1, column=1, sticky="e")
        ttk.Entry(frame, textvariable=low, width=9).grid(row=1, column=2, padx=(2, 5))
        ttk.Label(frame, text="max").grid(row=1, column=3, sticky="e")
        ttk.Entry(frame, textvariable=high, width=9).grid(row=1, column=4, padx=(2, 5))
        scale_label = ttk.Label(frame, text="log" if logarithmic else "linear", width=6)
        scale_label.grid(row=1, column=5)
        control = ParameterControl(parameter, frame, slider, value, low, high, scale_label)
        self.controls[parameter.name] = control
        self._set_slider_from_parameter(control)
        slider.configure(command=lambda position, name=parameter.name: self._slider_changed(name, position))
        for variable in (value, low, high):
            variable.trace_add("write", lambda *_args, name=parameter.name: self._entry_changed(name))

    def _bounds(self, control):
        low, high = float(control.low.get()), float(control.high.get())
        if high <= low:
            raise ValueError("upper bound must be greater than lower bound")
        if control.parameter.scale == "log" and (low <= 0 or high <= 0):
            raise ValueError("logarithmic parameters need positive bounds")
        return low, high

    def _set_slider_from_parameter(self, control):
        try:
            low, high = self._bounds(control)
            control.slider.set(value_to_slider(float(control.parameter.value), low, high,
                                               control.parameter.scale == "log"))
        except (ValueError, tk.TclError):
            pass

    def _slider_changed(self, name, position):
        control = self.controls[name]
        try:
            low, high = self._bounds(control)
            value = slider_to_value(position, low, high, control.parameter.scale == "log")
            control.parameter.set(value)
            control.value.set(f"{value:.6g}")
            self.include_manual.set(True)
            self._clear_caches()
        except ValueError:
            return

    def _entry_changed(self, name):
        control = self.controls[name]
        try:
            low, high = self._bounds(control)
            control.parameter.low, control.parameter.high = low, high
            value = float(control.value.get())
            control.parameter.set(value)
            control.slider.set(value_to_slider(value, low, high, control.parameter.scale == "log"))
            self.include_manual.set(True)
            self._clear_caches()
        except (ValueError, tk.TclError):
            # Incomplete numbers are normal while typing; preserve the last valid value.
            return

    def _manual_result(self):
        """Replay edited values while keeping the checkpoint's data selection/config."""
        checkpoint = self.active_checkpoint
        parameter_values = parameters.ModelParameters(self.working_parameters.get_config())
        model_run = run(copy.deepcopy(checkpoint.run_config), parameter_values, log_prediction=True)
        result = plot_checkpoint.plotted_checkpoint(
            "manual", checkpoint.source, checkpoint.metadata, checkpoint.run_config, model_run,
            comparison_name="manual values")
        return result

    def _clear_caches(self):
        """Discard stale replay/figure data after a source or value changes."""
        self._result_cache.clear()

    def _result_cache_key(self):
        selected = self.selected_checkpoints()
        if not selected and self.active_checkpoint is not None:
            selected = [self.active_checkpoint]
        selected_key = tuple(checkpoint.source for checkpoint in selected)
        manual_values = ()
        if self.include_manual.get() and self.working_parameters is not None:
            manual_values = tuple(sorted(
                (parameter.name, repr(parameter.value)) for parameter in self.working_parameters))
        return selected_key, manual_values

    def _comparison_results(self):
        cache_key = self._result_cache_key()
        cached = self._result_cache.get(cache_key)
        if cached is not None:
            return cached
        selected = self.selected_checkpoints()
        if not selected and self.active_checkpoint is not None:
            selected = [self.active_checkpoint]
        results = [plot_checkpoint.replay_saved_checkpoint(checkpoint) for checkpoint in selected]
        if self.include_manual.get() and self.active_checkpoint is not None:
            results.append(self._manual_result())
        self._result_cache[cache_key] = results
        return results

    def _build_figures(self):
        if self.dataset is None:
            self.dataset = data.get_data_before_day()
        results = self._comparison_results()
        if not results:
            raise ValueError("Select one or more checkpoints first.")
        names = plot_checkpoint.comparison_plot_names(results)
        if names != self.plot_names:
            current = self.current_plot.get()
            self.plot_names = names
            self.plot_picker.configure(values=names)
            self.current_plot.set(current if current in names else names[0])
        return results

    def _figures_for_results(self, results):
        """Return plot keys for a fixed set of replay results without drawing."""
        if self.dataset is None:
            self.dataset = data.get_data_before_day()
        return plot_checkpoint.comparison_plot_names(results)

    def _apply_y_scale(self, figure, log_y=None):
        log_y = self.log_y.get() if log_y is None else log_y
        for axis in figure.axes:
            if log_y:
                # The shared formatter safely refuses log scale for negative data.
                import plot
                plot.format_ax(axis, log_scale=True)
            else:
                axis.set_yscale("linear")

    def _choose_figure(self, name=None, log_y=None):
        results = self._build_figures()
        name = self.current_plot.get() if name is None else name
        if name not in self.plot_names:
            name = self.plot_names[0]
        if self.dataset is None:
            self.dataset = data.get_data_before_day()
        figure, _ = plot_checkpoint.comparison_figure(results, self.dataset, name)
        self._apply_y_scale(figure, log_y)
        figure.tight_layout(rect=(0, 0, 1, .92))
        return figure, name

    def refresh_main_plot(self):
        if self.active_checkpoint is None:
            return
        try:
            figure, name = self._choose_figure()
            if self._main_pane is not None:
                self._main_pane.destroy()
            self._main_pane = FigurePane(self.viewer, figure, f"Soil model — {name}")
        except Exception as error:
            traceback.print_exc()
            messagebox.showerror("Could not render plot", str(error))

    def detach_plot(self):
        if self.active_checkpoint is None:
            return
        try:
            # Replay once at detachment.  The separate window deliberately
            # retains these curves even if selection or sliders later change.
            frozen_results = self._comparison_results()
            frozen_names = self._figures_for_results(frozen_results)
            checkpoint = self.active_checkpoint
            identity = (f"sample {checkpoint.replica[0]}, replica {checkpoint.replica[1]}, "
                        f"{checkpoint.fit_mode}, {checkpoint.id}")
            if self.include_manual.get():
                identity += " + manual values"
            window = tk.Toplevel(self.root)
            window.geometry("900x700")
            window.title(f"Soil model comparison — {identity}")
            controls = ttk.Frame(window, padding=6)
            controls.pack(fill="x")
            ttk.Label(controls, text="Plot").pack(side="left")
            plot_name = tk.StringVar(value=self.current_plot.get())
            picker = ttk.Combobox(controls, state="readonly", textvariable=plot_name,
                                  values=frozen_names, width=38)
            picker.pack(side="left", padx=4)
            detached_log_y = tk.BooleanVar(value=self.log_y.get())
            ttk.Checkbutton(controls, text="Log y", variable=detached_log_y).pack(side="left")
            canvas_host = ttk.Frame(window)
            canvas_host.pack(fill="both", expand=True)
            state = {"pane": None}

            def redraw(*_):
                selected_name = plot_name.get()
                if selected_name not in frozen_names:
                    selected_name = frozen_names[0]
                    plot_name.set(selected_name)
                figure, _ = plot_checkpoint.comparison_figure(
                    frozen_results, self.dataset, selected_name)
                self._apply_y_scale(figure, detached_log_y.get())
                figure.tight_layout(rect=(0, 0, 1, .92))
                figure.suptitle(f"{selected_name} | {identity}", fontsize=10)
                if state["pane"] is not None:
                    state["pane"].destroy()
                state["pane"] = FigurePane(canvas_host, figure)

            picker.bind("<<ComboboxSelected>>", redraw)
            detached_log_y.trace_add("write", redraw)
            redraw()
            self._detached.append((window, state))
            window.protocol("WM_DELETE_WINDOW", lambda: self._close_detached(window, state))
        except Exception as error:
            traceback.print_exc()
            messagebox.showerror("Could not open comparison", str(error))

    def _close_detached(self, window, state):
        if state["pane"] is not None:
            state["pane"].destroy()
        self._detached = [(other_window, other_state) for other_window, other_state in self._detached
                          if other_window is not window]
        window.destroy()


def main():
    root = tk.Tk()
    ModelInspector(root)
    root.mainloop()


if __name__ == "__main__":
    main()
