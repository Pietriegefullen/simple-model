
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TKAgg")


import os
from datetime import datetime

from tkinter import *
from tkinter.ttk import Frame, Button

from demo import Model, View, Window, run_demo, MockModel, Slider

import matplotlib.pyplot as plt

from main import run_model, save_model
import pathways
import data
import plot
import USER_VARIABLES
import OPTIMIZATION_PARAMETERS

import numpy as np
import json

class FasterModel(Model):

    def __init__(self):
        self.specimen_index = '13760'
        self.site = 'all'
        self.measured_data = data.specimen_data(self.specimen_index, self.site)
        self.model_parameters = pathways.default_model_parameters(self.specimen_index, self.site)
        # self.figure_list = list()
        self.days = 4500
        self.pool_value_dict = run_model(self.model_parameters,
                                        np.arange(self.days),
                                        extended_output = None)

        same_plot = list()
        for k in self.pool_value_dict.keys():
            found = False
            if '_' in k:
                for i,s in enumerate(same_plot):
                    if '_' in s[0] and s[0].split('_')[-1] == k.split('_')[-1]:
                        same_plot[i].append(k)
                        found = True
                        break

            if not found:
                same_plot.append([k])

        self.plot_keys = same_plot

        super().__init__()

    def _get_views(self):
        return [PlotView(self, key_list) for key_list in self.plot_keys]

    def controls(self, container):

        self.controls = {
                        'run': Button(container,
                        text = 'run',
                        command = self.update_model),
                        'save': Button(container,
                        text = 'save',
                        command = self.save_parameters),
                        }
        self.controls.update({'days': Slider(container,
                                               value = self.days,
                                               name = 'days',
                                               low = 1,
                                               high = 4500)})
        
        guess_bounds = OPTIMIZATION_PARAMETERS.get_initial_guesses()
        for k, v in self.model_parameters.items():
            lower_bound = 0
            upper_bound = 10*v
            if k in guess_bounds:
                upper_bound = 5*guess_bounds[k][2]
            self.controls.update({k: Slider(container,
                                           value = v,
                                           name = k,
                                           low = lower_bound,
                                           high = upper_bound)})

        return self.controls

    def bindings(self):
        
        return {'<Return>':lambda e: self.update_model()}

    def update_model(self):
        self.update_parameters()
        self.pool_value_dict = run_model(self.model_parameters,
                                        np.arange(self.days),
                                        extended_output = None)
        self.notify()

    def update_parameters(self):
        for k, slider in self.controls.items():
            if k in self.model_parameters:
                self.model_parameters[k] = self.controls[k].get()

            elif k == 'days':
                self.days = self.controls[k].get()

    def save_parameters(self):
        save_model(self.specimen_index, self.site, self.model_parameters, prefix = 'slidertool')


class PlotView(View):
    def __init__(self,model, key_list):
        self.key_list = key_list
        super().__init__(model)

    def refresh(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        all_days = np.arange(self.model.days)
        title = self.key_list[0].split('_')[-1] if len(self.key_list) > 1 else self.key_list[0]
        for k in self.key_list:
            ax.plot(all_days,
                     self.model.pool_value_dict[k],
                     label = k.replace('_'+title, '') if len(self.key_list) > 1 else None)
        if len(self.key_list)>1:
            ax.legend()
        ax.set_title(title)

        if not self.model.measured_data is None and title in self.model.measured_data:
            ax.plot(self.model.measured_data['measured_time'],
                      self.model.measured_data[title],
                      'rx',
                      label = 'measured')

        ax.set_xlim(np.min(all_days), np.max(all_days))

        return fig


# TODO: group parameters, and/or use scroll bars.
# show more plots simultaneously
# select plot from list
# TODO: store parameters

if __name__ == '__main__':
    # model = MockModel()
    model = FasterModel()
    run_demo(model)
