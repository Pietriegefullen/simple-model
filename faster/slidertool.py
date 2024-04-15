
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

import data
import USER_VARIABLES
import numpy as np
import json
import model
import parameters

d = data.get_data_before_day()

class FasterModel(Model):

    def __init__(self):
        
        self.replica = d['13546']
        self.model = model.Model(model.get_pathways('simple'))
        self.model_parameters = parameters.default_model_parameters()
        self.model.parameters().set(self.model_parameters)
        self.pool_value_dict = self.model.predict(self.replica)

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
        return [PlotView(self, None)] + [PlotView(self, key_list) for key_list in self.plot_keys]

    def controls(self, container):

        self.controls = {
                        'run': Button(container,
                        text = 'run',
                        command = self.update_model),
                        'save': Button(container,
                        text = 'save',
                        command = self.save_parameters),
                        }
        
        #guess_bounds = OPTIMIZATION_PARAMETERS.get_initial_guesses()
        for p in self.model_parameters:
            if not p.is_variable(): continue
            k = p.name
            v = float(p)
            self.controls.update({k: Slider(container,
                                           value = float(v),
                                           name = k,
                                           low = p.lower(),
                                           high = p.upper(),
                                           log_scale = p.scale == 'log')})

        return self.controls

    def bindings(self):
        return {'<Return>':lambda e: self.update_model()}

    def update_model(self):
        self.update_parameters()
        self.pool_value_dict = self.model.predict(self.replica)
        self.notify()

    def update_parameters(self):
        parnames = [p.name for p in self.model_parameters]
        for k, slider in self.controls.items():
            if k in parnames:
                i = parnames.index(k)
                self.model_parameters[i].set(self.controls[k].get())
                
        self.model.parameters().set({p.name: float(p) for p in self.model_parameters})


    def save_parameters(self):
        save_model(self.specimen_index, self.site, self.model_parameters, prefix = 'slidertool')


class PlotView(View):
    def __init__(self,model, key_list):
        self.key_list = key_list
        super().__init__(model)

    def refresh(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        if self.key_list is None:
            plt.sca(ax)
            self.model.replica.plot()
            all_days, values = zip(*self.model.pool_value_dict['CO2'])
            ax.plot(all_days, values)
            all_days, values = zip(*self.model.pool_value_dict['CH4'])
            ax.plot(all_days, values)
            return fig

        title = self.key_list[0].split('_')[-1] if len(self.key_list) > 1 else self.key_list[0]
        for k in self.key_list:
            all_days, values = zip(*self.model.pool_value_dict[k])
            ax.plot(all_days,
                     values,
                     label = k.replace('_'+title, '') if len(self.key_list) > 1 else None)
        if len(self.key_list)>1:
            ax.legend()
        ax.set_title(title)
        
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
