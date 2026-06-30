
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
from optimizer import r2

d = data.get_data_before_day()

class FasterModel(Model):

    def __init__(self, replica):
        
        self.replica = d[replica]
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
                val = self.controls[k].get()
                self.model_parameters[i].set(val)
                
        self.model.parameters().set({p.name: float(p) for p in self.model_parameters})


    def save_parameters(self):
        save_model(self.specimen_index, self.site, self.model_parameters, prefix = 'slidertool')


class PlotView(View):
    def __init__(self,model, key_list):
        self.key_list = key_list
        super().__init__(model)

    def refresh(self):
        #fig = plt.figure()
        fig, (ax0,ax1) = plt.subplots(2,1)

        if self.key_list is None:
            plt.sca(ax0)
            self.model.replica.plot(newfigure = ax0, measurements = 'CO2')
            all_days, values = self.model.pool_value_dict['CO2']
            _,measured_co2 = self.model.replica.CO2()
            r2_co2 = r2(values, measured_co2, log = True)
            ax0.plot(all_days, values, label = f'CO2 (R² = {r2_co2:.3f})', color = 'b')
            ax0.legend()

            self.model.replica.plot(newfigure = ax1, measurements = 'CH4')
            all_days, values = self.model.pool_value_dict['CH4']
            _,measured_ch4 = self.model.replica.CH4()
            r2_ch4 = r2(values, measured_ch4, log = True)
            ax1.plot(all_days, values, label = f'CH4 (R² = {r2_ch4:.3f})', color = 'orange')
            ax1.legend()
            
            ax1.set_yscale('log')
            #ax1.set_ylim([1e-4, 1e0])
            fig.tight_layout()

            return fig

        title = self.key_list[0].split('_')[-1] if len(self.key_list) > 1 else self.key_list[0]
        for k in self.key_list:
            all_days, values = self.model.pool_value_dict[k]
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
    model = FasterModel('13544')
    run_demo(model)
