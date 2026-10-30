import matplotlib.pyplot as plt
import numpy as np

class ModelRun():
    def __init__(self):
        self._log = {}
        self._parameters = None
        self.rate = False

    def keys(self):
        return self._log.keys()

    def __eq__(self, other):
        return self._log == other._log

    def __getitem__(self, key):
        return self._log[key]

    def log_snap(self, name, t, value):
        if not name in self._log:
            ts = np.empty((0,))
            vs = np.empty((0,))
            self._log[name] = (ts, vs)
        ts, vs = self._log[name]
        ts = np.concatenate([ts, np.reshape(t,(1,))],
                            axis = 0)
        vs = np.concatenate([vs, np.reshape(value,(1,))],
                            axis = 0)
        self.log(name, ts, vs)

    def log_snap_2(self, name, t, value):
        if not name in self._log:
            self._log[name] = [[],[]]
        self._log[name][0].append(t)
        self._log[name][1].append(value)

    def log(self, name, ts, values):
        #if not name in self._log:
        #    self._log[name] = []
        #if 'v' in name:
        #    print('logging', name)
        self._log[name] = (ts, values)

    def reset(self, replica):
        self._log.clear()
        self.replica = replica

    def R2(self, gas_name, log_fit):
        t_pred, pred_val = self.__getitem__(gas_name)
        if gas_name == 'CO2':
            t_meas, meas_val = self.replica.CO2(rate = self.rate)
        elif gas_name == 'CH4':
            t_meas, meas_val = self.replica.CH4(rate = self.rate)
        else:
            raise Exception()

        t = np.intersect1d(np.round(t_pred,2), np.round(t_meas,2))
        idx_pred = np.squeeze([np.nonzero(np.round(t_pred,2) == np.round(_t,2))[0] for _t in t])
        idx_meas = np.squeeze([np.nonzero(np.round(t_meas,2) == np.round(_t,2))[0] for _t in t])
        r2_val = r2(pred_val[idx_pred], meas_val[idx_meas], log = log_fit)

        return r2_val

    def plot_factors(self, pathway):

        v_name = pathway + '_v'
        if not v_name in self._log:
            return

        fig, ax = plt.subplots()

        t, v = self._log[v_name]

        ax.plot(t,v, 'k-', linewidth = 1.)

        _pwy = pathway
        _pwy = _pwy.replace('Aceto', 'Ac')
        _pwy = _pwy.replace('Fermentation', 'Ferm')
        v_max = float(self._parameters[_pwy + '_v_max'])
        ax.plot(ax.get_xlim(), [v_max, v_max], 'k--', label = 'v_max', linewidth = 1.)

        ax.set_ylim([0,v_max])

        #ax.legend()

        ax2 = ax.twinx()

        all_factors = ['MM', 'inhib', 'thermodynamic_factor']
        factor_colors = ['r', 'b', 'g']

        for factor, col in zip(all_factors, factor_colors):
            name = pathway + '_' + factor
            if not name in self.keys():
                continue

            t, f = self._log[name]
            ax2.plot(t,f, col + '-', linewidth = .5, label = factor)
        ax2.set_ylim([0,1.0])
        #ax2.legend()
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        handles, labels = ax.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()

        fig.legend(handles + handles2, labels + labels2, loc = 'lower right')
        plt.title(v_name)

        return ax

    def plot(self, name = None,
             newfigure = True,
             log = False,
             label = None,
             save_target = None,
             xlim = None):
        if name is None:
            name = list(self._log.keys())

        if not isinstance(name, list):
            name = [name]

        ordered_names = name

        for i, n in enumerate(ordered_names):
            if 'MM' in n or 'inhib' in n or 'thermodynamic_factor' in n:
                pass #continue

            if n.endswith('_v'):
                pathway_name = n.replace('_v','')
                ax = self.plot_factors(pathway_name)
                n += '_factors'

            elif not n in self._log:
                print(n + ' not logged')
                continue

            else:
                if newfigure:
                    fig, ax = plt.subplots()
                else:
                    fig = plt.gcf()
                    ax = plt.gca()

                x, y = self._log[n]
                slabel = n
                if not label is None and not label == '':
                    slabel += ' ' + label

                mark = '-'
                #if 'R2' in self._log and n in self._log['R2']:
                 #   value = self._log['R2'][n]
                  #  label += ' ' + f'R² = {value:4.2f}'
                   # mark = 'x'
                if n == 'CH4':
                    ax = fig.axes
                    if len(name) > 1:
                        ax = ax[0]
                        ch4_ax = ax.twinx()
                    else:
                        ch4_ax = ax[0]
                    ch4_ax.plot(x, y, mark, label = slabel, color = 'orange')

                else:
                    ch4_ax = ax
                    ch4_ax.plot(x, y, mark, label = slabel)

                ax = ch4_ax
                ax.set_title(n)
                plt.legend(loc = 'lower right')

                if log:
                    ax.set_yscale('log')
                    ax.set_title(n + ' (log)')

                elif n in system.SYSTEM:
                    plt.yscale('log')
                    ax.set_ylim([0,ax.get_ylim()[1]])
                    pass

                elif 'inhib' in n or 'thermodynamic_factor' in n:
                    ax.set_ylim([0.,1])

                elif 'MM' in n:
                    ax.set_yscale('log')
                    ax.set_ylim([1e-5, 1e0])

                if not xlim is None:
                    ax.set_xlim(xlim)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)

            if not save_target is None:
                file_name = '_'.join([f'{i+1:02d}', n])
                plt.savefig(os.path.join(save_target, file_name + '.png'), dpi = 300)

    def __str__(self):
        run_string = 'Model run:'
        run_string += '\n' + '='*len(run_string) + '\n'
        run_string += '\n'.join([name + ' ' + str(self._log[name])
                          for name in sorted(self._log.keys())])
        ts, _ = self._log[next(iter(self._log.keys()))]
        t_start = np.min(ts)
        t_end = np.max(ts)
        run_string += '\n' + f't in [{t_start:.2f}, {t_end:.2f}]'
        return run_string

