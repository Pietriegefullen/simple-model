
import matplotlib.pyplot as plt

def all_pools(pool_value_dict, all_days, measured_data = None, show = True):

    del(pool_value_dict['weight'])
    del(pool_value_dict['pH'])
    del(pool_value_dict['water'])
    if 'HCO3' in pool_value_dict:
        del(pool_value_dict['HCO3'])


    # check keys:
    # if ending is same, use same subplot
    same_plot = list()
    for k in pool_value_dict.keys():
        found = False
        if '_' in k:
            for i,s in enumerate(same_plot):
                if '_' in s[0] and s[0].split('_')[-1] == k.split('_')[-1]:
                    same_plot[i].append(k)
                    found = True
                    break

        if not found:
            same_plot.append([k])


    plot_figures = list()
    for key_list in same_plot:
        fig = plt.figure()
        plot_figures.append(fig)
        title = key_list[0].split('_')[-1] if len(key_list) > 1 else key_list[0]
        for k in key_list:
            plt.plot(all_days,
                     pool_value_dict[k],
                     label = k.replace('_'+title, '') if len(key_list) > 1 else None)
            empty = [0] * 7500
            plt.plot(all_days, empty)
        if len(key_list)>1:
            plt.legend()
        plt.title(title)

        if not measured_data is None and title in measured_data:
            plt.plot(measured_data['measured_time'],
                      measured_data[title],
                      'rx',
                      label = 'measured')

    if show:
        plt.show()

    return fig

def plot_pool(name, values, time):

    plt.figure()
    plt.plot(time, values)
    plt.title(name)


def fit(days, CO2, pool_values, all_days):

    plt.figure()
    plt.plot(days, CO2, 'x')
    plt.plot(all_days, pool_values, '-')
    print('\007')
