
import matplotlib.pyplot as plt

def all_pools(pool_value_dict, all_days, specimen_index,  measured_data = None, show = True, plot_only = None):
#def all_pools(pool_value_dict, all_days,  measured_data = None, show = True):
    dont_plot = ['weight', 'pH', 'water', 'HCO3']
    for name in dont_plot:
        if name in pool_value_dict.keys():
            del pool_value_dict[name]
        
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

    plot_figures = dict()
    for key_list in same_plot:
        title = key_list[0].split('_')[-1] if len(key_list) > 1 else key_list[0]
        fig = plt.figure(title)
        for k in key_list:
            plt.plot(all_days,
                     pool_value_dict[k],
                     label = k.replace('_'+title, '') if len(key_list) > 1 else None)
            # fügt linie bei 0 ein
            #empty = [0] * 4500
            #plt.plot(all_days, empty)
        if len(key_list)>1:
            plt.legend()
        #plt.title(title)    
        
        plt.title(str(specimen_index) + ': ' + title)

        if not measured_data is None and title in measured_data:
            plt.plot(measured_data['measured_time'],
                      measured_data[title],
                      'rx',
                      label = 'measured')
            
        plot_figures[title] = fig

    if show:
        plt.show()

    return plot_figures

def plot_pool(name, values, time, style = '-'):

    plt.figure()
    plt.plot(time, values, style)
    plt.title(name)


def fit(days, CO2, pool_values, all_days):

    plt.figure()
    plt.plot(days, CO2, 'x')
    plt.plot(all_days, pool_values, '-')
    print('\007')
