
import matplotlib.pyplot as plt

def all_pools(pool_value_dict, all_days, measured_data = None):
    
    del(pool_value_dict['weight'])
    del(pool_value_dict['pH'])
    del(pool_value_dict['water'])
    del(pool_value_dict['HCO3'])

    for k, v in pool_value_dict.items():
        plot_pool(k,v,all_days)
        
        if not measured_data is None and k in measured_data:
            plt.plot(measured_data['measured_time'],
                     measured_data[k],
                     'rx')

    plt.show()

def plot_pool(name, values, time):

    plt.figure()
    plt.plot(time, values)
    plt.title(name)


def fit(days, CO2, pool_values, all_days):

    plt.figure()
    plt.plot(days, CO2, 'x')
    plt.plot(all_days, pool_values, '-')
