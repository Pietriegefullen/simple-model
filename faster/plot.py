
import matplotlib.pyplot as plt

def all_pools(pool_value_dict, all_days):

    for k, v in pool_value_dict.items():
        plot_pool(k,v,all_days)

    plt.show()

def plot_pool(name, values, time):

    plt.figure()
    plt.plot(time, values)
    plt.title(name)


def fit(days, CO2, pool_values, all_days):

    plt.figure()
    plt.plot(days, CO2, 'x')
    plt.plot(all_days, pool_values, '-')
