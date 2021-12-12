
from predict import run
import matplotlib.pyplot as plt

if __name__ == '__main__':

    m_gluc = 180                                                            # molar mass of glucose, g dw pro mol
    TOC =  4/100 #Realdata['Corg (%)']/100.0                                       # Knoblauchs Daten , g dw
    specimen_mass = 20.#Realdata['weight']                                      # g (Knoblauch proben)
    Cpool_init = (10**6)* specimen_mass * TOC / m_gluc

    model_parameters = {'M_Fe3_init':   0.2,
                        'C_init':       Cpool_init,
                        'Sensenmann':   0.002,
                        'Vmax_Fe3':     1.0}
    pool_value_dict = run(model_parameters)

    plt.plot(pool_value_dict['C'])
    plt.show()
