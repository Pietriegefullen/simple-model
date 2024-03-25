import matplotlib.pyplot as plt

import model
import pathways

# TODO: parallel coordinates plot for optimization runs
# TODO: throw unused pools out of the system!
#       -> build system as model is configured? add only as needed
#       -> substances and microbes
# TODO: shorter, cleaner print
# TODO: encapsulate optimization.

# TODO: rebuild compare_model functionality:
#       run simle or complex
#       save
#       load/plot
# TODO: check matlab data and compare. why different?
# TODO: global color scheme?
# TODO: compute measures of fit (Whose responsibility?)


def get_pathways(model_type):
    basic = [pathways.Hydrolysis,
             pathways.Fermentation,
             pathways.Hydro,
             pathways.Aceto]
    if model_type == 'complex':
        return basic + [pathways.Homo,
                        pathways.Fe3]
    elif model_type == 'simple':
        return basic
    else:
        raise NotImplementedError()

if __name__ == '__main__':
    import data
    d = data.get_data_before_carex()
        
    model_type = 'simple'
    chosen_pathways = get_pathways(model_type)
    model = model.Model(chosen_pathways)
    
    par = {
    "death_rate": 8.33e-05,
    #"Acetate": 1,
    #"temperature": 4.0,
    #"C": 2546.5533333333337,
    #"DOC": 50.93106666666667,
    #"pH": 3.95,
    #"weight": 11.82,
    #"water": 4.0,
    #"H2O": 222033.74024716797,
    #"M_Fe3": 0.15377556552732402,
    #"M_Ferm": 0.29386173195040044,
    #"M_Hydro": 0.43281383907529236,
    #"M_Homo": 0.2574949776769526,
    "Hydrolysis_v_max": 0.6166340111649226,
    "Ferm_v_max": 1.070276371909682,
    #"Fe3_v_max": 1.3804244562902253,
    #"Homo_v_max": 0.9318093189492231,
    "Hydro_v_max": 0.7064582813317815,
    "Ac_v_max": 0.4599047701351146,
    "Hydrolysis_Kmb": 288.99678942466437,
    "Aceto_Km_Ac": 145.99705636830586,
    #"Km_Homo_CO2": 376.8013896720814,
    #"Km_Homo_H2": 688.3240608121672,
    "Hydro_Km_CO2": 661.7562751340953,
    "Hydro_Km_H2": 497.8934720994153,
    #"Km_Fe3_Fe3": 173.34626557916957,
    #"Km_Fe3_Acetate": 637.4030208609411,
    "Ferm_Km": 160.15461453008587,
    "Ferm_inhibition": 4.643075236732733,
    #"Fe3": 81.99097055433658,
    #"M_Ac": 0.014042559314258995,
    "Ferm_CUE": 0.30944032735284144,
    #"CUE_Fe3": 0.012291327263939777,
    "Ac_CUE": 0.5724913805190271,
    #"Homo_CUE": 0.4988932700684491,
    "Hydro_CUE": 0.5054549655151662
}
 
    model.parameters().set('default')
    model.parameters().set(par)

    replica = d['13514']
    model.fit(replica)
    
    results = model.predict(replica, range(1501))
    #print(results['CO2'][-1])
    results.plot(['CH4', 'CO2'], newfigure = False)
    replica.plot()
    plt.show()
    
    
    
    