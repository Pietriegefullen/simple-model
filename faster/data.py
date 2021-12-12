import os
import scipy.io as sio
import pandas as pd
import numpy as np
import copy

from USER_VARIABLES import ROOT_DIRECTORY
import CONSTANTS

def model_parameters_from_data(specimen_index, site):

    data = specimen_data(specimen_index, site)

    TOC =  data['Corg (%)']/100.0                                          # Knoblauchs Daten , g dw
    specimen_mass = data['weight']                                         # g (Knoblauch proben)
    Cpool_init = (10**6)* specimen_mass * TOC / CONSTANTS.MOLAR_MASS_GLUCOSE

    model_parameters = {}

    model_parameters['C'] = float(Cpool_init)
    model_parameters['DOC'] = float(Cpool_init)*0.02                      # 0.02 ratio in song, 2% sind DOC laut den Christians.
    model_parameters['pH'] = data['pH (H2O)']#[0]
    model_parameters['weight'] = float(specimen_mass)
    model_parameters['water'] = float(data['water'])

    return model_parameters

def specimen_data(specimen_index, site):

    (superdata,
     replica_list,
     superdata_carex,
     superdata_Kuru,
     superdata_Sam,
     replica_list_Kuru,
     replica_list_Sam,
     superdata_2021_all,
     replica_list_superdata_2021_all,
     superdata_ohne_Fe3,
     Rep_ohne_Fe3,
     superdata_mit_Fe3,
     Rep_mit_Fe3) = load_matlab()

    if site == "S":
        data = superdata_Sam[replica_list_Sam[specimen_index]]
    elif site == "K":
        data = superdata_Kuru[replica_list_Kuru[specimen_index]]
    elif site == "all":
        data = superdata_2021_all[replica_list_superdata_2021_all[specimen_index]]

    return data

def load_matlab():

    """
     superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()

     superdata sind alle Datensätze VOR dem Carexexperiment,
     superdata_carex sind die Datensätze von Christian NACH dem Carex experiment ( eigentlich alt )
     superdata_all sind alle Datensätze von Christian vor und nach Carex (auch alt)
     Superdata_2021_all sind alle Datensätze vor und Nach Carex Zugabe von Knoblauch, darauf bauen auf:
     Superdata Kuru sind die Datensätze von Kurunak vor und nach Carex
     Superdata_sam sind die Datensätze von Samoylov vor und nach Carex
     superdata_ohne_Fe3 sind alle Datensätze deren exp  plot des CH4 vermuten lässt dass kein Fe3 vorhanden ist
     superdata_mit_Fe3 sind alle Datensätze deren exp  plot des CH4 vermuten lässt dass Fe3 vorhanden ist

      """


    print('load_data')
    mat_file = os.path.join(ROOT_DIRECTORY,'ActivityData_04062016')
    Data = sio.loadmat(mat_file, appendmat=True)

    anaerobic_indices = np.nonzero(Data['anaerobic'])[1] # die indizes im Matlab für alle anaeroben proben
    prob_array = Data['prob'] # das array mit den Probennummern
    measurement_time = Data['duration'] # das array mit allen Messtagen

    # erstelle die äußerste struktur, ein leeres dict für jede probe hat
    superdata = dict()
    for index_for_anaerobic_sample in anaerobic_indices:
        prob_name = prob_array[0,index_for_anaerobic_sample] #holt aus prob_array den probennamen an der indexstelle von index_for_anaerobic _sample.

        keys_so_far = [ k[:4]  for k in superdata.keys()] #schreibt die ersten vier zifFe3rn aller keys die wir schon haben in eine liste
        replica_count = keys_so_far.count(str(prob_name)) # zählt wie viele keys aus der liste die gleichen 4 anfangszifFe3rn haben wie die aktuelle probname

        replica_name = str(prob_name) + str(replica_count)
        superdata[replica_name] = dict()

        # erstellt das dict für die anaeroben proben
        superdata[replica_name]['measured_time'] = measurement_time[:,index_for_anaerobic_sample]
        superdata[replica_name]['CH4'] = Data['ch4'][:,index_for_anaerobic_sample]
        superdata[replica_name]['CO2'] = Data['co2'][:,index_for_anaerobic_sample]


    # fügt den Proben jeweils die relevanten Metadaten hinzu
    excel_file = os.path.join(ROOT_DIRECTORY,'Metadaten_all_86_cleanandneat.xlsx')
    Metadata = pd.read_excel(excel_file,engine = 'openpyxl')

    metadata_keys = [str(int(s)) for s in Metadata['Probe'] if not np.isnan(s)]
    for specimen_key, specimen_data in superdata.items():
        try:
            metadata_row = metadata_keys.index(specimen_key[:-1])
            for meta_key, meta_value in Metadata.items():
                if meta_key == 'Probe':
                    continue
                superdata[specimen_key][meta_key] = meta_value[metadata_row]

        except ValueError:
            print('No metadata for specimen', specimen_key)

    replica_list = list(superdata.keys())


    for key in superdata.keys():
        #ersetzt negative messwerte mit den jeweils vorherigen messwerten, und nan mit 0
        for index in range(len(superdata[key]['CH4'])):
            if superdata[key]['CH4'][index] < 0 :
                superdata[key]['CH4'][index] = superdata[key]['CH4'][index-1]
            if np.isnan(superdata[key]['CH4'][index]):
                superdata[key]['CH4'][index] = 0

            if superdata[key]['CO2'][index] < 0 :
                superdata[key]['CO2'][index] = superdata[key]['CO2'][index-1]
            if  np.isnan(superdata[key]['CO2'][index]):
                superdata[key]['CO2'][index] = 0

    # findet den ersten nan wert in measured_time und überträgt alles danach in superdata_carex
    superdata_carex = copy.deepcopy(superdata)
    superdata_all = copy.deepcopy(superdata)

    for key in superdata.keys():
        FirstNan = np.where(np.isnan(superdata[key]['measured_time']))[0][0]
        for column_name in ['CO2','CH4','measured_time']:
            superdata[key][column_name] = superdata[key][column_name][:FirstNan]
            superdata_carex[key][column_name] = superdata_carex[key][column_name][(FirstNan+1):]
              # key entry ab wann Carex zugegeben wurde
            superdata[key]['First_Carex_index'] =     len(superdata[key]['measured_time'])
            superdata[key]['Last_non_Carex_day'] =    max(superdata[key]['measured_time'])

            superdata_all[key][column_name] = superdata_all[key][column_name][:]


    # for key in superdata:
    #     plt.figure()
    #     plt.plot(superdata[key]['measured_time'], superdata[key]['CH4'], "r", label= "CH4")
    #     plt.plot(superdata[key]['measured_time'], superdata[key]['CO2'],"b",  label= "CO2")
    #     plt.title(str(superdata[key]['Probe'])  + "_____" + superdata[key]['Site'] + "_____" + superdata[key]['Location']+ "_____" + str(superdata[key]['depth']))





    # Ersetzen der alten Carexdaten mit den neuen Daten von Knoblauch bis 2021
    excel_file = os.path.join(ROOT_DIRECTORY, 'Carex_addition_2021_clean.xlsx')
    Carex_addition = pd.read_excel(excel_file,engine = 'openpyxl')
    Carex = np.array(Carex_addition)
    #Proben = np.unique(Carex[:,0])
   # Probenfree = np.unique(Proben[~np.isnan(Proben)])

   # vergleich der keys mit den Probennamen des array und ermitteln der Indizes
    for key in superdata_carex.keys():
        indices = []
        for row_index, row in enumerate(Carex):
            try:
                specimen_key = int(row[0])

            except ValueError:
                # Leere Zeilen verursachen NaN value als key, wird ignoriert
                continue

            if specimen_key == int(key):
                indices.append(row_index)

        if len(indices) == 0:
            print('No carex data for specimen', key)
        #hier werden die neuen Werte eingefügt
        superdata_carex[key]['measured_time'] = Carex[indices,1]
        superdata_carex[key]['CO2'] = Carex[indices,2]
        superdata_carex[key]['CH4'] = Carex[indices,3]

        superdata_carex[key]['measured_time']=[int(i) for i in superdata_carex[key]['measured_time']]# von float to int

    superdata_2021_all = copy.deepcopy(superdata)
    replica_list_superdata_2021_all = list(superdata_2021_all.keys())

    for key in  superdata_2021_all.keys():
        #anfügen der Carextage an die nicht carextage
          Time_to_append = [i + max(superdata_2021_all[key]['measured_time'])  for i in superdata_carex[key]['measured_time'][1:]]

          existing_time_values = superdata_2021_all[key]['measured_time']
          time_all_values = np.concatenate([existing_time_values,Time_to_append], axis = 0)
          superdata_2021_all[key]['measured_time'] = time_all_values

         #anfügen der Carex CH4 messwerte an die nicht carex CH4 Messwerte

          CH4_to_append = [i + superdata_2021_all[key]['CH4'][-1]  for i in superdata_carex[key]['CH4'][1:]]

          existing_CH4_values = superdata_2021_all[key]['CH4']
          CH4_all_values = np.concatenate([existing_CH4_values,CH4_to_append], axis = 0)
          superdata_2021_all[key]['CH4'] = CH4_all_values

          #anfügen der Carex CH4 messwerte an die nicht carex CH4 Messwerte

          CO2_to_append = [i + superdata_2021_all[key]['CO2'][-1]  for i in superdata_carex[key]['CO2'][1:]]

          existing_CO2_values = superdata_2021_all[key]['CO2']
          CO2_all_values = np.concatenate([existing_CO2_values,CO2_to_append], axis = 0)
          superdata_2021_all[key]['CO2'] = CO2_all_values


    #erstellt einen Datensatz nur für die Kurunak daten
    superdata_Kuru = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key]['Site']=='K':
            del superdata_Kuru[key]
    replica_list_Kuru = list(superdata_Kuru.keys())

    #erstellt einen datensatz nur für die Samoylov daten
    superdata_Sam = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key]['Site']=='S':
            del superdata_Sam[key]
    replica_list_Sam = list(superdata_Sam.keys())

     # erstellen der Datenreihen, die vermutlich ohne Fe3 Pool sind
    Rep_ohne_Fe3 = [13690, 13531,13530,13521,13520,13742, 13741, 13740,13732,13731, 13730,13721,13720]

    superdata_ohne_Fe3 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key] in Rep_ohne_Fe3:
            del superdata_ohne_Fe3[key]


     # erstellen der Datenreihen, die vermutlich MIT Fe3 Pool sind

    Rep_mit_Fe3 = [13510,13511,13512,13670,13671,13672,13691,13692,13700,13722, 13731, 13750,13751,13752]

    superdata_mit_Fe3 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not superdata_2021_all[key] in Rep_mit_Fe3:
            del superdata_mit_Fe3[key]

    #print(superdata_Kuru['13510']['First_Carex_index'])


    return superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3
