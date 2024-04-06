import os
import scipy.io as sio
import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
import traceback


from USER_VARIABLES import ROOT_DIRECTORY
import CONSTANTS

import loading

# TODO: check units (replica mass, water content, ...)
# TODO: check replica: water content must have correct OOM

knoblauch_data = None

def get_data_before_carex():
    global knoblauch_data
    if knoblauch_data is None:
        print('Loading Knoblauch data...')
        knoblauch_data = KnoblauchData()
        _ = [r.before_carex() for r in knoblauch_data.replicas()]
    return knoblauch_data

class KnoblauchData():
    def __init__(self):
        self.source_directory = 'D:\OneDrive - Universität Hamburg\daten und code'
        self.samples = []
        #superdata = load_matlab('superdata')
        
        print('loading incubation data')
        incubation_data = loading._load_raw_incubation(self.source_directory)
        print('loading ergaenzung data')
        ergaenzung_data = loading._load_raw_ergaenzung(self.source_directory)

        incubation_data.update(ergaenzung_data)
        
        print('loading metadata')
        metadata = loading.load_metadata(self.source_directory)
        
        print('building samples and replicas')
        for sample_name, sample_metadata in metadata.items():
            new_sample = Sample(sample_name)
            new_sample.site = sample_metadata['site']
            new_sample.TOC = sample_metadata['Corg']/100
            new_sample.pH = sample_metadata['pH']
            new_sample.depth = sample_metadata['depth']
            new_sample.origin = sample_metadata['origin']
            
            for replica_name, replica_data in incubation_data.items():
                if not replica_name.startswith(sample_name):
                    continue

                new_replica = Replica(replica_name[-1])
                new_replica.events = {event:day for (day, event) in replica_data['events']}
                new_replica.dry_weight = replica_data['dry weight']
                new_replica.incubation = {'days': np.reshape(replica_data['days'], (-1,)),
                                          'CO2': np.reshape(replica_data['CO2'], (-1,)),
                                          'CH4': np.reshape(replica_data['CH4'], (-1,))}
                new_replica.water_content = replica_data['water content']
                new_sample.add_replica(new_replica)

            self.add_sample(new_sample)

    def __getitem__(self, key):
        key = str(key).replace('/','')
        if len(key) == 4:
            return [s for s in self.samples if s.sample_name == key][0]
        elif len(key) == 5:
            return [r for r in self.replicas() if str(r) == key][0]
        raise Exception('Invalid sample or replica name.')

    def add_sample(self, sample):
        try:
            check_sample(sample)        
            self.samples.append(sample)
            
        except AssertionError as ex:
            print(f'Skipping sample {str(sample)}: {str(ex)}')
            print(traceback.format_exc())
    
    def replicas(self):
        return [r for s in self.samples for r in s.replicas]

    def plot_replicas(self):
        for r in self.replicas():
            plt.figure()
            r.plot()
        plt.show()
    
    def plot_samples(self):
        for s in self.samples:
            plt.figure()
            s.plot()
        plt.show()
    
    def __str__(self):
        total_replicas = sum([s.has_replicas() for s in self.samples])
        text = '\n'.join([str(s) for s in self.samples]) 
        text += f'\n {len(self.samples)} samples, {total_replicas} replicas'
        return text
    
class Sample():
    def __init__(self, sample_name):
        super().__init__()
        
        self.sample_name = sample_name
        self.replicas = []
        
        self.site = None
        self.origin = None
        self.depth = None
        self.pH = None
        self.TOC = None # as decimal, e.g. 3% is 0.03
        
        
    def add_replica(self, replica):
        replica.sample = self
        try:
            check_replica(replica)
            self.replicas.append(replica)
          
        except AssertionError as ex:
            print(f'Skipping replica {str(replica)}: {str(ex)}')
        
    def has_replicas(self):
        return len(self.replicas)
    
    def plot(self):
        marker = iter(['x','v','+', 's', 'o', '^'])
        for r in self.replicas:
            r.plot(marker = next(marker))
            
        plt.title(f'{self.sample_name} {self.site} ({self.origin})')
    
    def __str__(self):
        return f'{self.sample_name} {self.site} ({self.origin}) {len(self.replicas)} replicas'
    
    def leave_one_out_split(self):
        if len(self.replicas) <= 1:
            raise Exception('Cannot verify a fit to a sample with only one replica.')
            
        all_splits = []
        for i in range(len(self.replicas)):
            validation_replica = self.replicas[i]
            fit_replicas = [self.replicas[(i+k+1)%len(self.replicas)]
                            for k in range(len(self.replicas)-1)]
            all_splits.append({'fit': fit_replicas, 
                               'val': validation_replica})
        return all_splits
        
    
class Replica():
    def __init__(self, replica_number):

        self.sample = None
        self.replica_number = replica_number
        
        self.dry_weight = None # directly from Knoblauch, unit is g
        self.water_content = None # directly from Knoblauch, ml 
        self.temperature = CONSTANTS.SPECIMEN_TEMPERATURE
        
        self.events = {}
        self.incubation = {} # incubation data is from Knoblauch per g_dw
        
    def initial_C(self):
        # micro-mol per g dw
        return (10**6)*self.sample.TOC/CONSTANTS.MOLAR_MASS_GLUCOSE
        
    def initial_DOC(self):
        return self.initial_C()*.02
    
    def initial_H2O(self):
        # micro mol per g dw
        relative_water_content = self.water_content/self.dry_weight # g_H2O/g_dw
        return relative_water_content/CONSTANTS.MOLAR_MASS_H2O*1e6
    
    def CO2(self):
        return self.incubation['days'], self.incubation['CO2']
    
    def CH4(self):
        return self.incubation['days'], self.incubation['CH4']
    
    def carex(self):
        if 'carex' in self.events.keys():
            return self.events['carex']
        return None
    
    def __getitem__(self, key):
        if key == 'CO2':
            return self.incubation['CO2']
        elif key == 'CH4':
            return self.incubation['CH4']
        elif key == 'days':
            return self.incubation['days']
    
    def plot(self, events = True, marker = 'x'):
        plt.plot(*self.CO2(),'r' + marker, label = 'CO2')
        plt.plot(*self.CH4(),'b' + marker, label = 'CH4')
        plt.title(f'{str(self)} {self.sample.site} ({self.sample.origin})')
        plt.legend()
        plt.xlabel('day')
        plt.ylabel('gas')
        
        if not events:
            return
        
        ax = plt.gca()
        ylim = ax.get_ylim()
        for event, day in self.events.items():
            plt.plot([day, day], ylim, 'r-')
            plt.text(day-100, 0, event, rotation = 'vertical')
        
    def before_carex(self):
        carex_day = self.carex()
        if carex_day is None:
            print('Found no carex event in ', str(self))
            return 
        days_before = [d for d in self.incubation['days'] if d < carex_day]
        co2_before = self.incubation['CO2'][:len(days_before)]
        ch4_before = self.incubation['CH4'][:len(days_before)]
        self.incubation = {'days': days_before,
                           'CO2': co2_before,
                           'CH4': ch4_before}
        self.events = {event:day for event,day in self.events.items() 
                       if day <= carex_day}
            
    def __str__(self):
        return self.sample.sample_name + self.replica_number


def check_sample(sample):
    assert not sample.sample_name is None
    assert len(sample.sample_name) == 4
    assert all([c in '0123456789' for c in sample.sample_name])
    
    assert sample.site == 'Kurugnakh' or sample.site == 'Samoylov'
    assert sample.origin == 'cliff' or sample.origin == 'core'
    
    assert not sample.depth is None
    assert isinstance(sample.depth, (int, float))
        
    assert not sample.pH is None
    assert isinstance(sample.pH, (int, float))
    assert sample.pH >= 0 and sample.pH <= 14
    
    assert not sample.TOC is None
    assert isinstance(sample.TOC, (int, float))
    assert sample.TOC >= 0 and sample.TOC <= 1.0
    
def check_replica(replica):
    assert 'days' in replica.incubation
    assert 'CO2' in replica.incubation
    assert 'CH4' in replica.incubation
    
    assert all([isinstance(s, (int, float)) 
                for s in replica.incubation['days']]), replica.incubation['days']
    assert all([isinstance(s, (int, float)) 
                for s in replica.incubation['CO2']]), replica.incubation['CO2']
    assert all([isinstance(s, (int, float)) 
                for s in replica.incubation['CH4']]), replica.incubation['CH4']

    
    assert len(replica.incubation['CO2']) == len(replica.incubation['days'])
    assert len(replica.incubation['CH4']) == len(replica.incubation['days'])
    
    assert not replica.sample is None
    
    assert isinstance(replica.replica_number, (int, str))
    assert len(str(replica.replica_number)) == 1
    assert int(replica.replica_number) <=6 and int(replica.replica_number) > 0
    
    assert isinstance(replica.dry_weight, (int, float))
    print(str(replica), replica.dry_weight)
    assert replica.dry_weight > 0
    
    assert isinstance(replica.water_content, (int, float))
    assert replica.water_content > 0
    
    assert replica.temperature == 277.15
        





def specimen_sites(specimen_indices):    
    (replica_list_No_CH4,
     superdata_No_CH4_vor_Impfung, 
     superdata_No_CH4_nach_Impfung,
     superdata_after_No_CH4,
     superdata_bevor_No_CH4,
     superdata_No_CH4,
     superdata,
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
    
    all_sites = {}
    for specimen_index in specimen_indices:
        all_sites[specimen_index] = []
    
        if str(specimen_index) in superdata_2021_all:
            all_sites[specimen_index].append('all')
            
        if str(specimen_index) in superdata_No_CH4_vor_Impfung:
            all_sites[specimen_index].append('No-CH4')
            
        if str(specimen_index) in superdata_Kuru:
            all_sites[specimen_index].append('K')
    
        if str(specimen_index) in superdata_Sam:
            all_sites[specimen_index].append('S')

    return all_sites

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
    model_parameters['H2O'] = float(data['water'])/CONSTANTS.MOLAR_MASS_H2O*1e6

    return model_parameters

def specimen_data(specimen_index, site):

    (replica_list_No_CH4,
     superdata_No_CH4_vor_Impfung, 
     superdata_No_CH4_nach_Impfung,
     superdata_after_No_CH4,
     superdata_bevor_No_CH4,
     superdata_No_CH4,
     superdata,
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
        index = list(superdata_2021_all.keys()).index(specimen_index)
        data = superdata_2021_all[replica_list_superdata_2021_all[index]]
    elif site == "No-CH4":
        index = list(superdata_No_CH4_vor_Impfung.keys()).index(specimen_index)
        data = superdata_No_CH4_vor_Impfung[replica_list_No_CH4[index]]
        
    
    

    return data

def load_matlab(dataset = None):

    """
     # NOTE: order in returned tuple has changed.
     #       it is NO LONGER as listed here:
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


    #print('load_data')
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
        
    #___________________________________________________________________
    superdata_No_CH4 = copy.deepcopy(superdata)

    excel_file_No_CH4 = os.path.join(ROOT_DIRECTORY, 'Daten_Erganzung_temp.xlsx')
    No_CH4_addition = pd.read_excel(excel_file_No_CH4,engine = 'openpyxl')
    No_CH4 = np.array(No_CH4_addition)
    
    
    # finde unique keys
    no_ch4_keys = np.unique(No_CH4[:,0])
    

   # vergleich der keys und ermitteln der Indizes
    for key in no_ch4_keys:
        
        indices_for_key = np.nonzero(No_CH4[:,0]==key)[0]
    
        sample_name = str(int(key))
        superdata[sample_name] = {}
        superdata[sample_name]['measured_time'] = No_CH4[indices_for_key, 1]
        superdata[sample_name]['CO2'] = No_CH4[indices_for_key, 2]
        superdata[sample_name]['CH4'] = No_CH4[indices_for_key, 3]
    replica_list_superdata_2021_all = list(superdata.keys()) 
    replica_list_superdata_2021_all.extend(list(superdata_No_CH4.keys()))

   # print(superdata.keys())
   #_________________________________________________________________      


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
        
        all_nans_found = np.where(np.isnan(superdata[key]['measured_time']))[0]
        if all_nans_found.size == 0:
            continue
        FirstNan = all_nans_found[0]
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

        #if len(indices) == 0:
           # print('No carex data for specimen', key)
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
    
    
    Rep_No_CH4 = [str(k) for k in [13560, 13562, 13580, 13581, 13590, 13591, 
                                   13600, 13602, 13622, 13641]]
    superdata_No_CH4 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not key in Rep_No_CH4:
            del superdata_No_CH4[key]
    replica_list_No_CH4 = list (superdata_No_CH4.keys())        
            
            

    #die Daten bevor ich die Ergänzung mit No_CH4 vorgenommen hab
    superdata_bevor_No_CH4 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not key in [str(int(k)) for k in no_ch4_keys]:
            del superdata_bevor_No_CH4[key]
            
      
    #die No_CH4 daten 
    superdata_after_No_CH4 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if key in [str(int(k)) for k in no_ch4_keys]:
            del superdata_after_No_CH4[key]        
   
   
     # erstellen der Datenreihen, die vermutlich ohne Fe3 Pool sind
    Rep_ohne_Fe3 = [str(k) for k in 
                    [13690, 13531,13530,13521,13520,13742, 13741, 
                     13740,13732,13731, 13730,13721,13720]]
    

    superdata_ohne_Fe3 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not key in Rep_ohne_Fe3:
            del superdata_ohne_Fe3[key]


     # erstellen der Datenreihen, die vermutlich MIT Fe3 Pool sind
    Rep_mit_Fe3 = [str(k) for k in 
                    [13510,13511,13512,13670,13671,13672,13691,
                     13692,13700,13722, 13731, 13750,13751,13752]]

    superdata_mit_Fe3 = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not key in Rep_mit_Fe3:
            del superdata_mit_Fe3[key]
            
    #print('here we go',superdata_2021_all['13520']['measured_time']) 
      

    #print(superdata_Kuru['13510']['First_Carex_index'])
    

    # Alle CH4 Messwerte die nur 'Artefakte' sind werden auf 0 gesetzt. Grenzwert aus Knoblauch2018
    for i in superdata_2021_all.keys():
        superdata_2021_all[i]['CH4'] = [ 0 if x < 0.05 else x for x in superdata_2021_all[i]['CH4']]
        
        # [x for x in CH4_values]
        # [x for x in CH4_values if x < 0.05]#
        # [ (0 if x < 0.05 else x) for x in CH4_values]
        
        #plt.figure()
        #plt.plot(superdata_2021_all[i]['measured_time'] ,superdata_2021_all[i]['CH4'])
        #plt.title(i)
     #
    # Trennen der CH4 freien Proben in vor und nach der Impfung
    
   # I = next((index for index,value in enumerate(superdata_2021_all['13600']['CH4']) if value != 0), None)
    #print(I, 'index')
    
    # Trennen der daten in die proben die kein CH4 produziert haben ..
    # Vor Impfung sind die Daten wo CH4 0 ist
    #Nach Impfung sind die Daten voWO CH4 produziert wurde nach mikriobenzugabe
    superdata_No_CH4_vor_Impfung =  copy.deepcopy(superdata_2021_all)
    superdata_No_CH4_nach_Impfung = copy.deepcopy(superdata_2021_all)
    for key in superdata_2021_all.keys():
        if not key in Rep_No_CH4:
            del superdata_No_CH4_vor_Impfung[key]
            del superdata_No_CH4_nach_Impfung[key]
            
        else:
            #I = next((index for index,value in enumerate(superdata_2021_all[str(key)]['CH4']) if value != 0), None)
            I = np.argmax(np.array(superdata_2021_all[key]['CH4']) > 0)
            if np.max(superdata_2021_all[key]['CH4']) == 0:
                del superdata_No_CH4_nach_Impfung[key]
                continue
                
            for column, rows in superdata_2021_all[key].items():
                if isinstance(rows,list) or isinstance(rows, np.ndarray):
                    superdata_No_CH4_vor_Impfung[key][column] = rows[:I]
                    superdata_No_CH4_nach_Impfung[key][column] = rows[I:]
                    
            


    #superdata_No_CH4_nach_Impfung
    
    
    

    if not dataset is None:
        returned_datasets = []
        if not isinstance(dataset, list):
            dataset = [dataset]
        datasets = {'replica_list_no_CH4': replica_list_No_CH4,
                    'superdata_no_CH4_vor_impfung': superdata_No_CH4_vor_Impfung,
                    'superdata_no_CH4_nach_impfung': superdata_No_CH4_nach_Impfung,
                    'superdata_after_no_CH4': superdata_after_No_CH4,
                    'superdata_before_no_CH4': superdata_bevor_No_CH4,
                    'superdata_no_CH4': superdata_No_CH4,
                    'superdata': superdata,
                    'replica_list': replica_list,
                    'superdata_carex': superdata_carex,
                    'superdata_Kuru': superdata_Kuru,
                    'superdata_Sam': superdata_Sam,
                    'replica_list_Kuru': replica_list_Kuru,
                    'replica_list_Sam': replica_list_Sam,
                    'superdata_2021_all': superdata_2021_all,
                    'replica_list_superdata_2021_all': replica_list_superdata_2021_all,
                    'superdata_ohne_Fe3': superdata_ohne_Fe3,
                    'Rep_ohne_Fe3': Rep_ohne_Fe3
                    }
        for d in dataset:
            returned_datasets.append(datasets[d])
        
        if len(returned_datasets) == 1:
                return returned_datasets[0]
        return tuple(returned_datasets)
    
    """
    replica_list_No_CH4,
    superdata_No_CH4_vor_Impfung,
    superdata_No_CH4_nach_Impfung,
    
    superdata_after_No_CH4,
    superdata_bevor_No_CH4,
    superdata_No_CH4,
    
    superdata, 
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
    Rep_mit_Fe3
"""


    return replica_list_No_CH4,superdata_No_CH4_vor_Impfung, superdata_No_CH4_nach_Impfung, superdata_after_No_CH4,superdata_bevor_No_CH4, superdata_No_CH4, superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3


replica_list_No_CH4, superdata_No_CH4_vor_Impfung,superdata_No_CH4_nach_Impfung,superdata_after_No_CH4,superdata_bevor_No_CH4, superdata_No_CH4,  superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()


if __name__ == '__main__':
    d = get_data_before_carex()
    d.plot_samples()
    print(d)



