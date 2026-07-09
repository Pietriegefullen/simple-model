import os
import scipy.io as sio
import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
import traceback
import json

from USER_VARIABLES import ROOT_DIRECTORY
import CONSTANTS

import loading

DOC_per_TOC = 0.02

knoblauch_data = None

def save_to_json(d, target):
    cfg = d.get_config()
    target_file = os.path.join(target, 'KnoblauchData.json')
    if not os.path.isdir(target):
        os.makedirs(target)
    with open(target_file, 'w') as tf:
        json.dump(cfg, tf, indent = 4)

def get_data_before_carex():
    global knoblauch_data
    if knoblauch_data is None:
        knoblauch_data = load_knoblauch()
    _ = [r.before_day(r.carex()) for r in knoblauch_data.replicas()]
    return knoblauch_data

def load_knoblauch():
    print('Loading Knoblauch data...')
    d_file = os.path.join(ROOT_DIRECTORY, 'KnoblauchData.json')
    if os.path.isfile(d_file):
        with open(d_file, 'r') as df:
            d_dict = json.load(df)
        knoblauch_data = KnoblauchData(d_dict) # load from json file
        
    else:
        knoblauch_data = KnoblauchData() # load from Excel files
        save_to_json(knoblauch_data, ROOT_DIRECTORY)
    return knoblauch_data

def get_data_before_day():
    global knoblauch_data
    if knoblauch_data is None:
        knoblauch_data = load_knoblauch()
    _ = [r.before_day(r.last_day) for r in knoblauch_data.replicas()]
    return knoblauch_data
    

class KnoblauchData():
    def __init__(self, data_dict = None):
        if data_dict is None:
            samples = []

        else:
            samples = [Sample(**s) for s in data_dict['samples']]
            
        self.source_directory = ROOT_DIRECTORY
        
        self.last_days = {
                     '13514': '1274',
                     '13515': '1274',
                     '13516': '1526',#'(726,1526)',
                     
                     '13525': '1274',
                     '13526': '1274',
                     
                     '13534': '1274',
                     '13535': '1274',
                     
                     '13544': '1275',
                     '13546': '1275',
                     
                     '13554': '1309',
                     '13555': '1309',
                     
                     '13575': '1275',
                     '13576': '1275',
                     
                     '13584': '1308',
                     '13585': '1308',
                     
                     '13594': '1308',
                     '13595': '1308',
                     
                     '13604': '1308',
                     '13606': '1308',
                     
                     '13614': '1308',
                     '13616': '1308',
                     
                     '13624': '1308',
                     '13626': '1308',
                     
                     '13634': '1308',
                     '13635': '1308',
                     
                     '13654': '1274',
                     '13655': '1274',
                     '13656': '1274',
                     
                     '13665': '1274',
                     '13666': '1274',
                     
                     '13674': '1273',
                     '13675': '1273',
                     '13676': '1273',
                     
                     '13684': '1274',
                     '13685': '1274',
                     '13686': '1274',
                     
                     '13694': '1273',
                     '13695': '1273',
                     '13696': '1274',
                     
                     '13704': '1273',
                     '13706': '1274',
                     
                     '13724': '1260',
                     '13725': '1260',
                     '13726': '1260',
                     
                     '13734': '1260',
                     '13735': '1260',
                     '13736': '1260',
                     
                     '13744': '1260',
                     '13745': '1260',
                     '13746': '1260',
                     
                     '13754': '1260',
                     '13755': '1260',
                     '13756': '1260',
                     
                     '13764': '1261',
                     '13765': '1261',
                     '13766': '1261',
                     
                     '13774': '1260',
                     '13775': '1260',
                     '13776': '1260',
                     
                     '13784': '712',
                     '13785': '1260',
                     '13786': '1260',
                     
                     '13794': '1294',
                     '13795': '1294',
                     '13796': '1294',
                     
                     '13804': '1294',
                     '13805': '1294',
                     '13806': '1294',
                         }
    
        self.samples = []
        for s in samples:
            self.add_sample(s)

        if not self.samples:
            self.load()
        
    def load(self):
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
                new_replica.last_day = None
                if replica_name in self.last_days:
                    new_replica.last_day = float(self.last_days[replica_name])
                new_sample.add_replica(new_replica)

            self.add_sample(new_sample)

    def __getitem__(self, key):
        key = str(key).replace('/','')
        if len(key) == 4:
            return [s for s in self.samples if s.sample_name == key][0]
        elif len(key) == 5:
            return [r for r in self.replicas() if str(r) == key][0]
        raise Exception('Invalid sample or replica name.')

    def __contains__(self, key):
        try:
            self.__getitem__(key)
            return True
        except Exception:
            return False
        
    def add_sample(self, sample):
        try:
            check_sample(sample)        
            self.samples.append(sample)
            
        except AssertionError as ex:
            print(f'Skipping sample {str(sample)}: {str(ex)}')
    
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

    def get_config(self):
        cfg = {'samples': [s.get_config() for s in self.samples]}
        return cfg

class Sample():
    def __init__(self, sample_name, site = None, origin = None, depth = None, pH = None, TOC = None, replicas = None):
        super().__init__()
        
        self.sample_name = sample_name
        loaded_replicas = []
        if not replicas is None:
            for r in replicas:
                repl = Replica(**r)
                repl.sample = self
                loaded_replicas.append(repl)
        self.replicas = list(sorted(loaded_replicas,
                                     key = lambda r: r.replica_number))
        
        self.site = site
        self.origin = origin
        self.depth = depth
        self.pH = pH
        self.TOC = TOC # as decimal, e.g. 3% is 0.03
    
    def __contains__(self, replica_number):
        return any([str(r.replica_number) == str(replica_number) for r in self.replicas])
    
    def add_replica(self, replica):
        replica.sample = self
        try:
            check_replica(replica)
            self.replicas.append(replica)
          
        except AssertionError as ex:
            print(f'Skipping replica {str(replica)}: {str(ex)}')
            print(traceback.format_exc())
        
    def has_replicas(self):
        return len(self.replicas)
    
    def plot(self):
        marker = iter(['x','v','+', 's', 'o', '^'])
        for r in self.replicas:
            r.plot(marker = next(marker))
            
        plt.title(f'{self.sample_name} {self.site} ({self.origin})')
    
    def __str__(self):
        rep_names = ','.join([r.replica_number for r in self.replicas])
        return f'{self.sample_name} {self.site} ({self.origin}) {len(self.replicas)} replicas ({rep_names})'
    
    def leave_one_out_split(self):
        if len(self.replicas) <= 1:
            raise Exception('Cannot verify a fit to a sample with only one replica.')
            
        all_splits = []
        self.replicas = sorted(self.replicas, key = lambda r: int(r.replica_number))
        for i in range(len(self.replicas)):
            validation_replica = self.replicas[i]
            fit_replicas = [self.replicas[(i+k+1)%len(self.replicas)]
                            for k in range(len(self.replicas)-1)]
            all_splits.append({'fit': sorted(fit_replicas, key = lambda r: int(r.replica_number)), 
                               'val': validation_replica})
        return all_splits
        
    def get_config(self):
        cfg = {'sample_name': self.sample_name,
               'replicas': [r.get_config() for r in self.replicas],
               'site': self.site,
               'origin': self.origin,
               'depth': self.depth,
               'pH': self.pH,
               'TOC': self.TOC
               } 
        return cfg
    
class Replica():
    def __init__(self, replica_number,
                 dry_weight = None,
                 water_content = None,
                 temperature = None,
                 events = None,
                 incubation = None,
                 last_day = None):
        self.sample = None
        self.replica_number = replica_number
        
        self.dry_weight = dry_weight # directly from Knoblauch, unit is g
        self.water_content = water_content # directly from Knoblauch, ml 
        if temperature is None:
            temperature = CONSTANTS.SPECIMEN_TEMPERATURE
        self.temperature = temperature
        
        if events is None:
            events = {}
        self.events = events
        if incubation is None:
            incubation = {}
        self.incubation = {k:np.reshape(v, (-1,)) for k, v in incubation.items()} # incubation data is from Knoblauch per g_dw
        
        self.last_day = last_day
        
    def initial_TOC(self): # reines C (schwerverfügbar, nur Hydrolyse)
        # micro-mol per g dw
        return (10**6)*self.sample.TOC/CONSTANTS.MOLAR_MASS_C_g_mol

    def initial_DOC(self): # reines C, verfügbar (Fermentations)
        return self.initial_TOC()*DOC_per_TOC

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
        print('Found no carex event in ', str(self))

    def __getitem__(self, key):
        if key == 'CO2':
            return self.incubation['CO2']
        elif key == 'CH4':
            return self.incubation['CH4']
        elif key == 'days':
            return self.incubation['days']
    
    def plot(self, events = True, marker = 'x', log = False, newfigure = True, measurements = None, label = ''):
        if measurements is None:
            measurements = ['CO2', 'CH4']
        elif not isinstance(measurements, list):
            measurements = [measurements]
            
        if isinstance(newfigure, bool):
            if newfigure:
                fig, ax = plt.subplots()
            else:
                ax = plt.gca()
                fig = plt.gcf()
        else:
            ax = newfigure
        
        if not label == '':
            label = f' ({label})'
        
        if 'CO2' in measurements:
            ax.plot(*self.CO2(),'r' + marker, label = 'CO2' + label, color = 'b')
            ax.set_ylabel('gas (CO2)')
        
        if 'CH4' in measurements:
            ax.plot(*self.CH4(),'b' + marker, label = 'CH4' + label, color = 'orange')
            ax.set_yscale('log')
            ax.set_ylabel('gas (CH4)')

        ax.set_title(f'{str(self)} {self.sample.site} ({self.sample.origin})')
        
        handles, labels = [], []
        for axs in fig.axes:
            handles += axs.get_legend_handles_labels()[0]
            labels += axs.get_legend_handles_labels()[1]
        
        fig.legend(handles, labels)
        ax.set_xlabel('day')
        #ax.set_ylabel('gas')
        
        if not events:
            return
        
        if log and newfigure:
            ax.set_yscale('log')
        else:
            pass
            #max_ = max([np.max(self.CO2()[1]), np.max(self.CH4()[1])])
            #ax.set_ylim([0, max_])
        ylim = ax.get_ylim()
        for event, day in self.events.items():
            ax.plot([day, day], ylim, 'r-')
            ax.text(day-100, np.min(ylim), event, rotation = 'vertical')
        
    def before_day(self, last_day):
        if last_day is None:
            last_day = self.last_day
            if last_day is None:
                last_day = self.carex()
                if last_day is None:    
                    last_day = max(self.incubation['days'])
        days_before = [d for d in self.incubation['days'] if d < last_day]
        co2_before = self.incubation['CO2'][:len(days_before)]
        ch4_before = self.incubation['CH4'][:len(days_before)]
        self.incubation = {'days': np.reshape(days_before, (-1,)),
                           'CO2': np.reshape(co2_before,(-1,)),
                           'CH4': np.reshape(ch4_before, (-1,))}
        self.events = {event:day for event,day in self.events.items() 
                       if day < last_day}
            
    def __str__(self):
        return str(self.sample.sample_name) + str(self.replica_number)
    
    def get_config(self):
        cfg = {'replica_number': self.replica_number,
               'dry_weight': self.dry_weight,
               'water_content': self.water_content,
               'temperature': self.temperature,
               'events': self.events,
               'incubation': {k:v.tolist() for k, v in self.incubation.items()},
               'last_day': self.last_day
               }

        return cfg

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

    # make sure measurements are strictly increasing
    previous_day = None
    days = []
    co2_values = []
    ch4_values = []
    for d, co2, ch4 in zip(replica.incubation['days'], replica.incubation['CO2'], replica.incubation['CH4']):
        if any(['na' in s.lower() for s in [str(d), str(co2), str(ch4)]]):
            continue
        if not previous_day is None and d <= previous_day:
            continue
        days.append(d)
        co2_values.append(co2)
        ch4_values.append(ch4)
        prevous_day = d
    replica.incubation['days'] = np.reshape(days, (-1,))
    replica.incubation['CO2'] = np.reshape(co2_values, (-1,))
    replica.incubation['CH4'] = np.reshape(ch4_values, (-1,))

    assert len(replica.incubation['CO2']) == len(replica.incubation['days'])
    assert len(replica.incubation['CH4']) == len(replica.incubation['days'])

    assert all([isinstance(s, (int, float)) 
                for s in replica.incubation['days']]), replica.incubation['days']
    assert all([isinstance(s, (int, float)) 
                for s in replica.incubation['CO2']]), replica.incubation['CO2']
    assert all([isinstance(s, (int, float)) 
                for s in replica.incubation['CH4']]), replica.incubation['CH4']

 
    assert not replica.sample is None
    
    assert isinstance(replica.replica_number, (int, str))
    assert len(str(replica.replica_number)) == 1
    assert int(replica.replica_number) <=6 and int(replica.replica_number) > 0
    
    assert isinstance(replica.dry_weight, (int, float))
    assert replica.dry_weight > 0, f'{str(replica)} has dry weight {replica.dry_weight}'
    
    assert isinstance(replica.water_content, (int, float))
    assert replica.water_content > 0
    
    assert replica.temperature == 277.15

def save_data(d):
    d_file = os.path.join(ROOT_DIRECTORY, 'KnoblauchData.json')
    with open(d_file, 'w') as df:
        json.dump(d.get_config(), df, indent = 4) 

if __name__ == '__main__':
    d = get_data_before_day()
    
    for sample in d.samples:
        for replica in sample.replicas:
            t, ch4 = replica.CH4()
            print(replica, np.min(t), ch4[np.argmin(t)])
    1/0
    
    for s in d.samples:
        print(s)
        try:
            for r in s.leave_one_out_split():
                print('   '+', '.join([str(_t) for _t in r['fit']]))
        except:continue
    print('====')
    print(d)
    
