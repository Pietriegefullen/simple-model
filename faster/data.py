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

DOC_per_TOC = 6.

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
        self.source_directory = ROOT_DIRECTORY
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
        return self.initial_C()*DOC_per_TOC

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

if __name__ == '__main__':
    d = get_data_before_carex()
    print('====')
    print(d)


