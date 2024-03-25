import os

import numpy as np
import pandas as pd


def load_metadata(source_directory):
    raw_metadata = _load_raw_metadata(source_directory)
    
    use_columns = ['sample number',
                    'depth',    
                    'Corg (%)',
                    'pH (H2O)',
                    'Wassergehalt (%)',
                    'site',
                    'origin']
    raw_metadata = raw_metadata.loc[:,use_columns]        
    raw_metadata.loc[:,use_columns[:5]] = raw_metadata.loc[:,use_columns[:5]].astype(float)
    raw_metadata['sample number'] = raw_metadata['sample number'].astype(int).astype(str)
    
    metadata_dict = {}
    for _, row in raw_metadata.iterrows():
        metadata_dict[row['sample number']] = {'depth': row['depth'],
                                               'Corg': row['Corg (%)'],
                                               'pH': row['pH (H2O)'],
                                               'site': row['site'],
                                               'origin': row['origin']}
    
    return metadata_dict
    
    
def _load_raw_metadata(source_directory):
    metadata_file = os.path.join(source_directory, 'Metadaten_all_86.xlsx')
    raw_metadata = pd.read_excel(metadata_file, 
                                 engine = 'openpyxl',
                                 header = None)
    
    drop_rows = []
    headers = None
    site = {'Kurugnakh':[], 'Samoylov':[]}
    origin = {'cliff':[], 'core':[]}
    current_site = None
    current_origin = None
    for i, row in raw_metadata.iterrows():
        if all(row.isnull()):
            drop_rows.append(i)
            continue 
        if any(row.str.contains('from', na = False)):
            is_core = any(row.str.contains('core', case = False, na = False))
            is_cliff = any(row.str.contains('cliff', case = False, na = False))
            if is_core and not is_cliff:
                current_origin = 'core'
            elif is_cliff and not is_core:
                current_origin = 'cliff'
                    
            if any(row.str.contains('Kurugnakh', na = False)):
                current_site = 'Kurugnakh'
               
            elif any(row.str.contains('Samoylov', na = False)):
                current_site = 'Samoylov'
            
            drop_rows.append(i)
            continue
        
        if any(row.str.contains('Lab-number', na = False)):
            # header row
            if not headers:
                row = row.iloc[1:]
                row.iloc[1] = 'sample number'
                headers = row.array
            drop_rows.append(i)
            continue
        
        site[current_site].append(i)
        origin[current_origin].append(i)
    
    raw_metadata['site'] = np.nan
    for s, indices in site.items():
        raw_metadata.loc[indices, 'site'] = s
        
    raw_metadata['origin'] = np.nan
    for s, indices in origin.items():
        raw_metadata.loc[indices, 'origin'] = s

    raw_metadata = raw_metadata.drop(drop_rows, axis = 0)
    raw_metadata = raw_metadata.iloc[:,1:].rename(columns = {i+1:name for i, name in enumerate(headers)})
    
    return raw_metadata
                    
def _load_raw_incubation(source_directory):
    incubation_file = os.path.join(source_directory, 'Perm - alle.xlsx')
    raw_incubation = pd.read_excel(incubation_file, engine = 'openpyxl',
                                   sheet_name = 'incubation data')
    
    raw_constants = raw_incubation.iloc[:,:3]
    raw_incubation = raw_incubation.iloc[:,3:]
    
    raw_incubation.dropna(axis = 0, how = 'all', inplace = True)
    raw_incubation.dropna(axis = 1, how = 'all', inplace = True)
    
    headers = None
    current_sample = None
    samples = {}
    previous_row = None
    for i, row in raw_incubation.iterrows():
        if row['Unnamed: 3'] == 'Probe':
            previous_row = raw_incubation.loc[i-1, :]
            titles = [(p if not str(p) == 'nan' else r) 
                      for p, r in zip(previous_row.array, row.array)]
            raw_incubation.rename(columns = {('Unnamed: ' + str(k+3)):title
                                   for k, title in enumerate(titles)},
                                  inplace = True)
            break

    for i, row in raw_incubation.iterrows():
        if type(row['Probe']) == str and row['Probe'].startswith('09-'):
            sample_name = row['Probe'].replace('09-','').replace('/','')
            if sample_name == current_sample:
                day = row['duration days total']
                samples[current_sample]['events'].append((day, 'carex'))
                continue
            
            current_sample = sample_name
            # extract replica constants (weight wet sample)
            dry_weight = row['dry weight (g)']
            water_content = row['Water content (ml)']
            samples[current_sample] = {'days':[], 'CO2':[], 'CH4':[],
                                       'events': [], 
                                       'dry weight': dry_weight,
                                       'water content': water_content}
            continue
        
        if current_sample is None:
            continue
        
        co2_value = row['cummulative CO2 produced']
        ch4_value = row['CH4 produced']
        day = row['duration days total']
        if str(co2_value) == 'nan' or str(ch4_value) == 'nan' or str(day) == 'nan':
            print(current_sample, 'excluding measurements on day', day, 'CO2:', co2_value, 'CH4:', ch4_value)
        
        else:
            samples[current_sample]['days'].append(day)
            samples[current_sample]['CO2'].append(co2_value)
            samples[current_sample]['CH4'].append(ch4_value)
    
        if isinstance(day, (int, float)) and not type(row['Probe']) is float and  not str(row['Probe']).lower() == 'nan':
            event = (day,str(row['Probe']))
            samples[current_sample]['events'].append(event)
            
    return samples

def _load_raw_ergaenzung(source_directory):
    excel_file = os.path.join(source_directory,'Daten_Erganzung.xlsx')

    raw_ergaenzung = pd.read_excel(excel_file, engine = 'openpyxl',
                                   sheet_name = 'samples without priming anaerob',
                                   header = None)
    
    raw_ergaenzung_const = raw_ergaenzung.iloc[:,:3]
    raw_ergaenzung = raw_ergaenzung.iloc[:,3:]
    
    raw_ergaenzung.dropna(axis = 1, how = 'all', inplace = True)
    raw_ergaenzung.dropna(axis = 0, how = 'all', inplace = True)
    
    headers = None
    current_replica = None
    samples = {}
    previous_row = None
    for i, row in raw_ergaenzung.iterrows():
        if row[3] == 'Probe':
            previous_row = raw_ergaenzung.loc[i-1, :]
            titles = [(p if not str(p) == 'nan' else r) 
                      for p, r in zip(previous_row.array, row.array)]
            raw_ergaenzung.rename(columns = {k+3:title
                                   for k, title in enumerate(titles)},
                                  inplace = True)
            break

    for i, row in raw_ergaenzung.iterrows():
        if str(row['Probe']).startswith('09-'):
            replica_name = row['Probe'].replace('09-','').replace('/','')
            current_replica = replica_name
            # extract replica constants (weight wet sample)
            dry_weight = row['dry weight (g)']
            water_content = row['Water content (ml)']
            
            if str(dry_weight) == 'nan' and str(water_content) == 'nan':
                dry_weight = raw_ergaenzung.loc[i+1,:]['dry weight (g)']
                water_content = raw_ergaenzung.loc[i+1,:]['Water content (ml)']
            
            samples[current_replica] = {'days':[], 'CO2':[], 'CH4':[],
                                       'events': [], 
                                       'dry weight': dry_weight,
                                       'water content': water_content}
            continue
        
        
        if current_replica is None:
            continue
        
        co2_value = row['cummulative CO2 release']
        ch4_value = row['CH4 total']
        day = row['duration (d)']
        if str(co2_value) == 'nan' or str(ch4_value) == 'nan' or str(day) == 'nan':
            print(current_replica, 'excluding measurements on day', day, 'CO2:', co2_value, 'CH4:', ch4_value)
        else:
            samples[current_replica]['days'].append(day)
            samples[current_replica]['CO2'].append(co2_value)
            samples[current_replica]['CH4'].append(ch4_value)
    
        if isinstance(day, (int,float)) and not type(row['Probe']) is float and not str(row['Probe']).lower() == 'nan':
            event = (day,str(row['Probe']))
            samples[current_replica]['events'].append(event)
            
    return samples
