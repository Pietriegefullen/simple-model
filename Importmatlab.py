# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 11:04:59 2021

@author: Lara
"""
import scipy.io as sio
import os
import copy

os.chdir('C:/Users/Lara/Desktop/simple model')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_matlab():
    
    Data = sio.loadmat('ActivityData_04062016', appendmat=True)
    
    anaerobic_indices = np.nonzero(Data['anaerobic'])[1] # die indizes im Matlab für alle anaeroben proben
    prob_array = Data['prob'] # das array mit den Probennummern
    measurement_time = Data['duration'] # das array mit allen Messtagen
    
    # erstelle die äußerste struktur, ein leeres dict für jede probe hat
    superdata = dict()
    for index_for_anaerobic_sample in anaerobic_indices:
        prob_name = prob_array[0,index_for_anaerobic_sample] #holt aus prob_array den probennamen an der indexstelle von index_for_anaerobic _sample.
        
        keys_so_far = [ k[:4]  for k in superdata.keys()] #schreibt die ersten vier ziffern aller keys die wir schon haben in eine liste
        replica_count = keys_so_far.count(str(prob_name)) # zählt wie viele keys aus der liste die gleichen 4 anfangsziffern haben wie die aktuelle probname
        
        replica_name = str(prob_name) + str(replica_count)
        superdata[replica_name] = dict()
        
        superdata[replica_name]['measured_time'] = measurement_time[:,index_for_anaerobic_sample]
        superdata[replica_name]['CH4'] = Data['ch4'][:,index_for_anaerobic_sample]
        superdata[replica_name]['CO2'] = Data['co2'][:,index_for_anaerobic_sample]
        # analog für co2, ch4, ...
        
    Metadata = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Metadaten_all_86_cleanandneat.xlsx',engine = 'openpyxl')
    
    for meta_prob_names in Metadata['Probe']:
        for replica_name_iter in superdata.keys():
            if str(meta_prob_names) in replica_name_iter:
            #if str(meta_prob_names) == replica_name_iter[:4]:
                index = np.where(Metadata['Probe'] == meta_prob_names)[0]
                for key in Metadata.keys():
                    superdata[replica_name_iter][key] = np.array(Metadata[key][index])
     
    replica_list = list(superdata.keys())    
    
    for key in superdata.keys():
        for index in range(len(superdata[key]['CH4'])):
            if superdata[key]['CH4'][index] < 0 :
                superdata[key]['CH4'][index] = superdata[key]['CH4'][index-1]
            if np.isnan(superdata[key]['CH4'][index]):
                superdata[key]['CH4'][index] = 0    
    
            if superdata[key]['CO2'][index] < 0 :
                superdata[key]['CO2'][index] = superdata[key]['CO2'][index-1]         
            if  np.isnan(superdata[key]['CO2'][index]):
                superdata[key]['CO2'][index] = 0

             
    superdata_carex = copy.deepcopy(superdata)
    for key in superdata.keys():
        FirstNan = np.where(np.isnan(superdata[key]['measured_time']))[0][0]
        for column_name in ['CO2','CH4','measured_time']:
            superdata_carex[key][column_name] = superdata_carex[key][column_name][(FirstNan+1):]
            superdata[key][column_name] = superdata[key][column_name][:FirstNan]
             
 
    for key in superdata:
        plt.figure()  
        plt.plot(superdata[key]['measured_time'], superdata[key]['CH4'], "r", label= "CH4")       
        plt.plot(superdata[key]['measured_time'], superdata[key]['CO2'],"b",  label= "CO2")
        plt.title(str(superdata[key]['Probe'])  + "_____" + superdata[key]['Site'] + "_____" + superdata[key]['Location']+ "_____" + str(superdata[key]['depth']))
     
        
    superdata_Kuru = copy.deepcopy(superdata)
    for key in superdata.keys():
        if not superdata[key]['Site']=='K':
            del superdata_Kuru[key]        
    replica_list_Kuru = list(superdata_Kuru.keys()) 
        
    superdata_Sam = copy.deepcopy(superdata)
    for key in superdata.keys():
        if not superdata[key]['Site']=='S':
            del superdata_Sam[key]        
    replica_list_Sam = list(superdata_Sam.keys())   
    print(len(replica_list_Sam) , "Proben Sam")                 
    print(len(replica_list_Kuru) , "Proben Kuru")  
    print(len(replica_list) , "Proben insgesamt")          

    return superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam




if __name__ == '__main__':
    
    
   
    superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam = load_matlab()
