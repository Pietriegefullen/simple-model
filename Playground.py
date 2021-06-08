# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 08:45:18 2021

@author: Lara
"""
import pandas as pd
import os
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
#import numpy as np

#Importieren benutzter Funktionen und Daten
#from Realdataall import load_realdata
from Readalldata import load_realdata
from SimpleOUT import simplefun, simplesolve
from SimpleOUT import optifun
#from Babypascal import Mol_nach_Pa

#%%

plt.close('all')

#Realdata9 = load_realdata(24)

Data1and2and3and4and5and6and7and8and9 =[0,3,6,9,12,15,18,21,24]
#





for m in Data1and2and3and4and5and6and7and8and9:
    
    
    
    Realdata = load_realdata(m)
    
    t = Realdata[:,0]
    CH4 = Realdata[:,1]
    
    fig, axs = plt.subplots(2)

    #plt.figure()
    axs[0].plot(t,CH4, 'ro')
    #plt.show()
    
   # plt.figure()
    axs[1].plot(t,CH4, "r-")
    plt.yscale('log')
    plt.ylim([1e-8, 5e1]) 
    
    INDsmallremoved = np.where(CH4>1e-4)
    VALnewlistCH4 = CH4[INDsmallremoved]
    VALnewlistT = t[INDsmallremoved]
    #plt.plot(t[INDsmallremoved],VALnewlistCH4, 'red')

    INDlargeremoved = np.where(VALnewlistCH4 < (max(VALnewlistCH4)/2))

    axs[1].plot(VALnewlistT[INDlargeremoved],VALnewlistCH4[INDlargeremoved], 'blue')

    p = np.polyfit(VALnewlistT[INDlargeremoved], np.log(VALnewlistCH4[INDlargeremoved]), 1)
    
    #plt.plot(t[INDlargeremoved],CH4[INDlargeremoved], 'red')
    
    m = p[0]
    b = p[1]
    print(m,b)
    axs[1].plot(t, np.exp(m*t + b))

    plt.show()
    
#%%
from scipy.io import savemat
import scipy.io as sio
import copy

def load_matlab():
    
    Data = sio.loadmat('ActivityData_04062016', appendmat=True)
    
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
            superdata_carex[key][column_name] = superdata_carex[key][column_name][(FirstNan+1):]
            superdata[key][column_name] = superdata[key][column_name][:FirstNan]
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
    Carex_addition = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Carex_addition_2021_clean.xlsx',engine = 'openpyxl')
    Carex = np.array(Carex_addition)     
    #Proben = np.unique(Carex[:,0])     
   # Probenfree = np.unique(Proben[~np.isnan(Proben)])
    
   # vergleich der keys mit den Probennamen des array und ermitteln der Indizes
    for key in superdata_carex.keys():
        indices = []
        
        for i in range(len([Carex[:,0]][0])):
            if str(int([Carex[i,0]][0])) == str(key):
                indices.append(i)
        #hier werden die neuen Werte eingefügt
        superdata_carex[key]['measured_time'] = Carex[indices,1]
        superdata_carex[key]['CO2'] = Carex[indices,2]
        superdata_carex[key]['CH4'] = Carex[indices,3]
        
        superdata_carex[key]['measured_time']=[int(i) for i in superdata_carex[key]['measured_time']]# von float to int
        
    superdata_2021_all = copy.deepcopy(superdata)
    replica_list_superdata_2021_all = list(superdata_2021_all.keys()) 
    
    for key in  superdata_2021_all.keys():
        #anfügen der Carextage an die nicht carextage
          Time_to_append = [i + max(superdata_2021_all[key]['measured_time'] )  for i in superdata_carex[key]['measured_time']]

          existing_time_values = superdata_2021_all[key]['measured_time']
          time_all_values = np.concatenate([existing_time_values,Time_to_append], axis = 0)
          superdata_2021_all[key]['measured_time'] = time_all_values
          
         #anfügen der Carex CH4 messwerte an die nicht carex CH4 Messwerte
         
          CH4_to_append = [i + superdata_2021_all[key]['CH4'][-1]  for i in superdata_carex[key]['CH4']]

          existing_CH4_values = superdata_2021_all[key]['CH4']
          CH4_all_values = np.concatenate([existing_CH4_values,CH4_to_append], axis = 0)
          superdata_2021_all[key]['CH4'] = CH4_all_values
          
          #anfügen der Carex CH4 messwerte an die nicht carex CH4 Messwerte
         
          CO2_to_append = [i + superdata_2021_all[key]['CO2'][-1]  for i in superdata_carex[key]['CO2']]

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
    
    print(superdata_Kuru['13510']['First_Carex_index'])
    
                
    return superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3


#export of all files in matlab format
os.chdir('C:/Users/Lara/Desktop/Marius')
if __name__ == '__main__':
    
    superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()
    
    for probe in superdata_2021_all.keys():
        name = probe +'.mat'
        file = superdata_2021_all[probe]
        savemat(name, file)
        
        


# superdata sind alle Datensätze VOR dem Carexexperiment,
# superdata_carex sind die Datensätze von Christian NACH dem Carex experiment ( eigentlich alt )
# superdata_all sind alle Datensätze von Christian vor und nach Carex (auch alt)
# Superdata_2021_all sind alle Datensätze vor und Nach Carex Zugabe von Knoblauch, darauf bauen auf: 
# Superdata Kuru sind die Datensätze von Kurunak vor und nach Carex
# Superdata_sam sind die Datensätze von Samoylov vor und nach Carex
# superdata_ohne_Fe3 sind alle Datensätze deren exp  plot des CH4 vermuten lässt dass kein Fe3 vorhanden ist
# superdata_mit_Fe3 sind alle Datensätze deren exp  plot des CH4 vermuten lässt dass Fe3 vorhanden ist






















#%%

#for m in Data1and2and3and4and5and6and7and8and9:
#    
#    

#    Realdata = load_realdata(m)
#    
#    t = Realdata[:,0]
#    CH4 = Realdata[:,1]
#    
#
#    fig, axs = plt.subplots(2)
#    
#    axs[0].plot(t,CH4)
#    #plt.show()
#    
#    #plt.figure()
#    axs[1].plot(t,CH4)
#    plt.yscale('log')
#    plt.ylim([1e-8, 5e1])
#    
# 
#    
#    m, b = np.polyfit(t, CH4, 1)
#    print(m,b)
#    axs[1].plot(t, m*t + b)
#    plt.yscale('log')
#    plt.ylim([1e-8, 5e1])
#    plt.show()
#    
    
    
 


#%% 
      
for m in Data1and2and3and4and5and6and7and8and9:
    
    
    
    Realdata = load_realdata(m)
    
    t = Realdata[:,0]
    CO2 = Realdata[:,2]
    

    fig, axs = plt.subplots(2)
    
    axs[0].plot(t,CO2)
    #plt.show()
    
    #plt.figure()
    axs[1].plot(t,CO2)
    plt.yscale('log')
    plt.ylim([1e-8, 5e1])
    plt.show()
    
    
#%%

import pandas as pd
import os
df = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Gibbs.xlsx')
    
    
    
    
    
    