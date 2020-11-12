# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 10:03:03 2020

@author: Lara
"""


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_realdata():
    
    os.chdir('/Users/Lara/Desktop/simple model')
    
    
    df = pd.read_excel(r'C:/Users/Lara/Desktop/simple model/Trainingdatafull.xlsx')
    
    df = df.values
    
    df[:,3] = np.around(df[:,3]) # rundet die Tage auf ganze tage (trotzdem floats)
    
    df[:,4][np.where(df[:,4]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 
    
    #adding mock C data to realdata
    # =============================================================================
    # xneg = np.linspace(100,0,len(df[:,0]))
    # noiseCpool = np.random.uniform(-2,2,(int(len(df[:,0])),))
    # Cpoolval =  list(np.exp(0.05*xneg) + noiseCpool)
    # 
    # =============================================================================
    #Realdata = np.c_[df[:,0],Cpoolval, df[:,1],df[:,2]]
    Realdata  = df
    plt.plot(Realdata[:,4], 'ro',label="CH4")
    plt.plot(Realdata[:,5],'bo',label ="CO2")
    plt.legend()
    
    Realdata2 = df[:,3:6]
    
    return Realdata2


#%%

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt    

os.chdir('/Users/Lara/Desktop/simple model')


df = pd.read_excel(r'C:/Users/Lara/Desktop/simple model/Trainingdatafull.xlsx')

df = df.values

df[:,3] = np.around(df[:,3]) # rundet die Tage auf ganze tage (trotzdem floats)

df[:,4][np.where(df[:,4]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 

#adding mock C data to realdata
# =============================================================================
# xneg = np.linspace(100,0,len(df[:,0]))
# noiseCpool = np.random.uniform(-2,2,(int(len(df[:,0])),))
# Cpoolval =  list(np.exp(0.05*xneg) + noiseCpool)
# 
# =============================================================================
#Realdata = np.c_[df[:,0],Cpoolval, df[:,1],df[:,2]]
Realdata  = df
plt.plot(Realdata[:,4], 'ro',label="CH4")
plt.plot(Realdata[:,5],'bo',label ="CO2")
plt.legend()

plt.figure()
Realdata = df[:,3:6]
plt.plot(Realdata[:,1], 'ro',label="CH4")
plt.plot(Realdata[:,2],'bo',label ="CO2")
plt.legend()

