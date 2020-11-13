# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 10:03:03 2020

@author: Lara
"""


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_realdata(i):
    
    os.chdir('/Users/Lara/Desktop/simple model')
    
    
    df = pd.read_excel(r'C:/Users/Lara/Desktop/simple model/Trainingdatafull.xlsx')
    
    df = df.values
    
    
        
    df[:,i] = np.around(df[:,i]) # rundet die Tage auf ganze tage (trotzdem floats)
    
    df[:,i+1][np.where(df[:,i+1]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 
    
    #adding mock C data to realdata
    # =============================================================================
    # xneg = np.linspace(100,0,len(df[:,0]))
    # noiseCpool = np.random.uniform(-2,2,(int(len(df[:,0])),))
    # Cpoolval =  list(np.exp(0.05*xneg) + noiseCpool)
    # 
    # =============================================================================
    #Realdata = np.c_[df[:,0],Cpoolval, df[:,1],df[:,2]]
    Realdata  = df
    plt.plot(Realdata[:,i+1], 'ro',label="CH4")
    plt.plot(Realdata[:,i+2],'bo',label ="CO2")
    plt.legend()
    
    Realdata = df[:,i:i+3]
    
    return Realdata


#%%

#import os
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt    
#
#os.chdir('/Users/Lara/Desktop/simple model')
#
#
#df = pd.read_excel(r'C:/Users/Lara/Desktop/simple model/Trainingdatafull.xlsx')
#
#df = df.values
#
#for i in [0,3,6]:
#    
#    df[:,i] = np.around(df[:,i]) # rundet die Tage auf ganze tage (trotzdem floats)
#    
#    df[:,i+1][np.where(df[:,i+1]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 
#    
#    #adding mock C data to realdata
#    # =============================================================================
#    # xneg = np.linspace(100,0,len(df[:,0]))
#    # noiseCpool = np.random.uniform(-2,2,(int(len(df[:,0])),))
#    # Cpoolval =  list(np.exp(0.05*xneg) + noiseCpool)
#    # 
#    # =============================================================================
#    #Realdata = np.c_[df[:,0],Cpoolval, df[:,1],df[:,2]]
#    Realdata  = df
#    plt.figure()
#    plt.plot(Realdata[:,i+1], 'ro',label="CH4")
#    plt.plot(Realdata[:,i+2],'bo',label ="CO2")
#    plt.legend()
#    
#    Realdata2 = df[:,i:i+3]