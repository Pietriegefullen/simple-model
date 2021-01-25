# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:06:34 2021

@author: Lara
"""

#import os
import pandas as pd
import numpy as np


def load_realdata(m):
    
    #os.chdir('/Users/Lara/Desktop/simple model')
    
    df = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Trainingdatafull.xlsx')
    
    df = df.apply(pd.to_numeric, errors='coerce') # Macht " nicht zahlen" zu Nan
    df = df.dropna() #  löscht Zeilen mit NaN
    
    df = df.values
        
    df[:,m] = np.around(df[:,m]) # rundet die Tage auf ganze tage (trotzdem floats)
    
    df[:,m+1][np.where(df[:,m+1]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 
    df[:,m+2][np.where(df[:,m+2]<0)] = 0 # Messfehler mit neg. CO2 werten zu 0 
    
    #adding mock C data to realdata
    # =============================================================================
    # xneg = np.linspace(100,0,len(df[:,0]))
    # noiseCpool = np.random.uniform(-2,2,(int(len(df[:,0])),))
    # Cpoolval =  list(np.exp(0.05*xneg) + noiseCpool)
    #Realdata = np.c_[df[:,0],Cpoolval, df[:,1],df[:,2]]
    # =============================================================================
 
   
    
    Realdata = df[:,m:m+3]
    
    #print(Realdata[:10, :3])
    #Realdata = Realdata[Realdata[:,0]>200,:] # throw out all data before day 200
    
    return Realdata

#Realdata = load_realdata(9)