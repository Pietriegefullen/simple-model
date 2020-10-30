# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 09:36:26 2020

@author: Lara
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_realdata():
    
    os.chdir('/Users/Lara/Desktop/simple model')
    
    
    df = pd.read_excel(r'C:/Users/Lara/Desktop/simple model/Trainingdata.xlsx')
    
    df = df.values
    
    df[:,0] = np.around(df[:,0]) # rundet die Tage auf ganze tage (trotzdem floats)
    
    df[:,1][np.where(df[:,1]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 
    
    #adding mock C data to realdata
    # =============================================================================
    # xneg = np.linspace(100,0,len(df[:,0]))
    # noiseCpool = np.random.uniform(-2,2,(int(len(df[:,0])),))
    # Cpoolval =  list(np.exp(0.05*xneg) + noiseCpool)
    # 
    # =============================================================================
    #Realdata = np.c_[df[:,0],Cpoolval, df[:,1],df[:,2]]
    Realdata  = df
    plt.plot(Realdata[:,1], 'ro',label="CH4")
    plt.plot(Realdata[:,2],'bo',label ="CO2")
    plt.legend()
    
    return Realdata