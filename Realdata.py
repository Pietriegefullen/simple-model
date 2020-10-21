# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 09:36:26 2020

@author: Lara
"""
import os
import pandas as pd
import numpy as np

os.chdir('/Users/Lara/Desktop/simple model')


df = pd.read_excel(r'C:/Users/Lara/Desktop/simple model/Trainingdata.xlsx')

df = df.values

df[:,0] = np.around(df[:,0]) # rundet die Tage auf ganze tage (trotzdem floats)

df[:,1][np.where(df[:,1]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 

Realdata = df






