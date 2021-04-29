# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 11:04:59 2021

@author: Lara
"""
import scipy.io as sio
import os

os.chdir(C:\Users\Lara\Desktop\simple model)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


print("Current Working Directory " , os.getcwd())

test = sio.loadmat('ActivityData_04062016', appendmat=True)

uniqueValues = set(test.values())


    
Metadata = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Metadaten_all_86_clean.xlsx',
                       engine = 'openpyxl')


#%%

Data_dict_86 = {'measured_time':    data_array[:,day_column_index],
                'CH4':          data_array[:,CH4_column_index],
                'CO2':          data_array[:,CO2_column_index],
                'location':
                'Site':
                'Sample_number':
                'depth':
                'Canorg':
                'Corg':
                'C\N':
                'Ntot':
                'pH':
                'age':
                'd13C':
                'STD':
                '13C':
                'STD.1':
                'Wassergehalt': ],}


print(type(test))