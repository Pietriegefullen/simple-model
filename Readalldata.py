# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:06:34 2021

@author: Lara
"""

#import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#m = 0
def load_realdata(m):
    
    #os.chdir('/Users/Lara/Desktop/simple model')
    
    df = pd.read_excel(r'C:\Users\Lara\Desktop\simple model\Trainingdatafull.xlsx',
                       engine = 'openpyxl')

    df = df.apply(pd.to_numeric, errors='coerce') # Macht " nicht zahlen" zu Nan
    df = df.dropna() #  löscht Zeilen mit NaN
    
    df = df.values
        
    df[:,m] = np.around(df[:,m]) # rundet die Tage auf ganze tage (trotzdem floats)
    
    #df[:,m+1][np.where(df[:,m+1]<0)] = 0 # Messfehler mit neg. CH4 werten zu 0 
    #df[:,m+2][np.where(df[:,m+2]<0)] = 0 # Messfehler mit neg. CO2 werten zu 0 
    
    
    #CH4Listneg= []
    #CH4Listneg.append(df[:,m+1][np.where(df[:,m+1]<0)])  
    #indexList = []
    #indexList = list(df[:,m+1]).index(df[:,m+1][np.where(df[:,m+1]<0)])
    indexListCH4 = np.where(np.isin(df[:,m+1],df[:,m+1][np.where(df[:,m+1]<0)]))
    indexListCO2 = np.where(np.isin(df[:,m+1],df[:,m+1][np.where(df[:,m+1]<0)]))
    #print(indexListCH4[0])
    #print(len(indexListCH4[0]))
    #print(indexListCH4)
    
    for p in range(len(indexListCH4[0])):   
        df[:,m+1][indexListCH4[0][p]] = df[:,m+1][indexListCH4[0][p]-1]  # Messfehler mit neg. CH4 werten zum vorherigen wert
        
    for r in range(len(indexListCO2)):   
        df[:,m+1][indexListCO2[r]] = df[:,m+1][indexListCO2[r]-1]  # Messfehler mit neg. CO2 werten zum vorherigen wert
                            
    
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

if __name__=='__main__':
    Realdata = load_realdata(0)
    
    plt.figure()
    plt.plot(Realdata[:,0], Realdata[:,2])
    a = [5]* len(Realdata[:,0])
    plt.plot(Realdata[:,0], a)
    plt.show()  
    #%%
    plt.close('all')
    
    plt.figure(2)
    #126 tag bis  277
    plt.plot(Realdata[9:18,2])
    plt.show()  


#%%


GstAceto = (-50.72 + -586.77) - ( -369.31 + -237.15)
GstHydro= (0.25*-50.72 +0.75*-237.13) - (0.25*-586.77 +0 + 0.25*0)
GstIron = (9*0 + 2*-586.77+ 8*-78.9)- (-369.31+4*-237.13 + 8*-4.7)
GstHomo = (-396.46 +0+2*-237.13) - (4*0 +2*-394.36)
GstFerm = (2*-396.46+2*-394.36+4*0)- (-910 +2*-237.13)