import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expit

n = 101
x = np.linspace(0,100,n)
xneg = np.linspace(100,0,n)



noiseCO2 = np.random.uniform(-2,2,(n,))
noisemet = np.random.uniform(-2,2,(n,))
noiseCpool = np.random.uniform(-3,3,(n,))


transition = 50

CO2 = []

for i in range(0,101,1):
    if i < transition:
        CO2.append(20*expit(0.3*(i-3)) + noiseCO2[i])
    else:
        CO2.append(40*expit(0.1*(i-transition)) + noiseCO2[i])
        #CO2.append(0.01*i + CO2[transition-1] + 5*(i-transition) + 0*noiseCO2[i])
        #CO2.append(CO2[i-1] + noiseCO2[i])




Methane = []

for i in range(0,101,1):
    if i < transition:
        Methane.append(40*expit(0.3*(i-40)) + noisemet[i])
    else:
        Methane.append(80*expit(0.1*(i-transition)) + noisemet[i])
        #Methane.append(0.001*i + Methane[transition-1] +noisemet[i])


Cpool =  list(np.exp(0.05*xneg) +noiseCpool)



plt.clf()
plt.plot(CO2)
plt.plot(Methane)
plt.plot(Cpool)

aCpool = np.expand_dims(np.array(Cpool), axis = 1)
aMeth = np.expand_dims(np.array(Methane),axis = 1)
aCO2 = np.expand_dims(np.array(CO2), axis = 1)

Mockdata = np.concatenate((aCpool,aMeth,aCO2),axis = None)

#Mockdata = 
#import pandas as pd


#Mockdata_frame = pd.DataFrame(data=Mockdata, columns=["Cpool", "CH4", "CO2"])





