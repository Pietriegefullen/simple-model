#Generatin sample data to fit my model to

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expit

n = 101
x = np.linspace(0,100,n)
xneg = np.linspace(100,0,n) # n equal steps between 100 and 0


#creating noise to add to the data out of uniform distribution 
noiseCO2 = np.random.uniform(-2,2,(n,))
noisemet = np.random.uniform(-2,2,(n,))
noiseCpool = np.random.uniform(-3,3,(n,))

transition = 50

CO2 = []
 
# two expit distributions after each other. (i-3): the center of the expit dist
# is moved by 3 to the right. m* expit: m is the upper limit. 0.3*(i-3) stretches the values
for i in range(0,101,1):
    if i < transition:
        CO2.append(20*expit(0.3*(i-3)) + noiseCO2[i])
    else:
        CO2.append(40*expit(0.1*(i-transition)) + noiseCO2[i])
       



Methane = []

for i in range(0,101,1):
    if i < transition:
        Methane.append(40*expit(0.3*(i-40)) + noisemet[i])
    else:
        Methane.append(80*expit(0.1*(i-transition)) + noisemet[i])



Cpool =  list(np.exp(0.05*xneg) +noiseCpool)



plt.clf()
plt.plot(CO2)
plt.plot(Methane)
plt.plot(Cpool)

aCpool = np.expand_dims(np.array(Cpool), axis = 1)
aMeth = np.expand_dims(np.array(Methane),axis = 1)
aCO2 = np.expand_dims(np.array(CO2), axis = 1)

Mockdata = np.concatenate((aCpool,aMeth,aCO2),axis = None)


