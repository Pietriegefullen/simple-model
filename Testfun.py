# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 11:52:33 2021

@author: Lara
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def testfun2(y,x):
    k1 = .2
    k2 = .01
    return [k1, k2*y[0]]

def testfun(y,x):
    k = 0.2
    q = k #* y
    print(x)
    return q

#testfun(5)

xs = np.linspace(0,100,20)
y  = 1.0  # the initial condition
ys = odeint(testfun2, [0,0], xs)

ys = np.array(ys).flatten()

# Plot the numerical solution
plt.rcParams.update({'font.size': 14})  # increase the font size
plt.xlabel("x")
plt.ylabel("y")
plt.plot(xs, ys)

#%%

a =[1,2,3,4]
erstens,zweitens, drittens, viertens = a
print(erstens)