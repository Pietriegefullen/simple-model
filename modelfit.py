# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:45:44 2020

@author: Lara
"""
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from Mockdata import Mockdata, Mockdata_frame
from SimpleOUT import simplefun




xlist = list(range(0,101,1))
xdata = xlist
ydata = Mockdata[:,0]
optimal_parameters , _ = curve_fit(simplefun, xdata , ydata ,p0 = [0.001,0.001])


k_opt = optimal_parameters[0] # abbaugeschwindigkeit (Fressgeschwindigkeit)
h_opt = optimal_parameters[1] # Microbenwachstum

