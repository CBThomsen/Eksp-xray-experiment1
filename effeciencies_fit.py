# -*- coding: utf-8 -*-
"""
Created on Wed May  2 14:41:00 2018

@author: peter
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

plt.close('all')

effeciencies = np.array([3.1651874400744877e-10,
 2.702495146322019e-10,
 2.2454164406556931e-10,
 1.3519008349484144e-10,
 1.34693618642749e-10,
 1.5056754071442574e-10,
 9.056070801898412e-11,
 7.886386812318238e-11,
 5.541795870551268e-11])

energies =np.array([242.40493289838415,
 295.5934409936591,
 352.32950929188223,
 609.919154256204,
 768.9446623941775,
 934.9579062533103,
 1120.9168415851823,
 1238.725530507722,
 1765.103308261181])

def exponentiel(x, a, b):
    return a*np.exp(-b*x)

popt, pcov = curve_fit(exponentiel, energies, efficiencies, bounds = [[0, 0], [4.8*3.7*10**10, 1000000]])

plt.figure()
plt.plot(energies, efficiencies, 'rx')
plt.plot(energies, exponentiel(energies, *popt))

  


