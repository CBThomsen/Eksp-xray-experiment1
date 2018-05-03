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

efficiencies = np.array([3.1651874400744877e-10, 2.702495146322019e-10, 2.2454164406556931e-10, 1.3519008349484144e-10, 1.34693618642749e-10, 1.5056754071442574e-10, 9.056070801898412e-11, 7.886386812318238e-11, 5.541795870551268e-11])

efficiency_errors = np.array([1.5492372170649816e-10, 8.884372850751525e-11, 5.803657422025616e-11, 4.066156378591652e-11, 1.239246334907226e-10, 1.6729536442983763e-10, 5.8124081618633296e-11, 8.758693199566306e-11, 4.5021479115487303e-11])

energies = np.array([242.40493289838415, 295.5934409936591, 352.32950929188223, 609.919154256204, 768.9446623941775, 934.9579062533103, 1120.9168415851823, 1238.725530507722, 1765.103308261181])

energy_errors = np.array([0.91574251, 0.81807544, 0.77723467, 0.81136017, 1.76658178, 2.69299422, 1.34460184, 2.68546058, 1.78084622])

log_efficiencies = np.log(efficiencies)

yerr=efficiency_errors

xerr=energy_errors

def exponentiel(x, a, b):
    return np.exp(-a*x-b)

def fourDegPoly(x, c, d, e):
    return c*x**2+d*x+e

def linear(x, a, b):
    return a*x+b

popt, pcov = curve_fit(linear, energies, log_efficiencies)
poptexp, pcovexp = curve_fit(exponentiel, energies, efficiencies, bounds = [[abs(popt[0])-0.1*abs(popt[0]),abs(popt[1])-0.1*abs(popt[1])], [abs(popt[0])+0.1*abs(popt[0]),abs(popt[1])+0.1*abs(popt[1])]])

plt.figure()
plt.grid()
plt.xlabel('Energy [keV]')
plt.ylabel('Efficiency')
plt.title('Detector 1 Effeciency vs. Energy')
plt.errorbar(energies, efficiencies, yerr, xerr, fmt='rx')
plt.plot(np.linspace(5,2500), exponentiel(np.linspace(5,2500), *poptexp), 'k-')
plt.legend(labels = ['Data', 'Exp. Fit'])

  


