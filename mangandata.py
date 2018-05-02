# -*- coding: utf-8 -*-
"""
Created on Wed May  2 11:09:02 2018

@author: Jesper
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import energy_calibration as ec
import peakutils

plt.close('all')

datach, datacounts = ec.convert_list_to_hist('../eksp1_data/mn56_14400_ch001.txt')
backch, backcounts = ec.convert_list_to_hist('data_180417/cali_baggrund_live_1800_ch001.txt')
cal_coeffs, cal_coeffs_err = ec.get_energy_calibration()
a = cal_coeffs['ch001'][0]
b = cal_coeffs['ch001'][1]
datae = a * datach + b
backe = a * backch + b
realcounts = datacounts/14400-backcounts/1800
fitpar = np.linspace(1, max(datae), 10000)
# Names are the electron transitions, onezero is the transition from orbit 1 to orbit 0.
onezero, onezero0 = curve_fit(ec.gauss, datae, realcounts, bounds=([5, 840, 0],[6, 850, 5]))
threeone, threeone0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 1805, 0],[0.6, 1815, 5]))
fourone, fourone0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 2110, 0],[0.3, 2120, 5]))


plt.figure()
plt.plot(datae, realcounts, 'k-')
plt.plot(fitpar, ec.gauss(fitpar, *onezero), '--')
plt.plot(fitpar, ec.gauss(fitpar, *threeone), '--')
plt.plot(fitpar, ec.gauss(fitpar, *fourone), '--')
