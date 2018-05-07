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

datach, datacounts = ec.convert_list_to_hist('../eksp1_data/mn56_14400_ch000.txt')
backch, backcounts = ec.convert_list_to_hist('data_180417/cali_baggrund_live_1800_ch000.txt')
cal_coeffs, cal_coeffs_err = ec.get_energy_calibration()
a = cal_coeffs['ch000'][0]
b = cal_coeffs['ch000'][1]
datae = a * datach + b
backe = a * backch + b

# Define efficiency calibration function
epsilon = np.exp(-1.31385886*10**-3 * datae - 4.77730422)
realcounts = (datacounts/14400-backcounts/1800)/epsilon
fitpar = np.linspace(1, max(datae), 10000)

# Names are the electron transitions, onezero is the transition from orbit 1 to orbit 0.
onezero, onezero0 = curve_fit(ec.gauss, datae, realcounts, bounds=([1500, 842, 0],[2000, 852, 5]))
#fivetwo, fivetwo0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 1033, 0],[0.6, 1043, 5]))
#twoone, twoone0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 1233, 0],[0.6, 1243, 5]))
threeone, threeone0 = curve_fit(ec.gauss, datae, realcounts, bounds=([500, 1806, 0],[1000, 1816, 5]))
fourone, fourone0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 2108, 0],[700, 2118, 5]))
sixone, sixone0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 2518, 0],[100, 2528, 5]))
#sevenone, sevenone0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 2593, 0],[0.3*10**10, 2603, 5]))
threezero, threezero0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 2653, 0],[100, 2663, 5]))
fourzero, fourzero0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 2955, 0],[100, 2965, 5]))
#sixzero, sixzero0 = curve_fit(ec.gauss, datae, realcounts, bounds=([0, 3365, 0],[0.3*10**10, 3375, 5]))

plt.figure()
plt.xlabel('E [keV]')
plt.ylabel('Arbitrary unit')
plt.plot(datae, realcounts, 'k-')
plt.plot(fitpar, ec.gauss(fitpar, *onezero), '--')
#plt.plot(fitpar, ec.gauss(fitpar, *fivetwo), '--')
#plt.plot(fitpar, ec.gauss(fitpar, *twoone), '--')
plt.plot(fitpar, ec.gauss(fitpar, *threeone), '--')
plt.plot(fitpar, ec.gauss(fitpar, *fourone), '--')
plt.plot(fitpar, ec.gauss(fitpar, *sixone), '--')
#plt.plot(fitpar, ec.gauss(fitpar, *sevenone), '--')
plt.plot(fitpar, ec.gauss(fitpar, *threezero), '--')
plt.plot(fitpar, ec.gauss(fitpar, *fourzero), '--')
#plt.plot(fitpar, ec.gauss(fitpar, *sixzero), '--')

# Calculate brancing ratios DET ER FORKERT
total   = onezero[0]+threeone[0]+fourone[0]+sixone[0]+threezero[0]+fourzero[0]
totals  = np.sqrt(np.diag(onezero0))[0]+np.sqrt(np.diag(threeone0))[0]+np.sqrt(np.diag(fourone0))[0]+np.sqrt(np.diag(sixone0))[0]+np.sqrt(np.diag(threezero0))[0]+np.sqrt(np.diag(fourzero0))[0]
seven   = 0
sevens  = 0
six     = 100*sixone[0]/total
sixs    = np.sqrt(100**2*np.diag(sixone0)[0]**2/total**2 + 100**2*sixone[0]**2*totals**2/totals**4)
five    = 0
fives   = 0
four    = 100*(fourone[0]+fourzero[0])/total
fours   = np.sqrt(100**2*np.diag(fourone0)[0]**2/total**2 + 100**2*np.diag(fourzero0)[0]**2/total**2 + 100**2*totals**2*(np.diag(fourone0)[0]+np.diag(fourzero0)[0])**2/total**4)
three   = 100*(threeone[0]+threezero[0])/total
threes  = np.sqrt(100**2*np.diag(threeone0)[0]**2/total**2 + 100**2*np.diag(threezero0)[0]**2/total**2 + 100**2*totals**2*(np.diag(threeone0)[0]+np.diag(threezero0)[0])**2/total**4)
two     = 0
twos    = 0
one     = 100*(onezero[0]-(sixone[0]+fourone[0]+threeone[0]))/total
ones    = np.sqrt(np.diag(onezero0))[0]+np.sqrt(np.diag(sixone0))[0]+np.sqrt(np.diag(fourone0))[0]+np.sqrt(np.diag(threeone0))[0]

print('Branching ratios:')
print('First level: ' + str(one) + '+/-' + str(ones) + '%')
print('Second level: ' + str(two) + '+/-' + str(twos) + '%')
print('Third level: ' + str(three) + '+/-' + str(threes) + '%')
print('Fourth level: ' + str(four) + '+/-' + str(fours) + '%')
print('Fifth level: ' + str(five) + '+/-' + str(fives) + '%')
print('Sixth level: ' + str(six) + '+/-' + str(sixs) + '%')
print('Seventh level: ' + str(seven) + '+/-' + str(sevens) + '%')

# Calculate partial decay rates:
decay41 = fourone[0]/(fourone[0]+fourzero[0])