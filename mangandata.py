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

# Calculate brancing ratios
seven1  = 0
sevens1 = 0
six1    = sixone[0]
sixs1   = np.sqrt(np.diag(sixone0))[0]
five1   = 0
fives1  = 0
four1   = fourone[0] + fourzero[0]
fours1  = np.sqrt(np.diag(fourone0))[0] + np.sqrt(np.diag(fourzero0))[0]
three1  = threeone[0] + threezero[0]
threes1 = np.sqrt(np.diag(threeone0))[0] + np.sqrt(np.diag(threezero0))[0]
two1    = 0
twos1   = 0
one1    = onezero[0] - (sixone[0]+fourone[0]+threeone[0])
ones1   = np.sqrt(np.diag(onezero0))[0]+np.sqrt(np.diag(sixone0))[0]+np.sqrt(np.diag(fourone0))[0]+np.sqrt(np.diag(threeone0))[0]
total   = one1 + two1 + three1 + four1 + five1 + six1 + seven1
totals  = ones1 + twos1 + threes1 + fours1 + fives1 + sixs1 + sevens1

# Calculate percentages
seven2  = 100*seven1/total
sevens2 = np.sqrt((100**2/total**2)*sevens1**2 + (100**2*seven1**2/total**4)*totals**2)
six2    = 100*six1/total
sixs2   = np.sqrt((100**2/total**2)*sixs1**2 + (100**2*six1**2/total**4)*totals**2)
five2   = 100*five1/total
fives2  = np.sqrt((100**2/total**2)*fives1**2 + (100**2*five1**2/total**4)*totals**2)
four2   = 100*four1/total
fours2  = np.sqrt((100**2/total**2)*fours1**2 + (100**2*four1**2/total**4)*totals**2)
three2  = 100*three1/total
threes2 = np.sqrt((100**2/total**2)*threes1**2 + (100**2*three1**2/total**4)*totals**2)
two2    = 100*two1/total
twos2   = np.sqrt((100**2/total**2)*twos1**2 + (100**2*two1**2/total**4)*totals**2)
one2    = 100*one1/total
ones2   = np.sqrt((100**2/total**2)*ones1**2 + (100**2*one1**2/total**4)*totals**2)


print('Branching ratios:')
print('First level: ' + str(one2) + '+/-' + str(ones2) + '%')
print('Second level: ' + str(two2) + '+/-' + str(twos2) + '%')
print('Third level: ' + str(three2) + '+/-' + str(threes2) + '%')
print('Fourth level: ' + str(four2) + '+/-' + str(fours2) + '%')
print('Fifth level: ' + str(five2) + '+/-' + str(fives2) + '%')
print('Sixth level: ' + str(six2) + '+/-' + str(sixs2) + '%')
print('Seventh level: ' + str(seven2) + '+/-' + str(sevens2) + '%')

# Calculate partial decay rates:
decay41 = fourone[0]/(fourone[0]+fourzero[0])
decay40 = fourzero[0]/(fourone[0]+fourzero[0])
decay31 = threeone[0]/(threeone[0]+threezero[0])
decay30 = threezero[0]/(threeone[0]+threezero[0])