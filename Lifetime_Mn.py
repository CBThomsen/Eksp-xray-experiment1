
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 11:15:54 2018

@author: lasse
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
from energy_calibration import get_energy_calibration
from energy_calibration import convert_list_to_hist

#path = '/home/lasse/Documents/mn56_14400_ch000.txt'
path = '../Eksp1_HUGE_data/mn56_14400_ch000.txt'
bg_path = 'data_180417/cali_baggrund_live_1800_ch000.txt'
data = np.array(pd.read_csv(path, sep=' ', skiprows=5)) #Indlæs datafilen
bg_data = np.array(pd.read_csv(bg_path, sep=' ', skiprows=5)) #Indlæs datafilen

def eff_func(E):
    return np.exp(-1.31386054*10**(-3) * E - 4.77730422)

cal_coeffs, cal_err = get_energy_calibration()

def find_counts_per_second(data, enable_eff):
    #Første søjle er tid i 10⁻8 sekunder
    # anden søjle er counts
    dummy   = np.shape(data)

    j = 0
    counts_pr_sec = np.array([])
    time_step = 1e8
    b = time_step

    #for i in np.arange(100000): #1e5 kan køres relativt hurtigt. Meget mere, og så tager scriptet meget lang tid at køre
    for i in np.arange(dummy[0]):
        if data[i,0]>=b: #Når tiden summet op giver 1 sekund, lægges alle counts dertil sammen
            b = b + time_step

            counts = 0
            for count_index in range(j, i):
                if(enable_eff == True):
                    E = cal_coeffs['ch000'][0] * data[count_index, 1] + cal_coeffs['ch000'][1]
                    counts += 1 / eff_func(E)
                else:
                    counts += 1

            counts_pr_sec = np.append(counts_pr_sec,counts) #Lægger oveni array
            j = i #Omdefinerer j, så vi ikke tæller de samme to gange

    return counts_pr_sec

counts_pr_sec = find_counts_per_second(data, False)
tid = np.linspace(0, np.size(counts_pr_sec), np.size(counts_pr_sec))
bg_counts_pr_sec = find_counts_per_second(bg_data, False)
bg_counts_pr_sec = bg_counts_pr_sec[0:1800]

bg_tid = np.linspace(0, np.size(bg_counts_pr_sec), np.size(bg_counts_pr_sec))
bg_mean_activity = np.mean(bg_counts_pr_sec)

def activity_function(x,A,lam):
    return A*np.exp(-lam*x)


counts_pr_sec -= bg_mean_activity
counts_pr_sec = np.nan_to_num(counts_pr_sec)
popt,pcov = sp.curve_fit(activity_function, tid, counts_pr_sec, bounds=([0,0],[1e20,0.1]))

nydum2 = np.sqrt(np.diag(pcov))
sigma_A = nydum2[0]
sigma_lam = nydum2[1]
print('A =', popt[0], '+/-', sigma_A ,'   ','lambda =', popt[1], '+/-', sigma_lam)

plt.figure()
plt.plot(tid, counts_pr_sec, label='Data points')
plt.plot(tid, activity_function(tid,*popt), label='Fitted curve')
plt.plot(tid, activity_function(tid, popt[0], 7.46722 * 10**(-5)), label='Theoretical curve')
plt.xlabel('Time [s]')
plt.ylabel('Activity [Bq]')
plt.legend()
plt.yscale('log')

plt.show()