# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 12:36:03 2018

@author: peter
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import peakutils
from energy_calibration import get_energy_calibration
plt.close('all')

def gauss(x, a, x0,  sigma):
    return a*np.exp(-(x - x0)**2 / (2*sigma**2))

def convert_list_to_hist(path, plot=False):
    max_channel = 5000
    ch = np.linspace(0, max_channel, max_channel)
    data = np.array(pd.read_csv(path, sep=' ', skiprows=5))
    counts = np.zeros(max_channel)

    for i in range(0, len(data[:, 0])):
        if(data[i, 1] >= max_channel or data[i, 1] == 0):
            continue

        counts[int(data[i, 1])] += 1

    if(plot == True):
        plt.figure()
        plt.title(path)
        plt.plot(ch, counts, label="Data")

    return ch, counts

def find_peaks(ch, counts, plot=False):
    indexes = peakutils.indexes(counts, thres=0.085, min_dist=70)

    peaks = []
    sigma = []
    a = []
    peaks_error = []

    if(plot == True):
        plt.figure()
        plt.plot(ch, counts, label="Data")

    for peak_index in indexes:
        popt, pcov = curve_fit(gauss, ch, counts, bounds = [[0, peak_index-10, 0], [5000, peak_index+10, 4]])
        if(plot == True):
            plt.plot(ch, gauss(ch, *popt), label="Gauss fit")
            plt.legend()
            
        a.append(popt[0])
        peaks.append(popt[1])
        sigma.append(popt[2])
        peaks_error.append(np.sqrt(np.diag(pcov))[1])

    return peaks, peaks_error, sigma, a

def fitGauss(path, plot=False):
    ch, counts = convert_list_to_hist(path)

    if(plot == True):
        plt.figure()
        plt.plot(ch[0:500], counts[0:500], label="Data")

    popt, pcov = curve_fit(gauss, ch[0:240], counts[0:240], bounds = [[0, 100, 0], [50000, 240, 50]])
    if(plot == True):
        plt.plot(ch[0:500], gauss(ch[0:500], *popt), label="Gauss fit")
        plt.legend()

    return ch, counts, popt, pcov

def fractionToTotal(a, sigma):
    
    total_activity = 4.8*3.7*10**10
    # Bq
    theo_activities = np.array([0.0743, 0.193, 0.376, 0.461, 0.0494, 0.0303, 0.151, 0.0579, 0.154])*0.9998*total_activity
    
    theo_activity_errors = 0.01*theo_activities
    
    activity_errors = []
    
    activities = []
    for i in range(1,4):
        activities.append(np.sqrt(2)*a[i]*abs(sigma[i])*np.sqrt(np.pi)/1800)
        activity_errors.append(np.sqrt(np.sqrt(2)*a[i]*abs(sigma[i])*np.sqrt(np.pi)/1800))
    for i in range(5,11):
        activities.append(np.sqrt(2)*a[i]*abs(sigma[i])*np.sqrt(np.pi)/1800)
        activity_errors.append(np.sqrt(np.sqrt(2)*a[i]*abs(sigma[i])*np.sqrt(np.pi)/1800))
    
    efficiency_errors = []    
    efficiencies = []
    
    for i in range(0,9):
        efficiencies.append(activities[i]/theo_activities[i])
        efficiency_errors.append(np.sqrt((1/theo_activities[i]*activity_errors[i])**2+(activities[i]/(theo_activities[i]**2))**2*theo_activity_errors[i]**2))
        
    
    
    
    
    return efficiencies, efficiency_errors
    
    
    


ch_ra, counts_ra_mes, popt, pcov = fitGauss('../Eksp1_HUGE_data/ra226_3600_ch000.txt')
ch_bg, counts_bg, popt_bg, pcov_bg = fitGauss('data_180417/cali_baggrund_live_1800_ch000.txt')

counts_ra = counts_ra_mes * popt_bg[0]/popt[0] - counts_bg

plt.figure()
plt.plot(ch_ra, counts_ra)

peaks, peaks_error, sigma, a = find_peaks(ch_ra, counts_ra)
peaks = np.array(peaks)
peaks_error = np.array(peaks_error)

cal_coeffs, cal_coeffs_error = get_energy_calibration()

energy_peaks = cal_coeffs['ch000'][0]*peaks+cal_coeffs['ch000'][0]
energy_errors = cal_coeffs['ch000'][0]*peaks_error+cal_coeffs['ch000'][0]

efficiencies, efficiency_errors = fractionToTotal(a, sigma)

energy_peaks_ny = []
for i in energy_peaks[1:4]:
    energy_peaks_ny.append(i)
for i in energy_peaks[5:11]:
    energy_peaks_ny.append(i)

plt.figure()
plt.plot(energy_peaks_ny, efficiencies, 'rx')


#plt.figure()
#plt.plot(ch_ra, counts_ra_mes)
#
#plt.figure()
#plt.plot(ch_bg, counts_bg)









