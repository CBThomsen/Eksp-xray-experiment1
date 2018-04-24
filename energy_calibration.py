import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import peakutils

def gauss(x, a, x0,  sigma):
    return a*np.exp(-(x - x0)**2 / (2*sigma**2))

def lin_func(x, a, b):
    return x * a + b

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

def find_peaks(path, plot=False):
    ch, counts = convert_list_to_hist(path)
    counts[ch < 500] = 0
    indexes = peakutils.indexes(counts, thres=0.18, min_dist=50)

    peaks = []
    peaks_error = []

    if(plot == True):
        plt.figure()
        plt.title(path)
        plt.plot(ch, counts, label="Data")

    for peak_index in indexes:
        popt, pcov = curve_fit(gauss, ch[peak_index - 150:peak_index + 150], counts[peak_index - 150:peak_index + 150], bounds = [[0, peak_index-50, 0], [100000, peak_index+50, 50]])
        if(plot == True):
            plt.plot(ch[peak_index - 150:peak_index + 150], gauss(ch[peak_index - 150:peak_index + 150], *popt), label="Gauss fit")
            plt.legend()

        peaks.append(popt[1])
        peaks_error.append(np.sqrt(np.diag(pcov))[1])

    return peaks, peaks_error

def find_all_peaks():
    peaks_array = {'ch000': [], 'ch001':[], 'ch002':[]}
    peaks_error_array = {'ch000': [], 'ch001':[], 'ch002':[]}

    for mat in ['co60', 'na22', 'cs137']:
        for ch in ['ch000', 'ch001', 'ch002']:
            peaks, peaks_error = find_peaks('data_180417/cali_' + mat + '_' + ch + '.txt', False)

            for p in peaks:
                peaks_array[ch].append(p)
            for p_err in peaks_error:
                peaks_error_array[ch].append(p_err)

    return peaks_array

def get_energy_calibration():
    energy = [1173.2, 1332.5, 511, 1274.6, 661.7] # co60, na22, cs137
    peaks_array = find_all_peaks()
    channels = np.linspace(0, 3000, 3000)

    cal_coeffs = {'ch000': 0, 'ch001': 0, 'ch002': 0}
    cal_coeffs_err = {'ch000': 0, 'ch001': 0, 'ch002': 0}

    for ch in ['ch000', 'ch001', 'ch002']:
        plt.figure()

        popt, pcov = curve_fit(lin_func, peaks_array[ch], energy)
        plt.plot(channels, lin_func(channels, *popt), 'b-', label='Fit')
        plt.plot(peaks_array[ch], energy, 'rx', label='Data')

        plt.title('Energy calibration for channel ' + ch)
        plt.xlabel('Channel')
        plt.ylabel('Energy [keV]')
        plt.legend()

        plt.savefig('figurer/Energy calibration for channel ' + ch + '.pdf')

        cal_coeffs[ch] = popt
        cal_coeffs_err[ch] = np.sqrt(np.diag(pcov))

    return cal_coeffs, cal_coeffs_err

if __name__ == "__main__":
    cal_coeffs, cal_coeffs_err = get_energy_calibration()
    for ch in ['ch000', 'ch001', 'ch002']:
        print("Calibration result for channel", ch,":", cal_coeffs[ch], "error", cal_coeffs_err[ch])

    plt.show()