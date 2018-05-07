import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from energy_calibration import get_energy_calibration

def data_list_to_hist(data):
    max_channel = 2500
    ch = np.linspace(0, max_channel, max_channel)
    counts = np.zeros(max_channel)

    for i in range(0, len(data[:, 0])):
        if(data[i, 1] >= max_channel or data[i, 1] == 0):
            continue

        counts[int(data[i, 1])] += 1

    return ch, counts

path = '../Eksp1_HUGE_data/mn56_14400_ch000.txt'
bg_path = 'data_180417/cali_baggrund_live_1800_ch000.txt'
data = np.array(pd.read_csv(path, sep=' ', skiprows=5)) #Indlæs datafilen
bg_data = np.array(pd.read_csv(bg_path, sep=' ', skiprows=5)) #Indlæs datafilen

def eff_func(E):
    return np.exp(-1.31386054*10**(-3) * E - 4.77730422)

cal_coeffs, cal_err = get_energy_calibration()

#Only read first second
i = 0
time_cut = 1e8

while(data[i, 0] < time_cut):
    i += 1
one_sec_index = i

i = 0
while(bg_data[i, 0] < time_cut):
    i += 1
one_sec_index_bg = i

ch, counts = data_list_to_hist(data[0:one_sec_index, :])
ch_bg, counts_bg = data_list_to_hist(bg_data[0:one_sec_index_bg, :])
E = cal_coeffs['ch000'][0] * ch + cal_coeffs['ch000'][1]

counts /= eff_func(E)
#counts_bg /= eff_func(E)

print("Counts pr. first second", np.trapz(counts) - np.trapz(counts_bg))

plt.figure()
plt.plot(E, counts - counts_bg)
plt.show()
