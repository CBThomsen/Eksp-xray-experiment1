from energy_calibration import find_peaks
import matplotlib.pyplot as plt

find_peaks(r'./../Eksp1_HUGE_data/mn56_14400_ch000.txt', True)
find_peaks(r'./../Eksp1_HUGE_data/mn56_14400_ch001.txt', True)
find_peaks(r'./../Eksp1_HUGE_data/mn56_14400_ch002.txt', True)

plt.show()