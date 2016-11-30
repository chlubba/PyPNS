import matplotlib.pyplot as plt
import numpy as np

# now generate fitting filter
from scipy import signal

# b = signal.firwin(80, 0.5, window=('kaiser', 8))
# w, h = signal.freqz(1, (1,1))

from scipy.signal import butter, lfilter, freqz

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

fs = 20000

a, b = butter_lowpass(400, fs, order=1)

w, h = signal.freqz(a, b)

plt.plot(w/(2*np.pi)*fs, 20 * np.log10(abs(h)), 'b')
plt.ylabel('Amplitude [dB]')
plt.xlabel('Frequency [rad/sample]')
plt.grid()

plt.show()