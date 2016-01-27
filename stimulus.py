import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size(12)
import math

tstop = 3e1 #ms
dt = 0.0005  #ms
dur = 1e1 #ms
#cut_off = 0.995
freq = 0.1 #kHz
amp = 1 #nA
duty_cycle = 0.004
cut_off = math.sin((-duty_cycle+1.0/2)*math.pi)
t = np.linspace(0, dur, (tstop+2*dt)/dt, endpoint=True)
signal = amp*(cut_off < np.sin(2*math.pi*freq*t))-amp*(cut_off < np.sin(2*math.pi*freq*(t+duty_cycle/freq)))
sin = np.sin(2*math.pi*freq*t)
duty = float(np.count_nonzero(signal))/len(signal)
print 'Duty:' +str(duty)
plt.figure()
plt.plot(t,signal, 'r-', markersize=1)
plt.plot(t, sin, '-b', markersize =1)

 ### FFT part ###
"""plt.figure()
T = dt
N = int(float(dur)/dt)
t = np.linspace(0.0, T*N,N)
xf = np.linspace(0.0,1.0/(2.0*T),N/2)
freq = np.fft.fftfreq(N, d=dt)
yf  = np.fft.fft(signal)


plt.plot(freq[0:N/2], 2.0/N * np.abs(yf[0:N/2]))
filtered_yf = yf.copy()
filtered_yf [(xf < 0.1)] =0
filtered_yf [(xf > 30)] =0
plt.figure()
plt.plot(t, np.fft.ifft(yf))
plt.grid()
"""
plt.show()
