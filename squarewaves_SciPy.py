from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

timeRes = 0.005 # ms
stimDur  = 10 # ms
frequency = 1 # kHz
amplitude = 3.
dutyCycle = 0.05
delay = 5

t = np.arange(0, stimDur, timeRes)
stimulusSignal = amplitude * signal.square(2 * np.pi * frequency * t, duty=dutyCycle)

stimulusSignalDelayed = np.concatenate((np.zeros(delay/timeRes), stimulusSignal))
# tDelayed = np.concatenate((np.zeros(delay/timeRes), t))
# stimulusSignalMono = stimulusSignal/2 + amplitude/2

# t = np.linspace(0, 1, 500, endpoint=False)
plt.plot(stimulusSignalDelayed)
# plt.plot(t, stimulusSignalMono)
plt.show()