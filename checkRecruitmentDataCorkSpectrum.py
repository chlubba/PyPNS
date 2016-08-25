import numpy as np
import numpy.fft
import matplotlib.pylab as plt

from scipy import signal

import testWaveletDenoising as w

data = np.loadtxt('/home/carl/Dropbox/_Exchange/Project/SD_1ms_AllCurrents.txt')

# print np.shape(data)

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


t = data[:,0]
timeStep = t[1] - t[0]
Fs = 1/timeStep

windowTime = np.array([3.035, 3.18])
windowSampleIndices = windowTime*Fs
windowSampleIndices= windowSampleIndices.astype(int)

tWindowed = t[windowSampleIndices[0]:windowSampleIndices[1]]
sWindowed = data[windowSampleIndices[0]:windowSampleIndices[1],1]

sWindowedNorm = sWindowed/(max(sWindowed) - min(sWindowed))

# denoised data
denoisedVoltage = w.wden(data[:,1], level=12, threshold=1.5)

sWindowedDenoised = denoisedVoltage[windowSampleIndices[0]:windowSampleIndices[1]]
sWindowedDenoisedNorm = sWindowedDenoised/(max(sWindowedDenoised) - min(sWindowedDenoised))

# # load simulated CAP
# CAP = np.loadtxt('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00001/CAP_FEMRecCuff2D_recMech2/CAP_FEMRecCuff2D_recMech2.dat')
# print CAP.shape
#
# timeSim = CAP[:,0]
# signalSim = CAP[:,3]
# timeStepSim = (timeSim[1] - timeSim[0])*10**(-3)
#
# signalSimNorm = signalSim/(max(signalSim) - min(signalSim))

# load Siwoo's data
CAPSiwoo = np.loadtxt('/home/carl/PNPy/bundleFromSiwoo/bundle00000/CAP_RecordingMechanism_recMech1/CAP_RecordingMechanism_recMech1.dat')
timeSimSiwoo = CAPSiwoo[:,0]
signalSimSiwoo = CAPSiwoo[:,1]
timeStepSimSiwoo = (timeSimSiwoo[1] - timeSimSiwoo[0])*10**(-3)

# plt.figure()
# # plt.plot(timeSimSiwoo*10**-3, signalSimSiwoo/max(abs(signalSimSiwoo)))
# plt.plot(timeSim*10**-3, signalSim/max(abs(signalSim)))
# plt.plot(tWindowed - min(tWindowed), sWindowed/max(abs(sWindowed)))
# plt.show()


# plt.figure()
# plt.plot(tWindowed - min(tWindowed), sWindowed/max(abs(sWindowed)))
# plt.plot(timeSim*0.001, signalSim/max(abs(signalSim)))
# plt.show()

# plt.figure()
#
# plt.subplot(4,1,1)
# Pxx, freqs, bins, im = plt.specgram(sWindowed, Fs=Fs) # NFFT=NFFT, noverlap=900,
#
# plt.subplot(4,1,2)
# plt.plot(tWindowed - min(tWindowed), sWindowed)
#
# plt.subplot(4,1,3)
# Pxx, freqs, bins, im = plt.specgram(signalSim, Fs=1/timeStepSim) # NFFT=NFFT, noverlap=900,
#
# plt.subplot(4,1,4)
# plt.plot(timeSim, CAP[:,3])
# # plt.show()

plt.figure(),
sp = np.fft.fft(sWindowed)
spN = np.fft.fft(sWindowedNorm)
freq = np.fft.fftfreq(tWindowed.shape[-1], d=timeStep)

sortedInds = sorted(range(len(freq)), key=lambda k: freq[k])
absSp = np.abs(sp)
absSp = absSp[sortedInds]
freq = freq[sortedInds]
plt.semilogy(freq, smooth(absSp,100)/max(absSp), label='experimental recording')
# plt.plot(freq, abs(spN))
# plt.xlim([0,max(freq)])

# sp = np.fft.fft(sWindowedDenoised)
# spN = np.fft.fft(sWindowedDenoisedNorm)
# freq = np.fft.fftfreq(tWindowed.shape[-1], d=timeStep)
# plt.plot(freq, np.abs(sp)/max(np.abs(sp)))
# # plt.plot(freq, abs(spN))

# spSim = np.fft.fft(signalSim)
# spSimN = np.fft.fft(signalSimNorm)
# freqSim = np.fft.fftfreq(timeSim.shape[-1], d=timeStepSim)
# plt.plot(freqSim, np.abs(spSim)/max(np.abs(spSim)))
# # plt.plot(freqSim, abs(spSimN))

spSim = np.fft.fft(signalSimSiwoo)
# spSimN = np.fft.fft(signalSimNorm)
freqSim = np.fft.fftfreq(timeSimSiwoo.shape[-1], d=timeStepSimSiwoo)

sortedInds = sorted(range(len(freqSim)), key=lambda k: freqSim[k])
absSpSim = np.abs(spSim)
absSpSim = absSpSim[sortedInds]
freqSim = freqSim[sortedInds]

plt.semilogy(freqSim, smooth(absSpSim,100)/max(absSpSim), label='simulation')
plt.xlim([0,max(freq)]) # max(max(freqSim), max(freq))])
plt.grid()
plt.legend()
plt.title('FFT of simulated and recorded CAP')
plt.xlabel('frequency [Hz]')
plt.ylabel('smoothed Fourier coefficients of voltage, normed to max=1')
plt.show()

# plt.plot(timeSim, CAP[:,1])
# plt.plot(timeSim, CAP[:,2])
# plt.plot(timeSim, CAP[:,3])
# plt.show()


#
# plt.plot(data[:,0], data[:,1], color=[0.6,0.6,0.6])
# plt.plot(data[:,0], denoisedVoltage, color=[1,0,0], linewidth=2)
# plt.show()