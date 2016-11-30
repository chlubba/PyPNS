import numpy as np
import matplotlib.pylab as plt
import PyPN

import testWaveletDenoising as w

data = np.loadtxt('/media/carl/18D40D77D40D5900/Dropbox/_Exchange/Project/SD_1ms_AllCurrents.txt')

# print np.shape(data)

# denoise data
denoisedVoltage = w.wden(data[:,1], level=12, threshold=2)

# get bundle and then CAP
bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00001')

f, (ax1, ax2) = plt.subplots(2,1)

time, CAP = bundle.get_CAP_from_file(3)
ax1.plot(time, CAP[-1,:]*1000/2*0.3, label='homogeneous analytic', color=[0.4,0.7,0.4]) # , color=[0.3, 0.3, 1]) # [0.5,0.75,0.5]
time, CAP = bundle.get_CAP_from_file(2)
ax1.plot(time, CAP[-1,:]*1000, label='inhomogeneous FEM', linewidth=2, color='b') # , color=[0.1, 0.6, 0.1])
ax1.set_xlabel('time [ms]')
ax1.set_ylabel('voltage [$\mu$V]')
ax1.set_title('Simulated Recording')
ax1.legend()

ax2.plot(data[:,0], data[:,1], color=[0.6, 0.6, 0.6], label='raw recording')
ax2.plot(data[:,0], denoisedVoltage, color=[1,0,0], linewidth=2, label='wavelet denoised')
ax2.set_xlabel('time [s]')
ax2.set_ylabel('voltage [$\mu$V]')
ax2.set_title('Recording from Rat Vagus Nerve')
ax2.legend()

plt.tight_layout()
plt.show()