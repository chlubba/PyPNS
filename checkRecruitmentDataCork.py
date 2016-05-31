import numpy as np
import matplotlib.pylab as plt

import testWaveletDenoising as w

data = np.loadtxt('/home/carl/Dropbox/_Exchange/Project/SD_1ms_AllCurrents.txt')

# print np.shape(data)

# denoise data
denoisedVoltage = w.wden(data[:,1], level=12, threshold=8)

plt.plot(data[:,0], data[:,1])
plt.plot(data[:,0], denoisedVoltage, color=[1,0,0], linewidth=2)
plt.show()