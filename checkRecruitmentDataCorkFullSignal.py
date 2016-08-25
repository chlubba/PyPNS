import numpy as np
import matplotlib.pylab as plt

import testWaveletDenoising as w

data = np.loadtxt('/home/carl/Dropbox/_Exchange/Project/SD_1ms_AllCurrents.txt')

# print np.shape(data)

# denoise data
denoisedVoltage = w.wden(data[:,1], level=12, threshold=1.5)

distance = 70.
tStart = 3.0245 # 3.024 #

time = data[:,0]
tCut = (time[time>tStart] - tStart)*1000

# CVs = distance/tCut

vDenCut = denoisedVoltage[time>tStart]

# plt.plot(data[:,0], data[:,1], color=[0.6,0.6,0.6])
# plt.plot(data[:,0], denoisedVoltage, color=[1,0,0], linewidth=2)
f, ax = plt.subplots()
plt.plot(tCut, vDenCut)

timesTicks = np.arange(5,int(max(tCut)), 5)
tickLabelStrings = []
for tickTime in timesTicks:
    tickLabelStrings.append('%2.3f' % (distance / tickTime))

# ax.set_xticks(range(0,int(max(tCut)), 20))
# ax.set_xticklabels(tickLabelStrings)

plt.xticks(timesTicks, rotation='vertical')
ax.set_xticklabels(tickLabelStrings)
plt.xlabel('conduction velocity [m/s]')
plt.grid()
plt.show()