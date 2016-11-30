import numpy as np
import matplotlib.pylab as plt

import testWaveletDenoising as w

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

data = np.loadtxt('/media/carl/18D40D77D40D5900/Dropbox/_Exchange/Project/SD_1ms_AllCurrents.txt')

# print np.shape(data)

# denoise data
denoisedVoltage = w.wden(data[:,1], level=12, threshold=1.5)

distance = 70.
tStart = 3.0245 # 3.024 #

time = data[:,0]
tCut = (time[time>tStart] - tStart)*1000

# CVs = distance/tCut

vDenCut = denoisedVoltage[time>tStart]

vDenCutNoArt = vDenCut
vDenCutNoArt[:100] = np.zeros(100)
vDenCutAbsSmoo = smooth(np.abs(vDenCut), 50)

(f, axarr) = plt.subplots(2,1) # , sharex=True)
# plt.plot(tCut, vDenCutAbsSmoo)
# timesTicks = np.concatenate((np.arange(5, 30, 5), np.arange(26,int(max(tCut)), 10)))
timesTicks = np.arange(1 ,int(max(tCut)), 5)
plt.xticks(timesTicks, rotation=60, ha='right')
tickLabelStrings = []
for tickTime in timesTicks:
    if tickTime < 27:
        tickLabelStrings.append('%2.2f (%2.2f)' % (distance / tickTime, distance / tickTime / 5))
    else:
        tickLabelStrings.append('%2.2f (%2.2f)' % (distance / tickTime, (distance / tickTime / 2) ** 2))
plt.gca().set_xticklabels(tickLabelStrings)
axarr[1].set_xlabel('conduction velocity [m/s] (diameter [$\mu$m])')
# plt.show()
axarr[1].grid()

# yLims = plt.gca().axes.get_ylim()

import matplotlib.patches as patches
axarr[1].add_patch(
    patches.Rectangle(
        (distance/15, 0),   # (x,y)
        distance/3,          # width
        14,          # height
        alpha=0.1
    )
)

axarr[1].add_patch(
    patches.Rectangle(
        (distance/2.5, 0),   # (x,y)
        180 - distance/2.5,          # width
        14,          # height
        alpha=0.1,
        facecolor='g'
    )
)

fontB = {'family': 'serif',
        'color': 'b',
        'weight': 'normal',
        'size': 16,
        }

fontG = {'family': 'serif',
         'color': 'darkgreen',
         'weight': 'normal',
         'size': 16,
         }

plt.text(distance/15 + 1, 12, 'A$\delta$, B', fontdict=fontB)
plt.text(distance/2 + 1, 12, 'C', fontdict=fontG)

axarr[1].plot(tCut, vDenCutAbsSmoo, color='k')
axarr[0].set_ylabel('smoothed absolute value')
axarr[1].set_title('Intensity over conduction velocity and diameters')
# plt.show()

#
#
# tStep = 1
# tInSteps = np.arange(5, max(tCut), tStep)
# condVels = distance/(tInSteps[1:] + tStep/2)
# amplitudeAtCVs = []
# for i in range(1,len(tInSteps)):
#     amplitudeAtCVs.append(np.sum(vDenCutAbsSmoo[np.logical_and(tCut>tInSteps[i-1], tCut<tInSteps[i])]))
#
# plt.figure()
# plt.plot(condVels, amplitudeAtCVs)
# plt.show()
#
# # plt.plot(data[:,0], data[:,1], color=[0.6,0.6,0.6])
# # plt.plot(data[:,0], denoisedVoltage, color=[1,0,0], linewidth=2)
# f, ax = plt.subplots()
axarr[0].plot(tCut, vDenCut)
axarr[0].grid()
axarr[0].set_xlabel('time [ms]')
axarr[0].set_ylabel('$V_{ext}$ [$\mu$V]')
axarr[0].set_title('Recorded CAP from Rat Vagus')

#
# timesTicks = np.arange(5,int(max(tCut)), 5)
# tickLabelStrings = []
# for tickTime in timesTicks:
#     tickLabelStrings.append('%2.3f' % (distance / tickTime))

# ax.set_xticks(range(0,int(max(tCut)), 20))
# ax.set_xticklabels(tickLabelStrings)
#
# plt.xticks(timesTicks, rotation='vertical')
# ax.set_xticklabels(tickLabelStrings)
# plt.xlabel('conduction velocity [m/s]')
# plt.grid()
plt.show()