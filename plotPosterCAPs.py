import numpy as np
import matplotlib.pyplot as plt
import os

filenames = ['homogeneous2.npy', 'noCuff.npy', 'Cuff.npy']
legends = ['Homogeneous medium', 'No cuff', 'Cuff']

startIndex = 500
endIndex = 3500

scaling = True


def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size / 2:]


CAPs = []
for filename in filenames:
    # plt.figure()
    [t, CAP] = np.load(os.path.join('/media/carl/4ECC-1C44/PyPN/FEM_CAPs/forPoster', filename))
    CAPs.append(np.squeeze(CAP))

# plt.plot(t[startIndex:], np.squeeze(CAP)[startIndex:])
# print len(CAPs[0])
# print np.argmax(np.correlate(CAPs[0][startIndex:], CAPs[1][startIndex:], mode='full')) - len(CAPs[0][startIndex:])

meanCAPs = []
for CAP in CAPs:
    meanCAPs.append(np.mean(CAP[startIndex-200:startIndex]))

voltageUnitScaling = 1000000 # V to uV
for counter, CAP in  enumerate(CAPs):
    plt.plot(t[startIndex:endIndex], (np.squeeze(CAP)[startIndex:endIndex]-meanCAPs[counter])*voltageUnitScaling, label=legends[counter])
    plt.xlabel('time [ms]')
    plt.ylabel('voltage [$\mu$V]')

plt.legend(loc=3)
# plt.legend()
# plt.tight_layout()

if scaling:
    CAPsN = []
    for counter, CAP in enumerate(CAPs):
        maxCAP = np.max(CAP[startIndex:endIndex])
        minCAP = np.min(CAP[startIndex:endIndex])
        CAPsN.append((CAP - minCAP)/(maxCAP - minCAP))

    meanCAPsN = []
    for CAP in CAPsN:
        meanCAPsN.append(np.mean(CAP[startIndex - 200:startIndex]))

    plt.figure()
    for counter, CAP in enumerate(CAPsN):
        plt.plot(t[startIndex:endIndex], np.squeeze(CAP)[startIndex:endIndex] - meanCAPsN[counter], label=legends[counter])
        plt.xlabel('time [ms]')
        plt.ylabel('normalized voltage')
        plt.legend(loc=3)



plt.show()