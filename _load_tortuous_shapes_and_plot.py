import numpy as np
import matplotlib.pylab as plt
import cPickle as pickle
import os


def maxcor (signal1, signal2):
    return np.max(np.correlate(signal1, signal2, "full"))

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

startIndex = 30

saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/tortuous', 'toruousSignals.dict'), "rb" ))

signalArray = saveDict['signals']
RDCArray = saveDict['RDCs']
radiiArray = saveDict['radius']

RDCs = np.unique(RDCArray)
radii = np.unique(radiiArray)

maxcorMatrix = np.zeros([len(RDCs), len(radii)])

numberOfRuns = len(np.squeeze(np.where(np.logical_and(RDCArray == RDCs[0], radiiArray == radii[0]))))

radiusInd = 2
RDCInd = 1

f, axarr = plt.subplots(len(RDCs), len(radii))

for RDCInd, RDC in enumerate(RDCs):
    for radiusInd, radius in enumerate(radii):

        # radius = radii[radiusInd]
        # RDC = RDCs[RDCInd]

        axis = axarr[RDCInd][radiusInd]

        # plot all signals
        indices = np.squeeze(np.where(np.logical_and(RDCArray == RDC, radiiArray == radius)))
        for index in indices:
            signal1 = np.squeeze(signalArray[index])
            signal1 = signal1[startIndex:]
            signal1 = smooth(signal1,100)

            signal1Norm = signal1/np.linalg.norm(signal1)

            axis.plot(signal1)

        axis.set_title('radius ' + str(radius) + ', RDC ' + str(RDC))
# plt.tight_layout()
plt.show()
