import numpy as np
import matplotlib.pylab as plt
import cPickle as pickle
import os


def maxcor (signal1, signal2):
    return np.max(np.correlate(signal1, signal2, "full"))

startIndex = 30

saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/tortuous', 'toruousSignals.dict'), "rb" ))

signalArray = saveDict['signals']
RDCArray = saveDict['RDCs']
radiiArray = saveDict['radius']

RDCs = np.unique(RDCArray)
radii = np.unique(radiiArray)

amplitudeDiffMatrix = np.zeros([len(RDCs), len(radii)])
np.save('/media/carl/4ECC-1C44/PyPN/tortuous/amplitudeDiff', amplitudeDiffMatrix)

numberOfRuns = len(np.squeeze(np.where(np.logical_and(RDCArray == RDCs[0], radiiArray == radii[0]))))

for RDCInd, RDC in enumerate(RDCs):
    for radiusInd, radius in enumerate(radii):
        indices = np.squeeze(np.where(np.logical_and(RDCArray == RDC, radiiArray == radius)))


        vppDiffs = []
        otherIndices = indices[1:]
        for index in indices:
            print 'index:' + str(index)
            signal1 = np.squeeze(signalArray[index])
            signal1 = signal1[startIndex:]
            # signal1Norm = signal1/np.linalg.norm(signal1)

            for otherIndex in otherIndices:
                print 'otherIndex:' + str(otherIndex)
                signal2 = np.squeeze(signalArray[otherIndex])
                signal2 = signal2[startIndex:]
                # signal2Norm = signal2 / np.linalg.norm(signal2)

                vpp1 = max(signal1) - min(signal1)
                vpp2 = max(signal2) - min(signal2)

                vppDiffs.append(np.abs(vpp1 - vpp2))

            otherIndices  = otherIndices[1:]
            # print otherIndices

        amplitudeDiffMatrix[RDCInd, radiusInd] = np.std(vppDiffs)/np.mean(vppDiffs)

        #     plt.plot(signal)
        # plt.show()

        # print 'hm'

np.save('/media/carl/4ECC-1C44/PyPN/tortuous/amplitudeDiff', amplitudeDiffMatrix)
