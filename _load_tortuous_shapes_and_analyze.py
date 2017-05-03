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

# saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/tortuous', 'toruousSignals_50_300_moreRuns.dict'), "rb" ))
saveDict = pickle.load(open(os.path.join('tortuous', 'tortuousTest.dict'), "rb" ))

signalArray = saveDict['signals']
RDCArray = saveDict['RDCs']
radiiArray = saveDict['radius']

RDCs = np.unique(RDCArray)
radii = np.unique(radiiArray)

maxcorMatrix = np.zeros([len(RDCs), len(radii)])

numberOfRuns = len(np.squeeze(np.where(np.logical_and(RDCArray == RDCs[0], radiiArray == radii[0]))))


for RDCInd, RDC in enumerate(RDCs):
    print 'RDC ' + str(RDC)
    for radiusInd, radius in enumerate(radii):
        print 'radius ' + str(radius)
        indices = np.squeeze(np.where(np.logical_and(RDCArray == RDC, radiiArray == radius)))

        otherIndices = indices[1:]
        for index in indices:
            signal1 = np.squeeze(signalArray[index])
            signal1 = signal1[startIndex:]
            signal1 = smooth(signal1,100)

            signal1Norm = signal1/np.linalg.norm(signal1)

            for otherIndex in otherIndices:
                signal2 = np.squeeze(signalArray[otherIndex])
                signal2 = signal2[startIndex:]
                signal2 = smooth(signal2, 100)
                signal2Norm = signal2 / np.linalg.norm(signal2)

                maxcorMatrix[RDCInd, radiusInd] += maxcor(signal1Norm, signal2Norm)

            otherIndices  = otherIndices[1:]
            # print otherIndices

        maxcorMatrix[RDCInd, radiusInd] /= numberOfRuns / 2 * (numberOfRuns - 1)

        #     plt.plot(signal)
        # plt.show()

        # print 'hm'

np.save('/media/carl/4ECC-1C44/PyPN/tortuous/tortuous_50_300_maxcorr_box100_10_runs', maxcorMatrix)
