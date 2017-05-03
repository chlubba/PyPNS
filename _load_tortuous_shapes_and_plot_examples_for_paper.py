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

startIndex = 100

colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.), (240, 228, 66),(0, 114, 178), (213, 94, 0), (204, 121, 167)))/255

# saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/tortuous', 'toruousSignals_50_300_moreRuns.dict'), "rb" ))
saveDictMyel = pickle.load(open(os.path.join('tortuous', 'tortuousDTconstMyel2.dict'), "rb" ))
saveDictUnmyel = pickle.load(open(os.path.join('tortuous', 'tortuousDTconstUnmyel.dict'), "rb" ))

saveDicts = [saveDictUnmyel, saveDictMyel]
xlimits = [[[[21,24], [22.5, 27.5], [45,60]], [[21,24], [5,45], [13, 90]]],
           [[[0.5,1.5], [0.8,1.8], [2,3]], [[0, 3.5], [0, 3.5], [0,5]]]]
for fiberTypeInd, saveDict in enumerate(saveDicts):

    (f, axarr) = plt.subplots(2, 3, sharey='row')

    signalArray = saveDict['signals']
    RDCArray = saveDict['RDCs']
    radiiArray = saveDict['radius']

    RDCs = np.unique(RDCArray)
    radii = np.unique(radiiArray)

    maxcorMatrix = np.zeros([len(RDCs), len(radii)])

    numberOfRuns = len(np.squeeze(np.where(np.logical_and(RDCArray == RDCs[0], radiiArray == radii[0]))))


    for plotInd, RDCInd in enumerate((1,3,5)):
        RDC = RDCs[RDCInd]
        print 'RDC ' + str(RDC)
        for radiusInd, radius in enumerate(radii):
            print 'radius ' + str(radius)
            indices = np.squeeze(np.where(np.logical_and(RDCArray == RDC, radiiArray == radius)))

            otherIndices = indices[1:]
            for ii, index in enumerate(indices):
                # signal1 = np.squeeze(signalArray[index][0])
                # signal1 = signal1[startIndex:]

                time = np.arange(0,len(signalArray[index][1])*0.0025,0.0025)

                # axarr[0, RDCInd].plot(signalArray[index][0])
                axarr[0, plotInd].plot(time, signalArray[index][1], color=colors[np.mod(ii,8),:])
                axarr[1, plotInd].plot(time, signalArray[index][2], color=colors[np.mod(ii,8),:])

                axarr[0, plotInd].set_xlim(xlimits[fiberTypeInd][0][plotInd])
                axarr[1, plotInd].set_xlim(xlimits[fiberTypeInd][1][plotInd])

            #     signal1 = smooth(signal1,100)
            #
            #     signal1Norm = signal1/np.linalg.norm(signal1)
            #
            #     for otherIndex in otherIndices:
            #         signal2 = np.squeeze(signalArray[otherIndex])
            #         signal2 = signal2[startIndex:]
            #         signal2 = smooth(signal2, 100)
            #         signal2Norm = signal2 / np.linalg.norm(signal2)
            #
            #         maxcorMatrix[RDCInd, radiusInd] += maxcor(signal1Norm, signal2Norm)
            #
            #     otherIndices  = otherIndices[1:]
            #     # print otherIndices
            #
            # maxcorMatrix[RDCInd, radiusInd] /= numberOfRuns / 2 * (numberOfRuns - 1)



plt.show()
# np.save('/media/carl/4ECC-1C44/PyPN/tortuous/tortuous_50_300_maxcorr_box100_10_runs', maxcorMatrix)
