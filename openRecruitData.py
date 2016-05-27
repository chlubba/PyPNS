import numpy as np
import os
import cPickle as pickle
import matplotlib.pyplot as plt

analysisNames = ['unmyelinated,MONOPHASIC ,inverted-False, 50um distance, linear',
                 'unmyelinated, monophasic noninverted, 50um distance, linear',
                 'unmyelinated, biphasic inverted, 50um distance, linear',
                 'unmyelinated, biphasic noninverted, 50um distance, linear']

f, axarr = plt.subplots(2,2)

for k in range(4):

    analysisName = analysisNames[k]

    saveFolder = '/media/carl/4ECC-1C44/PyPN/recruitment'
    fileName = analysisName + '.dat'
    savePath = os.path.join(saveFolder, fileName)

    try:
        analysisDict = pickle.load( open(savePath, 'rb'))
    except:
        continue

    # print analysisDict['spikingTracker']

    spikeIndices = np.where(analysisDict['spikingTracker']==1)
    nonSpikeIndices = np.where(analysisDict['spikingTracker']==0)
    diameters = analysisDict['diams']
    amplitudes = analysisDict['amplitudes']

    for i in range(len(spikeIndices[0])):
        axarr[np.floor(k/2), k%2].scatter(diameters[spikeIndices[0][i]], amplitudes[spikeIndices[1][i]], color=[0,1,0])
    for i in range(len(nonSpikeIndices[0])):
        axarr[np.floor(k/2), k%2].scatter(diameters[nonSpikeIndices[0][i]], amplitudes[nonSpikeIndices[1][i]], color=[1,0,0])

    axarr[np.floor(k/2), k%2].set_xlabel('diameter [um]')
    axarr[np.floor(k/2), k%2].set_ylabel('amplitude [mA]')
    axarr[np.floor(k/2), k%2].set_title(analysisDict['name'])
plt.show()
# print analysisDict