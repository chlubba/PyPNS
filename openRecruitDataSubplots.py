import numpy as np
import os
import cPickle as pickle
import matplotlib.pyplot as plt

analysisName = 'unmyelinated, monophasic inverted, 50um distance, linear'

saveFolder = '/media/carl/4ECC-1C44/PyPN/recruitment'
fileName = analysisName + '.dat'
savePath = os.path.join(saveFolder, fileName)

analysisDict = pickle.load( open(savePath, 'rb'))

print analysisDict['spikingTracker']

spikeIndices = np.where(analysisDict['spikingTracker']==1)
nonSpikeIndices = np.where(analysisDict['spikingTracker']==0)
diameters = analysisDict['diams']
amplitudes = analysisDict['amplitudes']

for i in range(len(spikeIndices[0])):
    plt.scatter(diameters[spikeIndices[0][i]], amplitudes[spikeIndices[1][i]], color=[0,1,0])
for i in range(len(nonSpikeIndices[0])):
    plt.scatter(diameters[nonSpikeIndices[0][i]], amplitudes[nonSpikeIndices[1][i]], color=[1,0,0])

plt.xlabel('diameter [um]')
plt.ylabel('amplitude [mA]')
plt.title(analysisDict['name'])
plt.show()
# print analysisDict