import numpy as np
import os
import cPickle as pickle

# set path according to location of the exported field from COMSOL
# sourceFolder = '/media/carl/4ECC-1C44/ComsolData/cuffFiner'
# sourceFolder = '/media/carl/4ECC-1C44/ComsolData/stimulationField/noCuffStimulation'
# sourceFolder = '/media/carl/18D40D77D40D5900/COMSOL_data/oil_sigma_0.0000001_contact_xP0'
sourceFolder = '/media/carl/4ECC-1C44/ComsolData/oil_different_source_positions'
# sourceFolder = '/media/carl/4ECC-1C44/PyPN/voltageFieldDummyzVar'

# desired location and name of dictionary
destinationFolder = os.path.join(sourceFolder, 'numpy')

fieldDictArray = np.load(os.path.join(destinationFolder,'fieldDict.npy'))

fieldDict = fieldDictArray[()]


fieldImageOld = fieldDict['fieldImage']

xSteps = len(fieldDict['x'])
ySteps = len(fieldDict['y'])
zSteps = len(fieldDict['z'])
axonXSteps = len(fieldDict['axonX'])
try:
    axonYSteps = len(fieldDict['axonY'])
except:
    axonYSteps = 1
    pass
try:
    axonZSteps = len(fieldDict['axonZ'])
except:
    axonZSteps = 1
    pass

# correct the units of axon positions from um to m
fieldDict['axonX'] = np.array(fieldDict['axonX']).astype(float)
fieldDict['axonX'] = fieldDict['axonX'] / 1000000
try:
    fieldDict['axonY'] = np.array(fieldDict['axonY']).astype(float)
    fieldDict['axonY'] = fieldDict['axonY'] / 1000000
except:
    pass
try:
    fieldDict['axonZ'] = np.array(fieldDict['axonZ']).astype(float)
    fieldDict['axonZ'] = fieldDict['axonZ'] / 1000000
except:
    pass


fieldImageNew = np.zeros([xSteps, ySteps, zSteps, axonXSteps, axonYSteps, axonZSteps])

# this here does not care about the possible non existance of a axon z-coordinate and will crash then.
for axonXInd in range(axonXSteps):
    for axonYInd in range(axonYSteps):
        for axonZInd in range(axonZSteps):
            for xInd in range(xSteps):
                for yInd in range(ySteps):
                    for zInd in range(zSteps):
                        vIndexCalc = xInd + xSteps * (yInd + zInd * ySteps)
                        fieldImageNew[xInd, yInd, zInd, axonXInd, axonYInd, axonZInd] = fieldImageOld[xInd, yInd, zInd, axonXInd, axonZInd]


fieldDict['fieldImage'] = fieldImageNew
fieldDict['axonY'] = np.array([0])


# {'fieldImage': fieldImage,
#              'x': xValues,
#              'y': yValues,
#              'z': zValues,
#              'axonX': axonXs,
#              'axonZ': axonZs}


np.save(os.path.join(destinationFolder, 'fieldDictSymEnco.npy'), fieldDict)
