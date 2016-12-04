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

fieldDict['axonX'] = np.array(fieldDict['axonX'])/1000000
try:
    fieldDict['axonZ'] = np.array(fieldDict['axonZ'])/1000000
except:
    pass

# {'fieldImage': fieldImage,
#              'x': xValues,
#              'y': yValues,
#              'z': zValues,
#              'axonX': axonXs,
#              'axonZ': axonZs}


np.save(os.path.join(destinationFolder, 'fieldDictCor.npy'), fieldDict)
