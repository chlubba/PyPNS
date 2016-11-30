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

try:
    os.makedirs(destinationFolder)
except OSError:
    if not os.path.isdir(destinationFolder):
        raise

# parse folder and files, folders encoding xP
xPFolderNames = [o for o in os.listdir(sourceFolder) if os.path.isdir(os.path.join(sourceFolder,o))]
xPFolderNames.remove('numpy')
xPFolders = [os.path.join(sourceFolder,o) for o in xPFolderNames]

# axon positions from folder names
axonXs = [int(''.join(ele for ele in fn if ele.isdigit())) for fn in xPFolderNames]
axonXs.sort()

# number of axon positions
axonXSteps = len(axonXs)

def load_field(folder, axonXs):

    fields = []
    axonZSteps = 0
    for xPInd, xPFolderName in enumerate(xPFolders):

        # get file names
        filenamesZP = [f for f in sorted(os.listdir(xPFolderName)) if os.path.isfile(os.path.join(xPFolderName, f))]

        # axon (source) z-positions
        axonZs = [int(''.join(ele for ele in fn.split('_')[-1] if ele.isdigit())) for fn in filenamesZP]
        axonZSteps = len(filenamesZP)

        # load each field (different axon positions)
        fieldsZ = []
        for zPInd, filename in enumerate(filenamesZP):
            fieldsZ.append(np.loadtxt(os.path.join(xPFolderName, filename)))
            print 'loaded voltage field (z: ' + str(zPInd+1) + ' / ' + str(axonZSteps) + ', x: ' + str(xPInd+1) + ' / ' + str(axonXSteps) + ')'

        fields.append(fieldsZ)

    # get coordinates (should be equal for all field files, otherwise nothing works)
    x = fields[0][0][:, 0]
    y = fields[0][0][:, 1]
    z = fields[0][0][:, 2]

    # sort by coordinate values, x changing fastest, z slowest
    orderIndices = np.lexsort((x, y, z))
    x = x[orderIndices]
    y = y[orderIndices]
    z = z[orderIndices]

    # get coordinate values
    xValues = np.unique(x)
    yValues = np.unique(y)
    zValues = np.unique(z)

    # get number of steps
    xSteps = len(xValues)
    ySteps = len(yValues)
    zSteps = len(zValues)

    # voltages are different for each field
    voltages = []
    for i in range(axonXSteps):
        voltagesZ = []
        for j in range(axonZSteps):
            v = fields[i][j][:, 3]
            v = v[orderIndices]  # order voltages as well
            voltagesZ.append(v)
        voltages.append(voltagesZ)

    # transform data to 3D-field with integer indices replacing actual coordinate values
    fieldImage = np.zeros([xSteps, ySteps, zSteps, axonXSteps, axonZSteps])

    for axonXInd in range(axonXSteps):
        for axonZInd in range(axonZSteps):
            for xInd in range(xSteps):
                for yInd in range(ySteps):
                    for zInd in range(zSteps):
                        vIndexCalc = xInd + xSteps * (yInd + zInd * ySteps)
                        fieldImage[xInd, yInd, zInd, axonXInd, axonZInd] = voltages[axonXInd][axonZInd][vIndexCalc]

    fieldDict = {'fieldImage': fieldImage,
                 'x': xValues,
                 'y': yValues,
                 'z': zValues,
                 'axonX': axonXs,
                 'axonZ': axonZs}

    return fieldDict


fieldDict = load_field(sourceFolder, axonXs)
np.save(os.path.join(destinationFolder, 'fieldDict.npy'), fieldDict)
