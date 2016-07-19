import numpy as np
import os
import cPickle as pickle

def load_field(folder, axonXs):

        # get file names
        filenames = [f for f in sorted(os.listdir(folder)) if os.path.isfile(os.path.join(folder, f))]

        axonXSteps = len(axonXs)
        assert axonXSteps == len(filenames)

        # load each field (different axon positions)
        fields = []
        for filename in filenames:
            fields.append(np.loadtxt(os.path.join(folder, filename)))

        print 'loaded field'

        # get coordinates (should be equal for all field files, otherwise nothing works)
        x = fields[0][:, 0]
        y = fields[0][:, 1]
        z = fields[0][:, 2]

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
            v = fields[i][:, 3]
            v = v[orderIndices]  # order voltages as well
            voltages.append(v)

        # transform data to 3D-field with integer indices replacing actual coordinate values
        fieldImage = np.zeros([xSteps, ySteps, zSteps, axonXSteps])

        for axonXInd in range(axonXSteps):
            for xInd in range(xSteps):
                for yInd in range(ySteps):
                    for zInd in range(zSteps):
                        vIndexCalc = xInd + xSteps * (yInd + zInd * ySteps)
                        fieldImage[xInd, yInd, zInd, axonXInd] = voltages[axonXInd][vIndexCalc]

        fieldDict = {'fieldImage': fieldImage,
                     'x': xValues,
                     'y': yValues,
                     'z': zValues,
                     'axonX': axonXs}

        return fieldDict

folder = '/media/carl/4ECC-1C44/ComsolData/cuffFiner' # '/media/carl/4ECC-1C44/ComsolData/noCuffFiner/z0.03_1000,x0.0015_100,y_asym' # '/media/carl/4ECC-1C44/ComsolData/noCuffFiner/z0.003_100/Lagrange_Smoothing_Finer'
axonXs = [0, 180]

fieldDict = load_field(folder, axonXs)

# f = open('/media/carl/4ECC-1C44/ComsolData/noCuffFiner/z0.003_100/Lagrange_Smoothing_Finer/numpy/fieldDict.npy', 'wb')
np.save(os.path.join(folder, 'numpy', 'fieldDict.npy'), fieldDict)
# pickle.dump(fieldDict, open(os.path.join(folder, 'numpy', 'fieldDict.pickle'), 'wb'))
# f.close()

# loadedFieldArray = np.load(os.path.join(folder, 'numpy', 'fieldDict.npy'))
# loadedField = loadedFieldArray[0]
# print loadedField