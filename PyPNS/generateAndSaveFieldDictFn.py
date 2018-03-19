def genFieldDict(sourceFolder, axonXs):

    import numpy as np
    import os

    # desired location and name of dictionary
    destinationFolder = os.path.join(sourceFolder, 'numpy')

    try:
        os.makedirs(destinationFolder)
    except OSError:
        if not os.path.isdir(destinationFolder):
            raise

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


    fieldDict = load_field(sourceFolder, axonXs)
    np.save(os.path.join(destinationFolder, 'fieldDict.npy'), fieldDict)
