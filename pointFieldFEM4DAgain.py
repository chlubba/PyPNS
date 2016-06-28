import numpy as np
from scipy import ndimage
import os
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colorbar

# print fieldImage.shape
def getImageCoords(xValues, yValues, zValues, axonXValues, points):

    # assume equidistant original points

    xMin = min(xValues)
    xMax = max(xValues)
    xNum = len(xValues)

    yMin = min(yValues)
    yMax = max(yValues)
    yNum = len(yValues)

    zMin = min(zValues)
    zMax = max(zValues)
    zNum = len(zValues)

    axonXMin = min(axonXValues)
    axonXMax = max(axonXValues)
    axonXNum = len(axonXValues)

    points = np.array(points)

    if len(points.shape)>1:
        if points.shape[1] > 3:
            points = np.transpose(points)
        xCoords = np.add(points[:, 0], -xMin) / (xMax - xMin) * (xNum-1)
        yCoords = np.add(points[:, 1], -yMin)/ (yMax - yMin) * (yNum-1)
        zCoords = np.add(points[:, 2], -zMin) / (zMax - zMin) * (zNum-1)
        xAxonCoords = np.add(points[:, 3], -axonXMin) / (axonXMax - axonXMin) * (axonXNum-1)
    else:
        xCoords = (points[0] - xMin) / (xMax - xMin) * (xNum-1)
        yCoords = (points[1] - yMin) / (yMax - yMin) * (yNum-1)
        zCoords = (points[2] - zMin) / (zMax - zMin) * (zNum-1)
        xAxonCoords = (points[3] - axonXMin) / (axonXMax - axonXMin) * (axonXNum-1)

    return np.vstack([xCoords, yCoords, zCoords, xAxonCoords])

def interpolateFromImage(fieldDict, points, order=3):

    # first transform coordinates in points into position coordinates
    imageCoords = getImageCoords(fieldDict['x'], fieldDict['y'], fieldDict['z'], fieldDict['axonX'], points)

    # then with new coords to the interpolation
    return  ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

# print getImageCoords(xValues, yValues, zValues, np.array([[1,2,3], [1,2,3], [1,2,3]]))
#
# print ndimage.map_coordinates(fieldImage, [[0.01], [0.01], [0.01]], order=3)
#
# print 'hm'

# print ndimage.map_coordinates(a, [[0.5], [0.5]], order=3)

if __name__ == "__main__":



    folder = '/media/carl/4ECC-1C44/ComsolData/thickerEndoneurium'
    filename1 = 'xP_0.txt'
    filename2 = 'xP_180.txt'

    # these are the axon positions
    axonXs = [0, 180]
    axonXSteps = 2

    # load each field (different axon positions)
    fields = [[] for i in range(2)]
    fields[0] = np.loadtxt(os.path.join(folder, filename1))
    fields[1] = np.loadtxt(os.path.join(folder, filename2))

    print 'loaded field'

    # get coordinates (should be equal for all field files, otherwise nothing works)
    x = fields[0][:, 0]
    y = fields[0][:, 1]
    z = fields[0][:, 2]

    # sort by coordinate values, x changing fastest, z slowest
    orderIndices = np.lexsort((x,y,z))
    x=x[orderIndices]
    y=y[orderIndices]
    z=z[orderIndices]

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
        v = v[orderIndices] # order voltages as well
        voltages.append(v)

    # transform data to 3D-field with integer indices replacing actual coordinate values
    fieldImage = np.zeros([xSteps, ySteps, zSteps, axonXSteps])

    for axonXInd in range(axonXSteps):
        for xInd in range(xSteps):
            for yInd in range(ySteps):
                for zInd in range(zSteps):

                    vIndexCalc = xInd + xSteps*(yInd + zInd*ySteps)
                    fieldImage[xInd, yInd, zInd, axonXInd] = voltages[axonXInd][vIndexCalc]


    fieldDict = {'fieldImage': fieldImage,
                 'x': xValues,
                 'y': yValues,
                 'z': zValues,
                 'axonX': axonXs}

    # xPoints = 1000
    # xPlot = np.linspace(-0.005, 0.005, xPoints)
    # xAxonPlot = np.linspace(0, 180, 5) # np.array([0]) #
    #
    # # points = np.array([np.zeros(100), np.zeros(100), np.zeros(100), xAxonPlot])
    # for xAxonVal in xAxonPlot:
    #     points = np.array([xPlot, np.zeros(xPoints), np.ones(xPoints)*0.0, np.ones(xPoints)*xAxonVal])
    #
    #     for order in [1]: # range(1,4):
    #         values = interpolateFromImage(fieldDict, points, order=order)
    #         plt.semilogy(xPlot, values, label=str(xAxonVal)+' um')
    #         # plt.semilogy(xPlot, values, label='order '+str(order))
    #
    # plt.xlabel('x [um]')
    # plt.ylabel('V [V]')
    # plt.title('Voltage caused by current point source of of 1 nA')
    # plt.legend()
    # plt.show()

    # ------------------- 3D plot -----------------------

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # define points
    xPoints = 100
    xPlot = np.linspace(-0.005, 0.005, xPoints)
    yPlot = [0]
    zPoints = 30
    zPlot = np.linspace(0, 0.005, zPoints)

    xv, zv = np.meshgrid(xPlot, zPlot)

    xvLin = np.squeeze(xv.reshape([1, -1]))
    zvLin = np.squeeze(zv.reshape([1, -1]))

    points = np.array([xvLin, np.zeros(xPoints*zPoints), zvLin, np.ones(xPoints*zPoints) * 180])


    values = interpolateFromImage(fieldDict, points, order=1)
    # plt.semilogy(xPlot, values, label=str(180) + ' um')

    valuesLog = np.log10(values)

    # define colors
    jet = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=int(zPoints) - 1)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

    for zInd in range(zPoints):

        inds = np.add(np.arange(xPoints), zInd*xPoints)
        vValues = valuesLog[inds]

        zValue = zPlot[zInd]

        colorVal = scalarMap.to_rgba(int(zInd))
        ax.plot(xvLin[inds], np.ones(xPoints)*zValue, vValues, color=colorVal)

    ax.set_xlabel('x [um]')
    ax.set_ylabel('z [um]')
    ax.set_zlabel('log10(V [um])')

    # make room for colorbar
    fig.subplots_adjust(right=0.8)

    # add colorbar axis
    axColorbar = fig.add_axes([0.85, 0.15, 0.05, 0.7])

    cb1 = mpl.colorbar.ColorbarBase(axColorbar, cmap=jet,
                                    norm=cNorm,
                                    orientation='vertical')
    cb1.set_label('z [um]')


    # ax.set_zscale('log')
    # ax.legend()

    plt.show()

