import numpy as np
from numpy import inf, nan
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

    zCoords = np.abs(zCoords) # in the input FEM field, we only take one side of the z-value range due to symmetry

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

# ------------------ analytic model -----------------------

def distanceOneElectrode(x, y, electrodePosition):

    distance = np.sqrt((x - electrodePosition[0]) ** 2 + (y - electrodePosition[1]) ** 2)

    # voltage = I*resistanceFromInf
    return distance


def voltageAllElectrodes(x, y, electrodePositions, sigma=1., I=1, bipolar=False):

    """

    Args:
        x: x coordinate of interest
        y: y coordinate of interest
        electrodePositions: ...
        Ra: specific resistance in ohm*cm
        I: point current strength in mA
        bipolar: if True, every uneven electrodePosition will be weighted as cathode

    Returns:

        voltage in V

    """
    sumDistance = 0
    for i in range(electrodePositions.shape[0]):
        if bipolar:
            sign = (-1) ** i
        else:
            sign = 1
        electrodePosition = electrodePositions[i, :]
        sumDistance += sign / distanceOneElectrode(x, y, electrodePosition)

    constant = 1 / (4 * np.pi + sigma) * 0.001  # Ohm*m *0.001 = kOhm*um
    resistanceFromInf = constant * sumDistance  # kOhm*um/um = kOhm
    voltage = resistanceFromInf * I * 0.001 # kOhm*mA*0.001 = mV

    return voltage

def createFieldDict(folder, axonXs):

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


if __name__ == "__main__":



    folder1 = '/media/carl/4ECC-1C44/ComsolData/thickerEndoneurium'
    folder2 = '/media/carl/4ECC-1C44/ComsolData/cuff_0.5mm'

    fieldDict1 = createFieldDict(folder1, [0, 180])
    fieldDict2 = createFieldDict(folder2, [0, 180])

    # filename1 = 'xP_0.txt'
    # filename2 = 'xP_180.txt'

    # # these are the axon positions
    # axonXs = [0, 180]
    # axonXSteps = 2
    #
    # # load each field (different axon positions)
    # fields = [[] for i in range(2)]
    # fields[0] = np.loadtxt(os.path.join(folder, filename1))
    # fields[1] = np.loadtxt(os.path.join(folder, filename2))
    #
    # print 'loaded field'
    #
    # # get coordinates (should be equal for all field files, otherwise nothing works)
    # x = fields[0][:, 0]
    # y = fields[0][:, 1]
    # z = fields[0][:, 2]
    #
    # # sort by coordinate values, x changing fastest, z slowest
    # orderIndices = np.lexsort((x,y,z))
    # x=x[orderIndices]
    # y=y[orderIndices]
    # z=z[orderIndices]
    #
    # # get coordinate values
    # xValues = np.unique(x)
    # yValues = np.unique(y)
    # zValues = np.unique(z)
    #
    # # get number of steps
    # xSteps = len(xValues)
    # ySteps = len(yValues)
    # zSteps = len(zValues)
    #
    # # voltages are different for each field
    # voltages = []
    # for i in range(axonXSteps):
    #     v = fields[i][:, 3]
    #     v = v[orderIndices] # order voltages as well
    #     voltages.append(v)
    #
    # # transform data to 3D-field with integer indices replacing actual coordinate values
    # fieldImage = np.zeros([xSteps, ySteps, zSteps, axonXSteps])
    #
    # for axonXInd in range(axonXSteps):
    #     for xInd in range(xSteps):
    #         for yInd in range(ySteps):
    #             for zInd in range(zSteps):
    #
    #                 vIndexCalc = xInd + xSteps*(yInd + zInd*ySteps)
    #                 fieldImage[xInd, yInd, zInd, axonXInd] = voltages[axonXInd][vIndexCalc]
    #
    #
    # fieldDict = {'fieldImage': fieldImage,
    #              'x': xValues,
    #              'y': yValues,
    #              'z': zValues,
    #              'axonX': axonXs}



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

    #------------------- prepare data -------------------

    # define points
    xPoints = 100
    xMin = 0.0002
    xMax = 0.0003
    xPlot = np.linspace(xMin, xMax, xPoints)
    xTicks = np.arange(xMin, xMax+0.0001, 0.002)
    yPlot = [0]
    zPoints = 500
    zMin=0
    zMax=0.00005
    zPlot = np.linspace(zMin, zMax, zPoints)
    zTicks = np.arange(zMin, zMax+0.0001, 0.01)
    xAxon = 0

    xv, zv = np.meshgrid(xPlot, zPlot)

    xvLin = np.squeeze(xv.reshape([1, -1]))
    zvLin = np.squeeze(zv.reshape([1, -1]))

    points = np.array([xvLin, np.zeros(xPoints*zPoints), zvLin, np.ones(xPoints*zPoints) * xAxon])


    values1 = interpolateFromImage(fieldDict1, points, order=1)
    values2 = interpolateFromImage(fieldDict2, points, order=1)
    # plt.semilogy(xPlot, values, label=str(180) + ' um')

    with np.errstate(invalid='ignore'):
        valuesLog1 = np.log10(values1)
        valuesLog1[valuesLog1 == -inf] = nan

        valuesLog2 = np.log10(values2)
        valuesLog2[valuesLog2 == -inf] = nan

    # ---------------------------- plot voltage image -------------------
    fieldImage = fieldDict1['fieldImage']
    plt.imshow(fieldImage[:,0,:,0], interpolation='none')
    plt.show()

    # # ------------------------- comparison plot 1D --------------------
    #
    # # x
    # pointsComp = np.array([xPlot, np.zeros(np.shape(xPlot)), np.zeros(np.shape(xPlot))])
    #
    # valuesAna = voltageAllElectrodes(pointsComp[0,:], pointsComp[1,:], np.array([[0, 0], [0, 0]]), I=0.5*10**-6)
    # valuesInter1 = interpolateFromImage(fieldDict1, np.vstack([pointsComp, np.zeros(np.shape(xPlot))]), order=1)
    # valuesInter2 = interpolateFromImage(fieldDict2, np.vstack([pointsComp, np.zeros(np.shape(xPlot))]), order=1)
    #
    # f, (ax1, ax2) = plt.subplots(2,1)
    #
    # ax1.plot(xPlot, np.log10(valuesAna), label='analytical, $\sigma$=const')
    # ax1.plot(xPlot, np.log10(valuesInter1), label='FEM nerve only')
    # ax1.plot(xPlot, np.log10(valuesInter2), label='FEM cuff')
    # ax1.set_xlabel('x [um]')
    # # ax1.legend()
    # ax1.grid(True)
    # ax1.set_title('pointsource field over x-coordinate (y=0, z=0)')
    #
    # # z
    # # pointsComp = np.array([np.ones(np.shape(xPlot))*0.000, np.zeros(np.shape(xPlot)), xPlot])
    # zPoints = 500
    # zPlot = np.linspace(-0.005, 0.005, zPoints)
    # pointsComp = np.array([np.ones(np.shape(zPlot)) * 0.0002, np.zeros(np.shape(zPlot)), zPlot])
    #
    # valuesAna = voltageAllElectrodes(pointsComp[0, :], pointsComp[2, :], np.array([[0, 0], [0, 0]]), I=0.5*10**-6)
    # valuesInter1 = interpolateFromImage(fieldDict1, np.vstack([pointsComp, np.zeros(np.shape(zPlot))]), order=1)
    # valuesInter2 = interpolateFromImage(fieldDict2, np.vstack([pointsComp, np.zeros(np.shape(zPlot))]), order=1)
    #
    # ax2.plot(zPlot, np.log10(valuesAna), label='analytical, $\sigma$=const')
    # ax2.plot(zPlot, np.log10(valuesInter1), label='FEM nerve only')
    # ax2.plot(zPlot, np.log10(valuesInter2), label='FEM cuff')
    # ax2.set_xlabel('z [um]')
    # ax2.legend()
    # ax2.grid(True)
    # ax2.set_title('pointsource field over z-coordinate (x=200um, y=0)')

    # # ------------------------- comparison plot 2D --------------------
    #
    # # x
    # pointsComp = np.array([xPlot, np.zeros(np.shape(xPlot)), np.zeros(np.shape(xPlot))])
    #
    # valuesAna = np.log10(voltageAllElectrodes(xv, zv, np.array([[0, 0], [0, 0]])))
    # valuesInter1 = valuesLog1.reshape([zPoints, xPoints])
    # valuesInter2 = valuesLog2.reshape([zPoints, xPoints])
    #
    # # normalize field strengths
    # minValAna = np.nanmin(valuesAna)
    # maxValAna = np.nanmax(valuesAna)
    # valuesAnaNorm = (valuesAna - minValAna)/(maxValAna - minValAna)
    #
    # minValInter1 = np.nanmin(valuesInter1)
    # maxValInter1 = np.nanmax(valuesInter1)
    # valuesInterNorm1 = (valuesInter1 - minValInter1) / (maxValInter1 - minValInter1)
    #
    # minValInter2 = np.nanmin(valuesInter2)
    # maxValInter2 = np.nanmax(valuesInter2)
    # valuesInterNorm2 = (valuesInter2 - minValInter2) / (maxValInter2 - minValInter2)
    #
    # # minVal = min([minValAna, minValInter])
    # # maxVal = max([maxValAna, maxValInter])
    #
    # f, (ax1, ax2, ax3) = plt.subplots(1, 3)
    # cmap = plt.get_cmap('jet') # ('gist_stern')
    #
    #
    # im1 = ax1.imshow(valuesAnaNorm, cmap=cmap, interpolation='none', vmin=0, vmax=1)
    #
    # ax1.set_xlabel('x [m]')
    # ax1.set_ylabel('z [m]')
    # ax1.set_title('Homogeneous')
    # ttl = ax1.title
    # ttl.set_position([.5, 1.05])
    # # ax2.set_xticks([0,1,13,20]) # , ('la', 'li', 'lu', 'lo')
    #
    # # xTickNum  = 10
    # # xTickLabelStep = (max(xPlot) - min(xPlot))/xTickNum
    # # zTickPositions = np.array(ax2.get_xticks()*(xPoints-1))
    # # xTickLabels = np.round(xPlot[zTickPositions.astype(int)], 1)
    # # ax2.set_xticklabels(xTickLabels)
    # xTickPositions = (xTicks - np.min(xPlot))/(np.max(xPlot) - np.min(xPlot))*(xPoints-1)
    # ax1.set_xticks(xTickPositions)
    # ax1.set_xticklabels(xTicks, rotation='vertical')
    #
    # zTickPositions = (zTicks - np.min(zPlot))/(np.max(zPlot) - np.min(zPlot))*(zPoints-1)
    # ax1.set_yticks(zTickPositions)
    # ax1.set_yticklabels(zTicks)
    #
    #
    # # ax2.set_xticks( np.arange(5))
    # # ax2.set_xticklabels(('Tom', 'Dick', 'Harry', 'Sally', 'Sue') )
    #
    # ax2.imshow(valuesInterNorm1, cmap=cmap, interpolation='none', vmin=0, vmax=1) # , aspect='auto'
    # ax2.set_xlabel('x [m]')
    # # ax2.set_ylabel('z [um]')
    # ax2.set_title('Nerve in Saline')
    # ttl = ax2.title
    # ttl.set_position([.5, 1.05])
    #
    # xTickPositions = (xTicks - np.min(xPlot)) / (np.max(xPlot) - np.min(xPlot)) * (xPoints - 1)
    # ax2.set_xticks(xTickPositions)
    # ax2.set_xticklabels(xTicks, rotation='vertical')
    #
    # ax2.set_yticks([])
    #
    # ax3.imshow(valuesInterNorm2, cmap=cmap, interpolation='none', vmin=0, vmax=1) # , aspect='auto'
    # ax3.set_xlabel('x [m]')
    # # ax3.set_ylabel('z [um]')
    # ax3.set_title('Nerve with Cuff')
    # ttl = ax3.title
    # ttl.set_position([.5, 1.05])
    #
    # xTickPositions = (xTicks - np.min(xPlot)) / (np.max(xPlot) - np.min(xPlot)) * (xPoints - 1)
    # ax3.set_xticks(xTickPositions)
    # ax3.set_xticklabels(xTicks, rotation='vertical')
    #
    # ax3.set_yticks([])
    #
    # # plt.colorbar(im1)
    # # # make room for colorbar
    # # f.subplots_adjust(right=0.8)
    # #
    # # # add colorbar axis
    # # axColorbar = f.add_axes([0.85, 0.15, 0.05, 0.7])
    # # cNorm = colors.Normalize(vmin=0, vmax=1)
    # # cb1 = mpl.colorbar.ColorbarBase(axColorbar, cmap=cmap,
    # #                                 norm=cNorm,
    # #                                 orientation='vertical')
    # # cb1.set_label('normalized voltage')


    # # ------------------- 3D plot -----------------------
    #
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    #
    # # define colors
    # jet = plt.get_cmap('jet')
    # cNorm = colors.Normalize(vmin=0, vmax=int(zPoints) - 1)
    # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
    #
    # for zInd in range(zPoints):
    #
    #     inds = np.add(np.arange(xPoints), zInd*xPoints)
    #     vValues = valuesLog[inds]
    #
    #     zValue = zPlot[zInd]
    #
    #     colorVal = scalarMap.to_rgba(int(zInd))
    #     ax.plot(xvLin[inds], np.ones(xPoints)*zValue, vValues, color=colorVal)
    #
    # ax.set_xlabel('x [um]')
    # ax.set_ylabel('z [um]')
    # ax.set_zlabel('log10(V [um])')
    #
    # # make room for colorbar
    # fig.subplots_adjust(right=0.8)
    #
    # # add colorbar axis
    # axColorbar = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #
    # cb1 = mpl.colorbar.ColorbarBase(axColorbar, cmap=jet,
    #                                 norm=cNorm,
    #                                 orientation='vertical')
    # cb1.set_label('z [um]')




    # # ----------------- 3D contour plot ---------------
    #
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # X, Y, Z = xv, zv, np.reshape(valuesLog, (zPoints, xPoints))
    # ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
    # cset = ax.contour(X, Y, Z, zdir='z', offset=0.0, cmap=cm.coolwarm)
    # cset = ax.contour(X, Y, Z, zdir='x', offset=0.005, cmap=cm.coolwarm)
    # cset = ax.contour(X, Y, Z, zdir='y', offset=0, cmap=cm.coolwarm)
    #
    # ax.set_xlabel('X')
    # ax.set_xlim(-0.01, 0.01)
    # ax.set_ylabel('Z')
    # ax.set_ylim(-0.01, 0.01)
    # ax.set_zlabel('V')
    # ax.set_zlim(-10, 0)

    # # ---------------- 2D contour plot ---------------------
    # X, Y, Z = xv, zv, np.reshape(valuesLog, (zPoints, xPoints))
    # plt.figure()
    # CS = plt.contour(X, Y, Z)
    # # plt.clabel(CS, inline=1, fontsize=10)
    #
    # circle1 = plt.Circle((0, 0), .00019, color='b', fill=False)
    # circle2 = plt.Circle((0, 0), .00023, color='b', fill=False)
    # fig = plt.gcf()
    # fig.gca().add_artist(circle1)
    # fig.gca().add_artist(circle2)
    #
    # plt.grid()




    # plt.tight_layout()
    plt.show()