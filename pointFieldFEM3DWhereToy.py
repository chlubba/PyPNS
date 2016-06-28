import numpy as np
from scipy import ndimage
import os
import matplotlib.pyplot as plt

# print fieldImage.shape
def getImageCoords(xValues, yValues, zValues, points):

    # assume equidistant original points

    xMin = min(xValues)
    xMax = max(xValues)
    xNum = len(xValues)
    xStep = xValues[1] - xValues[0]

    yMin = min(yValues)
    yMax = max(yValues)
    yNum = len(yValues)
    yStep = yValues[1] - yValues[0]

    zMin = min(zValues)
    zMax = max(zValues)
    zNum = len(zValues)
    zStep = zValues[1] - zValues[0]

    points = np.array(points)

    if len(points.shape)>1:
        if points.shape[1] > 3:
            points = np.transpose(points)
        xCoords = np.add(points[:, 0], -xMin) / (xMax - xMin) * (xNum-1)
        yCoords = np.add(points[:, 1], -yMin)/ (yMax - yMin) * (yNum-1)
        zCoords = np.add(points[:, 2], -zMin) / (zMax - zMin) * (zNum-1)
    else:
        xCoords = (points[0] - xMin) / (xMax - xMin) * (xNum-1)
        yCoords = (points[1] - yMin) / (yMax - yMin) * (yNum-1)
        zCoords = (points[2] - zMin) / (zMax - zMin) * (zNum-1)

    return np.vstack([xCoords, yCoords, zCoords])

def interpolateFromImage(fieldDict, points):

    # first transform coordinates in points into position coordinates
    imageCoords = getImageCoords(fieldDict['x'], fieldDict['y'], fieldDict['z'], points)

    # then with new coords to the interpolation
    return  ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=3)

# print getImageCoords(xValues, yValues, zValues, np.array([[1,2,3], [1,2,3], [1,2,3]]))
#
# print ndimage.map_coordinates(fieldImage, [[0.01], [0.01], [0.01]], order=3)
#
# print 'hm'

# print ndimage.map_coordinates(a, [[0.5], [0.5]], order=3)

if __name__ == "__main__":

    z = [-1., 0., 1., 2.]
    y = [0., 1.]
    x = [0., 1.]

    field = np.array([   [0., 0.,-1.],
                [0., 0., 0.],
                [0., 0., 1.],
                [0., 0., 2.],
                [0., 1.,-1.],
                [0., 1., 0.],
                [0., 1., 1.],
                [0., 1., 2.],
                [1., 0.,-1.],
                [1., 0., 0.],
                [1., 0., 1.],
                [1., 0., 2.],
                [1., 1.,-1.],
                [1., 1., 0.],
                [1., 1., 1.],
                [1., 1., 2.],
                ])

    v = []
    for xInd in range(len(x)):
        for yInd in range(len(y)):
            for zInd in range(len(z)):
                v.append(x[xInd] + y[yInd] + z[zInd])
    v = np.array(v)
    # v = np.random.randn(16)

    # field = np.column_stack((fieldCoords, values))

    #
    # a = getImageCoords(x, y, z, point)
    # print a

    # get values
    xValues = np.unique(x)
    yValues = np.unique(y)
    zValues = np.unique(z)

    # get number of steps
    xSteps = len(xValues)
    ySteps = len(yValues)
    zSteps = len(zValues)
    # print xSteps, ySteps, zSteps

    # translate to 3D-array
    fieldImage = np.zeros([xSteps, ySteps, zSteps])
    for xInd in range(xSteps):
        for yInd in range(ySteps):
            for zInd in range(zSteps):
                xValue = xValues[xInd]
                yValue = yValues[yInd]
                zValue = zValues[zInd]

                vIndex = np.where(
                    np.logical_and(field[:, 0] == xValue, np.logical_and(field[:, 1] == yValue, field[:, 2] == zValue)))

                # vIndex = xInd + xSteps*(yInd + zInd*ySteps)
                fieldImage[xInd, yInd, zInd] = v[vIndex]

    fieldDict = {'fieldImage': fieldImage,
                 'x': xValues,
                 'y': yValues,
                 'z': zValues,}



    values = interpolateFromImage(fieldDict, field)

    print values - v


                # folder = '/media/carl/4ECC-1C44/ComsolData/thickerEndoneurium'
    # filename1 = 'xP_0.txt'
    # filename2 = 'xP_180.txt'
    #
    # field1 = np.loadtxt(os.path.join(folder, filename1))
    # # field2 = np.loadtxt(os.path.join(folder, filename2))
    #
    #
    # # print np.shape(field2)
    # # print field2.nbytes
    #
    # x = field1[:, 0]
    # y = field1[:, 1]
    # z = field1[:, 2]
    # v = field1[:, 3]
    #
    # # # sort by coordinate values, x changing fastest, z slowest
    # # orderIndices = np.lexsort((z,y,x))
    # # x=x[orderIndices]
    # # y=y[orderIndices]
    # # z=z[orderIndices]
    # # v=v[orderIndices]
    #
    # # transform data to 3D-field with integer indices replacing actual coordinate values
    #
    # # get values
    # xValues = np.unique(x)
    # yValues = np.unique(y)
    # zValues = np.unique(z)
    #
    # # get number of steps
    # xSteps = len(xValues)
    # ySteps = len(yValues)
    # zSteps = len(zValues)
    # # print xSteps, ySteps, zSteps
    #
    # # translate to 3D-array
    # fieldImage = np.zeros([xSteps, ySteps, zSteps])
    # for xInd in range(xSteps):
    #     for yInd in range(ySteps):
    #         for zInd in range(zSteps):
    #             xValue = xValues[xInd]
    #             yValue = xValues[yInd]
    #             zValue = xValues[zInd]
    #
    #             vIndex = np.where(
    #                 np.logical_and(field1[:, 1] == xValue, field1[:, 2] == yValue, field1[:, 2] == zValue))
    #
    #             # vIndex = xInd + xSteps*(yInd + zInd*ySteps)
    #             fieldImage[xInd, yInd, zInd] = v[vIndex]
    #
    # fieldDict = {'fieldImage': fieldImage,
    #              'x': xValues,
    #              'y': yValues,
    #              'z': zValues,}
    #
    # # ---------------- plot result of interpolation ---------------------------
    #
    # xPlot = np.linspace(-0.005, 0.005, 100)
    #
    # voltageOriginalInd = np.where(np.logical_and(field1[:,1] == min(abs(field1[:,1])), field1[:,2] == min(field1[:,2])))
    # voltageOriginal = np.squeeze(field1[voltageOriginalInd, 3])
    # xvOriginal = np.squeeze(field1[voltageOriginalInd, 0])
    #
    # plt.semilogy(xvOriginal, voltageOriginal)
    #
    # points = np.array([xPlot, np.zeros(100), np.zeros(100)])
    # values = interpolateFromImage(fieldDict, points)
    #
    # plt.figure()
    # plt.semilogy(xPlot, values)
    # plt.show()



