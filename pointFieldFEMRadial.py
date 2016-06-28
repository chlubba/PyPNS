import numpy as np
from scipy import ndimage
import os
import matplotlib.pyplot as plt

folder = '/media/carl/4ECC-1C44/ComsolData/thickerEndoneurium'
filename1 = 'xP_0.txt'
filename2 = 'xP_180.txt'

field1 = np.loadtxt(os.path.join(folder, filename1))
# field2 = np.loadtxt(os.path.join(folder, filename2))


# print np.shape(field2)
# print field2.nbytes

x = field1[:,0]
y = field1[:,1]
z = field1[:,2]
v = field1[:,3]

# transform into
rho = np.sqrt(x**2 + y**2)
phi = np.arctan2(y,x)

# sort by coordinate values, x changing fastest, z slowest
orderIndices = np.lexsort((z,phi,rho))
rho=rho[orderIndices]
phi=phi[orderIndices]
z=z[orderIndices]
v=v[orderIndices]

# transform data to 3D-field with integer indices replacing actual coordinate values

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
            vIndex = xInd + xSteps*(yInd + zInd*ySteps)
            fieldImage[xInd, yInd, zInd] = v[vIndex]

fieldDict = {'fieldImage': fieldImage,
             'x': xValues,
             'y': yValues,
             'z': yValues,}

# print fieldImage.shape
def getImageCoords(xValues, yValues, zValues, points):

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

    points = np.array(points)

    if min(points.shape)>1:
        if points.shape[1] > 3:
            points = np.transpose(points)
        xCoords = np.add(points[:, 0], -xMin) / (xMax - xMin) * xNum
        yCoords = np.add(points[:, 1], -yMin)/ (yMax - yMin) * yNum
        zCoords = np.add(points[:, 2], -zMin) / (zMax - zMin) * zNum
    else:
        xCoords = (points[0] - xMin) / (xMax - xMin) * xNum
        yCoords = (points[1] - yMin) / (yMax - yMin) * yNum
        zCoords = (points[2] - zMin) / (zMax - zMin) * zNum

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

xPlot = np.linspace(-0.005, 0.005, 100)

points = np.array([xPlot, np.zeros(100), np.zeros(100)])
values = interpolateFromImage(fieldDict, points)

plt.plot(xPlot, values)
plt.show()



