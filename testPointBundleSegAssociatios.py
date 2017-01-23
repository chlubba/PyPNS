import numpy as np
import PyPN
import PyPN.createGeometry as createGeometry
from PyPN.takeTime import takeTime

# set length of bundle and number of axons
lengthOfBundle = 40000

# bundle guide
segmentLengthAxon = 30
bundleGuide = np.array(PyPN.createGeometry.get_bundle_guide_straight_radius(lengthOfBundle, segmentLengthAxon))

# calculate the random axon coordinates
axonCoords = np.array(createGeometry.create_random_axon(bundleGuide, [0,0], segmentLengthAxon))
numAxonSegments = np.shape(axonCoords)[0]

bundleGuideLength = np.shape(bundleGuide)[0]
# print bundleGuideLength

bundleGuideStripped = bundleGuide[:,0:3] # no radius

# with takeTime('calculating squared distances'):
#     # variable to save squared distances between bundle guide segment centers and source positions
#     r2min = np.ones(numAxonSegments)*np.inf
#     closestSegInds = np.zeros(numAxonSegments)
#     for bundleGuideInd in range(bundleGuideLength-1):
#
#         bundleSegStart = bundleGuideStripped[bundleGuideInd,:]
#         bundleSegStop = bundleGuideStripped[bundleGuideInd+1, :]
#         bundleMiddle = (bundleSegStart + bundleSegStop) / 2
#
#
#         r2 = np.mean(np.square(axonCoords - bundleMiddle), axis=1)
#
#         compArray = r2 < r2min
#
#         closestSegInds[compArray] = bundleGuideInd
#         r2min[compArray] = r2[compArray]
#
#
# r2minseg = np.vstack([np.transpose(r2min), np.transpose(closestSegInds)])
# print 'bla'

def associatePointToBundleSegs(points, bundleGuide):

    # make sure orientation is fine
    pointsShape = np.shape(points)
    if pointsShape[0] == 3:
        points = np.transpose(points)
    numPoints = np.shape(points)[0]

    # cut away potential radius coordinate
    bundleGuide = bundleGuide[:,0:3]
    bundleGuideLength = np.shape(bundleGuide)[0]

    # variable to save squared distances between bundle guide segment centers and source positions
    r2min = np.ones(numPoints) * np.inf
    # indices of closest bundle segment for all points
    closestSegInds = np.zeros(numPoints)

    for bundleGuideInd in range(bundleGuideLength - 1):
        bundleSegStart = bundleGuideStripped[bundleGuideInd, :]
        bundleSegStop = bundleGuideStripped[bundleGuideInd + 1, :]
        bundleMiddle = (bundleSegStart + bundleSegStop) / 2

        r2 = np.sum(np.square(axonCoords - bundleMiddle), axis=1)

        # current distance smaller than current minimum?
        compArray = r2 < r2min

        closestSegInds[compArray] = bundleGuideInd
        r2min[compArray] = r2[compArray]

    return closestSegInds

with takeTime('calculating squared distances'):
    print associatePointToBundleSegs(axonCoords, bundleGuide)