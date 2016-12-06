import numpy as np
from scipy import ndimage
from scipy.interpolate import interp1d
import time
import createGeometry

def getImageCoords1D(physicalCoordinate1D, pointsOfInterest):

    if len(physicalCoordinate1D) > 1:

        # interpolate to obtain coordinate position
        physicalCoordinate1D.sort()
        coordInterp = interp1d(physicalCoordinate1D, range(len(physicalCoordinate1D)))

        pointsOfInterest = np.array(pointsOfInterest)

        # if only one half of coordinates was exported, is is to be assumed, that a symmetry in the geometry justifies that.
        # we therefore mirror here by taking the absolute value of the input coordinates
        if np.min(physicalCoordinate1D) >= 0:
            coords = coordInterp(np.abs(pointsOfInterest))
        else:
            coords = coordInterp(pointsOfInterest)

    else:
        # a single value only signifies that this coordinate does not interest. Only for source positions.
        coords = np.zeros(len(pointsOfInterest))

    return coords

def _getImageCoordsUnregXZ(fieldDict, points):

    xValues = fieldDict['x']
    yValues = fieldDict['y']
    zValues = fieldDict['z']
    axonXValues = fieldDict['axonX']
    axonZValues = fieldDict['axonZ']

    # interpolate to obtain coordinate position
    xCoordInterp = interp1d(xValues, range(len(xValues)))
    yCoordInterp = interp1d(yValues, range(len(yValues)))
    zCoordInterp = interp1d(zValues, range(len(zValues)))
    axonXCoordInterp = interp1d(axonXValues, range(len(axonXValues)))
    axonZCoordInterp = interp1d(axonZValues, range(len(axonZValues)))

    points = np.array(points)

    if len(points.shape) > 1:
        points = np.transpose(points)
        xCoords = xCoordInterp(points[:, 0])
        yCoords = yCoordInterp(points[:, 1])
        zCoords = zCoordInterp(points[:, 2])
        xAxonCoords = axonXCoordInterp(points[:, 3])
        zAxonCoords = axonZCoordInterp(points[:, 4])
    else:
        xCoords = xCoordInterp(points[0])
        yCoords = yCoordInterp(points[1])
        zCoords = zCoordInterp(points[2])
        xAxonCoords = axonXCoordInterp(points[3])
        zAxonCoords = axonZCoordInterp(points[4])

    # mirror coordinates where the symmetry allows it. Here only of y, NOT for z
    yCoords = np.abs(yCoords)

    return np.vstack([xCoords, yCoords, zCoords, xAxonCoords, zAxonCoords])

def _getImageCoords(fieldDict, points):

    xValues = fieldDict['x']
    yValues = fieldDict['y']
    zValues = fieldDict['z']
    axonXValues = fieldDict['axonX']

    ZGiven = False
    try:
        axonZValues = fieldDict['axonZ']
        ZGiven = True
    except:
        pass

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

    # todo: change from um to m! in generate dict script.
    axonXMin = min(axonXValues)
    axonXMax = max(axonXValues)
    axonXNum = len(axonXValues)

    if ZGiven:
        axonZMin = min(axonZValues)
        axonZMax = max(axonZValues)
        axonZNum = len(axonZValues)

    points = np.array(points)

    if len(points.shape) > 1:
        points = np.transpose(points)
        xCoords = np.add(points[:, 0], -xMin) / (xMax - xMin) * (xNum - 1)
        yCoords = np.add(points[:, 1], -yMin) / (yMax - yMin) * (yNum - 1)
        zCoords = np.add(points[:, 2], -zMin) / (zMax - zMin) * (zNum - 1)
        xAxonCoords = np.add(points[:, 3], -axonXMin) / (axonXMax - axonXMin) * (axonXNum - 1)
        if ZGiven:
            zAxonCoords = np.add(points[:, 4], -axonZMin) / (axonZMax - axonZMin) * (axonZNum - 1)
    else:
        xCoords = (points[0] - xMin) / (xMax - xMin) * (xNum - 1)
        yCoords = (points[1] - yMin) / (yMax - yMin) * (yNum - 1)
        zCoords = (points[2] - zMin) / (zMax - zMin) * (zNum - 1)
        xAxonCoords = (points[3] - axonXMin) / (axonXMax - axonXMin) * (axonXNum - 1)
        if ZGiven:
            zAxonCoords = np.add(points[4], -axonZMin) / (axonZMax - axonZMin) * (axonZNum - 1)

    if not ZGiven:
        zCoords = np.abs(zCoords)  # in the input FEM field, we only take one side of the z-value range thanks to symmetry

    yCoords = np.abs(yCoords)  # same for y-Coordinates

    if ZGiven:
        combinedCoords = np.vstack([xCoords, yCoords, zCoords, xAxonCoords, zAxonCoords])
    else:
        combinedCoords = np.vstack([xCoords, yCoords, zCoords, xAxonCoords])

    return combinedCoords



# def _interpolateFromImageGranular(fieldDict, points, order=3, zGiven=False):
#
#     fieldKeys = (('axonX', 'x'), ('axonY', 'y'), ('axonZ', 'z'))
#     receiverCoords = []
#     sourceCoords = []
#     coordInd = 0
#
#
#     print points
#
#     # process coordinate pairs (axon position and electrode position). If for a coordinate the axon position is only
#     # given in one direction (only positive coordinate values), for negative coordinate values the field of the positive
#     # position is taken. But then also the receiver position needs to be mirrored on the axis of symmetry.
#     for coordKeyPair in fieldKeys:
#
#         axonFieldKey = coordKeyPair[0]
#
#         signArray = []
#         try:
#             # exported (from FEM simulation) coordinate values for the axis of interest
#             physicalCoords = fieldDict[axonFieldKey]
#
#             # get the coordinates of interest
#             if len(points.shape) > 1:
#                 coordVals = points[coordInd, :]
#             else:
#                 coordVals = points[coordInd]
#
#             sourceCoords.append(_getImageCoords1D(physicalCoords, coordVals))
#
#             # if the exported coordinates for this dimension are only positive, this means the source position is taken
#             # as an absolute value. Therefore, if it actually was a negative value, the field needs to be mirrored along
#             # the symmetry axis. This mirroring is equivalent to a sign change of the receiver coordinate value.
#             if min(physicalCoords) >= 0:
#                 signArray = np.sign(coordVals)
#                 signArray[signArray == 0] = 1 # if no mirroring neccessary, don't.
#             else:
#                 signArray = np.ones(np.shape(physicalCoords))
#
#             coordInd += 1
#         except:
#             pass
#
#         # now process the receiver.
#         receiverFieldKey = coordKeyPair[1]
#         physicalCoords = fieldDict[receiverFieldKey]
#
#         # get the coordinates of interest
#         if len(points.shape) > 1:
#             coordVals = points[coordInd, :]
#         else:
#             coordVals = points[coordInd]
#
#         # mirror source positions along symmetry axis
#         if not np.size(signArray) == 0:
#             coordVals = np.multiply(coordVals, signArray)
#
#         # get the coordinated mapped to image 'pixel indices'
#         receiverCoords.append(_getImageCoords1D(physicalCoords, coordVals))
#
#         coordInd += 1
#
#     combinedCoords = np.vstack((np.array(receiverCoords), np.array(sourceCoords)))
#
#     print combinedCoords
#
#     imageCoords = np.array(combinedCoords)

# def _interpolateFromImageGranular(fieldDict, points, order=3, zGiven=False):
#
#     fieldKeys = ('x', 'y', 'z', 'axonX', 'axonY', 'axonZ')
#     combinedCoords = []
#     coordInd = 0
#     for fieldKey in fieldKeys:
#
#         try:
#             physicalCoords = fieldDict[fieldKey]
#         except:
#             continue
#
#         if len(points.shape) > 1:
#             combinedCoords.append(_getImageCoords1D(physicalCoords, points[coordInd,:]))
#         else:
#             combinedCoords.append(_getImageCoords1D(physicalCoords, points[coordInd]))
#
#         coordInd += 1
#
#     imageCoords = np.array(combinedCoords)
#
#     print imageCoords

#
#     # then with new coords to the interpolation
#     return ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

# def _interpolateFromImage(fieldDict, points, order=3, zGiven=False):
#
#     # first transform coordinates in points into position coordinates
#     # different function
#     if zGiven:
#         imageCoords = _getImageCoordsUnregXZ(fieldDict, points)
#     else:
#         imageCoords = _getImageCoords(fieldDict, points)
#
#     # then with new coords to the interpolation
#     return ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

def _interpolateFromImage(fieldDict, points, order=3, zGiven=False):

    # first transform coordinates in points into position coordinates
    # different function

    imageCoords = _getImageCoords(fieldDict, points)
    #
    # print imageCoords
    #
    # print '\n'

    # _interpolateFromImageGranular(fieldDict, points, order)



    # then with new coords to the interpolation
    return ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

# def _interpolateFromImage(fieldDict, points, order=3, zGiven=False):
#
#     # first transform coordinates in points into position coordinates
#     # different function
#
#     imageCoords = _getImageCoords(fieldDict, points)
#
#     # then with new coords to the interpolation
#     return ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

def compute_relative_positions_and_interpolate_symmetric_inhomogeneity(sourcePositions, sourceCurrents, receiverPositions, fieldDict, bundleGuide, receiverDisplacement=0, currentUnitFEM=-9, currentUnitSource=-9):

    """

    Structure as follows

    for bs in bundleSegments:
        for s in sources on this bundle segment:
            calculate distance of s from axon guide and distance along axon guide
        for r in all receivers calculate:
            for s in sources on this bundle segment:
                  calculate receiver position and interpolate potential based on spatial relation between source
                  and receiver and source and bundle

    Args:
        axon: needed for segment positions and membrane currents

    Returns:

    """


    nSourcePoints = np.shape(sourcePositions)[0]


    bundleCoords = bundleGuide[:, 0:3]

    # find the bundle guide segments that are closest to the electrode points
    # first calculate bundle segment centers
    bundleSegCenters = bundleCoords[1:,:] - bundleCoords[0:-1,:]

    # then calculate the distances between one electrode point and one bundle segment center
    closestSegInds = []
    receiverDistToOrigins = []
    receiverXDists = []
    for receiverPosition in receiverPositions:

        r2 = (bundleSegCenters[:,0] - receiverPosition[0]) ** 2 + (bundleSegCenters[:,1] - receiverPosition[1]) ** 2 + (bundleSegCenters[:,2] - receiverPosition[2]) ** 2
        r = np.sqrt(r2)

        closestSegInds.append(np.argmin(r))

        bundleSegStartPoint = bundleCoords[closestSegInds[-1]]

        # calculate projection and radius
        dir0 = bundleCoords[closestSegInds[-1]+1] - bundleSegStartPoint
        dir0norm = dir0/np.linalg.norm(dir0)

        # distance of receiver along bundle guide
        receiverPosPar = np.inner(receiverPosition - bundleCoords[closestSegInds[-1]], dir0norm)
        receiverDistToOrigin = createGeometry.length_from_coords(bundleCoords[:closestSegInds[-1]]) + receiverPosPar

        receiverDistToOrigins.append(receiverDistToOrigin)

        # normal from bundle guide to axon position
        receiverDir2D = receiverPosition - (bundleSegStartPoint + receiverPosPar * dir0norm)
        receiverXDist = np.linalg.norm(receiverDir2D)

        receiverXDists.append(receiverXDist)

        # if not sourceXDist == 0:
        #     sourceDir2DNorm = sourceDir2D / sourceXDist
        # else:
        #     sourceDir2DNorm = -1



    # first go through all bundle segments, find the associated source positions and calculate the needed
    # quantities for all sources
    receiverPotentials = np.zeros((receiverPositions.shape[0], np.shape(sourceCurrents)[1]))

    # # TODO: delete this again, only for testing
    # lastElectrodeSignal = []
    # f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    # jet = plt.get_cmap('jet')
    # cNorm = colors.Normalize(vmin=0, vmax=int(nSourcePoints) - 1)
    # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

    t0 = time.time()

    bundleSegInd = 1
    sourceInd = 0
    sourcesFinishedFlag = False
    lengthAlongBundle = 0
    while bundleSegInd < bundleCoords.shape[
        0] and not sourcesFinishedFlag:  # bundleSegInd points to endpoint of current bundle segment

        bundleSegEndPoint = bundleCoords[bundleSegInd, :]
        bundleSegStartPoint = bundleCoords[bundleSegInd - 1, :]

        # find the normal to the terminating surface of the bundle guide segment
        dir0 = bundleSegEndPoint - bundleSegStartPoint
        dir0norm = dir0 / np.linalg.norm(dir0)
        if bundleSegInd < bundleCoords.shape[0] - 1:
            dir1 = bundleCoords[bundleSegInd + 1, :] - bundleSegEndPoint
            n = (dir0 + dir1) / 2
        else:
            n = dir0

        # process axons associated with the current bundle segment
        # and store axon segment properties in these lists
        sourceDir2Ds = []
        sourceXDists = []
        sourceDistToOrigin = []

        sourcePosition = sourcePositions[sourceInd,:]  # why oh why, no do-while
        # print 'bundleSegEndPoint_ ' + str(bundleSegEndPoint)
        while np.inner(sourcePosition - bundleSegEndPoint, n) < 0:  # while in bundle guide segment
            # print 'sourcePosition ' + str(sourcePosition)

            sourcePosition = sourcePositions[sourceInd,:]
            sourcePosPar = np.inner(sourcePosition - bundleSegStartPoint,
                                  dir0norm)  # compontent parallel to bundle direction

            # normal from bundle guide to axon position
            sourceDir2D = sourcePosition - (bundleSegStartPoint + sourcePosPar * dir0norm)
            sourceXDist = np.linalg.norm(sourceDir2D)
            if not sourceXDist == 0:
                sourceDir2DNorm = sourceDir2D / sourceXDist
            else:
                sourceDir2DNorm = -1

            # save all computed values of axon into these lists, they are used when iterating over electrodes
            sourceXDists.append(sourceXDist)
            sourceDir2Ds.append(sourceDir2DNorm)
            sourceDistToOrigin.append(sourcePosPar + lengthAlongBundle)

            sourceInd += 1

            if sourceInd >= nSourcePoints:
                sourcesFinishedFlag = True
                break


        # now process receivers

        for sourceIndInner,sourceDir2D in enumerate(sourceDir2Ds):

            # TODO: x,y need to be set correctly, this is only a simplified implementation valid because no side
            # TODO: displacement of axon segments is modeled
            receiverX = receiverXDists
            receiverY = np.zeros(np.shape(receiverXDists))
            receiverZ = np.zeros(np.shape(receiverXDists))
            sourceXDist = sourceXDists[sourceIndInner]
            sourceZDist = np.multiply(receiverDistToOrigins,-1) + sourceDistToOrigin[sourceIndInner] + receiverDisplacement

            interpolationPoints = np.vstack([receiverX, receiverY, receiverZ, np.tile(sourceXDist, np.shape(receiverY)), sourceZDist])
            interpolationPoints = np.divide(interpolationPoints,
                                            1000000)  # from um to m TODO: numercal problem?

            # now interpolate from fieldImage
            receiverPotTempStatic = _interpolateFromImage(fieldDict, interpolationPoints, order=1)

            imemAxonSegInd = sourceInd - len(sourceDir2Ds) + sourceIndInner

            # scale potential-voltage-relation with current to obtain temporal signal
            # COMSOL gave V, we need mV, therefore multiply with 1000
            # also there can be a mismatch in current unit of the source, eliminate
            receiverPotTemp = np.outer(receiverPotTempStatic,
                                       sourceCurrents[imemAxonSegInd, :]
                                       * (10) ** (currentUnitSource - currentUnitFEM)) * 1000

        # # compontent parallel to bundle direction
        # bundleSegStartPointTiled = np.tile(bundleSegStartPoint, (receiverPositions.shape[0], 1))
        # receiverPosPar = np.inner(receiverPositions - bundleSegStartPointTiled, dir0norm)
        #
        # # electrode shifted to origin
        # receiverVector = receiverPositions - \
        #                  (bundleSegStartPointTiled + np.tile(dir0norm,(receiverPositions.shape[0], 1)) * receiverPosPar[:,np.newaxis])

        # for sourceLoc in range(len(sourceDir2Ds)):
        #     sourceDir2D = sourceDir2Ds[sourceLoc]
        #
        #     if isinstance(sourceDir2D, int):  # if the axon segment lies on the bundle middle exactly
        #         receiverX = np.ones(receiverPositions.shape[0]) * np.linalg.norm(receiverVector[0, :])
        #         receiverY = np.zeros(receiverPositions.shape[0])
        #     else:
        #         sourceDir2DTiled = np.tile(sourceDir2D, (receiverVector.shape[0], 1))
        #
        #         # electrode coordinates projected onto new base vectors
        #         receiverX = np.inner(receiverVector, sourceDir2D)
        #
        #         receiverYVec = receiverVector - sourceDir2DTiled * receiverX[:, np.newaxis]
        #         receiverY = np.linalg.norm(receiverYVec, axis=1)
        #
        #     sourceDistToOrigin = sourceZs[sourceLoc]
        #     receiverZ = np.add(receiverPosPar, -sourceZ)
        #
        #     sourceXDist = np.tile(sourceXDists[sourceLoc], (1, len(receiverZ)))
        #
        #     interpolationPoints = np.vstack([receiverX, receiverY, receiverZ, sourceXDist])
        #     interpolationPoints = np.divide(interpolationPoints,
        #                                     1000000)  # from um to m TODO: numercal problem?
        #
        #     # now interpolate from fieldImage
        #     receiverPotTempStatic = _interpolateFromImage(fieldDict, interpolationPoints, order=1)
        #
        #     imemAxonSegInd = sourceInd - len(sourceDir2Ds) + sourceLoc
        #
        #     # scale potential-voltage-relation with current to obtain temporal signal
        #     # COMSOL gave V, we need mV, therefore multiply with 1000
        #     # also there can be a mismatch in current unit of the source, eliminate
        #     receiverPotTemp = np.outer(receiverPotTempStatic,
        #                            sourceCurrents[imemAxonSegInd, :]
        #                            * (10)**(currentUnitSource-currentUnitFEM)) * 1000
        #
        #     # import matplotlib.pyplot as plt
        #     # plt.plot(receiverPotTemp.T)
        #     # plt.show()
        #
            # add contributions
            receiverPotentials = np.add(receiverPotentials, receiverPotTemp)

        # measure length along bundle to have the absolute position along it for all segments.
        lengthAlongBundle += np.linalg.norm(bundleSegEndPoint - bundleSegStartPoint)

        # go to next segment
        bundleSegInd += 1


    # import matplotlib.pyplot as plt
    # plt.plot(receiverPotentials.T, linewidth=2, color='r')
    # plt.show()

    return receiverPotentials

def compute_relative_positions_and_interpolate(sourcePositions, sourceCurrents, receiverPositions, fieldDict, bundleGuide, currentUnitFEM=-9, currentUnitSource=-9):

    """

    Structure as follows

    for bs in bundleSegments:
        for s in sources on this bundle segment:
            calculate distance of s from axon guide and distance along axon guide
        for r in all receivers calculate:
            for s in sources on this bundle segment:
                  calculate receiver position and interpolate potential based on spatial relation between source
                  and receiver and source and bundle

    Args:
        axon: needed for segment positions and membrane currents

    Returns:

    """

    nSourcePoints = np.shape(sourcePositions)[0]


    bundleCoords = bundleGuide[:, 0:3]

    # first go through all bundle segments, find the associated source positions and calculate the needed
    # quantities for all sources
    receiverPotentials = np.zeros((receiverPositions.shape[0], np.shape(sourceCurrents)[1]))

    # # TODO: delete this again, only for testing
    # lastElectrodeSignal = []
    # f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    # jet = plt.get_cmap('jet')
    # cNorm = colors.Normalize(vmin=0, vmax=int(nSourcePoints) - 1)
    # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

    t0 = time.time()

    bundleSegInd = 1
    sourceInd = 0
    sourcesFinishedFlag = False
    while bundleSegInd < bundleCoords.shape[
        0] and not sourcesFinishedFlag:  # bundleSegInd points to endpoint of current bundle segment

        bundleSegEndPoint = bundleCoords[bundleSegInd, :]
        bundleSegStartPoint = bundleCoords[bundleSegInd - 1, :]

        # find the normal to the terminating surface of the bundle guide segment
        dir0 = bundleSegEndPoint - bundleSegStartPoint
        dir0norm = dir0 / np.linalg.norm(dir0)
        if bundleSegInd < bundleCoords.shape[0] - 1:
            dir1 = bundleCoords[bundleSegInd + 1, :] - bundleSegEndPoint
            n = (dir0 + dir1) / 2
        else:
            n = dir0

        # # TODO: maybe this introduces problems?
        # randomAxonXDir = random_perpendicular_vectors(dir0norm)[0,:]

        # process axons associated with the current bundle segment
        # and store axon segment properties in these lists
        sourceDir2Ds = []
        sourceXDists = []
        sourceZs = []

        sourcePosition = sourcePositions[sourceInd,:]  # why oh why, no do-while
        # print 'bundleSegEndPoint_ ' + str(bundleSegEndPoint)
        while np.inner(sourcePosition - bundleSegEndPoint, n) < 0:  # while in bundle guide segment
            # print 'sourcePosition ' + str(sourcePosition)

            sourcePosition = sourcePositions[sourceInd,:]
            sourcePosPar = np.inner(sourcePosition - bundleSegStartPoint,
                                  dir0norm)  # compontent parallel to bundle direction

            # normal from bundle guide to axon position
            sourceDir2D = sourcePosition - (bundleSegStartPoint + sourcePosPar * dir0norm)
            sourceXDist = np.linalg.norm(sourceDir2D)
            if not sourceXDist == 0:
                sourceDir2DNorm = sourceDir2D / sourceXDist
            else:
                sourceDir2DNorm = -1

            # save all computed values of axon into these lists, they are used when iterating over electrodes
            sourceXDists.append(sourceXDist)
            sourceDir2Ds.append(sourceDir2DNorm)
            sourceZs.append(sourcePosPar)

            sourceInd += 1

            if sourceInd >= nSourcePoints:
                sourcesFinishedFlag = True
                break

        # now process receivers

        # compontent parallel to bundle direction
        bundleSegStartPointTiled = np.tile(bundleSegStartPoint, (receiverPositions.shape[0], 1))
        receiverPosPar = np.inner(receiverPositions - bundleSegStartPointTiled, dir0norm)

        # electrode shifted to origin
        receiverVector = receiverPositions - (bundleSegStartPointTiled + np.tile(dir0norm,
                                                                                  (
                                                                                      receiverPositions.shape[
                                                                                       0], 1)) * receiverPosPar[
                                                                                                 :,
                                                                                                 np.newaxis])

        for sourceLoc in range(len(sourceDir2Ds)):
            sourceDir2D = sourceDir2Ds[sourceLoc]

            # calculate spatial relation between source and receiver. First the angle. Projection of source position
            # onto bundle segment is new 'x-axis'. Vector perpendicular this projection and bundle guide segment is 2nd
            # base vector of the considered orthonormal base defined by the source position. Bundle can be turned
            # because of circular symmetry.
            if isinstance(sourceDir2D, int):  # if the axon segment lies on the bundle middle exactly
                receiverX = np.ones(receiverPositions.shape[0]) * np.linalg.norm(receiverVector[0, :])
                receiverY = np.zeros(receiverPositions.shape[0])
            else:
                sourceDir2DTiled = np.tile(sourceDir2D, (receiverVector.shape[0], 1))

                # electrode coordinates projected onto new base vectors
                receiverX = np.inner(receiverVector, sourceDir2D)

                receiverYVec = receiverVector - sourceDir2DTiled * receiverX[:, np.newaxis]
                receiverY = np.linalg.norm(receiverYVec, axis=1)

            # why are only sources on the same bundle segment as receiver considered?!
            sourceZ = sourceZs[sourceLoc]
            receiverZ = np.add(receiverPosPar, -sourceZ)

            sourceXDist = np.tile(sourceXDists[sourceLoc], (1, len(receiverZ)))

            interpolationPoints = np.vstack([receiverX, receiverY, receiverZ, sourceXDist])
            interpolationPoints = np.divide(interpolationPoints,
                                            1000000)  # from um to m TODO: numercal problem?

            # now interpolate from fieldImage
            receiverPotTempStatic = _interpolateFromImage(fieldDict, interpolationPoints, order=1)

            imemAxonSegInd = sourceInd - len(sourceDir2Ds) + sourceLoc

            # scale potential-voltage-relation with current to obtain temporal signal
            # COMSOL gave V, we need mV, therefore multiply with 1000
            # also there can be a mismatch in current unit of the source, eliminate
            receiverPotTemp = np.outer(receiverPotTempStatic,
                                   sourceCurrents[imemAxonSegInd, :]
                                   * (10)**(currentUnitSource-currentUnitFEM)) * 1000

            # import matplotlib.pyplot as plt
            # plt.plot(receiverPotTemp.T)
            # plt.show()

            # add contributions
            receiverPotentials = np.add(receiverPotentials, receiverPotTemp)

        bundleSegInd += 1

    # import matplotlib.pyplot as plt
    # plt.plot(receiverPotentials.T, linewidth=2, color='r')
    # plt.show()

    return receiverPotentials


def i_to_v_homogeneous(sourcePositions, sourceCurrents, receiverPositions, sigma=1., currentUnitSource=-9):
    """
    Idea and some implementation details from LFPy package

    Args:
        sourcePositions:
        sourceCurrents:
        receiverPositions:
        sigma:
        currentUnitSource:

    Returns:

    """

    # import matplotlib.pyplot as plt
    # plt.plot(sourceCurrents)
    # plt.show()

    nSourcePoints = np.shape(sourcePositions)[0]
    nReceiverPoints = np.shape(receiverPositions)[0]

    nTimePoints = len(sourceCurrents[:,0])

    receiverPotentials = []
    for rInd in range(nReceiverPoints):
        receiverPosition = receiverPositions[rInd,:]

        r2 = (sourcePositions[:,0] - receiverPosition[0]) ** 2 + (sourcePositions[:,1] - receiverPosition[1]) ** 2 + (sourcePositions[:,2] - receiverPosition[2]) ** 2
        r = np.sqrt(r2)

        receiverPotential = 1 / (4 * np.pi * sigma) * np.dot(sourceCurrents.T, 1 / r)

        receiverPotentials.append(receiverPotential)

    receiverPotentials = np.array(receiverPotentials)

    return receiverPotentials