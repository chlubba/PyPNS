import numpy as np
from scipy import ndimage
import time


def _getImageCoords(fieldDict, points):

    xValues = fieldDict['x']
    yValues = fieldDict['y']
    zValues = fieldDict['z']
    axonXValues = fieldDict['axonX']

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

    points = np.array(points)

    if len(points.shape) > 1:
        if points.shape[1] > 4:
            points = np.transpose(points)
        xCoords = np.add(points[:, 0], -xMin) / (xMax - xMin) * (xNum - 1)
        yCoords = np.add(points[:, 1], -yMin) / (yMax - yMin) * (yNum - 1)
        zCoords = np.add(points[:, 2], -zMin) / (zMax - zMin) * (zNum - 1)
        xAxonCoords = np.add(points[:, 3], -axonXMin) / (axonXMax - axonXMin) * (axonXNum - 1)
    else:
        xCoords = (points[0] - xMin) / (xMax - xMin) * (xNum - 1)
        yCoords = (points[1] - yMin) / (yMax - yMin) * (yNum - 1)
        zCoords = (points[2] - zMin) / (zMax - zMin) * (zNum - 1)
        xAxonCoords = (points[3] - axonXMin) / (axonXMax - axonXMin) * (axonXNum - 1)

    zCoords = np.abs(
        zCoords)  # in the input FEM field, we only take one side of the z-value range thanks to symmetry
    yCoords = np.abs(yCoords)  # same for y-Coordinates

    return np.vstack([xCoords, yCoords, zCoords, xAxonCoords])

def _interpolateFromImage(fieldDict, points, order=3):

    # first transform coordinates in points into position coordinates
    imageCoords = _getImageCoords(fieldDict, points)

    # then with new coords to the interpolation
    return ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

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
        while np.inner(sourcePosition - bundleSegEndPoint, n) < 0:  # while in bundle guide segment

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

            if isinstance(sourceDir2D, int):  # if the axon segment lies on the bundle middle exactly
                receiverX = np.ones(receiverPositions.shape[0]) * np.linalg.norm(receiverVector[0, :])
                receiverY = np.zeros(receiverPositions.shape[0])
            else:
                sourceDir2DTiled = np.tile(sourceDir2D, (receiverVector.shape[0], 1))

                # electrode coordinates projected onto new base vectors
                receiverX = np.inner(receiverVector, sourceDir2D)

                receiverYVec = receiverVector - sourceDir2DTiled * receiverX[:, np.newaxis]
                receiverY = np.linalg.norm(receiverYVec, axis=1)

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

            # add contributions
            receiverPotentials = np.add(receiverPotentials, receiverPotTemp)

        bundleSegInd += 1

    # import matplotlib.pyplot as plt
    # plt.plot(elecPotentials.T)
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