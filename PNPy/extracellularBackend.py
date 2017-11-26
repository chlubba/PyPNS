import numpy as np
from scipy import ndimage
import time
import createGeometry
# from takeTime import takeTime
# from scipy.interpolate import interp1d

# ================================================================================
# ====================== used by precomputedFEM to calculate =====================
# =================== image coordinates from spatial coordinates =================
# ================================================================================

def _getImageCoords(fieldDict, points):
    """
    This function transforms the coordinate values in points to positions on the field image. It assumes equidistant
    coordinate steps in fieldDict.

    Args:
        fieldDict: dictionary containing the coordinate values with their associated voltage values
        points: coordinate values this function transforms to image coordinates

    Returns: image coordinate values

    """

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

def interpolateFromImage(fieldDict, points, order=3, zGiven=False):

    # first transform coordinates in points into position coordinates
    # different function
    imageCoords = _getImageCoords(fieldDict, points)

    # then with new coords to the interpolation
    return ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

# ================================================================================
# ========= used by interpolator and precomputedFEM to calculate =================
# ======== the spatial relationships between points and the bundle ===============
# ================================================================================

def associatePointToBundleSegs(points, bundleGuide):

    # make sure orientation is fine
    pointsShape = np.shape(points)
    if pointsShape[0] == 3:
        points = np.transpose(points)
    numPoints = np.shape(points)[0]

    # cut away potential radius coordinate
    bundleGuideStripped = bundleGuide[:,0:3]
    bundleGuideLength = np.shape(bundleGuide)[0]

    # variable to save squared distances between bundle guide segment centers and source positions
    r2min = np.ones(numPoints) * np.inf
    # indices of closest bundle segment for all points
    closestSegInds = np.zeros(numPoints)

    for bundleGuideInd in range(bundleGuideLength - 1):
        bundleSegStart = bundleGuideStripped[bundleGuideInd, :]
        bundleSegStop = bundleGuideStripped[bundleGuideInd + 1, :]
        bundleMiddle = (bundleSegStart + bundleSegStop) / 2

        r2 = np.sum(np.square(points - bundleMiddle), axis=1)

        # current distance smaller than current minimum?
        compArray = r2 < r2min

        closestSegInds[compArray] = bundleGuideInd
        r2min[compArray] = r2[compArray]

    return closestSegInds

def rotationMatrixFromVectors(a, b):
    """np.dot(R,a) = b

    from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    user 'Kuba Ober'

    """

    if np.all(a == b):
        return np.diag(np.ones(3))
    else:

        G = np.array([[np.dot(a, b), -np.linalg.norm(np.cross(a, b)), 0],
                      [np.linalg.norm(np.cross(a, b)), np.dot(a, b), 0],
                      [0, 0, 1]])
        F = np.array(
            [a, (b - np.multiply(np.dot(a, b), a)) / np.linalg.norm(b - np.multiply(np.dot(a, b), a)), np.cross(b, a)])

        R = F.dot(G).dot(np.linalg.inv(F))

        return R

def spatialRelation(points, bundleGuide, segmentAssociations):
    """

    Args:
        points: points of interest (Nx3), can be electrode or axon segments
        bundleGuide: ... bundle guide.
        segmentAssociations: an array of length N (number of points), containing the bundle guide segment indices

    Returns:
        distPerpendiculars: distance of points from bundle guide
        lengthAlongs: distance along the bundle guide (from origin)
        angles: the angle between the perpendicular of the points towards the bundle guide segment and the y-axis (chosen arbitrarily)

    """

    # cut away potential radius coordinate
    bundleGuideStripped = bundleGuide[:, 0:3]
    bundleGuideLength = np.shape(bundleGuide)[0]

    # variables of interest that are calculated for all positions
    posPerpendiculars = np.zeros(np.shape(points))
    distPars = np.zeros(np.shape(points)[0])
    distPerpendiculars = np.zeros(np.shape(points)[0])
    lengthAlongs = np.zeros(np.shape(points)[0])
    angles = np.zeros(np.shape(points)[0])

    # loop over all bundle segments
    lengthAlongBundle = 0
    lastBundleDir = bundleGuideStripped[1, :] - bundleGuideStripped[0, :]
    lastBundleDirNorm = lastBundleDir / np.linalg.norm(lastBundleDir)
    overallR = np.diag(np.ones(3))
    for bundleSegInd in range(bundleGuideLength-1):

        bundleSegStart = bundleGuideStripped[bundleSegInd, :]
        bundleSegStop = bundleGuideStripped[bundleSegInd + 1, :]
        bundleSegLen = np.linalg.norm(bundleSegStop - bundleSegStart)
        bundleDirNorm = (bundleSegStop - bundleSegStart)/bundleSegLen

        # calulcate the rotation matrix between two following bundle segments
        R = rotationMatrixFromVectors(bundleDirNorm, lastBundleDirNorm)

        # overall rotation matrix from current segment to initial one
        overallR = overallR.dot(R)

        # look at points associated with current bundle segment
        pointIndicesCurrentSegment = np.where(segmentAssociations == bundleSegInd)[0]
        for pointInd in pointIndicesCurrentSegment:
            point = points[pointInd,:]

            # compontent parallel to bundle direction
            distPar = np.inner(point - bundleSegStart, bundleDirNorm)

            # normal from bundle guide to axon position
            posPerpendicular = point - (bundleSegStart + distPar * bundleDirNorm)
            distPerpendicular = np.linalg.norm(posPerpendicular)

            # save all computed values of axon into these lists
            distPerpendiculars[pointInd] = distPerpendicular
            posPerpendiculars[pointInd, :] = posPerpendicular
            distPars[pointInd] = distPar

            # distance from origin of bundle
            lengthAlongs[pointInd] = lengthAlongBundle + distPar

            # vector perpendicular onto bundle segment, turned to calculated the
            yzVec = np.dot(overallR, posPerpendicular)
            angle = np.arctan2(yzVec[1], yzVec[2])
            angles[pointInd] = angle

        lengthAlongBundle += bundleSegLen
        lastBundleDirNorm = bundleDirNorm

    return distPerpendiculars, lengthAlongs, angles

def compute_relative_positions_and_interpolate_fn_input(sourcePositions, sourceCurrents, receiverPositions, bundleGuide,
                                                positionToVoltageFn, currentUnitFEM=-9, currentUnitSource=-9):
    """

    Args:
        sourcePositions: position of current sources (Nx3) x,y,z
        sourceCurrents: source currents (NxM) with M time steps
        receiverPositions: positions of receivers (Nx3) x,y,z
        fieldDict: dictionary to interpolate voltage values with
        bundleGuide: bundle guide, needed to calculate spatial relation between source and recording positions
        currentUnitFEM: unit of current as used in the FEM simulation in powers of 10 (default -9)
        currentUnitSource: unit of current as used in the NEURON simulation in powers of 10 (default -9)

    Returns: (KxM) voltage matrix of K recording positions and M time steps

    """

    # precalculate the spatial relation between the bundle guide and the receivers
    segmentAssociationsRec = associatePointToBundleSegs(receiverPositions, bundleGuide)
    distPerpendicularRec, lengthAlongRec, anglesRec = spatialRelation(receiverPositions, bundleGuide,
                                                                      segmentAssociationsRec)

    # same for sources
    segmentAssociationsSource = associatePointToBundleSegs(sourcePositions, bundleGuide)
    distPerpendicularsSource, lengthAlongsSource, anglesSource = spatialRelation(sourcePositions, bundleGuide,
                                                                                 segmentAssociationsSource)
    # number of sources
    numSourcePos = np.shape(distPerpendicularsSource)[0]

    # matrix to save receiver potentials in
    receiverPots = np.array([]).reshape(0, np.shape(sourceCurrents)[1])

    # loop over all recording positions
    for recInd in range(np.shape(distPerpendicularRec)[0]):
        # distance to the bundle guide
        distPerpendicularRecTemp = distPerpendicularRec[recInd]
        # distance of recording along axon (important: this 'straightens' the bundle) No curves of the
        # bundle guide are considered. Axon tortuosity still plays a role.
        lengthAlongRecTemp = lengthAlongRec[recInd]
        # angle between y-axis and recording position perpendicular towards bundle guide
        # (bundle guide segment rotated to x-axis)
        angleRecTemp = anglesRec[recInd]

        # distance between source and recording positions
        distBetweenTemp = lengthAlongsSource - lengthAlongRecTemp

        # angle between them
        anglesTemp = anglesSource - angleRecTemp

        # calculate the interpolation points handed over to the fieldImage
        interpolationPoints = np.vstack(
            [np.cos(anglesTemp) * distPerpendicularRecTemp, np.sin(anglesTemp) * distPerpendicularRecTemp,
             distBetweenTemp, distPerpendicularsSource])
        interpolationPoints = np.divide(interpolationPoints, 1000000)  # from um to m

        # now interpolate from fieldImage
        receiverPotTempStatic = positionToVoltageFn(interpolationPoints)

        # clear to free memory
        interpolationPoints = None

        # scale potential-voltage-relation with current to obtain temporal signal
        # COMSOL gave V, we need mV, therefore multiply with 1000
        # also there can be a mismatch in current unit of the source, eliminate
        receiverPotTemp = np.sum(sourceCurrents * receiverPotTempStatic[:, np.newaxis], axis=0) \
                          * 10 ** (currentUnitSource - currentUnitFEM) * 1000

        receiverPots = np.vstack([receiverPots, receiverPotTemp])

    return receiverPots


