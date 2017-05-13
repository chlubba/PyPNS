from math import pi, tan, cos, sin, sqrt, floor
import numpy as np

def random_perpendicular_vectors(v):
    # adapted from: http://codereview.stackexchange.com/questions/43928/algorithm-to-get-an-arbitrary-perpendicular-vector
    # user: Jaime
    if v[0] == 0 and v[1] == 0:
        if v[2] == 0:
            raise ValueError('zero vector')
        # v is Vector(0, 0, v.z)
        v1 = np.array([0, 1, 0])

    v1 = np.array([-v[1], v[0], 0])
    v2 = np.cross(v,v1)

    randomAngle = np.random.uniform(0,2*pi,1)
    vRand1 = v1*cos(randomAngle) + v2*sin(randomAngle)
    vRand2 = v1*cos(randomAngle+pi/2) + v2*sin(randomAngle+pi/2)

    vRandNorm1 = vRand1/np.linalg.norm(vRand1)
    vRandNorm2 = vRand2/np.linalg.norm(vRand2)

    return np.row_stack((vRandNorm1,vRandNorm2))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    From: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    User: unutbu

    :param axis: rotation axis
    :param theta: rotation angle

    :return: rotation matrix
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis /= sqrt(np.dot(axis, axis))
    a = cos(theta/2.0)
    b, c, d = -axis*sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def create_random_axon(bundleCoords4D, axonCoords, segmentLengthAxon, randomDirectionComponent=0, distribution='normal'): # maximumAngle=pi / 10,
    """
    This function is used to generate the trajectory of axons. They can follow the bundle trajectory more or less loosely, as set by the parameter ``randomDirectionComponent``.

    :param bundleCoords4D: Coordinates of the bundle trajectory, 4th coordinate is radius
    :param axonCoords: initial position of the axon in the y-z-plane at x=0
    :param segmentLengthAxon: length of each straight step of the axon trajectory
    :param randomDirectionComponent: (0 - 1) Regulates the totuosity of the axon. If 0, axon follows bundle direction completely exactly. If 1, axon will only follow the bundle direction very loosely and get very long.

    :return: axon coordinates
    """
    bundleCoords = bundleCoords4D[:, :3]

    # take first two bundle segment start points as first two axon points
    pos1 = np.concatenate(([bundleCoords[0, 0]], axonCoords + bundleCoords[0, 1:3]))
    pos2 = np.concatenate(([bundleCoords[1, 0]], axonCoords + bundleCoords[1, 1:3]))

    coords = np.row_stack((pos1, pos2))

    # rhoMax = tan(maximumAngle) * segmentLengthAxon
    # rhoArray = np.zeros(5)

    bundleLength = np.shape(bundleCoords)[0]
    currentBundleSegment = 2
    while currentBundleSegment <= bundleLength - 1:

        # get last point
        lastPointAxon = coords[-1, :]
        # and vector between two last points
        lastAxonDirection = (coords[-1, :] - coords[-2, :])
        lastAxonDirectionNorm = lastAxonDirection / np.linalg.norm(lastAxonDirection)

        # current bundle direction
        bundleDirection = bundleCoords[currentBundleSegment, :] - bundleCoords[currentBundleSegment - 1, :]
        lastPointBundle = bundleCoords[currentBundleSegment - 1, :]
        bundleDirectionNorm = bundleDirection / np.linalg.norm(bundleDirection)

        # get orthogonal vector to current direction vector
        cp = np.inner(bundleDirectionNorm, lastPointAxon - lastPointBundle)
        if cp == 0:
            radiusVectorNorm = random_perpendicular_vectors(bundleDirection)[0,:] # [0, 0, 0]
            distance = 0
        else:
            radiusVector = -(lastPointAxon - (cp * bundleDirectionNorm + lastPointBundle))
            distance = np.linalg.norm(radiusVector)
            if not distance == 0:
                radiusVectorNorm = radiusVector / distance
            else:
                radiusVectorNorm = radiusVector

        # assure axon stays within bundle. If too far away -> next direction equals bundle direction
        bundleRadius = bundleCoords4D[currentBundleSegment, 3]

        # calculate the random component perpendicular to the bundle direction. If axon approaches bundle limit (radius)
        # the random component faces inwards, towards the bundle core.
        factorRadiusBias = min((max(0, distance / bundleRadius - 0.7)) * 6, 2.5)
        radialNorm = np.cross(radiusVectorNorm, bundleDirectionNorm)
        if distribution == 'uniform':
            randFactors = np.random.uniform(-1,1,2)
        elif distribution == 'normal':
            randFactors = np.random.normal(0, 1, 2)/3
        else:
            raise NameError('No valid distribution name given. Either ''uniform'' or ''normal''.')
        randomVectorToAdd = randFactors[0]*radialNorm + (randFactors[1] - factorRadiusBias)*np.array(radiusVectorNorm*(-1))
        if not np.sum(randomVectorToAdd) == 0:
            randomVectorToAddNorm = randomVectorToAdd/np.linalg.norm(randomVectorToAdd)
        else:
            randomVectorToAddNorm = randomVectorToAdd

        # direction of next axon segment is the sum of previous direction, bundle direction and random component.
        # The higher randomDirectionComponent, the more weight is put on the random vector and less weight on the
        # bundle direction.
        nextAxonDir = randomVectorToAddNorm*randomDirectionComponent + lastAxonDirectionNorm + bundleDirectionNorm*(1.1-randomDirectionComponent)
        nextAxonDirScaled = nextAxonDir*segmentLengthAxon
        nextPoint = lastPointAxon + nextAxonDirScaled

        # addpend to coordinate list
        coords = np.row_stack((coords, nextPoint))

        if np.inner(bundleDirection, (nextPoint - bundleCoords[currentBundleSegment, :])) > 0:
            currentBundleSegment += 1

    return coords

def length_from_coords(coords):
    """Calculates the length of a curve defined by the input coordinates.

    :param coords: coordinates

    :return: length curve

    """
    # do calculate that by summing over lenghts of segments, calculate the difference in coords between each consecutive
    # pair of segments
    dCoords = np.diff(coords,axis=0)

    # pythagoras
    radicand = np.sum(np.power(dCoords,2), axis=1)
    dL = np.sqrt(radicand)

    # sum over all segments
    return sum(dL)

def distance_along_bundle(bundleGuide, bundleLength, positionMax):
    # TODO: is this function useless?!
    #
    # A bundle is always a little longer than specified in order to complete myelinated axons that have large internodal
    # distances for higher diameters. This function finds the index of the bundle guide segment that
    #
    # Args:
    #     bundleGuide:
    #     bundleLength:
    #     positionMax:
    #
    # Returns:
    #
    #

    bundleGuide = bundleGuide[:, 0:3]

    # first find the bundle guide segment index that corresponds to the intendet bundle length (overlap for
    # myelinated axons gives longer bundle than specified by user)
    bundleLengthIndex = np.shape(bundleGuide)[0]-1
    bundleLengthTemp = length_from_coords(bundleGuide)
    while bundleLengthTemp > bundleLength:
        bundleLengthIndex -= 1
        bundleLengthTemp = length_from_coords(bundleGuide[:bundleLengthIndex])

    lastRecordedSegmentIndex = bundleLengthIndex*positionMax

    electrodeDistance = np.floor(length_from_coords(bundleGuide[:lastRecordedSegmentIndex]))

    return electrodeDistance


def circular_electrode(bundleGuide, positionAlongBundle, radius, numberOfPoles, poleDistance=None, numberOfPoints=20):
    """Calculate the set of electrode coordinates for a circular electrode. Can be used for ``RecordingMechanism`` or ``StimFieldQuasistatic``.

    :param bundleGuide: trajectory of the bundle
    :param positionAlongBundle: distance from origin of the bundle along the trajectory
    :param radius: radius of the ring
    :param numberOfPoles: number of poles (rings)
    :param poleDistance: distance between poles
    :param numberOfPoints: number of points per ring

    :return: 3D matrix of electrode coordinates

    """

    bundleGuide = bundleGuide[:, 0:3]

    if not numberOfPoles == 1:
        assert not poleDistance == None

    # first find the bundle guide segment index that corresponds to the intended bundle length (overlap for
    # myelinated axons gives longer bundle than specified by user)
    segmentIndex = np.shape(bundleGuide)[0]-1
    distanceTemp = length_from_coords(bundleGuide)
    while distanceTemp > positionAlongBundle:
        segmentIndex -= 1
        distanceTemp = length_from_coords(bundleGuide[:segmentIndex])

    # variable to save points of electrode
    # electrodePositions = np.squeeze(np.array([]).reshape(0, 3, numberOfPoles))
    electrodePositions = np.array([]).reshape(0, 3)

    # get the geometry of the segment, position and orientation.
    segmentStartingPos = bundleGuide[segmentIndex - 1, :]
    segmentEndPos = bundleGuide[segmentIndex, :]
    segmentMiddle = (segmentStartingPos + segmentEndPos) / 2

    segmentOrientation = segmentEndPos - segmentStartingPos
    segmentOrientation /= np.linalg.norm(segmentOrientation)

    electrodeProjectionOnBundleGuide = segmentStartingPos + (positionAlongBundle - distanceTemp)*segmentOrientation

    # get one random orthogonal vector
    orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0, :]

    # loop to generate one ring
    for j in range(numberOfPoints):
        # generate the coordinates for one ring for the first pole of the electrode
        pointPosition = np.dot(rotation_matrix(segmentOrientation, 2 * np.pi / numberOfPoints * j),
                               (orthogonalVector * radius)) + electrodeProjectionOnBundleGuide

        # append it to the list of coordinates for this pole
        electrodePositions = np.vstack([electrodePositions, pointPosition])

    # add axis for poles
    electrodePositions = np.expand_dims(electrodePositions, axis=2)

    # if the electrodes are bipolar
    for i in range(1,numberOfPoles):
        electrodePositionsPole = electrodePositions[:,:,0] + np.tile(segmentOrientation * poleDistance*i, (
        np.shape(electrodePositions)[0], 1))
        electrodePositionsPole = np.expand_dims(electrodePositionsPole, axis=2)
        electrodePositions = np.concatenate((electrodePositions, electrodePositionsPole), axis=2)

    return electrodePositions


def get_bundle_guide_corner(bundleLength, segmentLengthAxon, overlapLength=1000, lengthFactor=3):
    """Generate a bundle trajectory with a corner (only a demonstration).

    :param bundleLength: length of bundle
    :param segmentLengthAxon: length of straight axon segment. Important for this function as the bundle segment length needs to be longer than the axon segment length for ``create_random_axon`` to work properly.
    :param lengthFactor: (>2) factor (bundle segment length) / (axon segment length)
    :param overlapLength: additional length of the bundle trajectory to finish myelinated axons.

    :return: bundle trajectory coordinates. 3 coordinates, no radius

    """

    #length after bundle end. necessary for myelinated axons
    bundleLength += overlapLength

    segmentLengthBundle = segmentLengthAxon*lengthFactor

    numBundleGuideSteps = int(np.floor(bundleLength/segmentLengthBundle))

    cornerIndex = float(numBundleGuideSteps)/5
    turningPointIndex1 = int(np.floor(cornerIndex))
    turningPointIndex2 = (numBundleGuideSteps - turningPointIndex1)# + int(np.ceil(cornerIndex - turningPointIndex1))

    bundleCoords = np.zeros([numBundleGuideSteps, 3])
    bundleCoords[:,0] = range(0, numBundleGuideSteps*segmentLengthBundle, segmentLengthBundle)
    bundleCoords[:,1] = np.concatenate((np.zeros(turningPointIndex1),np.multiply(range(turningPointIndex2), segmentLengthBundle)))
    bundleCoords[:,2] = np.concatenate((np.zeros(turningPointIndex1),np.multiply(range(turningPointIndex2), segmentLengthBundle)))

    return bundleCoords

def get_bundle_guide_random(bundleLength, segmentLength = 200, overlapLength=1000):
    """Generate a random bundle trajectory.

    :param bundleLength: length of bundle
    :param segmentLength: length of straight bundle segment
    :param overlapLength: additional length of the bundle trajectory to finish myelinated axons.

    :return: coordinates of random bundle. 3 coordinates, no radius

    """

    bundleLength += overlapLength

    numBundleGuideSteps = int(np.floor(bundleLength/segmentLength))

    randomDeltaX = np.random.uniform(0,2,numBundleGuideSteps)
    randomDeltaYZ = np.random.uniform(-1,1,(numBundleGuideSteps,2))
    randomDelta = np.column_stack((randomDeltaX, randomDeltaYZ))

    for i in range(numBundleGuideSteps):
        randomDelta[i,:] = randomDelta[i,:]/np.linalg.norm(randomDelta[i,:])

    bundleGuide = np.cumsum(randomDelta,0)

    bundleGuideScaled = bundleGuide*segmentLength

    return bundleGuideScaled

def get_bundle_guide_straight(bundleLength, segmentLengthAxon, overlapLength=1000):
    """Generate a straight bundle. Default if no bundle guide is provided.

    :param bundleLength: length of bundle
    :param segmentLengthAxon: length of straight axon segment. Important for this function as the bundle segment length needs to be longer than the axon segment length for ``create_random_axon`` to work properly.
    :param overlapLength: additional length of the bundle trajectory to finish myelinated axons.

    :return: coordinates of a straight bundle. 3 coordinates, no radius

    """

    #length after bundle end. necessary for myelinated axons
    bundleLength += overlapLength

    segmentLengthBundle = segmentLengthAxon*3
    numBundleGuideSteps = int(np.floor(bundleLength/segmentLengthBundle))

    bundleCoords = np.zeros([numBundleGuideSteps, 3])
    bundleCoords[:,0] = range(0, numBundleGuideSteps*segmentLengthBundle, segmentLengthBundle)

    return bundleCoords

def get_bundle_guide_straight_radius(bundleLength, segmentLengthAxon, overlapLength=1000, radius=200):
    """Like ``get_bundle_guide_straight`` but with a radius as 4th coordinate.

    :param bundleLength: length of bundle
    :param segmentLengthAxon: length of straight axon segment. Important for this function as the bundle segment length needs to be longer than the axon segment length for ``create_random_axon`` to work properly.
    :param overlapLength: additional length of the bundle trajectory to finish myelinated axons.
    :param radius: radius of the nerve

    :return: coordinates of a straight bundle. Nx4, 4th coordinate for radius.

    """

    #length after bundle end. necessary for myelinated axons
    bundleLength += overlapLength

    segmentLengthBundle = segmentLengthAxon*3
    numBundleGuideSteps = int(np.floor(bundleLength/segmentLengthBundle))

    bundleCoords = np.zeros([numBundleGuideSteps, 4])
    bundleCoords[:,0] = range(0, numBundleGuideSteps*segmentLengthBundle, segmentLengthBundle)
    bundleCoords[:,-1] = np.ones(numBundleGuideSteps)*radius

    return bundleCoords

def get_bundle_guide_straight_2radii(bundleLength, segmentLengthAxon, overlapLength=1000, radii=(150, 150)):
    """Demo of a bundle with radius varying linearly along bundle length.

    :param bundleLength: length of bundle
    :param segmentLengthAxon: length of straight axon segment. Important for this function as the bundle segment length needs to be longer than the axon segment length for ``create_random_axon`` to work properly.
    :param overlapLength: additional length of the bundle trajectory to finish myelinated axons.
    :param radii: start and end radius of the nerve

    :return: bundle coordinates including radius

    """

    #length after bundle end. necessary for myelinated axons
    bundleLength += overlapLength

    segmentLengthBundle = segmentLengthAxon*3
    numBundleGuideSteps = int(np.floor(bundleLength/segmentLengthBundle))

    bundleCoords = np.zeros([numBundleGuideSteps, 4])
    bundleCoords[:,0] = range(0, numBundleGuideSteps*segmentLengthBundle, segmentLengthBundle)
    bundleCoords[:,-1] = np.linspace(radii[0], radii[1], numBundleGuideSteps)

    return bundleCoords

def get_bundle_guide_random_radius(bundleLength, segmentLength = 200, overlapLength=1000, radius=200):
    """Generate a random bundle trajectory.

    :param bundleLength: length of bundle
    :param segmentLength: length of straight bundle segment
    :param overlapLength: additional length of the bundle trajectory to finish myelinated axons.

    :return: coordinates of random bundle. 3 coordinates, no radius

    """

    bundleLength += overlapLength

    numBundleGuideSteps = int(np.floor(bundleLength/segmentLength))

    randomDeltaX = np.random.uniform(0,2,numBundleGuideSteps)
    randomDeltaYZ = np.random.uniform(-1,1,(numBundleGuideSteps,2))
    randomDelta = np.column_stack((randomDeltaX, randomDeltaYZ))

    for i in range(numBundleGuideSteps):
        randomDelta[i,:] = randomDelta[i,:]/np.linalg.norm(randomDelta[i,:])

    bundleGuide = np.cumsum(randomDelta,0)

    bundleGuideScaled = bundleGuide*segmentLength

    bundleGuideScaled = np.concatenate((bundleGuideScaled, np.expand_dims(np.ones(numBundleGuideSteps)*radius,axis=1)), axis=1)

    return bundleGuideScaled
