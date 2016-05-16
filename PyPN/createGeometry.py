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
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/sqrt(np.dot(axis, axis))
    a = cos(theta/2.0)
    b, c, d = -axis*sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def create_random_axon(bundleCoords, bundleRadius, axonCoords, segmentLengthAxon, maximumAngle = pi/10, randomDirectionComponent=0.1):

    pos1 = np.concatenate(([bundleCoords[0,0]], axonCoords+bundleCoords[0,1:3]))
    pos2 = np.concatenate(([bundleCoords[1,0]], axonCoords+bundleCoords[1,1:3]))

    coords = np.row_stack((pos1, pos2))

    rhoMax = tan(maximumAngle)*segmentLengthAxon
    rhoArray = np.zeros(5)

    bundleLength = np.shape(bundleCoords)[0]
    currentBundleSegment = 2
    while currentBundleSegment <= bundleLength-1:

        # get last point
        lastPointAxon = coords[-1,:]
        # and vector between two last points
        lastAxonDirection = (coords[-1,:] - coords[-2,:])
        lastAxonDirectionNorm = lastAxonDirection/np.linalg.norm(lastAxonDirection)

        # current bundle direction
        bundleDirection = bundleCoords[currentBundleSegment,:] - bundleCoords[currentBundleSegment-1,:]
        lastPointBundle = bundleCoords[currentBundleSegment-1,:]
        bundleDirectionNorm = bundleDirection/np.linalg.norm(bundleDirection)

        # get orthogonal vector to current direction vector
        cp = np.inner(bundleDirectionNorm, lastPointAxon-lastPointBundle)
        if cp == 0:
            radiusVectorNorm = [0, 0, 0]
            distance = 0
        else:
            radiusVector = -(lastPointAxon - (cp*bundleDirectionNorm + lastPointBundle))
            distance = np.linalg.norm(radiusVector)
            if not distance == 0:
                radiusVectorNorm = radiusVector/distance
            else:
                radiusVectorNorm = radiusVector

        # assure axon stays within bundle. If too far away -> next direction
        # equals bundle direction
        factorBundleDirection = min((max(0,distance/bundleRadius-0.7))*6,2.5)
        
        correctionVector = radiusVectorNorm + 0.1*bundleDirectionNorm
        correctionVector = correctionVector/np.linalg.norm(correctionVector)
        combinedDirection = lastAxonDirectionNorm + correctionVector*factorBundleDirection + 0.1*bundleDirection
        combinedDirectionNorm = combinedDirection/np.linalg.norm(combinedDirection)

        # get one random orthogonal vector to desired mean direction of next segment
        randomOrthogonalVectorNorm = random_perpendicular_vectors(combinedDirection)[0,:]

        # select a direction defined by cylindical coordinate rho
        rho = np.random.uniform(1)*rhoMax

        randomDirection = (1-randomDirectionComponent)*combinedDirectionNorm + randomDirectionComponent*randomOrthogonalVectorNorm*rho
        randomDirectionNorm = randomDirection/np.linalg.norm(randomDirection)
        nextDirectionScaled = randomDirectionNorm*segmentLengthAxon

        # add the direction to the last point to obtain the next point
        nextPoint = lastPointAxon + nextDirectionScaled

        # addpend to coordinate list
        coords = np.row_stack((coords,nextPoint))

        if np.inner(bundleDirection,(nextPoint-bundleCoords[currentBundleSegment,:])) > 0:
            currentBundleSegment = currentBundleSegment + 1

    return coords

def length_from_coords(coords):
    # get the length of the wanted axon geometry

    # do calculate that by summing over lenghts of segments, calculate the difference in coords between each consecutive
    # pair of segments
    dCoords = np.diff(coords,axis=0)

    # pythagoras
    radicand = np.sum(np.power(dCoords,2), axis=1)
    dL = np.sqrt(radicand)

    # sum over all segments
    return sum(dL)

def electrode_positions_bundle_guided(bundleGuide, bundleRadius, numberOfElectrodes, numberOfContacts, recElectrodePositions):

    electrodeRadius = bundleRadius*1.2 # allow 20% distance, not to intrude into nerve bundle

    numberOfPoles = len(recElectrodePositions)

    if numberOfPoles == 2:
        poleDistance = recElectrodePositions[1] - recElectrodePositions[0]
        bipolar = True
    elif numberOfPoles == 1:
        poleDistance = 0
        bipolar = False
    else:
        print 'Wrong number of recording poles.'
        return

    positionLastElectrode = recElectrodePositions[0]
    bundleGuideSegmentLength = np.linalg.norm(bundleGuide[1,:] - bundleGuide[0,:]) # assume segments are of equal length

    if numberOfElectrodes > 1:
        segmentIndices = np.linspace(1, positionLastElectrode/bundleGuideSegmentLength, numberOfElectrodes)
    else:
        segmentIndices = [positionLastElectrode/bundleGuideSegmentLength]


    for i in range(numberOfElectrodes):
        segmentNumber = int(segmentIndices[i])
        # segmentNumber = floor(np.shape(bundleGuide)[0]/numberOfElectrodes)*(i+1) - 1

        segmentStartingPos = bundleGuide[segmentNumber - 1,:]
        segmentEndPos = bundleGuide[segmentNumber,:]

        segmentMiddle = (segmentStartingPos + segmentEndPos)/2
        segmentOrientation = segmentEndPos - segmentStartingPos
        segmentOrientation = segmentOrientation/np.linalg.norm(segmentOrientation)

        # get one random orthogonal vector
        orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0,:]

        for j in range(numberOfContacts):
            electrodePositionPole1 = np.dot(rotation_matrix(segmentOrientation, 2*pi/numberOfContacts*j),(orthogonalVector*electrodeRadius)) + segmentMiddle
            if bipolar:
                electrodePositionPole2 = electrodePositionPole1 + segmentOrientation*poleDistance


            if j == 0: # how to get a truly empty array?!
                electrodePositionsPole1 = electrodePositionPole1
                if bipolar:
                    electrodePositionsPole2 = electrodePositionPole2
            else:
                electrodePositionsPole1 = np.row_stack((electrodePositionsPole1, electrodePositionPole1))
                if bipolar:
                    electrodePositionsPole2 = np.row_stack((electrodePositionsPole2, electrodePositionPole2))
        if bipolar:
            electrodePositionsAllPoles = np.row_stack((electrodePositionsPole1,electrodePositionsPole2))
        else:
            electrodePositionsAllPoles = electrodePositionsPole1

        if i == 0:
            allElectrodePositions = electrodePositionsAllPoles
        else:
            allElectrodePositions = np.row_stack((allElectrodePositions,electrodePositionsAllPoles))


    return allElectrodePositions

def get_bundle_guide_corner(bundleLength, segmentLengthAxon, overlapLength=1000):

    #length after bundle end. necessary for myelinated axons
    bundleLength = bundleLength + overlapLength

    segmentLengthBundle = segmentLengthAxon*3

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

    bundleLength = bundleLength + overlapLength

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

    #length after bundle end. necessary for myelinated axons
    bundleLength = bundleLength + overlapLength

    segmentLengthBundle = segmentLengthAxon*3
    numBundleGuideSteps = int(np.floor(bundleLength/segmentLengthBundle))

    bundleCoords = np.zeros([numBundleGuideSteps, 3])
    bundleCoords[:,0] = range(0, numBundleGuideSteps*segmentLengthBundle, segmentLengthBundle)

    return bundleCoords