from math import pi, tan, cos, sin, sqrt, floor
import numpy as np

def random_perpendicular_vectors(v):
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

        # if distance/bundleRadius > 1:
        #     print 'hm'

        # assure axon stays within bundle. If too far away -> next direction
        # equals bundle direction
        # factorAxonDirection = 1 - 1/(1+np.exp(-20*(distance/bundleRadius - 0.5)))
        # factorAxonDirection = 1 - 1/(1+np.exp(-20*(distance/bundleRadius - 0.7)))
        # factorBundleDirection = 1 - factorAxonDirection


        # factorBundleDirection = 1.5*1/(1+np.exp(-20*(distance/bundleRadius - 0.7)))
        factorBundleDirection = min((max(0,distance/bundleRadius-0.7))*6,2.5)
        #factorAxonDirection = min(0, 1 - factorBundleDirection)


        correctionVector = radiusVectorNorm + 0.1*bundleDirectionNorm
        correctionVector = correctionVector/np.linalg.norm(correctionVector)
        combinedDirection = lastAxonDirectionNorm + correctionVector*factorBundleDirection + 0.1*bundleDirection
        combinedDirectionNorm = combinedDirection/np.linalg.norm(combinedDirection)

        # get one random orthogonal vector to desired mean direction of next segment
        randomOrthogonalVectorNorm = random_perpendicular_vectors(combinedDirection)[0,:]

        # select a direction defined by cylindical coordinate rho
        rho = np.random.uniform(1)*rhoMax
        # rhoArray = rhoArray[1:-2]
        # rhoArray = np.concatenate((rhoArray, [rhoDrawn]))
        # rho = np.mean(rhoArray)
        # rho = np.random.uniform(1)*rhoMax

        randomDirection = (1-randomDirectionComponent)*combinedDirectionNorm + randomDirectionComponent*randomOrthogonalVectorNorm*rho
        # randomDirection = (1-randomDirectionComponent)*combinedDirectionNorm + factorAxonDirection*randomDirectionComponent*randomOrthogonalVectorNorm*rho
        randomDirectionNorm = randomDirection/np.linalg.norm(randomDirection)
        nextDirectionScaled = randomDirectionNorm*segmentLengthAxon

        # add the direction to the last point to obtain the next point
        nextPoint = lastPointAxon + nextDirectionScaled

        # addpend to coordinate list
        coords = np.row_stack((coords,nextPoint))

        if np.inner(bundleDirection,(nextPoint-bundleCoords[currentBundleSegment,:])) > 0:
            currentBundleSegment = currentBundleSegment + 1

    return coords

def lengthFromCoords(coords):
    # get the length of the wanted axon geometry

    # do calculate that by summing over lenghts of segments, calculate the difference in coords between each consecutive
    # pair of segments
    dCoords = np.diff(coords,axis=0)

    # pythagoras
    radicand = np.sum(np.power(dCoords,2), axis=1)
    dL = np.sqrt(radicand)

    # sum over all segments
    return sum(dL)

def electrodePositionsBundleGuided(bundleGuide, bundleRadius, numberOfElectrodes, numberOfContacts, recElectrodePositions):

    print '\nCaution, when setting the electrode positions along a non-stylized axon, the electrode position coordinates ' \
          'will be ignored.\n'

    electrodeRadius = bundleRadius*1.2

    numberOfPoles = len(recElectrodePositions)

    if numberOfPoles == 2:
        poleDistance = numberOfPoles[1] - numberOfPoles[0]
        bipolar = True
    elif numberOfPoles == 1:
        poleDistance = 0
        bipolar = False
    else:
        print 'Wrong number of recording poles.'
        return

    for i in range(numberOfElectrodes):
        segmentNumber = floor(np.shape(bundleGuide)[0]/numberOfElectrodes)*(i+1) - 1

        segmentStartingPos = bundleGuide[segmentNumber - 1,:]
        segmentEndPos = bundleGuide[segmentNumber,:]

        segmentMiddle = (segmentStartingPos + segmentEndPos)/2
        segmentOrientation = segmentStartingPos - segmentEndPos
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
            electrodPositionsAllPoles = np.row_stack((electrodePositionsPole1,electrodePositionsPole2))
        else:
            electrodPositionsAllPoles = electrodePositionsPole1

        if i == 0:
            allElectrodePositions = electrodPositionsAllPoles
        else:
            allElectrodePositions = np.row_stack((allElectrodePositions,electrodPositionsAllPoles))


    return allElectrodePositions

def getBundleGuideCorner(bundleLength, segmentLengthAxon):

    segmentLengthBundle = segmentLengthAxon*3

    numBundleGuideSteps = int(np.floor(bundleLength/segmentLengthBundle))

    halfIndex = float(numBundleGuideSteps)/2
    turningPointIndex1 = int(np.floor(halfIndex))
    turningPointIndex2 = turningPointIndex1 + int(np.ceil(halfIndex - turningPointIndex1))

    bundleCoords = np.zeros([numBundleGuideSteps, 3])
    bundleCoords[:,0] = range(0, numBundleGuideSteps*segmentLengthBundle, segmentLengthBundle)
    bundleCoords[:,1] = np.concatenate((np.zeros(turningPointIndex1),np.multiply(range(turningPointIndex2), segmentLengthBundle)))
    bundleCoords[:,2] = np.concatenate((np.zeros(turningPointIndex1),np.multiply(range(turningPointIndex2), segmentLengthBundle)))

    return bundleCoords