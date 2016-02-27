from math import pi, tan, cos, sin, sqrt, floor
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
# import scipy as sp

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

def create_random_axon(bundleCoords, bundleRadius, segmentLengthAxon, diameter = 1, maximumAngle = 0.314):

    coords = bundleCoords[0:2,:]

    rhoMax = tan(maximumAngle)*segmentLengthAxon
    rhoArray = np.zeros(5)

    currentBundleSegment = 2
    while currentBundleSegment <= bundleLength-1:

        # get last point
        lastPointAxon = coords[-1,:]
        # and vector between two last points
        lastAxonDirection = (coords[-1,:] - coords[-2,:])/segmentLengthAxon

        # current bundle direction
        bundleDirection = bundleCoords[currentBundleSegment,:] - bundleCoords[currentBundleSegment-1,:]
        lastPointBundle = bundleCoords[currentBundleSegment-1,:]

        # get orthogonal space to current direction vector

        cp = np.inner(bundleDirection, lastPointAxon-lastPointBundle)/np.linalg.norm(bundleDirection)
        if cp == 0:
            radiusVectorNorm = [0, 0, 0]
            distance = 0
        else:
            radiusVector = -(lastPointAxon - (cp*bundleDirection + lastPointBundle))
            distance = np.linalg.norm(radiusVector)
            radiusVectorNorm = radiusVector/distance


        # assure axon stays within bundle. If too far away -> next direction
        # equals bundle direction
        factorAxonDirection = 1 - 1/(1+np.exp(-20*(distance/bundleRadius - 0.5)))
        factorBundleDirection = 1 - factorAxonDirection


        correctionVector = radiusVectorNorm + 4*bundleDirection
        correctionVector = correctionVector/np.linalg.norm(correctionVector)
        combinedDirection = lastAxonDirection + correctionVector*factorBundleDirection + 0.1*bundleDirection

        # get one random orthogonal vector to desired mean direction of next segment
        randomOrthogonalVectorNorm = random_perpendicular_vectors(combinedDirection)[0,:]

        # select a direction defined by cylindical coordinate rho
        rho = np.random.uniform(1)*rhoMax
        # rhoArray = rhoArray[1:-2]
        # rhoArray = np.concatenate((rhoArray, [rhoDrawn]))
        # rho = np.mean(rhoArray)
        # rho = np.random.uniform(1)*rhoMax
        randomDirection = combinedDirection + randomOrthogonalVectorNorm*rho
        randomDirectionNorm = randomDirection/np.linalg.norm(randomDirection)
        nextDirectionScaled = randomDirectionNorm*segmentLengthAxon

        # add the direction to the last point to obtain the next point
        nextPoint = lastPointAxon + nextDirectionScaled

        # addpend to coordinate list
        coords = np.row_stack((coords,nextPoint))

        if np.inner(bundleDirection,(nextPoint-bundleCoords[currentBundleSegment,:])) > 0:
            currentBundleSegment = currentBundleSegment + 1

    return coords


bundleRadius = 10
bundleLength = 50
maximumAngle = pi/8
segmentLengthAxon = 1

numberOfAxons = 40

axonLength = bundleLength

bundleCoords = np.empty([bundleLength, 3])
bundleCoords[:,0] = range(bundleLength) 
bundleCoords[:,1] = np.concatenate((np.zeros(bundleLength/2),range(bundleLength/2)))
bundleCoords[:,2] = np.concatenate((np.zeros(bundleLength/2),range(bundleLength/2)))

# bundleCoords = importdata('bundleCoords.mat');

fig = plt.figure()
for axonIndex in range(numberOfAxons):

    print 'calculating axon ' + str(axonIndex)


    coords = create_random_axon(bundleCoords, bundleRadius, segmentLengthAxon)

    ax = fig.gca(projection='3d')
    ax.plot(coords[:,0], coords[:,1], coords[:,2], label='singleAxon')

electrodePositions = np.empty([1,3])

numberOfElectrodes = 10
numberOfContacts = 8


for i in range(numberOfElectrodes):
    segmentNumber = floor(np.shape(bundleCoords)[0]/numberOfElectrodes)*(i+1) - 1

    segmentStartingPos = bundleCoords[segmentNumber - 1,:]
    segmentEndPos = bundleCoords[segmentNumber,:]

    segmentMiddle = (segmentStartingPos + segmentEndPos)/2
    segmentOrientation = segmentStartingPos - segmentEndPos

    # get one random orthogonal vector
    orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0,:]

    for j in range(numberOfContacts):
        electrodePosition = np.dot(rotation_matrix(segmentOrientation, 2*pi/numberOfContacts*j),(orthogonalVector*bundleRadius)) + segmentMiddle

        if i == 0 and j == 0: # how to get a truly empty array?!
            electrodePositions = electrodePosition
        else:
            electrodePositions = np.row_stack((electrodePositions, electrodePosition))


ax.scatter(electrodePositions[:,0], electrodePositions[:,1], electrodePositions[:,2])
for i in range(np.shape(electrodePositions)[0]/numberOfContacts):
    indices = i*8 + np.concatenate(([8],range(1,9))) - 1
    ax.plot(electrodePositions[indices,0], electrodePositions[indices,1], electrodePositions[indices,2], color = '0.75', label='Electrodes')

ax.plot(bundleCoords[:,0], bundleCoords[:,1], bundleCoords[:,2], label='Bundle center', linewidth=5.0)
ax.legend()

plt.show()