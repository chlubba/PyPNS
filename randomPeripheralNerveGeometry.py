from math import pi, tan, cos, sin
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


bundleRadius = 10
bundleLength = 100
maximumAngle = pi/10
segmentLengthAxon = 1

numberOfAxons = 20

axonLength = bundleLength

bundleCoords = np.empty([bundleLength, 3])
bundleCoords[:,0] = range(bundleLength) 
bundleCoords[:,1] = np.concatenate((np.zeros(bundleLength/2),range(bundleLength/2)))
bundleCoords[:,2] = np.concatenate((np.zeros(bundleLength/2),range(bundleLength/2)))

# bundleCoords = importdata('bundleCoords.mat');

fig = plt.figure()
for axonIndex in range(numberOfAxons):

    print 'calculating axon ' + str(axonIndex)

    coords = bundleCoords[0:2,:]

    rhoMax = tan(maximumAngle)*segmentLengthAxon

    currentBundleSegment = 2
    cycleCount = 0
    while currentBundleSegment <= bundleLength-1:

        # get last point
        lastPointAxon = coords[-1,:]
        # and vector between two last points
        lastAxonDirection = coords[-1,:] - coords[-2,:]

        # current bundle direction
        bundleDirection = bundleCoords[currentBundleSegment,:] - bundleCoords[currentBundleSegment-1,:]

        # get orthogonal space to current direction vector

        orthogonalVectorsBundle = random_perpendicular_vectors(bundleDirection)

        # distance between extended axon segment and extendet bundle segment
        # project distance between axon segment endpoint and bundle segment
        # endpoint into orthogonal space
        # As bundle direction is irrelevant for this distance, only add axon
        # direction.
        distance = np.linalg.norm(np.dot(((lastPointAxon + lastAxonDirection) - bundleCoords[currentBundleSegment,:])
                                  ,np.transpose(orthogonalVectorsBundle)))

        # assure axon stays within bundle. If too far away -> next direction
        # equals bundle direction
        factorAxonDirection = 1 - 1/(1+np.exp(-20*(distance/bundleRadius - 0.5)))
        factorBundleDirection = 1 - factorAxonDirection

        combinedDirection = lastAxonDirection*factorAxonDirection + \
                            bundleDirection*factorBundleDirection


        # get one random orthogonal vector to desired mean direction of next segment
        randomOrthogonalVectorNorm = random_perpendicular_vectors(combinedDirection)[0,:]

        # select a direction defined by cylindical coordinate rho
        rho = np.random.uniform(1)*rhoMax
        randomDirection = combinedDirection + randomOrthogonalVectorNorm*rho
        randomDirectionNorm = randomDirection/np.linalg.norm(randomDirection)
        nextDirectionScaled = randomDirectionNorm*segmentLengthAxon

        # add the direction to the last point to obtain the next point
        nextPoint = lastPointAxon + nextDirectionScaled

        # addpend to coordinate list
        coords = np.row_stack((coords,nextPoint))

        if np.inner(bundleDirection,(nextPoint-bundleCoords[currentBundleSegment,:])) > 0:
            currentBundleSegment = currentBundleSegment + 1



    ax = fig.gca(projection='3d')
    ax.plot(coords[:,0], coords[:,1], coords[:,2], label='singleAxon')
#ax.plot(bundleCoords[:,0], bundleCoords[:,1], bundleCoords[:,2], label='Bundle center')
ax.legend()

plt.show()