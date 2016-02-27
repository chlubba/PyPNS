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

    coords = bundleCoords[0:2,:]

    rhoMax = tan(maximumAngle)*segmentLengthAxon
    rhoArray = np.zeros(5)

    currentBundleSegment = 2
    cycleCount = 0
    while currentBundleSegment <= bundleLength-1:

        # get last point
        lastPointAxon = coords[-1,:]
        # and vector between two last points
        lastAxonDirection = (coords[-1,:] - coords[-2,:])/segmentLengthAxon

        # current bundle direction
        bundleDirection = bundleCoords[currentBundleSegment,:] - bundleCoords[currentBundleSegment-1,:]
        lastPointBundle = bundleCoords[currentBundleSegment-1,:]

        # get orthogonal space to current direction vector

        orthogonalVectorsBundle = random_perpendicular_vectors(bundleDirection)

        cp = np.inner(bundleDirection, lastPointAxon-lastPointBundle)/np.linalg.norm(bundleDirection)
        if cp == 0:
            radiusVectorNorm = [0, 0, 0]
            distance = 0
        else:
            radiusVector = -(lastPointAxon - (cp*bundleDirection + lastPointBundle))
            distance = np.linalg.norm(radiusVector)
            radiusVectorNorm = radiusVector/distance


        # distance between extended axon segment and extendet bundle segment
        # project distance between axon segment endpoint and bundle segment
        # endpoint into orthogonal space
        # As bundle direction is irrelevant for this distance, only add axon
        # direction.
        # distance1 = np.linalg.norm(np.dot(((lastPointAxon + lastAxonDirection) - bundleCoords[currentBundleSegment,:])
        #                          ,np.transpose(orthogonalVectorsBundle)))

        # assure axon stays within bundle. If too far away -> next direction
        # equals bundle direction
        factorAxonDirection = 1 - 1/(1+np.exp(-20*(distance/bundleRadius - 0.5)))
        factorBundleDirection = 1 - factorAxonDirection

        # combinedDirection = lastAxonDirection*factorAxonDirection + \
        #                     bundleDirection*factorBundleDirection
        # if max((distance - bundleRadius),0) > 0:
        #     print 'ooh'
        #combinedDirection = lastAxonDirection + radiusVectorNorm*0.2*np.exp(0.01*(distance - bundleRadius)) - min(0,np.inner(lastAxonDirection, bundleDirection))*bundleDirection
        relativeDistance = (distance - 0.5*bundleRadius)
        if relativeDistance > 0:
            correctionFactor = 0.01*relativeDistance#**2
        else:
            correctionFactor = 0

        # if currentBundleSegment == 60:
        #     print 'wos'
        if distance > 15:
            print 'check.'

        correctionVector = radiusVectorNorm + 4*bundleDirection
        correctionVector = correctionVector/np.linalg.norm(correctionVector)
        #combinedDirection = lastAxonDirection + radiusVectorNorm*correctionFactor - min(0,np.inner(lastAxonDirection, bundleDirection))*bundleDirection
        #combinedDirection = lastAxonDirection + radiusVectorNorm*2*max(0,factorBundleDirection-0.5) - (min(-0.2,np.inner(lastAxonDirection, bundleDirection))+0.2)*bundleDirection
        combinedDirection = lastAxonDirection + correctionVector*factorBundleDirection + 0.1*bundleDirection

        # get one random orthogonal vector to desired mean direction of next segment
        randomOrthogonalVectorNorm = random_perpendicular_vectors(combinedDirection)[0,:]

        # select a direction defined by cylindical coordinate rho
        rhoDrawn = np.random.uniform(1)*rhoMax
        rhoArray = rhoArray[1:-2]
        rhoArray = np.concatenate((rhoArray, [rhoDrawn]))
        rho = np.mean(rhoArray)
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



    ax = fig.gca(projection='3d')
    ax.plot(coords[:,0], coords[:,1], coords[:,2], label='singleAxon')

#ax.plot(bundleCoords[:,0], bundleCoords[:,1], bundleCoords[:,2], label='Bundle center')
ax.legend()

plt.show()