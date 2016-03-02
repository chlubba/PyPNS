import neuron
from neuron import h
from math import pi
import numpy as np

# import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import createGeometry

h.load_file('noload.hoc')

def axonLengthFromCoords(coords):
    # get the length of the wanted axon geometry

    # do calculate that by summing over lenghts of segments, calculate the difference in coords between each consecutive
    # pair of segments
    dCoords = np.diff(coords,axis=0)

    # pythagoras
    radicand = np.sum(np.power(dCoords,2), axis=1)
    dL = np.sqrt(radicand)

    # sum over all segments
    return sum(dL)

def getSegmentPositions(axon):

    totnsegs = 0
    for seg in axon:
        totnsegs += 1

    areavec = np.zeros(totnsegs)
    diamvec = np.zeros(totnsegs)
    lengthvec = np.zeros(totnsegs)
    xstartvec = np.zeros(totnsegs)
    xendvec = np.zeros(totnsegs)
    ystartvec = np.zeros(totnsegs)
    yendvec = np.zeros(totnsegs)
    zstartvec = np.zeros(totnsegs)
    zendvec = np.zeros(totnsegs)

    counter = 0

    #loop over all segments
    sec = axon
    n3d = int(h.n3d())
    nseg = sec.nseg
    gsen2 = 1./2/nseg
    if n3d > 0:
        #create interpolation objects for the xyz pt3d info:
        L = np.zeros(n3d)
        x = np.zeros(n3d)
        y = np.zeros(n3d)
        z = np.zeros(n3d)
        for i in range(n3d):
            L[i] = h.arc3d(i)
            x[i] = h.x3d(i)
            y[i] = h.y3d(i)
            z[i] = h.z3d(i)
        #normalize as seg.x [0, 1]
        L /= sec.L

        #temporary store position of segment midpoints
        segx = np.zeros(nseg)
        for i, seg in enumerate(sec):
            segx[i] = seg.x

        #can't be >0 which may happen due to NEURON->Python float transfer:
        segx0 = (segx - gsen2).round(decimals=6)
        segx1 = (segx + gsen2).round(decimals=6)

        #fill vectors with interpolated coordinates of start and end points
        xstartvec[counter:counter+nseg] = np.interp(segx0, L, x)
        xendvec[counter:counter+nseg] = np.interp(segx1, L, x)

        ystartvec[counter:counter+nseg] = np.interp(segx0, L, y)
        yendvec[counter:counter+nseg] = np.interp(segx1, L, y)

        zstartvec[counter:counter+nseg] = np.interp(segx0, L, z)
        zendvec[counter:counter+nseg] = np.interp(segx1, L, z)

        #fill in values area, diam, length
        for i, seg in enumerate(sec):
            areavec[counter] = h.area(seg.x)
            diamvec[counter] = seg.diam
            lengthvec[counter] = sec.L/nseg
            counter += 1
    #set self attributes
    xstart = xstartvec
    ystart = ystartvec
    zstart = zstartvec
    xend = xendvec
    yend = yendvec
    zend = zendvec
    area = areavec
    diam = diamvec
    length = lengthvec

    junctionCoordsX = np.concatenate((np.array([xstart[0]]), xend))
    junctionCoordsY = np.concatenate((np.array([ystart[0]]), yend))
    junctionCoordsZ = np.concatenate((np.array([zstart[0]]), zend))


    '''Calculate midpoints of each segment'''
    xmid = .5*(xstart+xend).flatten()
    ymid = .5*(ystart+yend).flatten()
    zmid = .5*(zstart+zend).flatten()

    return xmid, ymid, zmid, junctionCoordsX, junctionCoordsY, junctionCoordsZ

def getPT3DCoords(axon):

    nn = int(h.n3d(sec=axon))
    xx = h.Vector(nn)
    yy = h.Vector(nn)
    zz = h.Vector(nn)
    length = h.Vector(nn)
    for ii in xrange(nn):
        xx.x[ii] = h.x3d(ii,sec=axon)
        yy.x[ii] = h.y3d(ii,sec=axon)
        zz.x[ii] = h.z3d(ii,sec=axon)
        length.x[ii] = h.arc3d(ii,sec=axon)

    xCoords = np.array(xx)
    yCoords = np.array(yy)
    zCoords = np.array(zz)

    retrievedCoordinates = np.column_stack((np.transpose(xCoords), np.transpose(yCoords), np.transpose(zCoords)))
    return retrievedCoordinates

# set neuron simulation parameters
h.celsius = 33 # set temperature in celsius
h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize(-100) #65 # initialize voltage state

# define axon
axon = h.Section(name='unmyelinated_axon')
# axon.nseg = 20
# axon.L = 100
# axon.diam = 1

# insert channels. Hodgkin-Huxley needed
axon.insert('hh')

# bundle properties
bundleRadius = 10
bundleLength = 50
maximumAngle = pi/8
segmentLengthAxon = 1
numberOfAxons = 40
axonLength = bundleLength
axonDiameter = 1

bundleCoords = np.empty([bundleLength, 3])
bundleCoords[:,0] = range(bundleLength)
bundleCoords[:,1] = np.concatenate((np.zeros(bundleLength/2),range(bundleLength/2)))
bundleCoords[:,2] = np.concatenate((np.zeros(bundleLength/2),range(bundleLength/2)))

coords = createGeometry.create_random_axon(bundleCoords, bundleRadius, segmentLengthAxon)



print 'length = '+str(axonLengthFromCoords(coords))

# bundleCoords = importdata('bundleCoords.mat');
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(coords[:,0], coords[:,1], coords[:,2], label='Intended Positions')

# set them in neuron
h.pt3dclear(sec=axon)

for i in range(np.shape(coords)[0]):
    h.pt3dadd(coords[i,0], coords[i,1], coords[i,2], axonDiameter, sec=axon)

h.define_shape()



retrievedPT3DCoords = getPT3DCoords(axon)

# fig = plt.figure()
#ax = fig.gca(projection='3d')
ax.plot(retrievedPT3DCoords[:,0], retrievedPT3DCoords[:,1], retrievedPT3DCoords[:,2], label='PT3D Positions')



# now look ad actual segment midpoints
xmid, ymid, zmid, junctionCoordsX, junctionCoordsY, junctionCoordsZ = getSegmentPositions(axon)
ax.plot(xmid, ymid, zmid, label='Segment Midpoints')
ax.plot(junctionCoordsX, junctionCoordsY, junctionCoordsZ, label='Segment Junctions')



plt.legend()
plt.show()


