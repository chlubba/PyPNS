import neuron
from neuron import h
from math import pi
import numpy as np

# import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import createGeometry

h.load_file('noload.hoc')

# set neuron simulation parameters
h.celsius = 33 # set temperature in celsius
h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize(-100) #65 # initialize voltage state

# define axon
axon = h.Section(name='unmyelinated_axon')
axon.nseg = 20
axon.L = 100
axon.diam = 1

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

# # bundleCoords = importdata('bundleCoords.mat');
#
# ax = fig.gca(projection='3d')
# ax.plot(coords[:,0], coords[:,1], coords[:,2], label='singleAxon')

# # set them in neuron
# h.pt3dclear(sec=axon)
#
# for i in range(np.shape(coords)[0]):
#     h.pt3dadd(coords[i,0], coords[i,1], coords[i,2], axonDiameter, sec=axon)

h.define_shape()

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

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(retrievedCoordinates[:,0], retrievedCoordinates[:,1], retrievedCoordinates[:,2], label='singleAxon')
plt.show()



