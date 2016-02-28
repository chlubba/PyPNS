import neuron
from neuron import h
# import neuron.gui
import numpy as np
from matplotlib import pyplot
h.load_file('noload.hoc')

# set neuron simulation parameters
h.celsius = 33 # set temperature in celsius
h.tstop = 3e2 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
#h.finitialize(-100) #65 # initialize voltage state


# define axon
axonTest = h.Section(name='unmyelinated_axon')
axonTest.nseg = 20
axonTest.L = 1000
axonTest.diam = 1


# insert channels. Hodgkin-Huxley needed
axonTest.insert('hh')

# here the problems arise

# extracellular
axonTest.insert('extracellular')

# configure ExpSyn synapse
stim = h.ExpSyn(1e-3, axonTest)
stim.e = 10
stim.i = 0.2
stim.tau = 3

# create random spike train


# configure input to synapse
vecStim = h.VecStim()
vec = h.Vector([5, 25])#[1, 2])
vecStim.play(vec)

# connect synapse and VecStim input
netCon = h.NetCon(vecStim, stim)
netCon.weight[0] = 1

#observe synapse current
iSynVec = h.Vector()
iSynVec.record(stim._ref_i)

# and observe axon potentials
vAxonVec0 = h.Vector()
vAxonVec0.record(axonTest(0)._ref_v)
vAxonVec1 = h.Vector()
vAxonVec1.record(axonTest(1)._ref_v)

# record time, too
timeVector = h.Vector()
timeVector.record(h._ref_t)

# run simulation
h.finitialize() #65 # initialize voltage state
h.run()

# convert recorded signals to numpy arrays
timeArray = np.array(timeVector)
iSynArr = np.array(iSynVec)
vAxonArr0 = np.array(vAxonVec0)
vAxonArr1 = np.array(vAxonVec1)

# plot synapse currents and axon potentials
pyplot.plot(timeVector,vAxonArr0, label='Axon potential at start-position')
pyplot.plot(timeVector,vAxonArr1, label='Axon potential at end-position')
pyplot.plot(timeVector,iSynArr, label='Synapse current')
pyplot.legend()
pyplot.title('Without extracellular mechanism')
pyplot.show()